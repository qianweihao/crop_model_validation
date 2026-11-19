"""Step 07 - calibrate cultivar parameters against the observed maize data.

Two stages, in this order (phenology first, otherwise the growth parameters
absorb timing errors):

  A. phenology   - Juvenile / FlagLeafToFlowering / FloweringToGrainFilling /
                   GrainFilling / Maturing thermal-time targets, fitted to the
                   observed 吐丝 and harvest dates of each variety.
  B. growth      - RUE and grain-set parameters, fitted to observed LAI,
                   biomass and yield across the plots of that variety.

Cultivars are injected with the mechanism verified in apsim_override.py
(a ``Models.Core.Folder`` + ``Models.PMF.Cultivar`` with ``.FixedValue`` paths).
Calibrated values land in analysis/out/calibrated_cultivars.json.
"""
from __future__ import annotations

import argparse
import json
import time
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_override import SimCase, run, run_many, stage_dates  # noqa: E402
from common import APSIM_DIR, DS, SIM  # noqa: E402

# ---------------------------------------------------------------------------
# grain parameters: upper bounds are capped so the optimiser cannot "buy" yield
# with an unrealistically large sink (previously MaxGrainsPerCob hit 1200 and
# produced HI = 0.70 against an observed 0.50).
# ---------------------------------------------------------------------------
YIELD_WEIGHT = 1.0       # dry-grain yield, relative error
BIOMASS_WEIGHT = 1.0     # above-ground biomass series, relative error
LAI_WEIGHT = 0.5         # LAI series, relative error
HI_WEIGHT = 8.0          # harvest index constrains sink size (see notes)

WORK = SIM / "calib_work"

# observed emergence / silking / harvest per variety (from 品种生育期_汇总 + harvest dates)
OBS_PHENO = {
    "乐农87":  dict(silking="2026-08-09", harvest="2026-10-08", leaf_stages={"4叶期": "2026-07-07", "12叶期": "2026-07-28"}),
    "MY73":    dict(silking="2026-08-08", harvest="2026-10-04", leaf_stages={"4叶期": "2026-07-05", "12叶期": "2026-07-26"}),
    "金华创826": dict(silking="2026-08-10", harvest="2026-10-04", leaf_stages={"4叶期": "2026-07-05", "12叶期": "2026-07-28"}),
}

PHENO_PARAMS = {
    "Phenology.Juvenile.Target.FixedValue": (80.0, 520.0),
    "Phenology.FlagLeafToFlowering.Target.FixedValue": (10.0, 220.0),
    "Phenology.FloweringToGrainFilling.Target.FixedValue": (40.0, 320.0),
    "Phenology.GrainFilling.Target.FixedValue": (250.0, 1000.0),
    "Phenology.Maturing.Target.FixedValue": (0.0, 80.0),
}

GROWTH_PARAMS = {
    # RUE must reach ~4.4: the field crop reaches 18 t/ha at LAI ~3.1, i.e.
    # thicker leaves and faster conversion than APSIM's maize defaults.
    "Leaf.Photosynthesis.RUE.FixedValue": (2.2, 4.8),
    "Grain.MaximumGrainsPerCob.FixedValue": (400.0, 950.0),
    "Grain.MaximumPotentialGrainSize.FixedValue": (0.15, 0.35),
    # canopy shape: where the widest leaf sits (lower -> smaller peak LAI)
    "Leaf.CohortParameters.MaxArea.LargestLeafPosition.RelativePositionOfLargestLeaf.FixedValue": (0.35, 0.70),
    # leaf structural fraction controls specific leaf weight (~68 g/m2 observed,
    # APSIM default gives ~30 g/m2, so leaves are modelled far too thin)
    "Leaf.CohortParameters.StructuralFraction.FixedValue": (0.30, 0.85),
    # senescence acceleration: with the default the canopy collapses by mid-August
    # (LAI 1.2 against an observed 3.1) - the crop stays green much longer
    "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue": (0.03, 1.20),
}


def plots_of(trial_id: str, variety: str) -> list[Path]:
    out = []
    for p in sorted(APSIM_DIR.glob(f"{trial_id}_*.apsimx")):
        case = SimCase.load(p)
        sow = next(c for c in case._field()["Children"]
                   if c.get("Name", "").startswith("Sow"))
        # variety is carried in the simulation name / source table; match by trial
        out.append(p)
    return out


def sim_case(path: Path, commands: dict) -> SimCase:
    return SimCase.load(path).set_commands(commands)


def eval_phenology(plot: Path, params: dict) -> dict:
    with tempfile.TemporaryDirectory() as td:
        df = run(sim_case(plot, params), td, "ph")
    return stage_dates(df)


def eval_phenology_batch(plot: Path, populations: list[dict]) -> list[dict]:
    """One parallel batch per DE generation (bounded by APSIM_WORKERS)."""
    jobs = [(sim_case(plot, p), f"ph{i:03d}") for i, p in enumerate(populations)]
    with tempfile.TemporaryDirectory() as td:
        tables = run_many(jobs, td)
    return [stage_dates(tables[tag]) for _, tag in jobs]


# --------------------------------------------------------------------------
def fit_phenology(variety: str, plot: Path, maxiter: int, popsize: int) -> dict:
    obs = OBS_PHENO[variety]
    sil = pd.Timestamp(obs["silking"])
    har = pd.Timestamp(obs["harvest"])

    def score(p: dict, s: dict) -> float:
        pen = 0.0
        if s["Flowering"] is None:
            pen += 500.0
        else:
            pen += (s["Flowering"].normalize() - sil).days ** 2
        if s["Maturity"] is None:
            pen += 500.0
        else:
            pen += (s["Maturity"].normalize() - har).days ** 2
        gf = p["Phenology.GrainFilling.Target.FixedValue"]
        if gf < 350:                      # keep grain filling physiological
            pen += (350 - gf) ** 2
        return pen

    def obj(x):
        p = {k: float(v) for k, v in zip(PHENO_PARAMS, x)}
        try:
            return score(p, eval_phenology(plot, p))
        except Exception:
            return 1e6

    bounds = list(PHENO_PARAMS.values())

    class _Batcher:
        """map-like callable for scipy: one parallel APSIM batch per generation."""
        def __call__(self, func, iterable):
            pops = [{k: float(v) for k, v in zip(PHENO_PARAMS, x)} for x in iterable]
            try:
                sims = eval_phenology_batch(plot, pops)
            except Exception as exc:
                print(f"  [batch failed: {type(exc).__name__}: {exc}]", flush=True)
                return [1e6] * len(pops)
            vals = [score(p, s) for p, s in zip(pops, sims)]
            _Batcher.gen += 1
            print(f"  gen {_Batcher.gen:3d}: n={len(pops)} best={min(vals):.1f} "
                  f"({time.time() - _Batcher.t0:.0f}s)", flush=True)
            return vals

        gen = 0
        t0 = time.time()

    res = differential_evolution(obj, bounds, maxiter=maxiter, popsize=popsize,
                                 tol=1e-6, seed=1, polish=False,
                                 init="sobol", updating="deferred",
                                 workers=_Batcher())
    best = {k: round(float(v), 1) for k, v in zip(PHENO_PARAMS, res.x)}
    s = eval_phenology(plot, best)
    return dict(params=best, score=float(res.fun), sim=s,
                obs=dict(silking=str(sil.date()), harvest=str(har.date())))


# --------------------------------------------------------------------------
def obs_targets(variety: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Observed series for calibration.

    Uses plot-scale LAI / biomass where a measured density exists.  For the
    Xinxiang irrigation trial only per-plant measurements were taken, so the
    series is expressed per plant (density_measured is NaN there); the
    simulation is then compared on the same per-plant basis.
    """
    ob = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    ob = ob[ob["variety"] == variety]

    # NOTE: do not require lai_plot here - the harvest date has no leaf-area
    # measurement but carries the peak biomass used for the harvest index.
    plot_scale = ob.dropna(subset=["density_measured", "biomass_kg_mu"])
    if len(plot_scale):
        plot_scale = plot_scale.assign(biomass_kg_ha=plot_scale["biomass_kg_mu"] * 15.0)
        plot_scale.attrs["scale"] = "plot"
    else:
        per_plant = ob.dropna(subset=["plant_la_cm2", "total_dw"])
        per_plant = per_plant.assign(
            # per-plant leaf area (cm2) and dry weight (g) stay on their own scale
            lai_plot=per_plant["plant_la_cm2"],
            biomass_kg_ha=per_plant["total_dw"])
        per_plant.attrs["scale"] = "plant"
        plot_scale = per_plant
    y = pd.read_csv(DS / "yield.csv")
    return plot_scale, y


def eval_growth(plot: Path, params: dict, obs: pd.DataFrame, t0: pd.Timestamp) -> dict:
    with tempfile.TemporaryDirectory() as td:
        df = run(sim_case(plot, params), td, "gr")
    return _interp(df, obs, t0)


def _interp(df: pd.DataFrame, obs: pd.DataFrame, t0: pd.Timestamp,
            scale: str = "plot", density: float | None = None) -> dict:
    """Interpolate the simulation onto the observation dates.

    scale='plant' converts LAI [m2/m2] -> leaf area per plant [cm2] and
    above-ground weight [kg/ha] -> [g/plant] using the planting density, so the
    Xinxiang irrigation trial (per-plant records only) can be used too.
    """
    sim = df.set_index("Clock.Today")
    xs = [(d - t0).days for d in sim.index]
    xq = [(d - t0).days for d in obs["date"]]
    lai = np.interp(xq, xs, sim["Maize.Leaf.LAI"])
    bio = np.interp(xq, xs, sim["Maize.AboveGround.Wt"] * 10.0)
    lai = np.where(np.isfinite(obs["lai_plot"].to_numpy(dtype=float)), lai, np.nan)
    if scale == "plant" and density:
        mu_m2 = 666.667
        lai = lai * mu_m2 / density * 1e4          # m2/m2 -> cm2/plant
        bio = bio / (density * 15.0) * 1000.0      # kg/ha -> g/plant
    return dict(lai=lai, bio=bio,
                yield_kg_ha=float(sim["Maize.Grain.Total.Wt"].max() * 10.0),
                ag_max=float(sim["Maize.AboveGround.Wt"].max() * 10.0))


def _growth_targets(variety: str, plots: list[Path]) -> tuple[dict, dict]:
    """obs LAI/biomass series per plot, plus the harvested yield per plot."""
    ob, yld = obs_targets(variety)
    ob = ob.assign(sid=ob["trial_id"] + "_" + ob["plot"])
    y = yld.assign(sid=yld["trial_id"] + "_" + yld["plot"])
    # observed yield is at 14% moisture; APSIM reports dry grain -> compare on dry
    ycol = "yield_dry_kg_ha" if "yield_dry_kg_ha" in y.columns else "yield_kg_ha"
    ymap = {r["sid"]: float(r[ycol]) for _, r in y.iterrows()
            if pd.notna(r.get(ycol))}
    series, yields = {}, {}
    for p in plots:
        sid = p.stem
        # keep every row with a biomass value: the harvest date (peak biomass)
        # has no leaf-area measurement, so LAI is masked per element below
        sub = ob[ob["sid"] == sid].sort_values("date").copy()
        sub = sub[sub["biomass_kg_ha"].notna()]
        if len(sub):
            sub.attrs["scale"] = ob.attrs.get("scale", "plot")
            series[sid] = sub
            if sid in ymap and ymap[sid] > 0:
                yields[sid] = ymap[sid]
    return series, yields


def fit_growth(variety: str, pheno: dict, plots: list[Path], maxiter: int,
               popsize: int, seed_json: str | None = None) -> dict:
    series, yields = _growth_targets(variety, plots)
    if not series:
        raise SystemExit(f"no observed LAI/biomass for {variety}")
    t0 = pd.Timestamp("2026-06-08")
    keys = list(GROWTH_PARAMS)

    # observed harvest index per plot: dry grain / peak above-ground dry biomass
    hi_rows = pd.read_csv(DS / "yield.csv")
    hi_rows = hi_rows.assign(sid=hi_rows["trial_id"] + "_" + hi_rows["plot"])
    obs_hi = {r["sid"]: float(r["yield_dry_kg_ha"]) for _, r in hi_rows.iterrows()
              if pd.notna(r.get("yield_dry_kg_ha"))}
    t_hi = {}
    for sid, t in series.items():
        ag = float(t["biomass_kg_ha"].max())
        if sid in obs_hi and ag > 0:
            t_hi[sid] = obs_hi[sid] / ag

    plots_df = pd.read_csv(DS / "plots.csv")
    dens = {(r["trial_id"] + "_" + r["plot"]): float(r["density_planted"])
            for _, r in plots_df.iterrows()}
    scale_of = {sid: ("plot" if "lai_plot" in t.columns
                      and t["lai_plot"].notna().any() else "plant")
                for sid, t in series.items()}
    # obs_targets() already put the right quantity in lai_plot / biomass_kg_ha
    scale_of = {sid: getattr(t, "attrs", {}).get("scale", "plot")
                for sid, t in series.items()}

    def score(g: dict, sims: dict) -> float:
        """Dimensionless relative error with an explicit harvest-index term.

        Without the HI term the optimiser satisfies yield by inflating the grain
        sink and simultaneously under-producing biomass, giving HI = 0.70 while
        the field data say 0.50.  HI ties the two together.
        """
        err = []
        for sid, t in series.items():
            if sid not in sims:
                return 1e6
            s = sims[sid]
            lai_o = t["lai_plot"].to_numpy(dtype=float)
            bio_o = t["biomass_kg_ha"].to_numpy(dtype=float)
            lai_s = np.asarray(s["lai"], dtype=float)
            bio_s = np.asarray(s["bio"], dtype=float)
            m_lai = np.isfinite(lai_o) & np.isfinite(lai_s)
            m_bio = np.isfinite(bio_o) & np.isfinite(bio_s)
            if m_lai.any() and lai_o[m_lai].mean() > 0:
                ref = lai_o[m_lai].mean()
                err += list(LAI_WEIGHT * ((lai_s[m_lai] - lai_o[m_lai]) / ref) ** 2)
            if m_bio.any() and bio_o[m_bio].mean() > 0:
                ref = bio_o[m_bio].mean()
                err += list(BIOMASS_WEIGHT * ((bio_s[m_bio] - bio_o[m_bio]) / ref) ** 2)

        for sid, yo in yields.items():
            if sid in sims and yo > 0:
                err.append(YIELD_WEIGHT * ((sims[sid]["yield_kg_ha"] - yo) / yo) ** 2)
                # harvest index: dry grain / above-ground dry biomass
                hi_o = t_hi.get(sid)
                ag = sims[sid].get("ag_max", np.nan)
                if hi_o is not None and np.isfinite(ag) and ag > 0:
                    hi_s = sims[sid]["yield_kg_ha"] / ag
                    err.append(HI_WEIGHT * (hi_s - hi_o) ** 2)
        return float(np.mean(err)) if err else 1e6

    class _Batcher:
        gen = 0
        t0 = time.time()

        def __call__(self, func, iterable):
            pops = [{k: float(v) for k, v in zip(keys, x)} for x in iterable]
            jobs, tags = [], []
            for i, g in enumerate(pops):
                params = {**pheno, **g}
                for sid in series:
                    tags.append((i, sid))
                    jobs.append((sim_case(APSIM_DIR / f"{sid}.apsimx", params),
                                 f"g{i:03d}_{sid}"))
            try:
                with tempfile.TemporaryDirectory() as td:
                    tables = run_many(jobs, td)
            except Exception as exc:
                print(f"  [batch failed: {type(exc).__name__}: {exc}]", flush=True)
                return [1e6] * len(pops)

            sims: list[dict] = [dict() for _ in pops]
            for (i, sid), (_, tag) in zip(tags, jobs):
                df = tables.get(tag)
                if df is None:          # APSIM rejected this parameter set
                    sims[i] = {}
                    break
                sims[i][sid] = _interp(df, series[sid], t0,
                                       scale=scale_of.get(sid, "plot"),
                                       density=dens.get(sid))
            vals = [score(g, sims[i]) for i, g in enumerate(pops)]
            _Batcher.gen += 1
            print(f"  gen {_Batcher.gen:3d}: n={len(pops)} best={min(vals):.4f} "
                  f"({time.time() - _Batcher.t0:.0f}s)", flush=True)
            return vals

    init = "sobol"
    if seed_json:
        guess = json.loads(Path(seed_json).read_text(encoding="utf-8"))
        base = np.array([float(guess[k]) for k in keys])
        lo = np.array([b[0] for b in GROWTH_PARAMS.values()])
        hi = np.array([b[1] for b in GROWTH_PARAMS.values()])
        rng = np.random.default_rng(7)
        n_pop = popsize * len(keys)
        init = np.clip(base + rng.normal(0, 0.05, (n_pop, len(keys))) * (hi - lo), lo, hi)
        init[0] = base
        print(f"  seeded from {seed_json}")

    res = differential_evolution(lambda x: 1e6, list(GROWTH_PARAMS.values()),
                                 maxiter=maxiter, popsize=popsize, tol=1e-4,
                                 seed=2, polish=False, init=init,
                                 workers=_Batcher())
    return dict(params={k: round(float(v), 3) for k, v in zip(keys, res.x)},
                score=float(res.fun))


def select_spanning_plots(mine: pd.DataFrame, k: int) -> list[Path]:
    """Pick k plots spread along the N gradient (and across trials) for calibration."""
    fert = pd.read_csv(DS / "fert.csv")
    n_tot = fert.groupby(["trial_id", "plot"])["n_kg_ha"].sum()
    m = mine.copy()
    m["n_kg_ha"] = [n_tot.get((r["trial_id"], r["plot"]), 0.0)
                    for _, r in m.iterrows()]
    m = m.sort_values(["trial_id", "n_kg_ha"]).reset_index(drop=True)
    if len(m) <= k:
        picked = m
    else:
        idx = np.unique(np.linspace(0, len(m) - 1, k).round().astype(int))
        picked = m.iloc[idx]
    return [APSIM_DIR / f"{r['trial_id']}_{r['plot']}.apsimx"
            for _, r in picked.iterrows()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variety", default=None)
    ap.add_argument("--stage", choices=["phenology", "growth", "both"], default="both")
    ap.add_argument("--maxiter", type=int, default=6)
    ap.add_argument("--popsize", type=int, default=6)
    ap.add_argument("--plots", type=int, default=3)
    ap.add_argument("--seed-json", default=None,
                    help="JSON of an initial parameter guess to seed the search")
    a = ap.parse_args()

    plots = pd.read_csv(DS / "plots.csv")
    out_path = SIM / "calibrated_cultivars.json"
    store = json.loads(out_path.read_text(encoding="utf-8")) if out_path.exists() else {}

    varieties = [a.variety] if a.variety else list(OBS_PHENO)
    for var in varieties:
        mine = plots[(plots["variety"] == var)
                     & plots["trial_id"].isin(["SD_SF", "XX_SF", "XX_IRR"])]
        if mine.empty:
            print(f"skip {var}: no simulable plots")
            continue
        ref = APSIM_DIR / f"{mine.iloc[0]['trial_id']}_{mine.iloc[0]['plot']}.apsimx"
        entry = store.setdefault(var, {})

        if a.stage in ("phenology", "both"):
            r = fit_phenology(var, ref, a.maxiter, a.popsize)
            entry["phenology"] = r["params"]
            print(f"[{var}] phenology score={r['score']:.1f}  sim={r['sim']}  obs={r['obs']}")

        if a.stage in ("growth", "both"):
            ph = entry.get("phenology")
            pheno_file = SIM / "phenology_calibrated.json"
            if pheno_file.exists():
                pc = json.loads(pheno_file.read_text(encoding="utf-8"))
                if var in pc:
                    ph = pc[var]["params"]
                    print(f"[{var}] using leaf-constrained phenology")
            if not ph:
                print(f"skip growth for {var}: run phenology first")
                continue
            ps = select_spanning_plots(mine, a.plots)
            g = fit_growth(var, ph, ps, a.maxiter, a.popsize, a.seed_json)
            entry["growth"] = g["params"]
            print(f"[{var}] growth    score={g['score']:.3f}  {g['params']}")

        store[var] = entry
        out_path.write_text(json.dumps(store, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
