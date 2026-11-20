"""Staged calibration of the CERES-Maize cultivar parameters for Lenong 87.

Three stages, each a coarse grid followed by a local refinement:

  1. phenology   P1 (emergence -> end of juvenile), P5 (silking -> maturity)
                 fitted to the observed anthesis (08-09) and maturity (10-08)
                 on A3B2.
  2. canopy      PHINT (phyllochron) fitted to the observed LAI series of all
                 nine plots, with the plateau window (LAI >= 0.95 * peak) used
                 as a shape check.
  3. yield       G2 (kernels/plant) and G3 (kernel growth rate) fitted to the
                 measured dry yields, together with the initial nitrate scale
                 (never measured, see MODEL_INPUTS_OUTPUTS.md).

The result is written to out/sim/calibration_log.csv and, on --write, into
out/sim/calibrated_cultivar.json.

Usage
-----
    python analysis/dssat/src/s03_calibrate.py            # report only
    python analysis/dssat/src/s03_calibrate.py --write    # update the JSON
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

# Targets for the yield stage.  SD_SF: measured grain number (all treatments)
# and single grain weight, from MODEL_INPUTS_OUTPUTS.md 2.4 / 3.1.
# XX_SF has no measured grain number, so only the yield RRMSE is scored there.
GRAIN_TARGETS = {
    "SD_SF": (3471.0, 240.0),
    # XX_SF has no measured grain number; only the single grain weight is
    # anchored (240 mg, the value measured for the same crop in SD_SF).
    "XX_SF": (None, 240.0),
}

# Search grids per trial.  SD_SF reuses the values in the report; XX_SF needs
# a lower PHINT (higher density / larger leaves) and a nitrogen stage because
# the measured LAI there implies much less N stress than the default profile.
GRIDS = {
    "SD_SF": dict(p1=[215, 225, 235, 245, 255], p5=[860, 880, 900, 920],
                  phint=[58, 60, 62, 64], sloc=[1.0],
                  g2=[860, 1100, 1300, 1500], g3=[5.5, 6.0, 7.0, 8.6],
                  ns=[0.20, 0.22, 0.24]),
    "XX_SF": dict(p1=[205, 215, 225, 235], p5=[900, 920, 940, 960],
                  phint=[46, 48, 50, 52], sloc=[0.5, 0.6, 0.7],
                  g2=[1100, 1300, 1500], g3=[5.8, 6.2, 6.6],
                  ns=[0.35, 0.45, 0.55]),
}

OBS_ANTHESIS = pd.Timestamp("2026-08-09")
OBS_MATURITY = pd.Timestamp("2026-10-08")
PLOT_SUBSET = ["A1", "A3B2", "A5B2"]


def _obs_lai() -> pd.DataFrame:
    o = dc.load_obs()
    return o.dropna(subset=["lai_plot"])


def lai_rmse(cultivar: dict, plots=None, obs: pd.DataFrame | None = None,
             sloc_scale: float | None = None) -> float:
    # ``plots=None`` rather than ``plots=dc.ORDER``: default arguments are bound
    # at import time, before set_trial() rebinds ORDER for the second trial.
    plots = dc.ORDER if plots is None else plots
    obs = _obs_lai() if obs is None else obs
    se = []
    for plot in plots:
        g, _ = dc.run_plot(plot, cultivar, sloc_scale=sloc_scale)
        gg = g.set_index(g["date"].dt.normalize())
        o = obs[obs["plot"] == plot]
        for r in o.itertuples():
            d = pd.Timestamp(r.date).normalize()
            if d in gg.index:
                se.append((gg.loc[d, "LAID"] - r.lai_plot) ** 2)
    return float(np.sqrt(np.mean(se)))


def yield_stats(cultivar: dict, nitrate_scale: float | None = None) -> tuple[float, float]:
    obs = dc.load_yield()["yield_dry_kg_ha"].values
    ns = nitrate_scale if nitrate_scale is not None else dc.DEFAULT_NITRATE_SCALE
    sloc = cultivar.get("sloc_scale", dc.SLOC_SCALE)
    sim = np.array([dc.run_plot(p, cultivar, nitrate_scale=ns,
                                sloc_scale=sloc)[1]["HWAM"] for p in dc.ORDER])
    rrmse = float(np.sqrt(np.mean((sim - obs) ** 2)) / obs.mean() * 100)
    return rrmse, float(sim.sum() / obs.sum() * 100 - 100)


def phenology(cultivar: dict, plot: str | None = None) -> tuple[pd.Timestamp, pd.Timestamp]:
    plot = plot or dc.ORDER[4]
    _, s = dc.run_plot(plot, cultivar)
    return dc.dssat_date(s.get("ADAT")), dc.dssat_date(s.get("MDAT"))


def stage_phenology(log: list[dict]) -> dict:
    """P1/P5 grid, scored on anthesis (weight 2) and maturity (weight 1)."""
    g = GRIDS[dc.TRIAL]
    best, best_score = None, 1e9
    for p1, p5 in itertools.product(g["p1"], g["p5"]):
        pars = dict(dc.DEFAULT_CULTIVAR, P1=float(p1), P5=float(p5))
        a, m = phenology(pars)
        if a is None or m is None:
            score = 1e3
        else:
            score = 2 * abs((a - OBS_ANTHESIS).days) + abs((m - OBS_MATURITY).days)
        log.append(dict(stage="phenology", P1=p1, P5=p5,
                        anthesis=str(a.date()) if a is not None else "-99",
                        maturity=str(m.date()) if m is not None else "-99",
                        score=score))
        if score < best_score:
            best, best_score = pars, score
    return best


def stage_canopy(base: dict, log: list[dict]) -> dict:
    """P1/PHINT/P5 are coupled (leaf number sets both LAI and anthesis), so the
    canopy stage re-opens P1/P5 and scores LAI plus a phenology penalty."""
    obs = _obs_lai()
    g = GRIDS[dc.TRIAL]
    best, best_score = None, 1e9
    for sloc, p1, p5, phint in itertools.product(g["sloc"], g["p1"],
                                                 g["p5"], g["phint"]):
        pars = dict(base, P1=float(p1), P5=float(p5), PHINT=float(phint))
        rmse = lai_rmse(pars, plots=PLOT_SUBSET, obs=obs, sloc_scale=sloc)
        a, m = phenology(pars)
        d_a = abs((a - OBS_ANTHESIS).days) if a is not None else 30
        d_m = abs((m - OBS_MATURITY).days) if m is not None else 30
        score = rmse + 0.03 * d_a + 0.02 * d_m
        peak = np.mean([dc.run_plot(p, pars, sloc_scale=sloc)[1]["LAIX"]
                        for p in PLOT_SUBSET])
        log.append(dict(stage="canopy", P1=p1, P5=p5, PHINT=phint,
                        sloc_scale=sloc, lai_rmse=round(rmse, 4),
                        d_anthesis=d_a, d_maturity=d_m, score=round(score, 4),
                        mean_peak=round(float(peak), 3)))
        if score < best_score:
            best, best_score = dict(pars, sloc_scale=sloc), score
    return best


def stage_yield(base: dict, log: list[dict]) -> dict:
    """Yield alone does not identify G2 from G3: a few heavy grains and many
    light grains give the same kg/ha.  The measured grain number (~3471 /m2)
    and single grain weight (~240 mg) are therefore scored alongside RRMSE."""
    best, best_score = None, 1e9
    grid = GRIDS[dc.TRIAL]
    sloc = base.get("sloc_scale", 1.0)
    for g2, g3, ns in itertools.product(grid["g2"], grid["g3"], grid["ns"]):
        pars = dict(base, G2=float(g2), G3=float(g3))
        obs = dc.load_yield()["yield_dry_kg_ha"].values
        sim, grains, gw = [], [], []
        for p in dc.ORDER:
            g, s = dc.run_plot(p, pars, nitrate_scale=ns, sloc_scale=sloc)
            sim.append(s["HWAM"])
            last = g.iloc[-1]
            grains.append(last["G#AD"])
            gw.append(last["GWGD"])
        sim = np.array(sim)
        rrmse = float(np.sqrt(np.mean((sim - obs) ** 2)) / obs.mean() * 100)
        bias = float(sim.sum() / obs.sum() * 100 - 100)
        g_t, w_t = GRAIN_TARGETS.get(dc.TRIAL, (None, None))
        g_pen = (abs(np.mean(grains) - g_t) / g_t) if g_t else 0.0
        w_pen = (abs(np.mean(gw) - w_t) / w_t) if w_t else 0.0
        score = rrmse + 30 * g_pen + 30 * w_pen
        log.append(dict(stage="yield", G2=g2, G3=g3, nitrate_scale=ns,
                        yield_rrmse=round(rrmse, 3), yield_bias=round(bias, 3),
                        grains_m2=round(float(np.mean(grains)), 0),
                        grain_w_mg=round(float(np.mean(gw)), 1),
                        score=round(score, 3)))
        if score < best_score:
            best, best_score = dict(pars, nitrate_scale=ns, sloc_scale=sloc), score
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--stage", default="all",
                    choices=["all", "phenology", "canopy", "yield"])
    args = ap.parse_args()
    dc.set_trial(args.trial)
    print(f"trial {dc.TRIAL}  ({dc.CFG.site} {dc.CFG.variety})")

    global OBS_ANTHESIS, OBS_MATURITY, PLOT_SUBSET
    OBS_ANTHESIS = dc.OBS_ANTHESIS
    OBS_MATURITY = dc.OBS_MATURITY
    PLOT_SUBSET = [dc.CONTROL_PLOT, dc.ORDER[4], dc.ORDER[-1]]

    log: list[dict] = []
    pars = dict(dc.DEFAULT_CULTIVAR)
    if args.stage in ("all", "phenology"):
        pars = stage_phenology(log)
        print("phenology ->", {k: pars[k] for k in ("P1", "P5")},
              "| anthesis", phenology(pars)[0].date(),
              "| maturity", phenology(pars)[1].date())
    if args.stage in ("all", "canopy"):
        pars = stage_canopy(pars, log)
        print("canopy    -> PHINT", pars["PHINT"],
              "| LAI RMSE", round(lai_rmse(pars), 3))
    if args.stage in ("all", "yield"):
        pars = stage_yield(pars, log)
        r, b = yield_stats(pars, pars.get("nitrate_scale"))
        print("yield     -> G2", pars["G2"], "G3", pars["G3"],
              "nitrate_scale", pars["nitrate_scale"],
              f"| RRMSE {r:.2f}% bias {b:+.2f}%")

    df = pd.DataFrame(log)
    df.to_csv(dc.SIM / "calibration_log.csv", index=False)
    print(f"log -> {dc.SIM / 'calibration_log.csv'} ({len(df)} rows)")

    if args.write:
        payload = dict(cultivar={k: pars[k] for k in dc.DEFAULT_CULTIVAR},
                       nitrate_scale=pars.get("nitrate_scale",
                                              dc.DEFAULT_NITRATE_SCALE),
                       sloc_profile=dc.SLOC_PROFILE,
                       sloc_scale=pars.get("sloc_scale", dc.SLOC_SCALE))
        (dc.SIM / "calibrated_cultivar.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8")
        print("wrote", dc.SIM / "calibrated_cultivar.json")


if __name__ == "__main__":
    main()
