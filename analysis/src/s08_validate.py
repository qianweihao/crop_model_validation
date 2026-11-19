"""Step 08 - apply the calibrated cultivars to every plot and validate.

Produces per-plot observed vs simulated tables and fit statistics.  Figures are
not drawn here: the final set comes from s12_figures.py (L1-L4) and
s20_conservative_figures.py (C1-C4).
Calibrated parameters come from analysis/out/sim/calibrated_cultivars.json.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_override import SimCase, run_many  # noqa: E402
from common import APSIM_DIR, DS, DRY_FACTOR, FIG, REPORTS, SIM  # noqa: E402

CAL_FILE = SIM / "calibrated_cultivars.json"
PHENO_FILE = SIM / "phenology_calibrated.json"


def load_params() -> dict[str, dict]:
    """Growth parameters from s07, phenology from the leaf-constrained s10 fit."""
    if not CAL_FILE.exists():
        raise SystemExit(f"run s07 first: {CAL_FILE} missing")
    raw = json.loads(CAL_FILE.read_text(encoding="utf-8"))
    phen = (json.loads(PHENO_FILE.read_text(encoding="utf-8"))
            if PHENO_FILE.exists() else {})
    out = {}
    for var, entry in raw.items():
        p = {}
        p.update(entry.get("phenology", {}))          # fallback
        if var in phen:
            p.update(phen[var]["params"])             # leaf-constrained wins
        p.update(entry.get("growth", {}))
        if p:
            out[var] = p
    return out


def run_all(params: dict[str, dict]) -> pd.DataFrame:
    plots = pd.read_csv(DS / "plots.csv")
    plots = plots[plots["trial_id"].isin(["SD_SF", "XX_SF", "XX_IRR"])]
    plots = plots[plots["variety"].isin(params)]

    jobs, meta = [], []
    for _, r in plots.iterrows():
        ident = f"{r['trial_id']}_{r['plot']}"
        src = APSIM_DIR / f"{ident}.apsimx"
        if not src.exists():
            continue
        jobs.append((SimCase.load(src).set_commands(params[r["variety"]]), ident))
        meta.append(dict(sim=ident, trial_id=r["trial_id"], plot=r["plot"],
                         variety=r["variety"]))
    with tempfile.TemporaryDirectory() as td:
        tables = run_many(jobs, td)

    rows = []
    for m in meta:
        df = tables[m["sim"]].copy()
        df["sim"] = m["sim"]
        for k in ("trial_id", "plot", "variety"):
            df[k] = m[k]
        rows.append(df)
    daily = pd.concat(rows, ignore_index=True)
    daily.to_csv(SIM / "sim_calibrated_daily.csv", index=False)
    return daily


def compare(daily: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    yld = pd.read_csv(DS / "yield.csv")

    ob = obs.dropna(subset=["density_measured"])
    ob = ob.assign(sid=ob["trial_id"] + "_" + ob["plot"],
                   biomass_kg_ha=ob["biomass_kg_mu"] * 15.0)

    rows = []
    for sid, g in daily.groupby("sim"):
        g = g.sort_values("Clock.Today").set_index("Clock.Today")
        t = ob[ob["sid"] == sid].sort_values("date")
        xs = [(d - pd.Timestamp("2026-06-08")).days for d in g.index]
        lai_pred = np.interp([(d - pd.Timestamp("2026-06-08")).days for d in t["date"]],
                             xs, g["Maize.Leaf.LAI"])
        bio_pred = np.interp([(d - pd.Timestamp("2026-06-08")).days for d in t["date"]],
                             xs, g["Maize.AboveGround.Wt"] * 10.0)
        yrow = yld[(yld["trial_id"] + "_" + yld["plot"]) == sid]
        obs_wet = float(yrow["yield_kg_ha"].iloc[0]) if len(yrow) else np.nan
        obs_dry = float(yrow["yield_dry_kg_ha"].iloc[0]) if len(yrow) else np.nan
        lai_o = t["lai_plot"].to_numpy(dtype=float)
        bio_o = t["biomass_kg_ha"].to_numpy(dtype=float)
        if len(t) and np.isfinite(lai_o).any():
            m = np.isfinite(lai_o) & np.isfinite(lai_pred)
            lai_rmse = float(np.sqrt(np.mean((lai_pred[m] - lai_o[m]) ** 2))) if m.any() else np.nan
        else:
            lai_rmse = np.nan
        if len(t) and np.isfinite(bio_o).any():
            m = np.isfinite(bio_o) & np.isfinite(bio_pred)
            bio_rmse = float(np.sqrt(np.mean((bio_pred[m] - bio_o[m]) ** 2))) if m.any() else np.nan
        else:
            bio_rmse = np.nan
        sim_y = float(g["Maize.Grain.Total.Wt"].max() * 10.0)     # dry
        sim_b = float(g["Maize.AboveGround.Wt"].max() * 10.0)    # dry
        rows.append(dict(
            sim=sid, trial_id=g["trial_id"].iloc[0], plot=g["plot"].iloc[0],
            variety=g["variety"].iloc[0],
            sim_yield=sim_y,
            obs_yield=obs_wet,
            obs_yield_dry=obs_dry,
            sim_yield_14=sim_y / DRY_FACTOR,
            sim_lai_max=float(g["Maize.Leaf.LAI"].max()),
            obs_lai_max=float(np.nanmax(lai_o)) if len(t) else np.nan,
            sim_bio_max=sim_b,
            obs_bio_max=float(np.nanmax(bio_o)) if len(t) else np.nan,
            sim_HI=sim_y / sim_b if sim_b > 0 else np.nan,
            obs_HI=obs_dry / float(np.nanmax(bio_o)) if len(t) and np.isfinite(bio_o).any() else np.nan,
            lai_rmse=lai_rmse,
            bio_rmse=bio_rmse,
        ))
    cmp = pd.DataFrame(rows)
    # all yield comparisons on a dry-weight basis (APSIM reports dry grain)
    cmp["yield_bias_pct"] = (cmp["sim_yield"] - cmp["obs_yield_dry"]) / cmp["obs_yield_dry"] * 100

    stats = []
    for (tid,), g in cmp.groupby(["trial_id"]):
        for var, gg in g.groupby("variety"):
            stats.append(dict(
                trial_id=tid, variety=var, n=len(gg),
                yield_rrmse=float(np.sqrt(np.mean(
                    (gg["sim_yield"] - gg["obs_yield_dry"]) ** 2))
                    / gg["obs_yield_dry"].mean() * 100),
                yield_bias=float(gg["yield_bias_pct"].mean()),
                lai_rmse=float(gg["lai_rmse"].mean()),
                sim_HI=float(gg["sim_HI"].mean()),
                obs_HI=float(gg["obs_HI"].mean()),
                bio_rrmse=float(np.sqrt(np.mean(
                    (gg["sim_bio_max"] - gg["obs_bio_max"]) ** 2))
                    / gg["obs_bio_max"].mean() * 100),
            ))
    return cmp, pd.DataFrame(stats)


def figures(daily: pd.DataFrame) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    made = []
    obs = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    obs = obs.assign(sid=obs["trial_id"] + "_" + obs["plot"])

    # 1: observed vs simulated LAI curves per trial
    trials = sorted(daily["trial_id"].unique())
    fig, axes = plt.subplots(1, len(trials), figsize=(5 * len(trials), 4.2), squeeze=False)
    for ax, tid in zip(axes[0], trials):
        for sid, g in daily[daily["trial_id"] == tid].groupby("sim"):
            ax.plot(g["Clock.Today"], g["Maize.Leaf.LAI"], lw=1.0, alpha=0.75)
            t = obs[obs["sid"] == sid]
            ax.plot(t["date"], t["lai_plot"], "o", ms=3, alpha=0.7)
        ax.set_title(f"{tid}  LAI (line=sim, dot=obs)", fontsize=10)
        ax.set_ylabel("LAI")
        ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    p = FIG / "08_lai_curves.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    made.append(p)

    # 2: simulated vs observed yield
    y = daily.groupby("sim").agg(sim=("Maize.Grain.Total.Wt", lambda s: s.max() * 10))
    yld = pd.read_csv(DS / "yield.csv")
    yld = yld.assign(sid=yld["trial_id"] + "_" + yld["plot"]).set_index("sid")
    y["obs"] = yld["yield_dry_kg_ha"].reindex(y.index)

    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    ax.scatter(y["obs"], y["sim"], s=42)
    lo = min(y["obs"].min(), y["sim"].min()) * 0.9
    hi = max(y["obs"].max(), y["sim"].max()) * 1.05
    ax.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax.set_xlabel("observed dry grain (kg/ha)")
    ax.set_ylabel("simulated dry grain (kg/ha)")
    ax.set_title("1:1 dry grain yield after calibration", fontsize=10)
    fig.tight_layout()
    p = FIG / "08_yield_1to1.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    made.append(p)
    return made


def main() -> None:
    params = load_params()
    print("cultivars with parameters:", list(params))
    daily = run_all(params)
    cmp, stats = compare(daily)
    cmp.to_csv(SIM / "validate_plots.csv", index=False)
    stats.to_csv(SIM / "validate_stats.csv", index=False)

    pd.set_option("display.width", 220)
    print("\n=== per-plot ===")
    print(cmp[["sim", "variety", "obs_yield_dry", "sim_yield", "yield_bias_pct",
               "obs_lai_max", "sim_lai_max", "lai_rmse",
               "obs_HI", "sim_HI"]].round(2).to_string(index=False))
    print("\n=== fit statistics ===")
    print(stats.round(2).to_string(index=False))

    print(f"\ntables -> {SIM / 'validate_plots.csv'}, {SIM / 'validate_stats.csv'}")


if __name__ == "__main__":
    main()
