"""Step 12 - the L series: the pre-optimisation baseline (sen = 0.33).

The C series (`s20_conservative_figures.py`) shows the adopted calibration, whose
leaf senescence acceleration is 0.20.  This module reruns the same nine plots with
the earlier value 0.33 and draws the identical four figures, so L1-L4 and C1-C4
sit side by side as a before / after pair.

Only one parameter differs between the two series:

    Leaf.CohortParameters.SenescenceAcceleration.One   0.33 (L) vs 0.20 (C)

Usage:  python analysis/src/s12_figures.py
"""
from __future__ import annotations

import io
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_override import SimCase, run_many  # noqa: E402
from common import DS, SIM  # noqa: E402
from s08_validate import compare  # noqa: E402
from s20_conservative_figures import (SEN, TRIAL, load, make_all,  # noqa: E402
                                      run_reference)

SEN_BASELINE = 0.33
PLOTS = ["A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"]


def run_baseline() -> pd.DataFrame:
    """Run the nine Shandong plots with the pre-optimisation senescence rate."""
    cal = json.load(io.open(SIM / "phenology_calibrated.json", encoding="utf-8"))["乐农87"]["params"]
    gro = json.load(io.open(SIM / "calibrated_cultivars.json", encoding="utf-8"))["乐农87"]["growth"]
    params = {**cal, **gro, SEN: SEN_BASELINE}

    plots = pd.read_csv(DS / "plots.csv")
    plots = plots[(plots["trial_id"] == TRIAL) & plots["plot"].isin(PLOTS)]

    jobs, meta = [], []
    for _, r in plots.iterrows():
        ident = f"{r['trial_id']}_{r['plot']}"
        jobs.append((SimCase.load(Path(SIM).parent / "apsim" / f"{ident}.apsimx")
                     .set_commands(params), ident))
        meta.append(dict(sim=ident, trial_id=r["trial_id"], plot=r["plot"],
                         variety=r["variety"]))
    with tempfile.TemporaryDirectory() as td:
        tables = run_many(jobs, td)

    rows = []
    for m in meta:
        df = tables[m["sim"]]
        if df is None:
            continue
        df = df.copy()
        df["sim"] = m["sim"]
        for k in ("trial_id", "plot", "variety"):
            df[k] = m[k]
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    daily = run_baseline()
    daily.to_csv(SIM / "sim_baseline_sen033_daily.csv", index=False)

    cmp, stats = compare(daily)
    sd = cmp[cmp["trial_id"] == TRIAL].set_index("plot")
    sd.to_csv(SIM / "validate_plots_baseline.csv")
    st = stats[stats["trial_id"] == TRIAL]
    st.to_csv(SIM / "validate_stats_baseline.csv", index=False)

    _sim_c, obs, _val_c = load()
    sim_pot = run_reference(SEN_BASELINE)
    print(f"sim rows={len(daily)}  plots={daily['sim'].nunique()}")
    for p in make_all("L", sim=daily[daily["trial_id"] == TRIAL], obs=obs, val=sd,
                      note="（sen=0.33，优化前）", sim_pot=sim_pot):
        print(" ", p)

    print("\n=== L 系列基线（sen=0.33） ===")
    print(st.round(2).to_string(index=False))
    print("\n=== C 系列优化后（sen=0.20） ===")
    v = pd.read_csv(SIM / "validate_plots.csv")
    c = v[v["trial_id"] == TRIAL]
    print(f"  产量 RRMSE {np.sqrt(np.mean((c.sim_yield - c.obs_yield_dry) ** 2)) / c.obs_yield_dry.mean() * 100:+.2f}%"
          f"   LAI RMSE {c.lai_rmse.mean():.3f}"
          f"   生物量 RRMSE {np.sqrt(np.mean((c.sim_bio_max - c.obs_bio_max) ** 2)) / c.obs_bio_max.mean() * 100:.2f}%")


if __name__ == "__main__":
    main()