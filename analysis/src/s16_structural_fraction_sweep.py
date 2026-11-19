"""Stage 1: restore the LAI curve after raising RUE.

RUE 2.70 fixes the -18% biomass bias but inflates LAI (bias +0.7, RMSE 0.87).
Raising Leaf.CohortParameters.StructuralFraction lifts specific leaf weight, so
the same leaf mass yields less area.  This sweep looks for the structural
fraction that brings LAI back, and reports what it does to biomass and yield.

Usage:  python analysis/src/s16_structural_fraction_sweep.py
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

CAL = json.load(io.open(SIM / "phenology_calibrated.json", encoding="utf-8"))["乐农87"]["params"]
GRO = json.load(io.open(SIM / "calibrated_cultivars.json", encoding="utf-8"))["乐农87"]["growth"]

SEN = "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue"
RUE = "Leaf.Photosynthesis.RUE.FixedValue"
SF = "Leaf.CohortParameters.StructuralFraction.FixedValue"
GPC = "Grain.MaximumGrainsPerCob.FixedValue"

PLOTS = ["A1", "A2B1", "A3B2", "A5B2"]
OBS_DATES = ["2026-07-03", "2026-07-11", "2026-07-23", "2026-08-04",
             "2026-08-14", "2026-08-21", "2026-08-29"]


def observed():
    d = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    d = d[d["trial_id"] == "SD_SF"]
    lai = d.pivot_table(index="date", columns="plot", values="lai_plot")
    v = pd.read_csv(SIM / "validate_plots.csv")
    v = v[v["trial_id"] == "SD_SF"].set_index("plot")
    return lai, v


def main() -> None:
    lai_obs, val = observed()
    grid = [
        {RUE: 2.70, SEN: 0.20, SF: sf, "_label": f"sf{sf:.2f}"}
        for sf in (0.827, 0.86, 0.89, 0.92, 0.95, 0.98)
    ]
    grid.append({RUE: 2.381, SEN: 0.20, SF: 0.827, "_label": "current"})
    cases = []
    for cfg in grid:
        params = {k: v for k, v in cfg.items() if k != "_label"}
        for plot in PLOTS:
            c = SimCase.load(Path(SIM).parent / "apsim" / f"SD_SF_{plot}.apsimx")
            cases.append((c.set_commands({**CAL, **GRO, **params}),
                          f"{cfg['_label']}__{plot}"))
    with tempfile.TemporaryDirectory() as td:
        res = run_many(cases, td)

    rows = []
    for cfg in grid:
        lab = cfg["_label"]
        errs, plat, bio_s, bio_o, yld_s, yld_o = [], [], [], [], [], []
        for plot in PLOTS:
            df = res.get(f"{lab}__{plot}")
            if df is None:
                continue
            df = df.copy()
            df["Clock.Today"] = pd.to_datetime(df["Clock.Today"])
            s = df.set_index("Clock.Today")["Maize.Leaf.LAI"]
            o = lai_obs[plot].dropna()
            for d in OBS_DATES:
                t = pd.Timestamp(d)
                if t in o.index and t in s.index:
                    errs.append(s.loc[t] - o.loc[t])
            p, q = s.get(pd.Timestamp("2026-08-14")), s.get(pd.Timestamp("2026-08-29"))
            if p and q:
                plat.append((q - p) / p * 100)
            bio_s.append(df["Maize.AboveGround.Wt"].max() * 10)
            yld_s.append(df["Maize.Grain.Total.Wt"].max() * 10)
            bio_o.append(val.loc[plot, "obs_bio_max"])
            yld_o.append(val.loc[plot, "obs_yield_dry"])
        if not errs:
            rows.append(dict(label=lab, note="全部失败"))
            continue
        e = np.array(errs)
        hi_s = np.mean(yld_s) / np.mean(bio_s)
        rows.append(dict(label=lab, lai_rmse=float(np.sqrt(np.mean(e ** 2))),
                         lai_bias=float(e.mean()),
                         plateau_pct=float(np.nanmean(plat)),
                         bio_sim=float(np.mean(bio_s)), bio_obs=float(np.mean(bio_o)),
                         bio_bias=float((np.mean(bio_s) / np.mean(bio_o) - 1) * 100),
                         yld_sim=float(np.mean(yld_s)), yld_obs=float(np.mean(yld_o)),
                         yld_bias=float((np.mean(yld_s) / np.mean(yld_o) - 1) * 100),
                         HI_sim=float(hi_s), HI_obs=0.486))
    tbl = pd.DataFrame(rows)
    pd.set_option("display.width", 240)
    print("实测：08-14→08-29 A1 -5.0%  A3B2 -1.1%   (4 小区均值)")
    print(tbl.to_string(index=False, float_format=lambda x: f"{x:8.2f}"))
    tbl.to_csv(SIM / "structural_fraction_sweep.csv", index=False)


if __name__ == "__main__":
    main()