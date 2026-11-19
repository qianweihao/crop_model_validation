"""Scan the leaf senescence rate against both the LAI curve and biomass.

Why: the observations run out on 08-29, so the only measured constraint on
post-peak decline is the 08-14 -> 08-29 window (A1 -5%, A3B2 -1%).  The
"50-70% seasonal decline" used earlier was an assumed target, not data.
Slowing senescence should both flatten the plateau and keep more radiation
intercepted, so it is the natural lever for the -19% biomass deficit.

Usage:  python analysis/src/s15_senescence_biomass_sweep.py
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

PLOTS = ["A1", "A2B1", "A3B2", "A5B2"]
OBS_DATES = ["2026-07-03", "2026-07-11", "2026-07-23", "2026-08-04",
             "2026-08-14", "2026-08-21", "2026-08-29"]


def observed():
    d = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    d = d[d["trial_id"] == "SD_SF"]
    lai = d.pivot_table(index="date", columns="plot", values="lai_plot")
    v = pd.read_csv(SIM / "validate_plots.csv")
    v = v[v["trial_id"] == "SD_SF"].set_index("plot")
    return lai, v["obs_bio_max"]


def main() -> None:
    lai_obs, bio_obs = observed()
    grid = [
        {SEN: v, RUE: r, "_label": f"sen{v:.2f}_rue{r:.2f}"}
        for v in (0.10, 0.15, 0.20, 0.25, 0.33)
        for r in (2.381, 2.70)
    ]
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
        errs, plat, bios, sim_bio, yields = [], [], [], [], []
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
            sim_bio.append(df["Maize.AboveGround.Wt"].max() * 10)
            bios.append(bio_obs.get(plot, np.nan))
            yields.append(df["Maize.Grain.Total.Wt"].max() * 10)
        e = np.array(errs)
        rows.append(dict(label=lab, lai_rmse=float(np.sqrt(np.mean(e ** 2))),
                         lai_bias=float(e.mean()),
                         plateau_pct=float(np.nanmean(plat)),
                         bio_sim=float(np.mean(sim_bio)),
                         bio_obs=float(np.nanmean(bios)),
                         bio_bias=float((np.mean(sim_bio) / np.nanmean(bios) - 1) * 100),
                         yield_sim=float(np.mean(yields))))
    tbl = pd.DataFrame(rows).sort_values("lai_rmse")
    pd.set_option("display.width", 200)
    print("实测：08-14→08-29 A1 -5.0%  A3B2 -1.1%")
    print(tbl.to_string(index=False, float_format=lambda x: f"{x:8.2f}"))
    tbl.to_csv(SIM / "senescence_biomass_sweep.csv", index=False)


if __name__ == "__main__":
    main()
