"""Stage 3: stay on the LAI-neutral line while raising RUE, then trim the grain sink.

Stage 2 showed that RUE and LargestLeafPosition can compensate each other: as RUE
goes up, pushing the largest leaf later keeps the canopy curve unchanged while
total dry matter rises.  This sweep walks along that line and, for each point,
also scales the grain sink so yield and hence HI can be brought back to the
observed values.

Usage:  python analysis/src/s18_rue_sink_grid.py
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
LLP = "Leaf.CohortParameters.MaxArea.LargestLeafPosition.RelativePositionOfLargestLeaf.FixedValue"
GPC = "Grain.MaximumGrainsPerCob.FixedValue"
GPS = "Grain.MaximumPotentialGrainSize.FixedValue"

PLOTS = ["A1", "A2B1", "A3B2", "A5B2"]
OBS_DATES = ["2026-07-03", "2026-07-11", "2026-07-23", "2026-08-04",
             "2026-08-14", "2026-08-21", "2026-08-29"]


def observed():
    d = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    d = d[d["trial_id"] == "SD_SF"]
    lai = d.pivot_table(index="date", columns="plot", values="lai_plot")
    v = pd.read_csv(SIM / "validate_plots.csv")
    return lai, v[v["trial_id"] == "SD_SF"].set_index("plot")


def main() -> None:
    lai_obs, val = observed()
    grid = []
    for rue, llp in [(2.70, 0.690), (2.85, 0.705), (3.00, 0.720)]:
        for gscale in (1.00, 0.85, 0.75):
            grid.append({RUE: rue, LLP: llp, SEN: 0.20,
                         GPC: GRO[GPC] * gscale, GPS: GRO[GPS],
                         "_label": f"r{rue:.2f}_p{llp:.3f}_g{gscale:.2f}"})
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
        errs, plat, bio_s, yld_s, bio_o, yld_o = [], [], [], [], [], []
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
            rows.append(dict(label=lab, lai_rmse=np.nan))
            continue
        e = np.array(errs)
        rows.append(dict(label=lab, lai_rmse=float(np.sqrt(np.mean(e ** 2))),
                         lai_bias=float(e.mean()),
                         plat=float(np.nanmean(plat)),
                         bio_bias=float((np.mean(bio_s) / np.mean(bio_o) - 1) * 100),
                         yld_bias=float((np.mean(yld_s) / np.mean(yld_o) - 1) * 100),
                         HI=float(np.mean(yld_s) / np.mean(bio_s))))
    tbl = pd.DataFrame(rows).dropna(subset=["lai_rmse"])
    pd.set_option("display.width", 200)
    print("目标：lai_rmse 尽量小、lai_bias≈0、bio_bias≈0、yld_bias≈0、HI≈0.486")
    print(tbl.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))
    tbl.to_csv(SIM / "rue_sink_grid.csv", index=False)


if __name__ == "__main__":
    main()