"""Stage 4: full 9-plot trade-off table for peak height and late decline.

Two facts drive this scan:

  * the simulated peak is 0.28 too low on the fertilised plots (A1 is exact),
    because leaf area does not respond to nitrogen in the model;
  * the late decline is too weak (season-end LAI 2.18 against 1.0-1.5 expected),
    because a slow senescence rate was needed to reproduce the plateau.

Lowering RelativePositionOfLargestLeaf raises the canopy peak, so this grid crosses
it with the senescence rate and reports every metric that matters, on all nine
plots, using the same LAI RMSE definition as s08.

Usage:  python analysis/src/s22_peak_decline_grid.py
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
LLP = "Leaf.CohortParameters.MaxArea.LargestLeafPosition.RelativePositionOfLargestLeaf.FixedValue"

PLOTS = ["A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"]
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
    for sen in (0.20, 0.35, 0.50):
        for llp in (0.65, 0.61, 0.57):
            grid.append({SEN: sen, LLP: llp, "_label": f"s{sen:.2f}_p{llp:.2f}"})
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
        errs, oc, sc, plat, end, decl, ysq, yb, pk = [], [], [], [], [], [], [], [], []
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
            oc.append(o.max())
            sc.append(s.max())
            pk.append(s.idxmax().dayofyear)
            p, q = s.get(pd.Timestamp("2026-08-14")), s.get(pd.Timestamp("2026-08-29"))
            if p and q:
                plat.append((q - p) / p * 100)
            e = s.loc[s.index >= "2026-10-01"].mean()
            end.append(e)
            decl.append((e / s.max() - 1) * 100)
            y = df["Maize.Grain.Total.Wt"].max() * 10
            ysq.append((y - val.loc[plot, "obs_yield_dry"]) ** 2)
            yb.append(y / val.loc[plot, "obs_yield_dry"] - 1)
        e_arr = np.array(errs)
        rows.append(dict(
            label=lab,
            LAI_RMSE=float(np.sqrt(np.mean(e_arr ** 2))),
            peak_sim=float(np.mean(sc)), peak_obs=float(np.mean(oc)),
            peak_bias=float(np.mean(np.array(sc) - np.array(oc))),
            peak_day=int(np.mean(pk)),
            plateau=float(np.nanmean(plat)),
            end_LAI=float(np.mean(end)),
            decline=float(np.mean(decl)),
            yield_RRMSE=float(np.sqrt(np.mean(ysq)) / val["obs_yield_dry"].mean() * 100),
            yield_bias=float(np.mean(yb) * 100)))
    tbl = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print("实测：峰值 3.31、平台期 -3.2%；收获期 LAI 应为 1.0-1.5（衰减 55-65%）\n")
    print(tbl.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))
    tbl.to_csv(SIM / "peak_decline_grid.csv", index=False)


if __name__ == "__main__":
    main()