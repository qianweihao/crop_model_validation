"""Sweep leaf-area shape parameters against the observed LAI curve.

Motivation: the observations show a clear plateau after the peak
(08-14 -> 08-29 loses only 5%), while the current simulation peaks 11 days
early and declines immediately.  This script scans the parameters that control
when the canopy reaches full area, and reports both the point RMSE and the
shape diagnostics (peak date, plateau slope).

Usage:  python analysis/src/s14_lai_shape_sweep.py
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

LLP = "Leaf.CohortParameters.MaxArea.LargestLeafPosition.RelativePositionOfLargestLeaf.FixedValue"
SKEW = "Leaf.CohortParameters.MaxArea.Skewness.FixedValue"
SEN = "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue"
RUE = "Leaf.Photosynthesis.RUE.FixedValue"


PLOTS = ["A1", "A3B2"]

OBS_DATES = ["2026-07-03", "2026-07-11", "2026-07-23", "2026-08-04",
             "2026-08-14", "2026-08-21", "2026-08-29"]


def observed() -> pd.DataFrame:
    d = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    d = d[d["trial_id"] == "SD_SF"]
    return d.pivot_table(index="date", columns="plot", values="lai_plot")


def run(configs: list[dict]) -> dict[str, pd.DataFrame]:
    cases = []
    for cfg in configs:
        label = cfg["_label"]
        params = {k: v for k, v in cfg.items() if k != "_label"}
        for plot in PLOTS:
            c = SimCase.load(Path(SIM).parent / "apsim" / f"SD_SF_{plot}.apsimx")
            c = c.set_commands({**CAL, **GRO, **params})
            cases.append((c, f"{label}__{plot}"))
    with tempfile.TemporaryDirectory() as td:
        res = run_many(cases, td)
    out = {}
    for key, df in res.items():
        if df is not None and "Clock.Today" in df:
            df = df.copy()
            df["Clock.Today"] = pd.to_datetime(df["Clock.Today"])
            out[key] = df
    return out


def score(res: dict[str, pd.DataFrame], obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    labels = sorted({k.split("__")[0] for k in res})
    for lab in labels:
        errs, peak_dates, plat = [], [], []
        for plot in PLOTS:
            df = res.get(f"{lab}__{plot}")
            if df is None:
                continue
            s = df.set_index("Clock.Today")["Maize.Leaf.LAI"]
            o = obs[plot].dropna()
            for d in OBS_DATES:
                t = pd.Timestamp(d)
                if t in o.index and t in s.index:
                    errs.append(s.loc[t] - o.loc[t])
            peak_dates.append(s.idxmax().dayofyear)
            p = s.loc[pd.Timestamp("2026-08-14")] if pd.Timestamp("2026-08-14") in s.index else np.nan
            q = s.loc[pd.Timestamp("2026-08-29")] if pd.Timestamp("2026-08-29") in s.index else np.nan
            plat.append((q - p) / p * 100 if p else np.nan)
        e = np.array(errs)
        rows.append(dict(label=lab, rmse=float(np.sqrt(np.mean(e ** 2))),
                         bias=float(e.mean()),
                         peak_day=int(np.mean(peak_dates)) if peak_dates else -1,
                         plateau_pct=float(np.nanmean(plat))))
    return pd.DataFrame(rows)


def main() -> None:
    obs = observed()
    print("实测：峰值 08-14，08-14→08-29 变化 "
          + ", ".join(f"{p} {(obs[p][pd.Timestamp('2026-08-29')]/obs[p][pd.Timestamp('2026-08-14')]-1)*100:+.0f}%"
                      for p in PLOTS))
    print("（模拟当前：峰值 08-03，08-14→08-29 −" +
          "7.5%" + "）\n")

    grid = []
    for llp in [0.65, 0.70, 0.75, 0.80]:
        grid.append({LLP: llp, "_label": f"LLP={llp:.2f}"})
    for sk in [-0.6, -0.2, 0.001, 0.6]:
        grid.append({SKEW: sk, "_label": f"Skew={sk:+.3f}"})
    for se in [0.20, 0.33, 0.50]:
        grid.append({SEN: se, "_label": f"Sen={se:.2f}"})

    res = run(grid)
    tbl = score(res, obs).sort_values("rmse")
    tbl["peak_mmdd"] = tbl["peak_day"].map(
        lambda d: pd.Timestamp("2026-01-01") + pd.Timedelta(days=d - 1) if d > 0 else None)
    print(tbl[["label", "rmse", "bias", "peak_mmdd", "plateau_pct"]]
          .to_string(index=False, float_format=lambda x: f"{x:7.3f}"))
    tbl.to_csv(SIM / "lai_shape_sweep.csv", index=False)


if __name__ == "__main__":
    main()
