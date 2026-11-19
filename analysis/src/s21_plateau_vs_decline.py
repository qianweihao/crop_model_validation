"""Can a real plateau and a steep late decline be had at the same time?

The conservative scheme (sen=0.20) reproduces the measured 08-14 -> 08-29 window
but leaves the canopy far too green at harvest (LAI 2.2, decline 29%).  Raising
the senescence rate steepens the tail but destroys the plateau.  This script
quantifies that trade-off and tests an alternative, physically motivated lever:
the September drought in Dezhou (6.4 mm of rain) should accelerate senescence on
its own if the profile is not held at field capacity all season.

Two experiments:

  A. senescence scan          - sen in 0.20..0.90, water as now
  B. initial water scan       - sen fixed at 0.20, initial PAW in 100..40 %

Usage:  python analysis/src/s21_plateau_vs_decline.py
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

PLOTS = ["A1", "A2B1", "A3B2", "A5B2"]
OBS_DATES = ["2026-07-03", "2026-07-11", "2026-07-23", "2026-08-04",
             "2026-08-14", "2026-08-21", "2026-08-29"]


def observed():
    d = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    d = d[d["trial_id"] == "SD_SF"]
    lai = d.pivot_table(index="date", columns="plot", values="lai_plot")
    v = pd.read_csv(SIM / "validate_plots.csv")
    return lai, v[v["trial_id"] == "SD_SF"].set_index("plot")


def set_initial_paw(case: SimCase, frac: float) -> SimCase:
    fld = case._field()
    soil = next(c for c in fld["Children"] if c.get("Name") == "Soil")
    nodes = {c["Name"]: c for c in soil["Children"]}
    ll, dul = nodes["Physical"]["LL15"], nodes["Physical"]["DUL"]
    nodes["Water"]["InitialValues"] = [l + frac * (d - l) for l, d in zip(ll, dul)]
    return case


def build(grid, tag_prefix):
    cases = []
    for cfg in grid:
        lab = cfg["_label"]
        for plot in PLOTS:
            c = SimCase.load(Path(SIM).parent / "apsim" / f"SD_SF_{plot}.apsimx")
            c = c.set_commands({**CAL, **GRO, **{k: v for k, v in cfg.items()
                                                 if k != "_label" and k != "_paw"}})
            if "_paw" in cfg:
                c = set_initial_paw(c, cfg["_paw"])
            cases.append((c, f"{lab}__{plot}"))
    return cases


def evaluate(res, grid, lai_obs, val):
    rows = []
    for cfg in grid:
        lab = cfg["_label"]
        errs, plat, peak_d, bio, yld = [], [], [], [], []
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
            peak_d.append(s.idxmax().dayofyear)
            bio.append(df["Maize.AboveGround.Wt"].max() * 10)
            yld.append(df["Maize.Grain.Total.Wt"].max() * 10)
            end = s.loc[s.index >= "2026-10-01"].mean()
            pk = s.max()
            rows.append(dict(label=lab, plot=plot, peak=float(pk), peak_day=int(s.idxmax().dayofyear),
                             end=float(end), decline=float((end / pk - 1) * 100),
                             plateau=float((q - p) / p * 100) if p and q else np.nan,
                             err=float(np.mean([s.loc[pd.Timestamp(d)] - o.loc[pd.Timestamp(d)]
                                                for d in OBS_DATES
                                                if pd.Timestamp(d) in o.index and pd.Timestamp(d) in s.index]))
                             if any(pd.Timestamp(d) in o.index for d in OBS_DATES) else np.nan))
    per = pd.DataFrame(rows)
    agg = per.groupby("label").agg(
        LAI_RMSE=("err", lambda x: float(np.sqrt(np.nanmean(np.square(x))))),
        peak_mean=("peak", "mean"), peak_day=("peak_day", "mean"),
        end_mean=("end", "mean"), decline=("decline", "mean"),
        plateau=("plateau", "mean")).reset_index()
    return agg, per


def main() -> None:
    lai_obs, val = observed()
    print("目标：峰值 ≈ 3.31（实测）；平台期 ≈ -3.2%；季末 LAI ≈ 1.0-1.5（衰减 55-65%）\n")

    grid_a = [{SEN: v, "_label": f"A_sen{v:.2f}"} for v in (0.20, 0.35, 0.50, 0.70, 0.90)]
    with tempfile.TemporaryDirectory() as td:
        res_a = run_many(build(grid_a, "A"), td)
    agg_a, per_a = evaluate(res_a, grid_a, lai_obs, val)
    print("=== A. 衰老速率扫描（水分维持现状）===")
    print(agg_a.to_string(index=False, float_format=lambda x: f"{x:8.2f}"))

    grid_b = [{SEN: 0.20, "_paw": f, "_label": f"B_paw{int(f*100)}"} for f in (1.00, 0.75, 0.55, 0.40)]
    with tempfile.TemporaryDirectory() as td:
        res_b = run_many(build(grid_b, "B"), td)
    agg_b, per_b = evaluate(res_b, grid_b, lai_obs, val)
    print("\n=== B. 初始含水量扫描（sen 固定 0.20）===")
    print(agg_b.to_string(index=False, float_format=lambda x: f"{x:8.2f}"))

    out = pd.concat([agg_a.assign(exp="A"), agg_b.assign(exp="B")], ignore_index=True)
    out.to_csv(SIM / "plateau_vs_decline.csv", index=False)
    per = pd.concat([per_a.assign(exp="A"), per_b.assign(exp="B")], ignore_index=True)
    per.to_csv(SIM / "plateau_vs_decline_per_plot.csv", index=False)
    print("\n累计 ->", SIM / "plateau_vs_decline.csv")


if __name__ == "__main__":
    main()