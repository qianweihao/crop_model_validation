"""Quantify the post-peak LAI plateau and test a senescence fix.

The observations show the canopy reaching its peak on 08-14 and then holding
almost flat until 08-29 (A1 -5.0%, A3B2 -1.1%).  The current simulation declines
from the peak immediately (A1 -8.1%).  This script reruns all nine plots with a
slower senescence rate, scores the curve against the observations, and writes a
figure that shows the observed plateau next to both simulations.

Usage:  python analysis/src/s19_lai_plateau.py
"""
from __future__ import annotations

import io
import json
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_override import SimCase, run_many  # noqa: E402
from common import DS, FIG, SIM  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["SimHei", "Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 130, "savefig.dpi": 150,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    "legend.fontsize": 8, "axes.labelsize": 10,
})

CAL = json.load(io.open(SIM / "phenology_calibrated.json", encoding="utf-8"))["乐农87"]["params"]
GRO = json.load(io.open(SIM / "calibrated_cultivars.json", encoding="utf-8"))["乐农87"]["growth"]
SEN = "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue"

PLOTS = ["A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"]
OBS_DATES = ["2026-07-03", "2026-07-11", "2026-07-23", "2026-08-04",
             "2026-08-14", "2026-08-21", "2026-08-29"]
CANDIDATES = [("sen0.33 (当前)", 0.33), ("sen0.20", 0.20), ("sen0.15", 0.15), ("sen0.10", 0.10)]


def load_obs():
    d = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    d = d[d["trial_id"] == "SD_SF"]
    lai = d.pivot_table(index="date", columns="plot", values="lai_plot")
    v = pd.read_csv(SIM / "validate_plots.csv")
    return lai, v[v["trial_id"] == "SD_SF"].set_index("plot")


def run_all() -> dict[str, pd.DataFrame]:
    cases = []
    for label, sen in CANDIDATES:
        tag = f"s{int(sen*100):03d}"
        for plot in PLOTS:
            c = SimCase.load(Path(SIM).parent / "apsim" / f"SD_SF_{plot}.apsimx")
            cases.append((c.set_commands({**CAL, **GRO, SEN: sen}), f"{tag}__{plot}"))
    with tempfile.TemporaryDirectory() as td:
        res = run_many(cases, td)
    out = {}
    for k, df in res.items():
        if df is not None:
            df = df.copy()
            df["Clock.Today"] = pd.to_datetime(df["Clock.Today"])
            out[k] = df
    return out


def score(res, lai_obs, val) -> pd.DataFrame:
    rows = []
    for tag, label in [(f"s{int(sen*100):03d}", lab) for lab, sen in CANDIDATES]:
        errs, plat, rrmse, bias, bio_b, yld_b, his, tail = [], [], [], [], [], [], [], []
        for plot in PLOTS:
            df = res.get(f"{tag}__{plot}")
            if df is None:
                continue
            s = df.set_index("Clock.Today")["Maize.Leaf.LAI"]
            o = lai_obs[plot].dropna()
            for d in OBS_DATES:
                t = pd.Timestamp(d)
                if t in o.index and t in s.index:
                    errs.append(s.loc[t] - o.loc[t])
            p, q = s.get(pd.Timestamp("2026-08-14")), s.get(pd.Timestamp("2026-08-29"))
            if p and q:
                plat.append((q - p) / p * 100)
            bio = df["Maize.AboveGround.Wt"].max() * 10
            yld = df["Maize.Grain.Total.Wt"].max() * 10
            rrmse.append((yld - val.loc[plot, "obs_yield_dry"]) ** 2)
            bias.append((yld / val.loc[plot, "obs_yield_dry"] - 1) * 100)
            bio_b.append((bio / val.loc[plot, "obs_bio_max"] - 1) * 100)
            yld_b.append(yld)
            his.append(yld / bio)
            end = s.loc[s.index >= "2026-10-01"].mean()
            tail.append(100 * (end / s.max() - 1))
        e = np.array(errs)
        rows.append(dict(
            config=label, LAI_RMSE=float(np.sqrt(np.mean(e ** 2))),
            LAI_bias=float(e.mean()), plateau_pct=float(np.nanmean(plat)),
            yield_RRMSE=float(np.sqrt(np.mean(rrmse)) / val["obs_yield_dry"].mean() * 100),
            yield_bias=float(np.mean(bias)), bio_bias=float(np.mean(bio_b)),
            HI_sim=float(np.mean(his)), HI_obs=0.486,
            season_decline_pct=float(np.nanmean(tail))))
    return pd.DataFrame(rows)


def figure(res, lai_obs) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, plot in zip(axes, ["A1", "A3B2", "A5B2"]):
        o = lai_obs[plot].dropna()
        ax.plot(o.index, o.values, "o", ms=7, color="black", zorder=5,
                label="实测")
        for tag, ls, col, lab in (("s033", "--", "#1f77b4", "当前 sen=0.33"),
                                  ("s020", "-", "#d62728", "sen=0.20"),
                                  ("s010", ":", "#2ca02c", "sen=0.10")):
            df = res.get(f"{tag}__{plot}")
            if df is None:
                continue
            s = df.set_index("Clock.Today")["Maize.Leaf.LAI"]
            s = s["2026-06-25":"2026-10-10"]
            ax.plot(s.index, s.values, ls, lw=1.9, color=col, label=lab)
        ax.axvspan(pd.Timestamp("2026-08-14"), pd.Timestamp("2026-08-29"),
                   color="orange", alpha=0.12)
        ax.set_title(f"{plot}")
        ax.set_xlabel("日期 (2026)")
        ax.set_ylabel("LAI")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
        ax.set_ylim(0, 4.2)
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle("LAI 峰值后的平台期：实测（黑点）vs 模拟（橙带=实测平台区间）",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    p = FIG / "P1_LAI平台期.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    lai_obs, val = load_obs()
    print("实测平台期（08-14 → 08-29）：")
    for p in PLOTS:
        o = lai_obs[p].dropna()
        a, b = o.get(pd.Timestamp("2026-08-14")), o.get(pd.Timestamp("2026-08-29"))
        if a and b:
            print(f"  {p:5s} {a:.2f} → {b:.2f}  ({(b/a-1)*100:+.1f}%)")
    print()
    res = run_all()
    tbl = score(res, lai_obs, val)
    pd.set_option("display.width", 200)
    print(tbl.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))
    tbl.to_csv(SIM / "lai_plateau_compare.csv", index=False)
    print("\n图:", figure(res, lai_obs))


if __name__ == "__main__":
    main()