"""Figure: what the senescence rate can and cannot buy.

Three candidate rates are drawn against the observations so the trade-off between
the measured 08-14 -> 08-29 plateau and a realistic end-of-season collapse is
visible at a glance.

Usage:  python analysis/src/s23_decline_figure.py
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
    "axes.labelsize": 10, "axes.titlesize": 11, "legend.fontsize": 8,
})

CAL = json.load(io.open(SIM / "phenology_calibrated.json", encoding="utf-8"))["乐农87"]["params"]
GRO = json.load(io.open(SIM / "calibrated_cultivars.json", encoding="utf-8"))["乐农87"]["growth"]
SEN = "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue"

SENS = [(0.20, "#1f77b4", "sen=0.20（当前，保守）"),
        (0.35, "#2ca02c", "sen=0.35（折中）"),
        (0.50, "#d62728", "sen=0.50（后期真实）")]
PLOTS = ["A1", "A3B2", "A5B2"]
PLATEAU = (pd.Timestamp("2026-08-14"), pd.Timestamp("2026-08-29"))


def main() -> None:
    obs = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    obs = obs[obs["trial_id"] == "SD_SF"]
    lai_obs = obs.pivot_table(index="date", columns="plot", values="lai_plot")

    cases = []
    for sen, _, _ in SENS:
        for plot in PLOTS:
            c = SimCase.load(Path(SIM).parent / "apsim" / f"SD_SF_{plot}.apsimx")
            cases.append((c.set_commands({**CAL, **GRO, SEN: sen}),
                          f"s{int(sen*100):03d}__{plot}"))
    with tempfile.TemporaryDirectory() as td:
        res = run_many(cases, td)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    stats = {}
    for ax, plot in zip(axes[:2], PLOTS[:2]):
        o = lai_obs[plot].dropna()
        ax.plot(o.index, o.values, "o", ms=7.5, color="black", zorder=6, label="实测")
        for sen, col, lab in SENS:
            df = res.get(f"s{int(sen*100):03d}__{plot}")
            if df is None:
                continue
            df = df.copy()
            df["Clock.Today"] = pd.to_datetime(df["Clock.Today"])
            s = df.set_index("Clock.Today")["Maize.Leaf.LAI"]["2026-06-20":"2026-10-12"]
            ax.plot(s.index, s.values, lw=2, color=col, label=lab)
            e = s.loc["2026-10-01":].mean()
            stats.setdefault(plot, {})[sen] = ((e / s.max() - 1) * 100, e)
        ax.axvspan(*PLATEAU, color="orange", alpha=0.13, lw=0)
        ax.set_title(f"{plot}")
        ax.set_xlabel("日期 (2026)")
        ax.set_ylabel("LAI")
        ax.set_ylim(0, 4.1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
    axes[0].legend(frameon=False, loc="upper left")

    # right panel: the trade-off
    ax = axes[2]
    tbl = pd.DataFrame([
        dict(sen=0.20, plateau=-5.2, decline=-29.3),
        dict(sen=0.35, plateau=-8.5, decline=-47.3),
        dict(sen=0.50, plateau=-11.3, decline=-59.7)])
    ax.plot(tbl["plateau"], tbl["decline"], "o-", lw=1.6, ms=9, color="grey", zorder=2)
    for (sen, col, lab), (_, r) in zip(SENS, tbl.iterrows()):
        ax.scatter(r["plateau"], r["decline"], s=150, color=col, zorder=3,
                   edgecolor="k", lw=0.6)
        ax.annotate(f"sen={sen:.2f}", (r["plateau"], r["decline"]), fontsize=8.5,
                    xytext=(8, 5), textcoords="offset points")
    ax.axvline(-3.2, color="crimson", ls="--", lw=1.3, label="实测平台期 -3.2%")
    ax.axhspan(-65, -55, color="green", alpha=0.13, lw=0)
    ax.text(-3.0, -60, "典型玉米季末衰减 55-65%", fontsize=8, color="darkgreen")
    ax.set_xlabel("模拟平台期变化 (%)")
    ax.set_ylabel("模拟季末衰减 (%)")
    ax.set_title("两个目标不在同一条曲线上")
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("D1  衰老速率决定一切：平台期与季末陡降不可兼得", fontsize=13, y=1.01)
    fig.tight_layout()
    p = FIG / "D1_平台期与衰减取舍.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)

    print("季末衰减 / 季末 LAI：")
    for plot, d in stats.items():
        print(f"  {plot}: " + "  ".join(f"sen{s:.2f} {v[0]:+.0f}% (LAI {v[1]:.2f})"
                                        for s, v in sorted(d.items())))
    print("\n图:", p)


if __name__ == "__main__":
    main()