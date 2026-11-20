"""Leaf-number check: simulated L#SD against the recorded "N 叶期" dates.

The phenology table records the date at which each leaf stage was reached
(4叶期, 6叶期, ...).  CERES-Maize writes the appeared leaf number as `L#SD`,
so the two can be compared directly.

Outputs
-------
out/<trial>/diag/leaf_number.csv
fig/<trial>/diag/H_叶片数.png

Usage
-----
    python analysis/dssat/src/s12_leafnumber.py --trial SD_MY73
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402


def _diag(trial: str) -> Path:
    d = dc.BASE / "out" / trial / "diag"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fig(trial: str) -> Path:
    d = dc.BASE / "fig" / trial / "diag"
    d.mkdir(parents=True, exist_ok=True)
    return d


def obs_leaf_stages() -> pd.DataFrame:
    ph = pd.read_csv(dc.DATASET / "phenology.csv", parse_dates=["date"])
    site = str(dc.CFG.site)
    sub = ph[ph["site"] == site]
    if sub.empty:                      # 山东潍坊 -> filed under 山东德州
        sub = ph[ph["site"].astype(str).str.startswith(site[:2])]
    if sub.empty:                      # 河南新乡 -> filed under 新乡
        sub = ph[ph["site"].astype(str).str.endswith(site[-2:])]
    ph = sub[sub["variety"] == dc.CFG.variety]
    rows = []
    for r in ph.itertuples():
        m = re.match(r"\s*(\d+)\s*叶", str(r.stage))
        if m:
            rows.append(dict(stage=str(r.stage).strip(),
                             leaf=int(m.group(1)), date=pd.Timestamp(r.date)))
    return pd.DataFrame(rows).sort_values(["leaf", "date"])


def sim_reach(daily: pd.DataFrame, leaf: int) -> pd.Series:
    """First date each plot reaches `leaf` appeared leaves."""
    out = {}
    for plot, g in daily.groupby("plot"):
        g = g.sort_values("date")
        hit = g[g["L#SD"] >= leaf]
        out[plot] = hit["date"].iloc[0] if len(hit) else pd.NaT
    return pd.Series(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_MY73", choices=list(dc.TRIALS))
    args = ap.parse_args()
    dc.set_trial(args.trial)

    daily = pd.read_csv(dc.SIM / "dssat_daily.csv", parse_dates=["date"])
    obs = obs_leaf_stages()
    if obs.empty:
        raise SystemExit(f"{dc.TRIAL}: no leaf-stage records for "
                         f"{dc.CFG.site} {dc.CFG.variety}")

    # observed: one date per leaf stage (some stages repeat, keep the first)
    obs = obs.drop_duplicates("leaf", keep="first").set_index("leaf")
    rows = []
    for leaf, r in obs.iterrows():
        d_sim = sim_reach(daily, int(leaf))
        rows.append(dict(leaf=int(leaf), obs_date=r["date"],
                         obs_stage=r["stage"],
                         sim_date=d_sim.median(),
                         sim_min=d_sim.min(), sim_max=d_sim.max(),
                         d_days=(d_sim.median() - r["date"]).days
                         if pd.notna(d_sim.median()) else np.nan))
    tab = pd.DataFrame(rows)
    tab.to_csv(_diag(dc.TRIAL) / "leaf_number.csv", index=False)
    print(f"{dc.TRIAL}  {dc.CFG.site} {dc.CFG.variety}")
    print(tab[["leaf", "obs_date", "sim_date", "d_days"]].to_string(index=False))
    print(f"\n  平均误差 {tab.d_days.mean():+.1f} d，"
          f"最大 |误差| {tab.d_days.abs().max():.0f} d")

    # final leaf number
    final = daily.groupby("plot")["L#SD"].max()
    print(f"  最终叶数 模拟 {final.median():.0f}（{final.min():.0f}–{final.max():.0f}）"
          f"，实测记录到 {int(tab.leaf.max())} 叶")

    _figure(daily, obs, tab)


def _figure(daily: pd.DataFrame, obs: pd.DataFrame, tab: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": 0.25, "grid.linestyle": "--",
                         "axes.labelsize": 10, "axes.titlesize": 11})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2),
                             gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axes[0]
    for plot, g in daily.groupby("plot"):
        ax.plot(g["date"], g["L#SD"], lw=1.2, alpha=0.75,
                label="模拟" if plot == dc.ORDER[0] else None)
    ax.plot(obs["date"], obs.index, "o", ms=8, color="k", zorder=5,
            label="实测（N 叶期）")
    for leaf, r in obs.iterrows():
        ax.annotate(f"{leaf}叶", (r["date"], leaf), xytext=(4, -11),
                    textcoords="offset points", fontsize=7.5, color="grey")
    ax.set_xlabel("日期 (2026)")
    ax.set_ylabel("叶片数")
    ax.set_title("A  出叶进程：模拟 L#SD vs 实测 N 叶期")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    ax.plot([-6, 6], [-6, 6], color="grey", lw=1)
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.scatter(tab["d_days"], tab["leaf"], s=80, c="#c0392b",
               edgecolor="k", lw=0.5, zorder=3)
    for r in tab.itertuples():
        ax.annotate(f"{r.leaf}叶", (r.d_days, r.leaf), xytext=(5, -3),
                    textcoords="offset points", fontsize=7.5)
    ax.set_xlabel("模拟 - 实测 (天)")
    ax.set_ylabel("叶序")
    ax.set_title(f"B  到达各叶位的日期差（平均 {tab.d_days.mean():+.1f} d）")
    ax.set_xlim(-8, 8)

    fig.suptitle(f"H｜{dc.DISPLAY} {dc.CFG.variety}｜叶片数模拟 vs 实测", fontsize=12)
    fig.tight_layout()
    p = _fig(dc.TRIAL) / "H_叶片数.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  figure ->", p)


if __name__ == "__main__":
    main()
