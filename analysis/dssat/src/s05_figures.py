"""The four validation figures.

  D1_LAI.png      3x3 panels, observed LAI vs simulated curve, plateau shaded
  D2_生物量.png    3x3 panels, above-ground biomass
  D3_产量.png      per-plot bars + 1:1 scatter
  D4_氮响应.png    yield / biomass / harvest index against nitrogen rate

Usage
-----
    python analysis/dssat/src/s05_figures.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["SimHei", "Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 130, "savefig.dpi": 150,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    "axes.labelsize": 9.5, "axes.titlesize": 10, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
})
SIMCOL = "#c0392b"
OBSCOL = "black"
PLATEAU_COL = "#f6c85f"


def load():
    daily = pd.read_csv(dc.SIM / "dssat_daily.csv", parse_dates=["date"])
    val = pd.read_csv(dc.SIM / "validate_plots.csv").set_index("plot")
    obs = dc.load_obs()
    return daily, obs, val


def _grid():
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 11.2), sharex=True, sharey=True)
    return fig, axes


def _finish(fig, axes, title, ylabel):
    for ax in axes[-1]:
        ax.set_xlabel("日期 (2026)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    for ax in axes.ravel():
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle(title, fontsize=12.5, y=0.995)


def fig_lai(daily, obs, val) -> Path:
    fig, axes = _grid()
    ymax = max(4.4, float(daily["LAID"].max()) * 1.12)
    for ax in axes.ravel():
        ax.set_ylim(0, ymax)
    for ax, plot in zip(axes.ravel(), dc.ORDER):
        g = daily[daily["plot"] == plot].sort_values("date")
        o = obs[obs["plot"] == plot].dropna(subset=["lai_plot"])
        p = dc.plateau_metrics(g["date"], g["LAID"])
        ax.axvspan(p["start"], p["end"], color=PLATEAU_COL, alpha=0.35, lw=0,
                   label="平台期 (≥95% 峰值)")
        ax.plot(o["date"], o["lai_plot"], "o", ms=6.5, color=OBSCOL,
                zorder=5, label="实测")
        if o.empty:
            ax.text(0.5, 0.06, "无 LAI 观测", transform=ax.transAxes,
                    ha="center", fontsize=8.5, color="grey")
        ax.plot(g["date"], g["LAID"], lw=1.9, color=SIMCOL, label="模拟")
        ax.set_title(
            f"{plot}  N={dc.N_RATE[plot]}  "
            f"RMSE={val.loc[plot, 'lai_rmse']:.2f}  "
            f"平台 {p['span_days']}d ({p['change_pct']:+.1f}%)", fontsize=9)
    axes[0, 0].legend(frameon=False, loc="upper left", fontsize=7)
    _finish(fig, axes,
            f"{dc.PREFIX}1｜{dc.DISPLAY} {dc.CFG.variety}｜{dc.SCEN_LABEL}"
            f"　逐小区 LAI：实测 vs 模拟（阴影为模拟平台期）", "LAI")
    fig.tight_layout()
    p = dc.FIG / f"{dc.PREFIX}1_LAI.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_biomass(daily, obs, val) -> Path:
    fig, axes = _grid()
    # observed maxima reach ~25 t/ha in XX_SF, so the axis has to follow the data
    bmax = max(float(daily["CWAD"].max()), float(obs["bio_kg_ha"].max()))
    for ax in axes.ravel():
        ax.set_ylim(0, bmax * 1.08)
    for ax, plot in zip(axes.ravel(), dc.ORDER):
        g = daily[daily["plot"] == plot].sort_values("date")
        o = obs[(obs["plot"] == plot) & obs["bio_kg_ha"].notna()]
        good = o[~o["organ_split_suspect"].fillna(False).astype(bool)]
        bad = o[o["organ_split_suspect"].fillna(False).astype(bool)]
        ax.plot(good["date"], good["bio_kg_ha"], "o", ms=6.5, color=OBSCOL,
                zorder=5, label="实测")
        if len(bad):
            ax.plot(bad["date"], bad["bio_kg_ha"], "^", ms=7, mfc="none",
                    mec="grey", mew=1.2, zorder=5, label="实测（器官分配可疑）")
        ax.plot(g["date"], g["CWAD"], lw=1.9, color=SIMCOL, label="模拟")
        ax.set_title(f"{plot}  N={dc.N_RATE[plot]}  "
                     f"偏差 {val.loc[plot, 'bio_bias']:+.0f}%", fontsize=9)
    axes[0, 0].legend(frameon=False, loc="upper left", fontsize=7)
    _finish(fig, axes,
            f"{dc.PREFIX}2｜{dc.DISPLAY} {dc.CFG.variety}｜{dc.SCEN_LABEL}"
            f"　逐小区地上部生物量：实测 vs 模拟", "生物量 (kg/ha)")
    fig.tight_layout()
    p = dc.FIG / f"{dc.PREFIX}2_生物量.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_yield(val) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.0),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ob = val["obs_yield"].reindex(dc.ORDER).values
    sb = val["sim_yield"].reindex(dc.ORDER).values

    ax = axes[0]
    x = np.arange(len(dc.ORDER))
    ax.bar(x - 0.2, ob, 0.4, label="实测", color="#7f8c8d", edgecolor="k", lw=0.4)
    ax.bar(x + 0.2, sb, 0.4, label="模拟", color=SIMCOL, edgecolor="k", lw=0.4)
    for i, (a, b) in enumerate(zip(ob, sb)):
        ax.text(i, max(a, b) + 180, f"{(b / a - 1) * 100:+.0f}%",
                ha="center", fontsize=7.5)
    ax.set_xticks(x, [f"{p}\nN={dc.N_RATE[p]}" for p in dc.ORDER], fontsize=7.5)
    ax.set_ylabel("产量 (kg/ha, 干重)")
    ax.set_title("逐小区产量")
    ax.legend(frameon=False)
    ax.set_ylim(0, max(ob.max(), sb.max()) * 1.16)

    ax = axes[1]
    ax.scatter(ob, sb, s=85, c=SIMCOL, edgecolor="k", lw=0.5, zorder=3)
    lim = [min(ob.min(), sb.min()) * 0.9, max(ob.max(), sb.max()) * 1.06]
    ax.plot(lim, lim, "k--", lw=1.2)
    for p, a, b in zip(dc.ORDER, ob, sb):
        ax.annotate(p, (a, b), fontsize=7.5, xytext=(6, 3),
                    textcoords="offset points")
    rrmse = np.sqrt(np.mean((sb - ob) ** 2)) / ob.mean() * 100
    ax.text(0.04, 0.95,
            f"RRMSE = {rrmse:.2f}%\n偏差 = {(sb.sum() / ob.sum() - 1) * 100:+.2f}%",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(fc="white", ec="grey", alpha=0.9))
    ax.set_xlabel("实测产量 (kg/ha)")
    ax.set_ylabel("模拟产量 (kg/ha)")
    ax.set_title("1:1")

    fig.suptitle(f"{dc.PREFIX}3｜{dc.DISPLAY} {dc.CFG.variety}｜{dc.SCEN_LABEL}"
                 f"　逐小区产量：实测 vs 模拟", fontsize=12.5, y=1.01)
    fig.tight_layout()
    p = dc.FIG / f"{dc.PREFIX}3_产量.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_nitrogen(val) -> Path:
    n = val["n_rate"].reindex(dc.ORDER).values.astype(float)
    panels = [
        ("产量 (kg/ha, 干重)", "obs_yield", "sim_yield", "产量对氮的响应"),
        ("地上部生物量 (kg/ha)", "obs_bio", "sim_bio", "生物量对氮的响应"),
        ("收获指数 HI", "obs_hi", "sim_hi", "收获指数"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (ylabel, oc, sc, title) in zip(axes, panels):
        o = val[oc].reindex(dc.ORDER).values
        s = val[sc].reindex(dc.ORDER).values
        ax.plot(n, o, "o", ms=8, color=OBSCOL, label="实测", zorder=4)
        ax.plot(n, s, "s", ms=7, color=SIMCOL, label="模拟", zorder=4)
        grp_o = pd.Series(o).groupby(n).mean()
        grp_s = pd.Series(s).groupby(n).mean()
        ax.plot(grp_o.index, grp_o.values, "-", color=OBSCOL, lw=1.4, alpha=0.7)
        ax.plot(grp_s.index, grp_s.values, "--", color=SIMCOL, lw=1.4, alpha=0.7)
        ax.set_xlabel("施氮量 (kg N/ha)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
        lo = min(np.nanmin(o), np.nanmin(s))
        hi = max(np.nanmax(o), np.nanmax(s))
        pad = (hi - lo) * 0.18 or 1.0
        ax.set_ylim(lo - pad, hi + pad)
    d_o = (val["obs_yield"].reindex(dc.ORDER).iloc[-1]
           / val["obs_yield"].reindex(dc.ORDER).iloc[0] - 1) * 100
    d_s = (val["sim_yield"].reindex(dc.ORDER).iloc[-1]
           / val["sim_yield"].reindex(dc.ORDER).iloc[0] - 1) * 100
    axes[0].text(0.04, 0.94, f"A5B2/A1 产量比\n实测 {d_o:+.0f}%  模拟 {d_s:+.0f}%",
                 transform=axes[0].transAxes, va="top", fontsize=8.5,
                 bbox=dict(fc="white", ec="grey", alpha=0.9))
    fig.suptitle(f"{dc.PREFIX}4｜{dc.DISPLAY} {dc.CFG.variety}｜{dc.SCEN_LABEL}"
                 f"　氮素响应：模拟 vs 实测", fontsize=12.5, y=1.02)
    fig.tight_layout()
    p = dc.FIG / f"{dc.PREFIX}4_氮响应.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    ap.add_argument("--scenario", default="base", choices=list(dc.SCENARIOS))
    args = ap.parse_args()
    dc.set_trial(args.trial)
    dc.set_scenario(args.scenario)
    print(f"trial {dc.TRIAL}  ({dc.CFG.site} {dc.CFG.variety})  "
          f"scenario {dc.SCENARIO} ({dc.SCEN_LABEL})")

    daily, obs, val = load()
    made = [fig_lai(daily, obs, val), fig_biomass(daily, obs, val),
            fig_yield(val)]
    if val["n_rate"].nunique() >= 3:
        made.append(fig_nitrogen(val))
    else:
        print("  skipped D4: only one nitrogen level in this trial")
    for p in made:
        print(" ", p, f"({p.stat().st_size / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
