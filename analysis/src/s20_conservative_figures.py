"""Comparison figures, produced twice under two names.

The same four drawings are emitted as the C series and as the L series, so either
naming can be cited from the reports without the content drifting apart:

  x1  LAI              nine per-plot panels, observed vs simulated
  x2  地上部生物量       nine per-plot panels, observed vs simulated
  x3  产量              per-plot bars with bias, plus a 1:1 scatter
  x4  氮响应             yield / biomass / HI against nitrogen rate

Usage
-----
    python analysis/src/s20_conservative_figures.py   # C series
    python analysis/src/s12_figures.py               # L series
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DS, FIG, SIM  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["SimHei", "Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 130, "savefig.dpi": 150,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    "axes.labelsize": 9.5, "axes.titlesize": 10, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
})

TRIAL = "SD_SF"
ORDER = ["A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"]
N_RATE = dict(zip(ORDER, [0, 139, 139, 177, 177, 215, 215, 254, 254]))
SIMCOL = "#c0392b"
OBSCOL = "black"
PLANTS_PER_MU = 667.0        # density is plants per mu (1 mu = 667 m2)


# --------------------------------------------------------------------------
def load():
    sim = pd.read_csv(SIM / "sim_calibrated_daily.csv", parse_dates=["Clock.Today"])
    sim = sim[sim["trial_id"] == TRIAL]
    obs = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    obs = obs[obs["trial_id"] == TRIAL].copy()
    obs["bio_kg_ha"] = obs["total_dw"] * obs["density_measured"] / PLANTS_PER_MU * 10.0
    val = pd.read_csv(SIM / "validate_plots.csv")
    return sim, obs, val[val["trial_id"] == TRIAL].set_index("plot")


def _series(sim, plot) -> pd.DataFrame:
    return sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")


def _grid():
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 11.2), sharex=True, sharey=True)
    return fig, axes


def _finish(fig, axes, title, ylabel) -> None:
    for ax in axes[-1]:
        ax.set_xlabel("日期 (2026)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    for ax in axes.ravel():
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle(title, fontsize=13.5, y=0.995)


# --------------------------------------------------------------------------
def fig_lai(sim, obs, val, prefix: str, note: str = "") -> Path:
    fig, axes = _grid()
    for ax in axes.ravel():
        ax.set_ylim(0, 4.2)
    for ax, plot in zip(axes.ravel(), ORDER):
        o = obs[(obs["plot"] == plot)].dropna(subset=["lai_plot"])
        ax.plot(o["date"], o["lai_plot"], "o", ms=6.5, color=OBSCOL, zorder=5, label="实测")
        g = _series(sim, plot)
        ax.plot(g["Clock.Today"], g["Maize.Leaf.LAI"], lw=1.9, color=SIMCOL, label="模拟")
        ax.set_title(f"{plot}  N={N_RATE[plot]}  LAI RMSE={val.loc[plot, 'lai_rmse']:.2f}",
                     fontsize=9.5)
    axes[0, 0].legend(frameon=False, loc="upper left")
    _finish(fig, axes, f"{prefix}1  逐小区 LAI：实测 vs 模拟{note}", "LAI")
    fig.tight_layout()
    p = FIG / f"{prefix}1_LAI.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_biomass(sim, obs, prefix: str, note: str = "") -> Path:
    fig, axes = _grid()
    for ax in axes.ravel():
        ax.set_ylim(0, 20000)
    for ax, plot in zip(axes.ravel(), ORDER):
        o = obs[(obs["plot"] == plot) & obs["bio_kg_ha"].notna()]
        good = o[~o["organ_split_suspect"].fillna(False).astype(bool)]
        bad = o[o["organ_split_suspect"].fillna(False).astype(bool)]
        ax.plot(good["date"], good["bio_kg_ha"], "o", ms=6.5, color=OBSCOL,
                zorder=5, label="实测")
        if len(bad):
            ax.plot(bad["date"], bad["bio_kg_ha"], "^", ms=7, mfc="none",
                    mec="grey", mew=1.2, zorder=5, label="实测（器官分配可疑）")
        g = _series(sim, plot)
        ax.plot(g["Clock.Today"], g["Maize.AboveGround.Wt"] * 10, lw=1.9,
                color=SIMCOL, label="模拟")
        ax.set_title(f"{plot}  N={N_RATE[plot]}", fontsize=9.5)
    axes[0, 0].legend(frameon=False, loc="upper left", fontsize=7)
    _finish(fig, axes, f"{prefix}2  逐小区地上部生物量：实测 vs 模拟{note}", "生物量 (kg/ha)")
    fig.tight_layout()
    p = FIG / f"{prefix}2_生物量.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_yield(val, prefix: str, note: str = "") -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.0),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ob = val["obs_yield_dry"].reindex(ORDER).values
    sb = val["sim_yield"].reindex(ORDER).values

    ax = axes[0]
    x = np.arange(len(ORDER))
    ax.bar(x - 0.2, ob, 0.4, label="实测", color="#7f8c8d", edgecolor="k", lw=0.4)
    ax.bar(x + 0.2, sb, 0.4, label="模拟", color=SIMCOL, edgecolor="k", lw=0.4)
    for i, (a, b) in enumerate(zip(ob, sb)):
        ax.text(i, max(a, b) + 180, f"{(b / a - 1) * 100:+.0f}%", ha="center", fontsize=7.5)
    ax.set_xticks(x, [f"{p}\nN={N_RATE[p]}" for p in ORDER], fontsize=7.5)
    ax.set_ylabel("产量 (kg/ha, 干重)")
    ax.set_title("逐小区产量")
    ax.legend(frameon=False)
    ax.set_ylim(0, max(ob.max(), sb.max()) * 1.16)

    ax = axes[1]
    ax.scatter(ob, sb, s=85, c=SIMCOL, edgecolor="k", lw=0.5, zorder=3)
    lim = [min(ob.min(), sb.min()) * 0.9, max(ob.max(), sb.max()) * 1.06]
    ax.plot(lim, lim, "k--", lw=1.2)
    for p, a, b in zip(ORDER, ob, sb):
        ax.annotate(p, (a, b), fontsize=7.5, xytext=(6, 3), textcoords="offset points")
    rrmse = np.sqrt(np.mean((sb - ob) ** 2)) / ob.mean() * 100
    ax.text(0.04, 0.95, f"RRMSE = {rrmse:.2f}%\n偏差 = {(sb.sum() / ob.sum() - 1) * 100:+.2f}%",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(fc="white", ec="grey", alpha=0.9))
    ax.set_xlabel("实测产量 (kg/ha)")
    ax.set_ylabel("模拟产量 (kg/ha)")
    ax.set_title("1:1")

    fig.suptitle(f"{prefix}3  逐小区产量：实测 vs 模拟{note}", fontsize=13.5, y=1.01)
    fig.tight_layout()
    p = FIG / f"{prefix}3_产量.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_nitrogen_response(val, prefix: str, note: str = "") -> Path:
    """Yield, biomass and harvest index against nitrogen rate (the old L3 panel)."""
    n = np.array([N_RATE[p] for p in ORDER], dtype=float)
    panels = [
        ("产量 (kg/ha, 干重)", "obs_yield_dry", "sim_yield", "产量对氮的响应"),
        ("地上部生物量 (kg/ha)", "obs_bio_max", "sim_bio_max", "生物量对氮的响应"),
        ("收获指数 HI", "obs_HI", "sim_HI", "收获指数"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (ylabel, oc, sc, title) in zip(axes, panels):
        o = val[oc].reindex(ORDER).values
        s = val[sc].reindex(ORDER).values
        ax.plot(n, o, "o-", ms=8, lw=1.8, color=OBSCOL, label="实测", zorder=4)
        ax.plot(n, s, "s--", ms=7, lw=1.8, color=SIMCOL, label="模拟", zorder=4)
        ax.set_xlabel("施氮量 (kg N/ha)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
        lo = min(np.nanmin(o), np.nanmin(s))
        hi = max(np.nanmax(o), np.nanmax(s))
        pad = (hi - lo) * 0.18 or 1.0
        ax.set_ylim(lo - pad, hi + pad)
    d_o = (val["obs_yield_dry"].reindex(ORDER).iloc[-1] / val["obs_yield_dry"].reindex(ORDER).iloc[0] - 1) * 100
    d_s = (val["sim_yield"].reindex(ORDER).iloc[-1] / val["sim_yield"].reindex(ORDER).iloc[0] - 1) * 100
    axes[0].text(0.04, 0.94, f"0→254 kg N/ha\n实测 {d_o:+.0f}%  模拟 {d_s:+.0f}%",
                 transform=axes[0].transAxes, va="top", fontsize=8.5,
                 bbox=dict(fc="white", ec="grey", alpha=0.9))
    fig.suptitle(f"{prefix}4  氮素响应：模拟 vs 实测{note}", fontsize=13.5, y=1.02)
    fig.tight_layout()
    p = FIG / f"{prefix}4_氮响应.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def make_all(prefix: str, sim=None, obs=None, val=None, note: str = "") -> list[Path]:
    """Draw the four figures.  Data can be injected so a second parameter set
    (the pre-optimisation baseline) can reuse exactly the same code path."""
    if sim is None or obs is None or val is None:
        sim, obs, val = load()
    return [fig_lai(sim, obs, val, prefix, note), fig_biomass(sim, obs, prefix, note),
            fig_yield(val, prefix, note), fig_nitrogen_response(val, prefix, note)]


def main() -> None:
    sim, obs, val = load()
    print(f"sim rows={len(sim)}  plots={sim['sim'].nunique()}")
    for p in make_all("C", note="（sen=0.20，优化后）"):
        print(" ", p)
    o, s = val["obs_yield_dry"], val["sim_yield"]
    print(f"\n产量 RRMSE {np.sqrt(np.mean((s - o) ** 2)) / o.mean() * 100:.2f}%  "
          f"偏差 {(s.sum() / o.sum() - 1) * 100:+.2f}%  "
          f"LAI RMSE {val['lai_rmse'].mean():.3f}")


if __name__ == "__main__":
    main()