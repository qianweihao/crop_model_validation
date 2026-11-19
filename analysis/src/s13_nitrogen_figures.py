"""Step 13 - figures for the nitrogen-stress response in the Shandong trial.

FN (the photosynthesis nitrogen factor) is a constant 1.0 in APSIM's maize
model, so it cannot be used to show the nitrogen response.  The stress instead
acts on post-anthesis grain filling: the zero-N control sets the same grain
number as the fertilised plots but fills each grain less.  These panels show
that pathway using the variables that actually differ.

  N1  individual grain weight through grain filling (all 9 treatments)
  N2  grain nitrogen accumulation
  N3  above-ground nitrogen uptake
  N4  soil nitrate in the top layer
  N5  summary: yield / grain N / HI against nitrogen rate
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
    "axes.unicode_minus": False,
    "figure.dpi": 130,
    "savefig.dpi": 150,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

TRIAL = "SD_SF"
ORDER = ["A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"]
CMAP = plt.get_cmap("viridis")


def colour(plot: str) -> tuple:
    return CMAP(ORDER.index(plot) / (len(ORDER) - 1))


def load():
    sim = pd.read_csv(SIM / "sim_calibrated_daily.csv", parse_dates=["Clock.Today"])
    sim = sim[sim["trial_id"] == TRIAL].copy()
    fert = pd.read_csv(DS / "fert.csv")
    n = fert[fert.trial_id == TRIAL].groupby("plot")["n_kg_ha"].sum()
    return sim, n


def per_grain_mg(g: pd.DataFrame) -> pd.Series:
    n = g["Maize.Grain.Number"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(n > 0, g["Maize.Grain.Total.Wt"].to_numpy() / n * 1000.0, np.nan)


def fig_grain_weight(sim: pd.DataFrame, n: pd.Series) -> Path:
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for plot in ORDER:
        g = sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")
        y = per_grain_mg(g)
        m = np.isfinite(y)
        ax.plot(g["Clock.Today"][m], y[m], lw=1.8, color=colour(plot),
                label=f"{plot} (N={n.get(plot, 0):.0f})")
    ax.set_xlabel("日期 (2026)")
    ax.set_ylabel("单粒重 (mg)")
    ax.set_title("灌浆期单粒重：氮效应的直接体现")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.tick_params(axis="x", rotation=45)
    ax.legend(ncol=3, frameon=False, loc="upper left")
    ax.axvspan(pd.Timestamp("2026-09-10"), pd.Timestamp("2026-10-01"),
               color="orange", alpha=0.10)
    ax.annotate("分化主要发生在此区间", xy=(pd.Timestamp("2026-09-18"), 150),
                xytext=(pd.Timestamp("2026-08-22"), 90), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", lw=0.8, color="grey"))
    fig.tight_layout()
    p = FIG / "N1_单粒重.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_grain_n(sim: pd.DataFrame, n: pd.Series) -> Path:
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for plot in ORDER:
        g = sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")
        v = g["Maize.Grain.N"] * 10.0            # g/m2 -> kg/ha
        ax.plot(g["Clock.Today"], v, lw=1.8, color=colour(plot),
                label=f"{plot} (N={n.get(plot, 0):.0f})")
    ax.set_xlabel("日期 (2026)")
    ax.set_ylabel("籽粒氮累积 (kg N/ha)")
    ax.set_title("籽粒氮累积")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.tick_params(axis="x", rotation=45)
    ax.legend(ncol=3, frameon=False, loc="upper left")
    fig.tight_layout()
    p = FIG / "N2_籽粒氮累积.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_n_uptake(sim: pd.DataFrame, n: pd.Series) -> Path:
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    for plot in ORDER:
        g = sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")
        ax.plot(g["Clock.Today"], g["Maize.AboveGround.N"] * 10.0, lw=1.8,
                color=colour(plot), label=f"{plot} (N={n.get(plot, 0):.0f})")
    ax.set_xlabel("日期 (2026)")
    ax.set_ylabel("地上部吸氮量 (kg N/ha)")
    ax.set_title("作物吸氮量")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.tick_params(axis="x", rotation=45)
    ax.legend(ncol=3, frameon=False, loc="upper left")
    fig.tight_layout()
    p = FIG / "N3_地上部吸氮.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_soil_no3(sim: pd.DataFrame, n: pd.Series) -> Path:
    cols = [c for c in sim.columns if c.startswith("Soil.Nutrient.NO3.kgha")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, (label, subset) in zip(
            axes, [("表层 0–30 cm", cols[:2]), ("剖面合计 0–180 cm", cols)]):
        for plot in ORDER:
            g = sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")
            ax.plot(g["Clock.Today"], g[subset].sum(axis=1), lw=1.6,
                    color=colour(plot), label=f"{plot} (N={n.get(plot, 0):.0f})")
        ax.set_xlabel("日期 (2026)")
        ax.set_ylabel("NO3-N (kg/ha)")
        ax.set_title(label)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
    axes[0].legend(ncol=2, frameon=False, fontsize=6.5)
    fig.suptitle("土壤硝态氮：施肥的残留差异", fontsize=13, y=1.02)
    fig.tight_layout()
    p = FIG / "N4_土壤硝态氮.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_summary(sim: pd.DataFrame, n: pd.Series) -> Path:
    rows = []
    for plot in ORDER:
        g = sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")
        i = g["Maize.Grain.Total.Wt"].idxmax()
        r = g.loc[i]
        gn = r["Maize.Grain.Number"]
        rows.append(dict(plot=plot, N=n.get(plot, 0.0),
                         yield_kg=float(r["Maize.Grain.Total.Wt"] * 10),
                         grain_n=float(r["Maize.Grain.N"] * 10),
                         per_grain=float(r["Maize.Grain.Total.Wt"] / gn * 1000) if gn > 0 else np.nan,
                         ag_n=float(g["Maize.AboveGround.N"].max() * 10)))
    d = pd.DataFrame(rows).sort_values("N")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    ax = axes[0]
    ax.plot(d["N"], d["yield_kg"], "o-", lw=1.8, ms=7, color="#1f77b4")
    for _, r in d.iterrows():
        ax.annotate(r["plot"], (r["N"], r["yield_kg"]), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("施氮量 (kg N/ha)")
    ax.set_ylabel("产量 (kg/ha, 干重)")
    ax.set_title("产量")

    ax = axes[1]
    ax.plot(d["N"], d["per_grain"], "s-", lw=1.8, ms=7, color="#d62728")
    for _, r in d.iterrows():
        ax.annotate(r["plot"], (r["N"], r["per_grain"]), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("施氮量 (kg N/ha)")
    ax.set_ylabel("单粒重 (mg)")
    ax.set_title("单粒重")

    ax = axes[2]
    ax.plot(d["N"], d["grain_n"], "^-", lw=1.8, ms=7, color="#2ca02c", label="籽粒氮")
    ax.plot(d["N"], d["ag_n"], "v--", lw=1.6, ms=6, color="#ff7f0e", label="地上部吸氮")
    ax.set_xlabel("施氮量 (kg N/ha)")
    ax.set_ylabel("kg N/ha")
    ax.set_title("氮累积")
    ax.legend(frameon=False)

    fig.suptitle("山东德州 乐农87：氮响应的三个层面", fontsize=13, y=1.02)
    fig.tight_layout()
    p = FIG / "N5_氮响应汇总.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p, d


def main() -> None:
    sim, n = load()
    made = [fig_grain_weight(sim, n), fig_grain_n(sim, n),
            fig_n_uptake(sim, n), fig_soil_no3(sim, n)]
    p5, d = fig_summary(sim, n)
    made.append(p5)
    d.to_csv(SIM / "nitrogen_response.csv", index=False)
    print("氮响应汇总：")
    print(d.round(1).to_string(index=False))
    print()
    for p in made:
        print(" ", p)


if __name__ == "__main__":
    main()
