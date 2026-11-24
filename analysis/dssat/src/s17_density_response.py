"""SD_SF：密度响应曲线（用于生产应用时的"密度敏感度"评估）。

三个施氮水平 x 7 个密度，看 CERES 的产量/LAI/单穗重如何响应。

Outputs
-------
out/SD_SF/diag/density_response.csv
fig/SD_SF/diag/L_密度响应.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

DENS = (3000, 3750, 4500, 5250, 6000, 6750, 7500)
PLOTS = ("A1", "A2B1", "A5B2")          # 0 / 139 / 254 kg N/ha


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="all", choices=["all", "table", "figure"])
    args = ap.parse_args()
    dc.set_trial("SD_SF")
    mu = dc.PLANTS_PER_MU

    rows = []
    for plot in PLOTS:
        for d_mu in DENS:
            dens = d_mu / mu
            g, s = dc.run_plot(plot, density=dens)
            m = dc.plot_metrics(g, s, dens)
            rows.append(dict(
                plot=plot, N=dc.N_RATE[plot], dens_mu=d_mu,
                yield_mu=m["yield_kg"] / 15,            # 干重 kg/亩
                bio_mu=m["bio_final"] / 15,
                lai_peak=m["lai_peak"],
                w_ear=m["yield_kg"] / (dens * 10000) * 1000,   # g/穗
                grains=m["grains_per_ear"], hi=m["hi"]))
    d = pd.DataFrame(rows)
    d.to_csv(dc.OUT / "diag" / "density_response.csv", index=False)

    if args.section in ("all", "table"):
        for p, g in d.groupby("plot", sort=False):
            base = g.loc[g.dens_mu == 3750, "yield_mu"].iloc[0]
            g = g.assign(rel=g.yield_mu / base)
            print(f"--- {p}  N {g.N.iloc[0]:.0f} kg/ha ---")
            print(g[["dens_mu", "yield_mu", "rel", "lai_peak", "w_ear", "hi"]]
                  .round(2).to_string(index=False))
    if args.section in ("all", "figure"):
        _figure(d)


def _figure(d: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": 0.25, "grid.linestyle": "--",
                         "legend.fontsize": 8.5})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    cols = ["#2980b9", "#e67e22", "#c0392b"]
    for (p, g), c in zip(d.groupby("plot", sort=False), cols):
        lab = f"{p}（{g.N.iloc[0]:.0f} kg N/ha）"
        axes[0].plot(g.dens_mu, g.yield_mu, "o-", color=c, lw=1.8, ms=5, label=lab)
        axes[1].plot(g.dens_mu, g.lai_peak, "o-", color=c, lw=1.8, ms=5, label=lab)
        axes[2].plot(g.dens_mu, g.w_ear, "o-", color=c, lw=1.8, ms=5, label=lab)
    axes[0].set_ylabel("产量 kg/亩（干重）")
    axes[0].set_title("A  产量-密度")
    axes[1].set_ylabel("峰值 LAI")
    axes[1].set_title("B  LAI-密度")
    axes[2].set_ylabel("单穗粒重 g/穗")
    axes[2].set_title("C  单穗重-密度")
    for ax in axes:
        ax.set_xlabel("密度 株/亩")
        ax.axvspan(3757, 4675, color="grey", alpha=0.12)
        ax.axvline(5640, color="k", ls=":", lw=1.0)
        ax.legend(frameon=False, loc="best")
    axes[0].annotate("灰带=标定密度\n虚线=播种密度", xy=(0.02, 0.02),
                     xycoords="axes fraction", fontsize=8, va="bottom")
    fig.suptitle("L｜CERES 密度响应：3757→5640 株/亩 仅 +9%，低氮下反向", fontsize=12)
    fig.tight_layout()
    p = dc.BASE / "fig" / "SD_SF" / "diag" / "L_密度响应.png"
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  figure ->", p)


if __name__ == "__main__":
    main()
