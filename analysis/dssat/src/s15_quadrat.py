"""德州：小区实测产量 vs 10 米双行实测 vs 两种密度下的模拟（kg/亩）。

原始记录（`产量记录` sheet，山东德州段）有三套产量：
  测产密度（株/亩）
  小区实测产量（kg/亩）        —— 样方实收，对应实际留苗
  10 米双行实测（kg/亩）       —— 按 5500 株/亩 折算的单穗产量
模拟给出两套：实测留苗密度、播种密度（5640 株/亩）。

Outputs
-------
out/SD_SF/diag/yield_density_compare.csv
fig/SD_SF/diag/K_产量对比_kg亩.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

# 原始记录：data/.../玉米试验数据分类整理_施肥灌溉记录_产量.xlsx -> 产量记录
# (测产密度, 小区实测产量 kg/亩, 10 米双行实测 kg/亩)
RAW = {
    "A1":   (4315, 471.968453683525, 522.6),
    "A2B1": (4568, 587.553686878378, 608.3),
    "A2B2": (4675, 592.244959010204, 612.4),
    "A3B1": (4626, 645.934961958904, 698.1),
    "A3B2": (4636, 671.916833625,    714.4),
    "A4B1": (4003, 606.547709485714, 775.7),
    "A4B2": (4026, 637.5259335,      812.4),
    "A5B1": (3903, 637.70948119708,  861.4),
    "A5B2": (3757, 670.1657473125,   934.9),
}


def rrmse(s, o):
    s, o = np.asarray(s, float), np.asarray(o, float)
    return float(np.sqrt(((s - o) ** 2).mean()) / o.mean() * 100)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.parse_args()
    dc.set_trial("SD_SF")
    meas = pd.read_csv(dc.SIM / "dssat_summary.csv").set_index("plot")
    nom = pd.read_csv(dc.OUT / "sim_nominal" / "dssat_summary.csv").set_index("plot")

    rows = []
    for plot, (dens, plot_y, quad_y) in RAW.items():
        rows.append(dict(
            plot=plot, dens_plot=dens,
            obs_plot=plot_y, obs_quad=quad_y,
            # HWAM 是干重 kg/ha；/15 得干重 kg/亩，/0.86 折回 14% 含水口径
            sim_meas=meas.loc[plot, "HWAM"] / 15 / 0.86,
            sim_nom=nom.loc[plot, "HWAM"] / 15 / 0.86))
    d = pd.DataFrame(rows)
    d["b_meas_vs_plot"] = (d.sim_meas / d.obs_plot - 1) * 100
    d["b_nom_vs_plot"] = (d.sim_nom / d.obs_plot - 1) * 100
    d["b_nom_vs_quad"] = (d.sim_nom / d.obs_quad - 1) * 100
    d.to_csv(dc.OUT / "diag" / "yield_density_compare.csv", index=False)

    print(d.round(1).to_string(index=False))
    print()
    for tag, s, o in [
        ("模拟@实测密度 vs 小区实测", d.sim_meas, d.obs_plot),
        ("模拟@播种密度 vs 小区实测", d.sim_nom, d.obs_plot),
        ("模拟@播种密度 vs 10米双行", d.sim_nom, d.obs_quad),
        ("模拟@实测密度 vs 10米双行", d.sim_meas, d.obs_quad),
    ]:
        print(f"  {tag:26s} RRMSE {rrmse(s, o):5.1f}%  "
              f"偏差 {np.sum(s) / np.sum(o) * 100 - 100:+5.1f}%")
    _figure(d)


def _figure(d: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": 0.25, "grid.linestyle": "--",
                         "legend.fontsize": 8.5})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6),
                             gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    x = np.arange(len(d))
    ax.bar(x - 0.3, d.obs_plot, 0.2, label="小区实测（实际留苗）",
           color="#7f8c8d", edgecolor="k", lw=0.4)
    ax.bar(x - 0.1, d.sim_meas, 0.2, label="模拟（实测留苗密度）",
           color="#c0392b", edgecolor="k", lw=0.4)
    ax.bar(x + 0.1, d.obs_quad, 0.2, label="10 米双行实测（折 5500 株/亩）",
           color="#95a5a6", edgecolor="k", lw=0.4)
    ax.bar(x + 0.3, d.sim_nom, 0.2, label="模拟（播种密度 5640 株/亩）",
           color="#e67e22", edgecolor="k", lw=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.plot}\n{r.dens_plot:.0f}" for r in d.itertuples()],
                       fontsize=8)
    ax.set_ylabel("产量 kg/亩（14% 含水）")
    ax.set_title("A  德州逐小区：实测 vs 模拟（两种密度）")
    ax.legend(frameon=False, ncol=2, loc="upper left")

    ax = axes[1]
    ax.scatter(d.obs_plot, d.sim_meas, s=80, c="#c0392b", edgecolor="k",
               lw=0.5, label="实测密度 vs 小区实测")
    ax.scatter(d.obs_quad, d.sim_nom, s=80, c="#e67e22", edgecolor="k",
               lw=0.5, marker="s", label="播种密度 vs 10 米双行")
    lim = [0, max(d.obs_plot.max(), d.obs_quad.max(),
                  d.sim_nom.max()) * 1.1]
    ax.plot(lim, lim, "k--", lw=1.1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("实测 kg/亩")
    ax.set_ylabel("模拟 kg/亩")
    ax.set_title("B  1:1")
    ax.legend(frameon=False, loc="upper left")

    fig.suptitle("K｜德州 乐农87｜产量对比：kg/亩", fontsize=12)
    fig.tight_layout()
    p = dc.BASE / "fig" / "SD_SF" / "diag" / "K_产量对比_kg亩.png"
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  figure ->", p)


if __name__ == "__main__":
    main()
