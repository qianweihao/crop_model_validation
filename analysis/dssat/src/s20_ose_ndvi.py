"""OSE（观测系统实验）：NDVI 同化能否把 SLPF 估准、把产量预测改善？

设计
----
真值    : SLPF_true ∈ {1.00, 1.20, 1.35} 各跑一次，得到"真实"LAI 轨迹
观测    : NDVI_obs = NDVI(LAI_true) + N(0, sigma)，sigma=0.03（Sentinel-2 量级）
代价函数: 在 **NDVI 空间** 比较（避免 LAI 反演的非线性偏差）
同化    : SLPF 网格上最小二乘，等价于变分参数估计
评估    : 20 组噪声实现；SLPF 估计偏差/RMSE；同化前后产量误差
窗口    : 早窗 / 全季 / 晚窗（检验 NDVI 饱和的影响）

Outputs: out/SD_SF/diag/ose_ndvi.csv, fig/SD_SF/diag/M_OSE_NDVI同化.png
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc
from DSSATTools import DSSAT
from DSSATTools.base.sections import TabularSubsection

SLPF_GRID = np.round(np.arange(0.85, 1.56, 0.05), 2)
SLPF_TRUE = (1.00, 1.20, 1.35)
SLPF_PRIOR = 1.15
DENS_MU, N_RATE, NOISE = 5640, 254, 0.03
NREAL = 20
WINDOWS = {
    "早窗 06-28~07-23": ["2026-06-28", "2026-07-03", "2026-07-08",
                        "2026-07-13", "2026-07-18", "2026-07-23"],
    "全季 06-28~08-22": ["2026-06-28", "2026-07-03", "2026-07-08", "2026-07-13",
                        "2026-07-18", "2026-07-23", "2026-07-28", "2026-08-02",
                        "2026-08-07", "2026-08-12", "2026-08-17", "2026-08-22"],
    "晚窗 08-02~08-22": ["2026-08-02", "2026-08-07", "2026-08-12",
                        "2026-08-17", "2026-08-22"],
}


def lai_to_ndvi(lai):
    return 0.15 + 0.75 * (1 - np.exp(-0.5 * np.asarray(lai, float)))


def run(slpf, dens_mu=DENS_MU, n_rate=N_RATE):
    dc.set_trial("SD_SF")
    with dc.patch_species(**dc.CFG.species):
        soil = dc.build_soil("A5B2"); soil.SLPF = slpf
        wth = dc.build_weather(); crop = dc.build_crop(None)
        man = dc.build_management("A5B2", soil, nitrogen=True,
                                  nitrate_scale=dc.CFG.nitrate_scale,
                                  density=dens_mu / dc.PLANTS_PER_MU)
        t = man.fertilizers["table"].copy(); cur = float(t["FAMN"].sum())
        t["FAMN"] = (t["FAMN"] * n_rate / cur).round(2)
        man.fertilizers["table"] = TabularSubsection(t)
        d = DSSAT(); d.run(soil, wth, crop, man, verbose=False)
        g = d.output["PlantGro"].sort_index()
        y = float(g.iloc[-1]["GWAD"]) / 15 / 0.86
        lai = g["LAID"].astype(float); d.close()
    return y, lai


def main() -> None:
    ap = argparse.ArgumentParser(); ap.parse_args()
    rng = np.random.default_rng(2026)
    cache = {float(s): run(float(s)) for s in SLPF_GRID}
    obs = {}
    for tag, dates in WINDOWS.items():
        idx = [d for d in dates if pd.Timestamp(d) in cache[SLPF_TRUE[0]][1].index]
        obs[tag] = [pd.Timestamp(d) for d in idx]
    rows = []
    for t in SLPF_TRUE:
        lai_t = cache[t][1]
        for tag, idx in obs.items():
            ndvi_true = lai_to_ndvi(lai_t.loc[idx].values)
            for rep in range(NREAL):
                ndvi_obs = ndvi_true + rng.normal(0, NOISE, len(idx))
                errs = {s: float(np.sum((lai_to_ndvi(cache[s][1].loc[idx].values)
                                        - ndvi_obs) ** 2)) for s in cache}
                best = min(errs, key=errs.get)
                rows.append(dict(真值=t, 窗口=tag, 实现=rep, 估计=best,
                                 真值产量=cache[t][0], 先验产量=cache[SLPF_PRIOR][0],
                                 同化产量=cache[best][0]))
    d = pd.DataFrame(rows)
    d["先验误差"] = (d.先验产量 / d.真值产量 - 1) * 100
    d["同化误差"] = (d.同化产量 / d.真值产量 - 1) * 100
    d.to_csv(dc.OUT / "diag" / "ose_ndvi.csv", index=False)

    print("== SLPF 反演精度（20 组噪声实现）==")
    g = d.groupby(["窗口", "真值"]).agg(估计均值=("估计", "mean"), 估计std=("估计", "std")).reset_index()
    g["偏差"] = g.估计均值 - g.真值
    print(g.round(3).pivot(index="真值", columns="窗口", values=["估计均值", "偏差"]).to_string())
    print()
    print("== 产量误差（绝对值平均）==")
    s = d.groupby("窗口").agg(先验=("先验误差", lambda x: x.abs().mean()),
                              同化=("同化误差", lambda x: x.abs().mean()),
                              同化p90=("同化误差", lambda x: np.percentile(np.abs(x), 90)))
    print(s.round(2).to_string())
    s.to_csv(dc.OUT / "diag" / "ose_ndvi_summary.csv")
    _fig(d, cache, obs)


def _fig(d, cache, obs):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": .25, "grid.linestyle": "--", "legend.fontsize": 8})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    lai = cache[1.35][1]
    ax[0].plot(lai.index, lai.values, "k-", lw=1.8, label="LAI (SLPF 1.35)")
    for tag, c, mk in [("早窗 06-28~07-23", "#27ae60", "o"),
                       ("晚窗 08-02~08-22", "#7f8c8d", "s")]:
        e = obs[tag]
        ax[0].scatter(e, lai.loc[e].values, c=c, marker=mk, zorder=3, label=tag)
    ax[0].axhspan(4, 7, color="grey", alpha=.15)
    ax[0].text(pd.Timestamp("2026-07-25"), 4.3, "NDVI 饱和区", fontsize=8)
    ax[0].set_ylabel("LAI"); ax[0].set_title("A  观测窗口 vs LAI")
    ax[0].legend(frameon=False)
    for i, (tag, c) in enumerate(zip(obs, ["#27ae60", "#2980b9", "#7f8c8d"])):
        gg = d[d.窗口 == tag]
        ax[1].scatter(gg.真值, gg.估计, s=22, alpha=.5, c=c, label=tag)
    ax[1].plot([.85, 1.55], [.85, 1.55], "k--", lw=1)
    ax[1].set_xlabel("真值 SLPF"); ax[1].set_ylabel("估计 SLPF")
    ax[1].set_title("B  SLPF 反演（20 组噪声）"); ax[1].legend(frameon=False, loc="upper left")
    x = np.arange(len(SLPF_TRUE)); w = .26
    for i, (tag, c) in enumerate(zip(obs, ["#27ae60", "#2980b9", "#7f8c8d"])):
        gg = d[d.窗口 == tag].groupby("真值")
        ax[2].bar(x + (i - 1) * w, gg.同化误差.apply(lambda s: s.abs().mean()).values,
                  w, label=tag, color=c, edgecolor="k", lw=.4)
    prior = d.groupby("真值").先验误差.apply(lambda s: s.abs().mean()).values
    ax[2].plot(x, prior, "r^", ms=9, label="先验（不同化）")
    ax[2].set_xticks(x); ax[2].set_xticklabels([f"{v:.2f}" for v in SLPF_TRUE])
    ax[2].set_xlabel("真值 SLPF"); ax[2].set_ylabel("产量误差绝对值 %")
    ax[2].set_title("C  同化前后产量误差"); ax[2].legend(frameon=False)
    fig.suptitle("M｜OSE：NDVI 同化估 SLPF（5640 株/亩、254 kg N/ha、NDVI 噪声 σ=0.03）", fontsize=12)
    fig.tight_layout()
    p = dc.BASE / "fig" / "SD_SF" / "diag" / "M_OSE_NDVI同化.png"
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=200, bbox_inches="tight"); plt.close(fig)
    print("  figure ->", p)


if __name__ == "__main__":
    main()
