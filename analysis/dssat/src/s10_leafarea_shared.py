"""Tune the per-plant leaf area (IXIM AX) and re-fit a shared MY73 cultivar.

1. leaf area : CERES has no leaf-area parameter, IXIM has cultivar-level AX/ALL.
   Scan AX per trial, pick the value that matches the observed peak per-plant
   leaf area, and draw observed vs simulated for every plot.
2. shared cultivar : one MY73 parameter set for both sites, scored on yield and
   above-ground biomass (phenology may be off by up to 4 days).

Outputs
-------
out/<trial>/diag/per_plant_la.csv         Q1 per plot: obs / CERES / IXIM
out/<trial>/diag/ixim_ax_grid.csv         Q1 AX scan
fig/<trial>/diag/F_单株叶面积.png
out/XX_SF/diag/shared_cultivar.csv        Q2 joint grid
fig/XX_SF/diag/G_共用品种参数.png

Usage
-----
    python analysis/dssat/src/s10_leafarea_shared.py --section all
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

MU = dc.PLANTS_PER_MU
TRIALS_LA = ["XX_SF", "SD_SF"]        # both have LAI observations
AX_GRID = [850.0, 950.0, 1050.0, 1150.0, 1250.0]
ALL_FIXED = 900.0


def _diag(trial: str) -> Path:
    d = dc.BASE / "out" / trial / "diag"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _fig(trial: str) -> Path:
    d = dc.BASE / "fig" / trial / "diag"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _dens(plot: str) -> float:
    pl = dc.load_plots().loc[plot]
    v = pl["density_measured"]
    if pd.isna(v):
        v = pl["density_planted"]
    return float(v) / MU


def obs_per_plant_la(trial: str) -> pd.Series:
    o = pd.read_csv(dc.DATASET / "obs_plant.csv", parse_dates=["date"])
    o = o[o.trial_id == trial]
    out = {}
    for plot in dc.load_plots().index:
        lai = o.loc[o["plot"] == plot, "lai_plot"].dropna()
        out[plot] = float(lai.max()) / _dens(plot) if len(lai) else np.nan
    return pd.Series(out)


# --------------------------------------------------------------------------
def leaf_area() -> None:
    for trial in TRIALS_LA:
        dc.set_trial(trial)
        obs = obs_per_plant_la(trial)

        dc.set_model("CERES")
        ceres = {}
        for plot in dc.ORDER:
            g, s = dc.run_plot(plot)
            m = dc.plot_metrics(g, s, _dens(plot))
            ceres[plot] = m["la_per_plant"]

        grid, best = [], None
        for ax in AX_GRID:
            dc.set_model("IXIM")
            vals = {}
            for plot in dc.ORDER:
                cv = dict(dc.DEFAULT_CULTIVAR, AX=ax, ALL=ALL_FIXED)
                g, s = dc.run_plot(plot, cv)
                m = dc.plot_metrics(g, s, _dens(plot))
                vals[plot] = m["la_per_plant"]
            rmse = float(np.sqrt(np.nanmean(
                [(vals[p] - obs[p]) ** 2 for p in dc.ORDER])))
            bias = float(np.nanmean([vals[p] / obs[p] - 1 for p in dc.ORDER])) * 100
            grid.append(dict(trial=trial, AX=ax, ALL=ALL_FIXED,
                             la_rmse=rmse, la_bias_pct=bias,
                             **{f"la_{p}": vals[p] for p in dc.ORDER}))
            if best is None or rmse < best[1]:
                best = (ax, rmse, vals)
        pd.DataFrame(grid).to_csv(_diag(trial) / "ixim_ax_grid.csv", index=False)
        ax, rmse, ixim = best
        print(f"  {trial}: best AX={ax:.0f}  LA RMSE={rmse:.3f} m2/plant")

        rows = [dict(plot=p, density=_dens(p), obs=obs[p], ceres=ceres[p],
                     ixim=ixim[p], ax=ax) for p in dc.ORDER]
        d = pd.DataFrame(rows)
        d.to_csv(_diag(trial) / "per_plant_la.csv", index=False)
        print(d.round(3).to_string(index=False))
        _fig_leafarea(trial, d, ax)
    dc.set_model("CERES")


def _fig_leafarea(trial: str, d: pd.DataFrame, ax_val: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": 0.25, "grid.linestyle": "--",
                         "axes.labelsize": 10, "axes.titlesize": 11,
                         "legend.fontsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4),
                             gridspec_kw={"width_ratios": [1.55, 1]})
    ax = axes[0]
    x = np.arange(len(d))
    ax.bar(x - 0.27, d.obs, 0.27, label="实测", color="#7f8c8d",
           edgecolor="k", lw=0.4)
    ax.bar(x, d.ceres, 0.27, label="CERES-Maize", color="#c0392b",
           edgecolor="k", lw=0.4)
    ax.bar(x + 0.27, d.ixim, 0.27, label=f"IXIM（AX={ax_val:.0f}）",
           color="#2980b9", edgecolor="k", lw=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.plot}\n{dc.N_RATE[r.plot]}N" for r in d.itertuples()],
                       fontsize=8)
    ax.set_ylabel("峰值单株叶面积 (m2/株)")
    ax.set_title("A  逐处理：实测 vs CERES vs IXIM")
    ax.legend(frameon=False, ncol=3, loc="upper left")

    ax = axes[1]
    lim = [0, max(d.obs.max(), d.ceres.max(), d.ixim.max()) * 1.12]
    ax.plot(lim, lim, "k--", lw=1.1)
    ax.scatter(d.obs, d.ceres, s=70, c="#c0392b", edgecolor="k", lw=0.4,
               label=f"CERES (RMSE {np.sqrt(np.mean((d.ceres-d.obs)**2)):.3f})")
    ax.scatter(d.obs, d.ixim, s=70, c="#2980b9", edgecolor="k", lw=0.4,
               label=f"IXIM (RMSE {np.sqrt(np.mean((d.ixim-d.obs)**2)):.3f})")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("实测 单株叶面积 (m2/株)")
    ax.set_ylabel("模拟 单株叶面积 (m2/株)")
    ax.set_title("B  1:1")
    ax.legend(frameon=False, loc="lower right")

    fig.suptitle(f"F｜{dc.DISPLAY} {dc.CFG.variety}｜峰值单株叶面积"
                 f"（IXIM 品种参数 AX 可调，CERES 无此参数）", fontsize=12)
    fig.tight_layout()
    p = _fig(trial) / "F_单株叶面积.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  figure ->", p)


# --------------------------------------------------------------------------
def _yield_obs(trial: str) -> pd.Series:
    y = pd.read_csv(dc.DATASET / "yield.csv")
    return y[y.trial_id == trial].set_index("plot")["yield_dry_kg_ha"]


def shared_cultivar() -> None:
    """Joint (P5, G3) fit for MY73 at both sites: yield + biomass first."""
    grid = []
    for p5 in (940.0, 960.0, 980.0):
        for g3 in (5.6, 5.9, 6.2):
            rec = dict(P5=p5, G3=g3)
            ok = True
            for trial in ("XX_SF", "SD_MY73"):
                dc.set_trial(trial)
                dc.set_model("CERES")
                cv = dict(dc.DEFAULT_CULTIVAR, P5=p5, G3=g3)
                obs_y = _yield_obs(trial)
                obs_b = pd.read_csv(dc.DATASET / "obs_plant.csv",
                                    parse_dates=["date"])
                obs_b = obs_b[obs_b.trial_id == trial]
                sy, sb, ob, oy, ant, mat = [], [], [], [], [], []
                for p in dc.ORDER:
                    g, s = dc.run_plot(p, cv)
                    m = dc.plot_metrics(g, s, _dens(p))
                    sy.append(m["yield_kg"])
                    sb.append(m["bio_final"])
                    oy.append(float(obs_y[p]))
                    ob.append(float((obs_b.loc[obs_b["plot"] == p, "total_dw"]
                                     * _dens(p) * 10).max()))
                    ant.append(m["anthesis"])
                    mat.append(m["maturity"])
                sy, sb, oy, ob = map(np.array, (sy, sb, oy, ob))
                rec[f"{trial}_y_rrmse"] = float(
                    np.sqrt(np.mean((sy - oy) ** 2)) / oy.mean() * 100)
                rec[f"{trial}_b_bias"] = float(sb.sum() / ob.sum() * 100 - 100)
                rec[f"{trial}_ant"] = str(pd.Series(ant).median().date())
                rec[f"{trial}_mat"] = str(pd.Series(mat).median().date())
            obs_ant = {"XX_SF": "2026-08-08", "SD_MY73": "2026-08-08"}
            obs_mat = {"XX_SF": "2026-10-05", "SD_MY73": "2026-10-08"}
            for trial in ("XX_SF", "SD_MY73"):
                da = (pd.Timestamp(rec[f"{trial}_ant"])
                      - pd.Timestamp(obs_ant[trial])).days
                dm = (pd.Timestamp(rec[f"{trial}_mat"])
                      - pd.Timestamp(obs_mat[trial])).days
                rec[f"{trial}_d_ant"] = da
                rec[f"{trial}_d_mat"] = dm
                if abs(da) > 4 or abs(dm) > 4:
                    ok = False
            rec["pheno_ok"] = ok
            # the user's priority: yield and biomass first, phenology <= 4 d
            rec["score"] = (rec["XX_SF_y_rrmse"] + rec["SD_MY73_y_rrmse"]
                            + abs(rec["XX_SF_b_bias"]) + abs(rec["SD_MY73_b_bias"]))
            grid.append(rec)
    d = pd.DataFrame(grid).sort_values(["pheno_ok", "score"],
                                       ascending=[False, True])
    d.to_csv(_diag("XX_SF") / "shared_cultivar.csv", index=False)
    print(d.round(2).to_string(index=False))
    best = d.iloc[0]
    _fig_shared(best)


def _fig_shared(best: pd.Series) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": 0.25, "grid.linestyle": "--"})
    p5, g3 = float(best["P5"]), float(best["G3"])
    cv = dict(dc.DEFAULT_CULTIVAR, P5=p5, G3=g3)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    for ax, trial in zip(axes, ("XX_SF", "SD_MY73")):
        dc.set_trial(trial)
        dc.set_model("CERES")
        obs_y = _yield_obs(trial)
        obs_b = pd.read_csv(dc.DATASET / "obs_plant.csv", parse_dates=["date"])
        obs_b = obs_b[obs_b.trial_id == trial]
        sy, sb, oy, ob, names = [], [], [], [], []
        for p in dc.ORDER:
            g, s = dc.run_plot(p, cv)
            m = dc.plot_metrics(g, s, _dens(p))
            sy.append(m["yield_kg"])
            sb.append(m["bio_final"])
            oy.append(float(obs_y[p]))
            ob.append(float((obs_b.loc[obs_b["plot"] == p, "total_dw"]
                             * _dens(p) * 10).max()))
            names.append(p)
        x = np.arange(len(names))
        ax.bar(x - 0.3, oy, 0.2, label="实测产量", color="#7f8c8d", edgecolor="k", lw=0.4)
        ax.bar(x - 0.1, sy, 0.2, label="模拟产量", color="#c0392b", edgecolor="k", lw=0.4)
        ax.bar(x + 0.1, ob, 0.2, label="实测生物量", color="#95a5a6", edgecolor="k", lw=0.4)
        ax.bar(x + 0.3, sb, 0.2, label="模拟生物量", color="#e67e22", edgecolor="k", lw=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels(names, fontsize=8, rotation=45)
        ax.set_ylabel("kg/ha (干重)")
        ax.set_title(f"{dc.DISPLAY} {dc.CFG.variety}　产量 RRMSE "
                     f"{best[f'{trial}_y_rrmse']:.1f}%　生物量偏差 "
                     f"{best[f'{trial}_b_bias']:+.1f}%")
        ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.suptitle(f"G｜MY73 两地共用一套品种参数（P5 {p5:.0f} / G3 {g3:.2f}，"
                 f"物候误差 ≤4 d）", fontsize=12)
    fig.tight_layout()
    p = _fig("XX_SF") / "G_共用品种参数.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  figure ->", p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="all",
                    choices=["all", "leafarea", "shared"])
    args = ap.parse_args()
    if args.section in ("all", "leafarea"):
        print("\n[1] per-plant leaf area: CERES vs IXIM(AX)")
        leaf_area()
    if args.section in ("all", "shared"):
        print("\n[2] shared MY73 cultivar: yield + biomass first")
        shared_cultivar()


if __name__ == "__main__":
    main()
