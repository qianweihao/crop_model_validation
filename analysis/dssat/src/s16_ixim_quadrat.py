"""SD_SF：IXIM vs CERES 对"10 米双行"产量的表现。

两套观测：
  obs_plot  小区实测产量（与实际留苗对应）
  obs_quad  10 米双行实测 = 单穗重 x 5500 株/亩（与留苗无关，密度归一化）

两种模拟密度：
  base     实测留苗 3757-4675
  nominal  播种密度 5640

问：换上 IXIM 能否减小与 obs_quad 的误差？

Outputs
-------
out/SD_SF/diag/ixim_quadrat.csv
out/SD_SF/diag/ixim_quadrat_grid.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402
from s15_quadrat import RAW, rrmse  # noqa: E402

QUAD_DENS = 5500.0
NOM_DENS = 5640.0


def _run(model: str, scen: str, overrides: dict | None = None) -> pd.DataFrame:
    """Per-plot yield (kg/亩, 14%), final biomass (kg/ha dry) and LAI RMSE."""
    dc.set_model(model)
    dc.set_scenario(scen)
    cv = dict(dc.DEFAULT_CULTIVAR)
    cv.update(overrides or {})
    obs = dc.load_obs()
    rows = {}
    for p in dc.ORDER:
        g, s = dc.run_plot(p, cv)
        g = g.set_index("date")
        o = obs[obs["plot"] == p]
        d = [(g.loc[pd.Timestamp(r.date).normalize(), "LAID"], r.lai_plot)
             for r in o.itertuples() if pd.Timestamp(r.date).normalize() in g.index]
        lai_rmse = (float(np.sqrt(np.mean([(a - b) ** 2 for a, b in d])))
                    if d else np.nan)
        last = g.sort_index().iloc[-1]
        rows[p] = dict(yield_ha=s["HWAM"], bio_ha=float(last["CWAD"]),
                       lai_rmse=lai_rmse)
    dc.set_scenario("base")
    dc.set_model("CERES")
    return pd.DataFrame(rows).T


def _table(meas: pd.DataFrame, nom: pd.DataFrame, model: str) -> dict:
    plot = np.array([RAW[p][1] for p in meas.index])
    quad = np.array([RAW[p][2] for p in meas.index])
    dens = np.array([RAW[p][0] for p in meas.index])
    sm = meas["yield_ha"].values / 15 / 0.86
    sn = nom["yield_ha"].values / 15 / 0.86
    bm = meas["bio_ha"].values
    return dict(
        model=model,
        rrmse_plot_meas=rrmse(sm, plot), bias_plot_meas=sm.sum() / plot.sum() * 100 - 100,
        rrmse_plot_nom=rrmse(sn, plot),  bias_plot_nom=sn.sum() / plot.sum() * 100 - 100,
        rrmse_quad_meas=rrmse(sm, quad), bias_quad_meas=sm.sum() / quad.sum() * 100 - 100,
        rrmse_quad_nom=rrmse(sn, quad),  bias_quad_nom=sn.sum() / quad.sum() * 100 - 100,
        w_meas_mean=(sm / dens * 1000).mean(),
        w_nom_mean=(sn / NOM_DENS * 1000).mean(),
        lai_rmse=float(meas["lai_rmse"].mean()),
        bio_bias=bm.mean() / _obs_bio().mean() * 100 - 100,
    )


def _obs_bio() -> pd.Series:
    o = dc.load_obs()
    return o.groupby("plot")["bio_kg_ha"].max().reindex(dc.ORDER)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="all", choices=["all", "naive", "grid"])
    args = ap.parse_args()
    diag = dc.BASE / "out" / "SD_SF" / "diag"

    if args.section in ("all", "naive"):
        rows = []
        for model, ov in [("CERES", {}),
                          ("IXIM", dict(AX=650.0, ALL=800.0))]:
            rows.append(_table(_run(model, "base", ov), _run(model, "nominal", ov), model))
        d = pd.DataFrame(rows)
        d.to_csv(diag / "ixim_quadrat.csv", index=False)
        print("== 直接套各自已标定参数 ==")
        print(d[["model", "rrmse_plot_meas", "bias_plot_meas", "rrmse_quad_meas",
                 "bias_quad_meas", "rrmse_quad_nom", "bias_quad_nom", "w_nom_mean",
                 "lai_rmse", "bio_bias"]].round(2).to_string(index=False))
        print("\n  单穗重 g/穗 : 田间(小区实测口径) "
              f"{np.mean([RAW[p][1] * 1000 / RAW[p][0] for p in RAW]):.1f} ; "
              f"田间(10米双行口径) {np.mean([RAW[p][2] * 1000 / QUAD_DENS for p in RAW]):.1f}")

    if args.section in ("all", "grid"):
        rows = []
        for ax in (650.0, 750.0, 850.0):          # MZIXM048.CUL: AX 200-850
            for g2 in (900.0, 1100.0):            # G2 上限 1100 kernels/ear
                for g3 in (6.2, 7.0, 8.0, 9.0):   # G3 4-20 mg/核/天
                    ov = dict(AX=ax, ALL=800.0, G2=g2, G3=g3)
                    r = _table(_run("IXIM", "base", ov), _run("IXIM", "nominal", ov), "IXIM")
                    r.update(AX=ax, G2=g2, G3=g3)
                    rows.append(r)
        d = pd.DataFrame(rows)
        d["score"] = d.rrmse_plot_meas + d.bias_plot_meas.abs()   # 只用小区实测标定
        d = d.sort_values("score")
        d.to_csv(diag / "ixim_quadrat_grid.csv", index=False)
        print("\n== IXIM 网格（按小区实测排序，前 8）==")
        cols = ["AX", "G2", "G3", "rrmse_plot_meas", "bias_plot_meas",
                "rrmse_quad_nom", "bias_quad_nom", "w_nom_mean", "lai_rmse",
                "bio_bias"]
        print(d.head(8)[cols].round(2).to_string(index=False))
        best = d.iloc[0]
        print(f"\n最优 IXIM: AX {best.AX:.0f} / G2 {best.G2:.0f} / G3 {best.G3:.1f}"
              f"  ->  小区实测 RRMSE {best.rrmse_plot_meas:.2f}%"
              f" ; 10米双行 RRMSE {best.rrmse_quad_nom:.2f}%"
              f" (CERES 18.0%)")
    dc.set_model("CERES")


if __name__ == "__main__":
    main()
