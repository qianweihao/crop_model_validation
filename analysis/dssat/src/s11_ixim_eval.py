"""How well does IXIM-Maize do on yield and biomass, compared with CERES?

Two evaluations:

  naive  : IXIM run with the CERES-calibrated cultivar record (AX 950 / ALL 900)
  fitted : IXIM with its own G2 / G3, searched on XX_SF and then applied to
           the other two trials

Outputs
-------
out/<trial>/diag/ixim_eval.csv
out/<trial>/diag/ixim_grid.csv
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
AX = 950.0          # from the leaf-area calibration (s10)
ALL = 900.0


def _diag(trial: str) -> Path:
    d = dc.BASE / "out" / trial / "diag"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _dens(plot: str) -> float:
    pl = dc.load_plots().loc[plot]
    v = pl["density_measured"]
    if pd.isna(v):
        v = pl["density_planted"]
    return float(v) / MU


def _obs(trial: str):
    y = pd.read_csv(dc.DATASET / "yield.csv")
    y = y[y.trial_id == trial].set_index("plot")["yield_dry_kg_ha"]
    o = pd.read_csv(dc.DATASET / "obs_plant.csv", parse_dates=["date"])
    o = o[o.trial_id == trial]
    return y, o


def evaluate(trial: str, model: str, overrides: dict | None = None) -> dict:
    """overrides are applied *after* set_trial, so the trial's own calibrated
    cultivar record is the starting point."""
    dc.set_trial(trial)
    dc.set_model(model)
    cultivar = dict(dc.DEFAULT_CULTIVAR)
    cultivar.update(overrides or {})
    y, o = _obs(trial)
    sy, sb, oy, ob, ant, mat = [], [], [], [], [], []
    for p in dc.ORDER:
        g, s = dc.run_plot(p, cultivar)
        m = dc.plot_metrics(g, s, _dens(p))
        sy.append(m["yield_kg"])
        sb.append(m["bio_final"])
        oy.append(float(y[p]))
        ob.append(float((o.loc[o["plot"] == p, "total_dw"] * _dens(p) * 10).max()))
        ant.append(m["anthesis"])
        mat.append(m["maturity"])
    sy, sb, oy, ob = map(np.array, (sy, sb, oy, ob))
    return dict(
        trial=trial, model=model, n=len(sy),
        y_rrmse=float(np.sqrt(np.mean((sy - oy) ** 2)) / oy.mean() * 100),
        y_bias=float(sy.sum() / oy.sum() * 100 - 100),
        b_rrmse=float(np.sqrt(np.mean((sb - ob) ** 2)) / ob.mean() * 100),
        b_bias=float(sb.sum() / ob.sum() * 100 - 100),
        y_sim=float(sy.mean()), y_obs=float(oy.mean()),
        b_sim=float(sb.mean()), b_obs=float(ob.mean()),
        ant=str(pd.Series(ant).median().date()),
        mat=str(pd.Series(mat).median().date()),
    )


def grid(trial: str) -> pd.DataFrame:
    rows = []
    for g2 in (900.0, 1100.0, 1300.0):
        for g3 in (6.0, 8.0, 10.0):
            cv = dict(dc.DEFAULT_CULTIVAR, AX=AX, ALL=ALL, G2=g2, G3=g3)
            r = evaluate(trial, "IXIM", dict(AX=AX, ALL=ALL, G2=g2, G3=g3))
            r.update(G2=g2, G3=g3)
            rows.append(r)
    d = pd.DataFrame(rows)
    d["score"] = d["y_rrmse"] + d["b_bias"].abs()
    d = d.sort_values("score")
    d.to_csv(_diag(trial) / "ixim_grid.csv", index=False)
    print(d[["G2", "G3", "y_rrmse", "y_bias", "b_bias", "y_sim", "b_sim",
             "score"]].round(2).to_string(index=False))
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="all",
                    choices=["all", "naive", "grid"])
    args = ap.parse_args()

    if args.section in ("all", "naive"):
        print("\n== 直接套 CERES 品种参数（AX 950 / ALL 900）==")
        rows = []
        for trial in ("XX_SF", "SD_SF", "SD_MY73"):
            rows.append(evaluate(trial, "CERES"))
            rows.append(evaluate(trial, "IXIM", dict(AX=AX, ALL=ALL)))
        d = pd.DataFrame(rows)
        d.to_csv(_diag("XX_SF") / "ixim_eval.csv", index=False)
        print(d[["trial", "model", "y_rrmse", "y_bias", "b_rrmse", "b_bias",
                 "y_sim", "y_obs", "b_sim", "b_obs"]].round(1).to_string(index=False))

    if args.section in ("all", "grid"):
        print("\n== IXIM 自己的 G2 x G3 网格（XX_SF, AX 950 / ALL 900）==")
        g = grid("XX_SF")
        best = g.iloc[0]
        cv = dict(AX=AX, ALL=ALL, G2=float(best.G2), G3=float(best.G3))
        print(f"\n== 用最优 IXIM 参数 (G2 {best.G2:.0f} / G3 {best.G3:.0f}) "
              f"评估三个试验 ==")
        rows = [evaluate(t, "IXIM", cv) for t in ("XX_SF", "SD_SF", "SD_MY73")]
        d = pd.DataFrame(rows)
        d.to_csv(_diag("XX_SF") / "ixim_eval_fitted.csv", index=False)
        print(d[["trial", "model", "y_rrmse", "y_bias", "b_rrmse", "b_bias",
                 "y_sim", "y_obs", "b_sim", "b_obs"]].round(1).to_string(index=False))
    dc.set_model("CERES")


if __name__ == "__main__":
    main()
