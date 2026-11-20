"""Can the leaf-appearance timing be fixed with PHINT, and does it cost?

PHINT is the phyllochron (GDD between successive leaf tips).  Raising it delays
leaf appearance but also reduces the number of leaves that fit into P1, so P1
has to grow with it.  This script scans (P1, PHINT) and reports, for each pair,
the leaf-stage timing error and the effect on LAI / yield / biomass / phenology.

Outputs
-------
out/<trial>/diag/phint_scan.csv
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

COMBOS = {
    "XX_SF": [(235, 50), (240, 51), (245, 52), (250, 53), (255, 54)],
    "SD_SF": [(250, 51), (260, 53), (265, 55), (235, 60)],
}

# Raw field book (展开叶 / 可见叶), used as the reference for `L#SD`.
# 新乡 MY73: 原始记录/2026新乡玉米试验数据汇总.xlsx -> "MY73试验田"
# 德州 乐农87: 原始记录/2026年德州玉米试验数据汇总.xlsx -> "试验田_乐农87"
RAW_VISIBLE = {
    "XX_SF": [("2026-07-03", 3, 5), ("2026-07-08", 5, 7), ("2026-07-11", 6, 8),
              ("2026-07-16", 7, 10), ("2026-07-18", 8, 11), ("2026-07-22", 10, 12),
              ("2026-07-24", 11, 13), ("2026-07-26", 12, 14), ("2026-07-29", 14, 15),
              ("2026-08-01", 15, 17), ("2026-08-05", 19, 20), ("2026-08-08", 20, 21)],
    "SD_SF": [("2026-07-07", 4, 7), ("2026-07-10", 5, 9), ("2026-07-13", 6, 10),
              ("2026-07-16", 7, 11), ("2026-07-19", 8, 12), ("2026-07-22", 9, 13),
              ("2026-07-25", 10, 15), ("2026-07-28", 12, 17), ("2026-07-31", 14, 18),
              ("2026-08-03", 15, 19), ("2026-08-04", 16, 19), ("2026-08-05", 17, 20)],
}


def leaf_obs() -> pd.Series:
    ph = pd.read_csv(dc.DATASET / "phenology.csv", parse_dates=["date"])
    site = str(dc.CFG.site)
    sub = ph[ph["site"] == site]
    if sub.empty:
        sub = ph[ph["site"].astype(str).str.startswith(site[:2])]
    if sub.empty:
        sub = ph[ph["site"].astype(str).str.endswith(site[-2:])]
    sub = sub[sub["variety"] == dc.CFG.variety]
    rows = []
    for r in sub.itertuples():
        m = re.match(r"\s*(\d+)\s*叶", str(r.stage))
        if m:
            rows.append((int(m.group(1)), pd.Timestamp(r.date)))
    d = pd.DataFrame(rows, columns=["leaf", "date"]).drop_duplicates("leaf")
    return d.set_index("leaf")["date"].sort_index()


def _dens(plot: str) -> float:
    pl = dc.load_plots().loc[plot]
    v = pl["density_measured"]
    if pd.isna(v):
        v = pl["density_planted"]
    return float(v) / dc.PLANTS_PER_MU


def _med_date(values) -> str:
    s = pd.Series(values).dropna()
    if s.empty:
        return "NaT"
    med = s.median()
    return str(med.date()) if hasattr(med, "date") else "NaT"


def run_combo(trial: str, p1: float, phint: float) -> dict:
    dc.set_trial(trial)
    dc.set_model("CERES")
    cv = dict(dc.DEFAULT_CULTIVAR, P1=p1, PHINT=phint)
    y = pd.read_csv(dc.DATASET / "yield.csv")
    y = y[y.trial_id == trial].set_index("plot")["yield_dry_kg_ha"]
    o = pd.read_csv(dc.DATASET / "obs_plant.csv", parse_dates=["date"])
    o = o[o.trial_id == trial]

    sy, sb, oy, ob, lai_p, ant, mat = [], [], [], [], [], [], []
    reach = {}
    for p in dc.ORDER:
        dens = _dens(p)
        g, s = dc.run_plot(p, cv)
        m = dc.plot_metrics(g, s, dens)
        sy.append(m["yield_kg"])
        sb.append(m["bio_final"])
        oy.append(float(y[p]))
        ob.append(float((o.loc[o["plot"] == p, "total_dw"] * dens * 10).max()))
        gg = g.set_index(g["date"].dt.normalize())
        for r in o[(o["plot"] == p) & o.lai_plot.notna()].itertuples():
            k = pd.Timestamp(r.date).normalize()
            if k in gg.index:
                lai_p.append((gg.loc[k, "LAID"], r.lai_plot))
        ant.append(m["anthesis"])
        mat.append(m["maturity"])
        for lf in range(2, 24):
            hit = g[g["L#SD"] >= lf]
            if len(hit):
                reach.setdefault(lf, []).append(hit["date"].iloc[0])

    errs, errs_exp = [], []
    for dt_s, exp, vis in RAW_VISIBLE[trial]:
        dt = pd.Timestamp(dt_s)
        for target, bucket in ((vis, errs), (exp, errs_exp)):
            if target in reach:
                med = pd.Series(reach[target]).median()
                if pd.notna(med):
                    bucket.append((med - dt).days)
    sy, sb, oy, ob = map(np.array, (sy, sb, oy, ob))
    lai_err = np.array([a - b for a, b in lai_p])
    return dict(
        trial=trial, P1=p1, PHINT=phint,
        leaf_err=np.mean(errs) if errs else np.nan,
        leaf_err_expanded=np.mean(errs_exp) if errs_exp else np.nan,
        leaf_n=int(max(reach)) if reach else 0,
        y_rrmse=float(np.sqrt(((sy - oy) ** 2).mean()) / oy.mean() * 100),
        y_bias=float(sy.sum() / oy.sum() * 100 - 100),
        b_bias=float(sb.sum() / ob.sum() * 100 - 100),
        lai_rmse=float(np.sqrt(np.mean(lai_err ** 2))),
        anthesis=_med_date(ant), maturity=_med_date(mat),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="all", choices=["all", "XX_SF", "SD_SF"])
    args = ap.parse_args()
    trials = ["XX_SF", "SD_SF"] if args.trial == "all" else [args.trial]
    rows = []
    for trial in trials:
        dc.set_trial(trial)
        print(f"\n== {trial}  {dc.CFG.variety}  "
              f"(当前 P1 {dc.DEFAULT_CULTIVAR['P1']:.0f} / PHINT "
              f"{dc.DEFAULT_CULTIVAR['PHINT']:.0f}) ==")
        for p1, ph in COMBOS[trial]:
            r = run_combo(trial, float(p1), float(ph))
            rows.append(r)
            print(f"  P1 {p1:3d} PHINT {ph:2d} | 叶位误差 {r['leaf_err']:+5.1f} d | "
                  f"最终叶数 {r['leaf_n']:2d} | 抽雄 {r['anthesis']} 成熟 {r['maturity']} | "
                  f"LAI RMSE {r['lai_rmse']:.3f} | 产量 RRMSE {r['y_rrmse']:5.2f}% | "
                  f"生物量 {r['b_bias']:+5.1f}%")
    out = dc.BASE / "out" / "XX_SF" / "diag" / "phint_scan.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("\n->", out)


if __name__ == "__main__":
    main()
