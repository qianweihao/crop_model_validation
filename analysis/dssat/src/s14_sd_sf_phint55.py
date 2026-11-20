"""SD_SF (乐农87) with PHINT 55: can yield and biomass errors be balanced?

PHINT 55 + P1 265 fixes the leaf number (20, matching the field book) and pulls
biomass from -10.9% to -3.7%, but pushes yield RRMSE from 3.70% to 10.11%.
This scans P1 x G3 on top of PHINT 55 and ranks the candidates by the mean of
|yield bias| and |biomass bias|, keeping an eye on LAI and phenology.

Outputs
-------
out/SD_SF/diag/phint55_grid.csv
fig/SD_SF/diag/I_PHINT55.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

PHINT = 55.0
P1S = [250.0, 265.0, 280.0]
G3S = [5.2, 5.6, 6.0]
BASE = dict(P1=235.0, PHINT=60.0)


def _dens(plot: str) -> float:
    pl = dc.load_plots().loc[plot]
    v = pl["density_measured"]
    if pd.isna(v):
        v = pl["density_planted"]
    return float(v) / dc.PLANTS_PER_MU


def evaluate(cv: dict) -> dict:
    y = pd.read_csv(dc.DATASET / "yield.csv")
    y = y[y.trial_id == dc.TRIAL].set_index("plot")["yield_dry_kg_ha"]
    o = pd.read_csv(dc.DATASET / "obs_plant.csv", parse_dates=["date"])
    o = o[o.trial_id == dc.TRIAL]
    sy, sb, oy, ob, lai, ant, mat = [], [], [], [], [], [], []
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
                lai.append(gg.loc[k, "LAID"] - r.lai_plot)
        ant.append(m["anthesis"])
        mat.append(m["maturity"])
    sy, sb, oy, ob = map(np.array, (sy, sb, oy, ob))
    y_bias = float(sy.sum() / oy.sum() * 100 - 100)
    b_bias = float(sb.sum() / ob.sum() * 100 - 100)
    return dict(
        P1=cv["P1"], PHINT=cv["PHINT"], G3=cv["G3"],
        y_rrmse=float(np.sqrt(((sy - oy) ** 2).mean()) / oy.mean() * 100),
        y_bias=y_bias, b_bias=b_bias,
        balance=(abs(y_bias) + abs(b_bias)) / 2,
        lai_rmse=float(np.sqrt(np.mean(np.array(lai) ** 2))),
        anthesis=str(pd.Series(ant).median().date()),
        maturity=str(pd.Series(mat).median().date()),
        leaves=int(np.median([g["L#SD"].max() for _, g in
                              [(p, None) for p in []]]) if False else 0),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    args = ap.parse_args()
    dc.set_trial("SD_SF")
    dc.set_model("CERES")

    rows = []
    base = dict(dc.DEFAULT_CULTIVAR, **BASE)
    r = evaluate(base)
    r["tag"] = "当前 P1 235 / PHINT 60 / G3 6.0"
    rows.append(r)

    for p1 in P1S:
        for g3 in G3S:
            cv = dict(dc.DEFAULT_CULTIVAR, P1=p1, PHINT=PHINT, G3=g3)
            r = evaluate(cv)
            r["tag"] = f"PHINT 55 / P1 {p1:.0f} / G3 {g3}"
            rows.append(r)

    d = pd.DataFrame(rows).sort_values("balance")
    cols = ["tag", "y_rrmse", "y_bias", "b_bias", "balance", "lai_rmse",
            "anthesis", "maturity"]
    print(d[cols].round(2).to_string(index=False))
    out = dc.BASE / "out" / "SD_SF" / "diag" / "phint55_grid.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    d[cols].to_csv(out, index=False)
    print("\n->", out)
    _figure(d)


def _figure(d: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.sans-serif": ["SimHei", "Microsoft YaHei"],
                         "axes.unicode_minus": False, "axes.grid": True,
                         "grid.alpha": 0.25, "grid.linestyle": "--"})
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    cur = d[d.tag.str.startswith("当前")]
    cand = d[~d.tag.str.startswith("当前")]
    ax.scatter(cand.y_bias, cand.b_bias, s=90, c="#2980b9", edgecolor="k",
               lw=0.5, zorder=3, label="PHINT 55 候选")
    ax.scatter(cur.y_bias, cur.b_bias, s=150, marker="*", c="#c0392b",
               edgecolor="k", lw=0.6, zorder=4, label="当前 (PHINT 60)")
    best = d.iloc[0]
    ax.scatter([best.y_bias], [best.b_bias], s=190, marker="o", facecolor="none",
               edgecolor="#e67e22", lw=2.4, zorder=5,
               label=f"最优：{best.tag}")
    for r in cand.itertuples():
        ax.annotate(f"{r.tag.split('/')[1].strip()}/{r.tag.split('/')[2].strip()}",
                    (r.y_bias, r.b_bias), xytext=(6, 4),
                    textcoords="offset points", fontsize=7.5)
    ax.axhline(0, color="grey", lw=0.8)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("产量偏差 (%)")
    ax.set_ylabel("生物量偏差 (%)")
    ax.set_title("I｜德州 乐农87｜产量 vs 生物量 偏差（越靠原点越好）")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    p = dc.BASE / "fig" / "SD_SF" / "diag" / "I_PHINT55.png"
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("figure ->", p)


if __name__ == "__main__":
    main()
