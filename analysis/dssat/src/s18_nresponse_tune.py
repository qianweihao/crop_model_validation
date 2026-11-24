"""用田间实际氮饱和点约束模型，主动牺牲整体精度换取 N 响应。

目标饱和点（由小区实测产量线性-平台拟合得到）：
    德州 SD_SF  177 kg N/ha
    新乡 XX_SF  150 kg N/ha

可调参数：
    CTCNP1       临界氮浓度上限（SPE，控制作物需氮量）
    nitrate_scale 初始土壤 NO3 折算系数（控制土壤供氮）

Outputs: out/<trial>/diag/nresponse_tune.csv
"""
from __future__ import annotations
import argparse, sys
from functools import lru_cache
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc
from DSSATTools import DSSAT
from DSSATTools.base.sections import TabularSubsection

TARGET = {"SD_SF": (177.0, 4500), "XX_SF": (150.0, 5200)}
SWEEP = (60, 100, 140, 180, 220, 260)


def _run(trial, plot, dens_mu, n_target=None, nitrate_scale=None):
    dc.set_trial(trial)
    dens = dens_mu / dc.PLANTS_PER_MU
    soil = dc.build_soil(plot)
    wth = dc.build_weather()
    crop = dc.build_crop(None)
    man = dc.build_management(plot, soil, nitrogen=True,
                              nitrate_scale=nitrate_scale, density=dens)
    if n_target is not None:
        t = man.fertilizers["table"].copy()
        cur = float(t["FAMN"].sum())
        t["FAMN"] = (t["FAMN"] * n_target / cur).round(2)
        man.fertilizers["table"] = TabularSubsection(t)
    d = DSSAT(); d.run(soil, wth, crop, man, verbose=False)
    g = d.output["PlantGro"].sort_index(); last = g.iloc[-1]
    out = (float(last["GWAD"]) / 15, float(last["CWAD"]) / 15)
    d.close()
    return out


def _plateau_n(N, y):
    N, y = np.asarray(N, float), np.asarray(y, float)
    best = None
    for n0 in np.arange(40, 320, 2.0):
        X = np.column_stack([np.ones_like(N), np.minimum(N, n0)])
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        sse = ((y - X @ b) ** 2).sum()
        if best is None or sse < best[0]:
            best = (sse, n0)
    return best[1]


def evaluate(trial, c1, ns):
    tgt, dens = TARGET[trial]
    dc.set_trial(trial)
    yy = pd.read_csv(dc.DATASET / "yield.csv")
    yy = yy[yy.trial_id == trial].set_index("plot")
    pl = dc.load_plots(); ob = dc.load_obs()
    sim, obs, bs, bo = [], [], [], []
    for p in dc.ORDER:
        y, b = _run(trial, p, float(pl.loc[p, "density_measured"]),
                    nitrate_scale=ns)
        sim.append(y * 15); obs.append(float(yy.loc[p, "yield_dry_kg_ha"]))
        bs.append(b * 15); bo.append(float(ob[ob["plot"] == p]["bio_kg_ha"].max()))
    sim, obs, bs, bo = map(np.array, (sim, obs, bs, bo))
    sweep_y = []
    for n in SWEEP:
        y, _ = _run(trial, "A5B2" if trial == "SD_SF" else "A4B2",
                    dens, n_target=n, nitrate_scale=ns)
        sweep_y.append(y)
    n0 = _plateau_n(SWEEP, sweep_y)
    return dict(trial=trial, CTCNP1=c1, nitrate=ns,
                yield_rrmse=round(float(np.sqrt(((sim - obs) ** 2).mean())
                                       / obs.mean() * 100), 2),
                yield_bias=round(float(sim.sum() / obs.sum() * 100 - 100), 2),
                bio_bias=round(float(bs.sum() / bo.sum() * 100 - 100), 2),
                N0_model=n0, N0_target=tgt,
                N_gap=round(n0 - tgt, 0),
                resp_139_253=round(float(sweep_y[SWEEP.index(140)]) and
                                   (sweep_y[-1] / sweep_y[SWEEP.index(140)] - 1) * 100, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="both")
    ap.add_argument("--grid", default="a")
    args = ap.parse_args()
    trials = ["SD_SF", "XX_SF"] if args.trial == "both" else [args.trial]
    rows = []
    c1s = ((1.52, 1.70, 1.85, 2.00) if args.grid == "a"
           else (1.52, 1.60, 1.70, 1.80))
    nss = ((0.24, 0.40, 0.55) if args.grid == "a" else (0.45, 0.55, 0.65))
    for c1 in c1s:
        for ns in nss:
            for tr in trials:
                sp = dict(dc.TRIALS[tr].species or {})   # 试验自带物种覆盖
                with dc.patch_species(CTCNP1=c1, **sp):
                    r = evaluate(tr, c1, ns)
                rows.append(r)
                print(f"  {tr} CTCNP1={c1} nitrate={ns}: "
                      f"RRMSE {r['yield_rrmse']}%  N0 {r['N0_model']:.0f}"
                      f" (目标 {r['N0_target']:.0f})  139→253 {r['resp_139_253']}%")
    d = pd.DataFrame(rows)
    for tr, g in d.groupby("trial"):
        out = dc.OUT / "diag" / f"nresponse_tune_{tr}.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        g.to_csv(out, index=False)
    d.to_csv(dc.OUT / "diag" / "nresponse_tune.csv", index=False)
    print("\n", d.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
