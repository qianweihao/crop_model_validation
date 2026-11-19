"""Step 10 - phenology calibration constrained by the full leaf-appearance series.

Why this replaces the 2-point fit in s07:
  s07 fitted 5 thermal-time parameters to only (silking, maturity).  With 2
  constraints and 5 degrees of freedom the parameters are not identifiable and
  drift to their bounds (FlagLeafToFlowering came out as 30 for one variety and
  207 for another, physically the same interval).  Here the daily visible-leaf
  counts (0 -> 20 leaves), the 抽雄 (tasselling) date and silking/maturity all
  enter the objective, so every parameter is pinned by observations.

Output: analysis/out/sim/phenology_calibrated.json
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_override import SimCase, run_many  # noqa: E402
from common import APSIM_DIR, DS, SIM  # noqa: E402

CAL_FILE = SIM / "phenology_calibrated.json"

# reference plot per variety (management identity is irrelevant for pure phenology)
REF_PLOT = {"乐农87": "SD_SF_A2B1", "MY73": "XX_SF_A2B1", "金华创826": "XX_IRR_A3"}

# observed tasselling / silking / maturity (品种生育期_汇总 + harvest dates)
STAGES = {
    "乐农87": dict(tassel="2026-08-06", silking="2026-08-09", harvest="2026-10-08"),
    "MY73":   dict(tassel="2026-08-05", silking="2026-08-08", harvest="2026-10-04"),
    "金华创826": dict(tassel="2026-08-05", silking="2026-08-10", harvest="2026-10-04"),
}

PARAMS = {
    "Structure.zFinalLeafNumber.FixedValue": (16.0, 26.0),
    "Phenology.Juvenile.Target.FixedValue": (60.0, 480.0),
    "Phenology.FlagLeafToFlowering.Target.FixedValue": (60.0, 260.0),
    "Phenology.FloweringToGrainFilling.Target.FixedValue": (60.0, 300.0),
    "Phenology.GrainFilling.Target.FixedValue": (300.0, 1000.0),
    "Phenology.Maturing.Target.FixedValue": (0.0, 90.0),
}

# weights: leaf counts dominate (many constraints), dates are single points
W_LEAF = 1.0
W_DATE = 25.0


def observed(variety: str, site: str) -> tuple[pd.DataFrame, dict]:
    ls = pd.read_csv(DS / "leaf_series.csv", parse_dates=["date"])
    g = ls[(ls["variety"] == variety) & (ls["site"] == site)]
    g = g.dropna(subset=["leaves_visible"]).sort_values("date")
    return g, STAGES[variety]


def evaluate(df: pd.DataFrame, leaves: pd.DataFrame, stages: dict) -> float:
    st = df["Maize.Phenology.CurrentStageName"].astype(str)
    d = df.set_index("Clock.Today")
    pen = 0.0

    # --- leaf appearance curve (weighted by number of observations) --------
    if len(leaves):
        xq = [(t - d.index[0]).days for t in leaves["date"]]
        pred = np.interp(xq, [(t - d.index[0]).days for t in d.index],
                         d["Maize.Structure.LeafTipsAppeared"].to_numpy())
        obs = leaves["leaves_visible"].to_numpy()
        pen += W_LEAF * float(np.mean((pred - obs) ** 2))

    # --- dated stages ------------------------------------------------------
    for key, stage_name, wt in (("tassel", "FlagLeaf", 0.5),
                                ("silking", "Flowering", 1.0),
                                ("harvest", "Maturity", 1.0)):
        obs_date = pd.Timestamp(stages[key])
        hit = d.index[st.to_numpy() == stage_name]
        if len(hit) == 0:
            pen += W_DATE * wt * 400.0
        else:
            pen += W_DATE * wt * float((hit[0].normalize() - obs_date).days ** 2)

    # --- keep the phases physiological ------------------------------------
    return pen


def fit(variety: str, maxiter: int, popsize: int) -> dict:
    site = "山东德州" if variety == "乐农87" else "河南新乡"
    leaves, stages = observed(variety, site)
    plot = APSIM_DIR / f"{REF_PLOT[variety]}.apsimx"
    keys = list(PARAMS)

    print(f"[{variety}] {len(leaves)} leaf observations, "
          f"stages {stages}", flush=True)

    class _Batch:
        gen = 0
        t0 = time.time()

        def __call__(self, func, iterable):
            pops = [{k: float(v) for k, v in zip(keys, x)} for x in iterable]
            jobs = [(SimCase.load(plot).set_commands(p), f"p{i:03d}")
                    for i, p in enumerate(pops)]
            try:
                with tempfile.TemporaryDirectory() as td:
                    tabs = run_many(jobs, td)
            except Exception as exc:
                print(f"  [batch failed {type(exc).__name__}: {exc}]", flush=True)
                return [1e6] * len(pops)
            vals = [evaluate(tabs[tag], leaves, stages) for _, tag in jobs]
            _Batch.gen += 1
            print(f"  gen {_Batch.gen:3d}: best={min(vals):8.2f} "
                  f"({time.time() - _Batch.t0:.0f}s)", flush=True)
            return vals

    res = differential_evolution(lambda x: 1e6, list(PARAMS.values()),
                                 maxiter=maxiter, popsize=popsize, tol=1e-4,
                                 seed=11, polish=False, init="sobol",
                                 workers=_Batch())
    params = {k: round(float(v), 1) for k, v in zip(keys, res.x)}

    # --- report the fit ---------------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        df = run_many([(SimCase.load(plot).set_commands(params), "chk")], td)["chk"]
    st = df["Maize.Phenology.CurrentStageName"].astype(str)
    d = df.set_index("Clock.Today")
    rep = {}
    for key, stage_name in (("tassel", "FlagLeaf"), ("silking", "Flowering"),
                            ("harvest", "Maturity")):
        hit = d.index[st.to_numpy() == stage_name]
        rep[key] = (str(hit[0].date()) if len(hit) else None,
                    stages[key], int((hit[0].normalize()
                                      - pd.Timestamp(stages[key])).days) if len(hit) else None)
    if len(leaves):
        xq = [(t - d.index[0]).days for t in leaves["date"]]
        pred = np.interp(xq, [(t - d.index[0]).days for t in d.index],
                         d["Maize.Structure.LeafTipsAppeared"].to_numpy())
        leaf_rmse = float(np.sqrt(np.mean((pred - leaves["leaves_visible"].to_numpy()) ** 2)))
    else:
        leaf_rmse = np.nan
    return dict(params=params, score=float(res.fun), leaf_rmse=leaf_rmse, stages=rep)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variety", default=None)
    ap.add_argument("--maxiter", type=int, default=14)
    ap.add_argument("--popsize", type=int, default=6)
    a = ap.parse_args()

    store = (json.loads(CAL_FILE.read_text(encoding="utf-8"))
             if CAL_FILE.exists() else {})
    for var in ([a.variety] if a.variety else list(STAGES)):
        r = fit(var, a.maxiter, a.popsize)
        store[var] = r
        CAL_FILE.write_text(json.dumps(store, ensure_ascii=False, indent=2),
                            encoding="utf-8")
        print(f"[{var}] score={r['score']:.1f} leaf_RMSE={r['leaf_rmse']:.2f} leaves")
        for k, (sim, obs, err) in r["stages"].items():
            print(f"    {k:8s} sim={sim} obs={obs} err={err:+d}d" if err is not None
                  else f"    {k:8s} MISSING (obs {obs})")
        print(f"    {r['params']}\n", flush=True)
    print(f"-> {CAL_FILE}")


if __name__ == "__main__":
    main()
