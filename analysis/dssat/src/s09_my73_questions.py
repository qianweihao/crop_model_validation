"""Three follow-up questions about the two MY73 trials.

Q1  SD_MY73 has no stand count - how much would the density actually matter?
Q2  can the per-plant leaf area gap be tuned with IXIM-Maize (AX / ALL)?
Q3  can one cultivar parameter set describe MY73 at both sites?

Outputs (out/<trial>/diag/)
    density_sensitivity_<trial>.csv
    ixim_ax_sweep.csv
    my73_shared_cultivar.csv
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


def _diag(trial: str | None = None) -> Path:
    """Always write under the trial the question belongs to."""
    d = dc.BASE / "out" / (trial or dc.TRIAL) / "diag"
    d.mkdir(parents=True, exist_ok=True)
    return d


def yield_rrmse(sim, obs):
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    return float(np.sqrt(np.mean((sim - obs) ** 2)) / obs.mean() * 100)


def obs_yield(trial: str) -> pd.Series:
    y = pd.read_csv(dc.DATASET / "yield.csv")
    y = y[y.trial_id == trial].set_index("plot")
    return y["yield_dry_kg_ha"]


def run_all(plots, cultivar, density=None):
    out = {}
    for p in plots:
        g, s = dc.run_plot(p, cultivar, density=density)
        m = dc.plot_metrics(g, s, density or _dens(p))
        out[p] = m
    return out


def _dens(plot):
    pl = dc.load_plots().loc[plot]
    v = pl["density_measured"]
    if pd.isna(v):
        v = pl["density_planted"]
    return float(v) / MU


# --------------------------------------------------------------------------
def q1_density_sensitivity() -> pd.DataFrame:
    """SD_MY73 has no stand count: show what a +/- 500 plants/mu would do."""
    dc.set_trial("SD_MY73")
    dc.set_model("CERES")
    plot = dc.ORDER[0]
    rows = []
    for per_mu in (5000, 5500, 6000, 6500, 7000):
        dens = per_mu / MU
        g, s = dc.run_plot(plot, density=dens)
        m = dc.plot_metrics(g, s, dens)
        rows.append(dict(plot=plot, plants_per_mu=per_mu, density_m2=dens,
                         lai_peak=m["lai_peak"], la_per_plant=m["la_per_plant"],
                         bio=m["bio_final"], yield_kg=m["yield_kg"], hi=m["hi"],
                         grains_per_ear=m["grains_per_ear"]))
    d = pd.DataFrame(rows)
    d.to_csv(_diag() / "density_sensitivity_SD_MY73.csv", index=False)
    print(d.round(3).to_string(index=False))
    return d


# --------------------------------------------------------------------------
def q2_ixim() -> pd.DataFrame:
    """IXIM has cultivar-level AX (largest leaf area) and ALL (longevity)."""
    rows = []
    for trial, plot in (("XX_SF", "A2B2"), ("SD_SF", "A5B2")):
        dc.set_trial(trial)
        dens = _dens(plot)
        # CERES reference
        dc.set_model("CERES")
        g, s = dc.run_plot(plot)
        m = dc.plot_metrics(g, s, dens)
        rows.append(dict(trial=trial, plot=plot, model="CERES", AX=np.nan,
                         ALL=np.nan, lai_peak=m["lai_peak"],
                         la_per_plant=m["la_per_plant"], bio=m["bio_final"],
                         yield_kg=m["yield_kg"]))
        dc.set_model("IXIM")
        for ax, allv in ((650, 800), (850, 800), (950, 800), (950, 900)):
            cv = dict(dc.DEFAULT_CULTIVAR, AX=float(ax), ALL=float(allv))
            g, s = dc.run_plot(plot, cv)
            m = dc.plot_metrics(g, s, dens)
            rows.append(dict(trial=trial, plot=plot, model="IXIM", AX=ax,
                             ALL=allv, lai_peak=m["lai_peak"],
                             la_per_plant=m["la_per_plant"], bio=m["bio_final"],
                             yield_kg=m["yield_kg"]))
    dc.set_model("CERES")
    d = pd.DataFrame(rows)
    d.to_csv(_diag("XX_SF") / "ixim_ax_sweep.csv", index=False)
    print(d.round(3).to_string(index=False))
    return d


# --------------------------------------------------------------------------
def q3_shared_cultivar() -> pd.DataFrame:
    """One set of cultivar coefficients for MY73 at both sites."""
    rows = []
    for p5 in (900.0, 940.0, 980.0):
        for trial in ("XX_SF", "SD_MY73"):
            dc.set_trial(trial)
            dc.set_model("CERES")
            cv = dict(dc.DEFAULT_CULTIVAR, P5=p5)
            obs = obs_yield(trial)
            res = run_all(dc.ORDER, cv)
            sim = np.array([res[p]["yield_kg"] for p in dc.ORDER])
            o = obs.reindex(dc.ORDER).values
            ant = pd.Series([res[p]["anthesis"] for p in dc.ORDER])
            mat = pd.Series([res[p]["maturity"] for p in dc.ORDER])
            rows.append(dict(
                P5=p5, trial=trial, n_plots=len(dc.ORDER),
                yield_rrmse=yield_rrmse(sim, o),
                yield_bias=(sim.sum() / o.sum() - 1) * 100,
                anthesis_med=ant.median().date(),
                maturity_med=mat.median().date(),
                lai_peak=float(np.mean([res[p]["lai_peak"] for p in dc.ORDER])),
                grains_per_ear=float(np.mean([res[p]["grains_per_ear"]
                                              for p in dc.ORDER]))))
    d = pd.DataFrame(rows)
    d["obs_anthesis"] = d["trial"].map({"XX_SF": "2026-08-08",
                                        "SD_MY73": "2026-08-08"})
    d["obs_maturity"] = d["trial"].map({"XX_SF": "2026-10-05 (推定)",
                                        "SD_MY73": "2026-10-08"})
    d["d_anthesis"] = [(pd.Timestamp(a) - pd.Timestamp(b)).days
                       for a, b in zip(d.anthesis_med, d.obs_anthesis)]
    d.to_csv(_diag("XX_SF") / "my73_shared_cultivar.csv", index=False)
    print(d.round(2).to_string(index=False))
    return d


def q3_grid() -> pd.DataFrame:
    """Joint (P5, G3) search: one cultivar set scored on both trials."""
    rows = []
    for p5 in (940.0, 960.0, 980.0):
        for g3 in (5.4, 5.8, 6.2):
            rec = dict(P5=p5, G3=g3)
            for trial in ("XX_SF", "SD_MY73"):
                dc.set_trial(trial)
                dc.set_model("CERES")
                cv = dict(dc.DEFAULT_CULTIVAR, P5=p5, G3=g3)
                obs = obs_yield(trial)
                res = run_all(dc.ORDER, cv)
                sim = np.array([res[p]["yield_kg"] for p in dc.ORDER])
                o = obs.reindex(dc.ORDER).values
                rec[f"{trial}_rrmse"] = yield_rrmse(sim, o)
                rec[f"{trial}_bias"] = (sim.sum() / o.sum() - 1) * 100
                rec[f"{trial}_mat"] = pd.Series(
                    [res[p]["maturity"] for p in dc.ORDER]).median().date()
            rec["score"] = np.mean([rec["XX_SF_rrmse"], rec["SD_MY73_rrmse"]])
            rows.append(rec)
    d = pd.DataFrame(rows).sort_values("score")
    d.to_csv(_diag("XX_SF") / "my73_shared_cultivar_grid.csv", index=False)
    print(d.round(2).to_string(index=False))
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="all",
                    choices=["all", "q1", "q2", "q3", "q3grid"])
    args = ap.parse_args()
    if args.section in ("all", "q1"):
        print("\n[Q1] SD_MY73 density sensitivity (no stand count in the data)")
        q1_density_sensitivity()
    if args.section in ("all", "q2"):
        print("\n[Q2] IXIM-Maize cultivar-level leaf area (AX) vs CERES")
        q2_ixim()
    if args.section in ("all", "q3"):
        print("\n[Q3] one MY73 cultivar for both sites")
        q3_shared_cultivar()
    if args.section in ("q3grid",):
        print("\n[Q3] joint P5 x G3 grid")
        q3_grid()


if __name__ == "__main__":
    main()
