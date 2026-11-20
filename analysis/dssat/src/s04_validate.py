"""Validation statistics for the DSSAT run.

Reads out/sim/dssat_daily.csv + dssat_summary.csv and the standardised
observations, then writes

  out/sim/validate_plots.csv    one row per plot
  out/sim/validate_stats.csv    one row per metric
  out/sim/lai_plateau.csv       observed vs simulated plateau summary

Usage
-----
    python analysis/dssat/src/s04_validate.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

PLATEAU_THRESHOLD = 0.95


def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    daily = pd.read_csv(dc.SIM / "dssat_daily.csv", parse_dates=["date"])
    summ = pd.read_csv(dc.SIM / "dssat_summary.csv").set_index("plot")
    obs = dc.load_obs()
    yld = dc.load_yield()
    return daily, summ, obs, yld


def per_plot(daily, summ, obs, yld) -> pd.DataFrame:
    rows = []
    for plot in dc.ORDER:
        g = daily[daily["plot"] == plot].sort_values("date")
        o = obs[obs["plot"] == plot].dropna(subset=["lai_plot"])
        gg = g.set_index(g["date"].dt.normalize())
        matched = [(gg.loc[pd.Timestamp(r.date).normalize(), "LAID"], r.lai_plot)
                   for r in o.itertuples()
                   if pd.Timestamp(r.date).normalize() in gg.index]
        sim_lai = np.array([m[0] for m in matched])
        obs_lai = np.array([m[1] for m in matched])

        sim_p = dc.plateau_metrics(g["date"], g["LAID"], PLATEAU_THRESHOLD)
        obs_p = dc.plateau_metrics(o["date"], o["lai_plot"], PLATEAU_THRESHOLD)

        o_bio = obs.loc[obs["plot"] == plot, "bio_kg_ha"].dropna()
        sim_bio = float(g["CWAD"].max())
        obs_bio = float(o_bio.max()) if len(o_bio) else np.nan
        obs_y = float(yld.loc[plot, "yield_dry_kg_ha"])
        sim_y = float(summ.loc[plot, "HWAM"])

        rows.append(dict(
            plot=plot, n_rate=dc.N_RATE[plot],
            density=float(summ.loc[plot, "density"]),
            obs_yield=obs_y, sim_yield=sim_y,
            yield_bias=(sim_y / obs_y - 1) * 100,
            obs_bio=obs_bio, sim_bio=sim_bio,
            bio_bias=(sim_bio / obs_bio - 1) * 100 if obs_bio else np.nan,
            obs_hi=obs_y / obs_bio if obs_bio else np.nan,
            sim_hi=sim_y / sim_bio if sim_bio else np.nan,
            lai_rmse=(float(np.sqrt(np.mean((sim_lai - obs_lai) ** 2)))
                      if len(matched) else np.nan),
            lai_n=len(matched),
            lai_peak_obs=float(obs_lai.max()) if len(obs_lai) else np.nan,
            lai_peak_sim=float(g["LAID"].max()),
            plateau_days_sim=sim_p["span_days"],
            plateau_change_sim=sim_p["change_pct"],
            plateau_days_obs=obs_p["span_days"],
            plateau_change_obs=obs_p["change_pct"],
            lai_end=float(g["LAID"].iloc[-1]),
            lai_decline_pct=(g["LAID"].iloc[-1] / g["LAID"].max() - 1) * 100,
            anthesis=dc.dssat_date(summ.loc[plot, "ADAT"]),
            maturity=dc.dssat_date(summ.loc[plot, "MDAT"]),
            nuptake=summ.loc[plot, "NUCM"], nmineralised=summ.loc[plot, "NMINC"],
            n_applied=summ.loc[plot, "NICM"],
        ))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    ap.add_argument("--scenario", default="base", choices=list(dc.SCENARIOS))
    args = ap.parse_args()
    dc.set_trial(args.trial)
    dc.set_scenario(args.scenario)
    print(f"trial {dc.TRIAL}  ({dc.CFG.site} {dc.CFG.variety})  "
          f"scenario {dc.SCENARIO} ({dc.SCEN_LABEL})")

    daily, summ, obs, yld = load()
    val = per_plot(daily, summ, obs, yld)
    val.to_csv(dc.SIM / "validate_plots.csv", index=False)

    o, s = val["obs_yield"].values, val["sim_yield"].values
    ob, sb = val["obs_bio"].values, val["sim_bio"].values
    stats = pd.DataFrame([
        dict(metric="yield_rrmse_pct",
             value=float(np.sqrt(np.mean((s - o) ** 2)) / o.mean() * 100)),
        dict(metric="yield_bias_pct", value=float(s.sum() / o.sum() * 100 - 100)),
        dict(metric="lai_rmse", value=float(np.nanmean(val["lai_rmse"]))
             if val["lai_rmse"].notna().any() else np.nan),
        dict(metric="biomass_rrmse_pct",
             value=float(np.sqrt(np.mean((sb - ob) ** 2)) / ob.mean() * 100)),
        dict(metric="biomass_bias_pct", value=float(sb.sum() / ob.sum() * 100 - 100)),
        dict(metric="hi_obs", value=float(val["obs_hi"].mean())),
        dict(metric="hi_sim", value=float(val["sim_hi"].mean())),
        dict(metric="plateau_days_sim", value=float(val["plateau_days_sim"].mean())),
        dict(metric="plateau_days_obs", value=float(val["plateau_days_obs"].mean())),
        dict(metric="plateau_change_sim_pct",
             value=float(val["plateau_change_sim"].mean())),
        dict(metric="plateau_change_obs_pct",
             value=float(val["plateau_change_obs"].mean())),
        dict(metric="anthesis_median", value=str(val["anthesis"].median().date())),
        dict(metric="maturity_median", value=str(val["maturity"].median().date())),
    ])
    stats.to_csv(dc.SIM / "validate_stats.csv", index=False)

    plat = val[["plot", "n_rate", "lai_peak_obs", "lai_peak_sim",
                "plateau_days_obs", "plateau_days_sim",
                "plateau_change_obs", "plateau_change_sim",
                "lai_end", "lai_decline_pct"]].copy()
    plat.to_csv(dc.SIM / "lai_plateau.csv", index=False)

    show = val[["plot", "n_rate", "obs_yield", "sim_yield", "yield_bias",
                "lai_rmse", "lai_peak_obs", "lai_peak_sim",
                "plateau_days_sim", "plateau_change_sim"]]
    print(show.round(2).to_string(index=False))
    print()
    print(stats.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
