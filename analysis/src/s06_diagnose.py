"""Step 06 - diagnose where the observed N / water gradients are lost.

Compares simulated vs observed yield, biomass and LAI per plot, and quantifies
how much of each gradient the model reproduces (slope of yield on total N).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DS, REPORTS, SIM  # noqa: E402

FOCUS = ["SD_SF", "XX_SF", "XX_IRR"]


def load_sim() -> pd.DataFrame:
    w = pd.read_csv(SIM / "sim_daily_wide.csv", parse_dates=["Clock.Today"])
    for base in ("Soil.Nutrient.NO3.kgha", "Soil.SoilWater.PAW",
                 "Soil.SoilWater.ESW", "Soil.SoilWater.SWmm"):
        cols = [c for c in w.columns if c.startswith(base + "(")]
        if cols:
            w[base + ".total"] = w[cols].sum(axis=1)
    w["trial_id"] = w["sim"].str.rsplit("_", n=1).str[0]
    w["plot"] = w["sim"].str.rsplit("_", n=1).str[1]
    return w


def summary(sim: pd.DataFrame, obs_plant: pd.DataFrame, yld: pd.DataFrame,
            fert: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n_tot = fert.groupby(["trial_id", "plot"])["n_kg_ha"].sum()
    for (tid, plot), g in sim.groupby(["trial_id", "plot"]):
        g = g.sort_values("Clock.Today").reset_index(drop=True)
        st = g["Maize.Phenology.CurrentStageName"].astype(str)
        fl = g.loc[st.eq("Flowering"), "Clock.Today"]
        mt = g.loc[st.eq("Maturity"), "Clock.Today"]
        rows.append(dict(
            trial_id=tid, plot=plot,
            sim_yield=g["Maize.Grain.Total.Wt"].max() * 10.0,
            sim_biomass=g["Maize.AboveGround.Wt"].max() * 10.0,
            sim_lai_max=g["Maize.Leaf.LAI"].max(),
            sim_flower_date=fl.iloc[0].date() if len(fl) else None,
            sim_mature_date=mt.iloc[0].date() if len(mt) else None,
            n_kg_ha=n_tot.get((tid, plot), np.nan),
            days_FN_below_09=int((g["Maize.Leaf.Photosynthesis.FN"] < 0.9).sum()),
            min_FN=float(g["Maize.Leaf.Photosynthesis.FN"].min()),
            days_FW_below_09=int((g["Maize.Leaf.Photosynthesis.FW"] < 0.9).sum()),
            min_FW=float(g["Maize.Leaf.Photosynthesis.FW"].min()),
            min_NO3_kgha=float(g["Soil.Nutrient.NO3.kgha.total"].min()),
            min_ESW_mm=float(g["Soil.SoilWater.ESW.total"].min()),
            lai_on_flower=g.loc[st.eq("Flowering"), "Maize.Leaf.LAI"].max()
            if st.eq("Flowering").any() else np.nan,
        ))
    cmp = pd.DataFrame(rows)

    y = yld[["trial_id", "plot", "yield_measured"]].copy()
    y["obs_yield"] = y["yield_measured"] * 15.0
    cmp = cmp.merge(y[["trial_id", "plot", "obs_yield"]], on=["trial_id", "plot"], how="left")

    ob = obs_plant.dropna(subset=["biomass_kg_mu"])
    pk = ob.assign(v=ob["biomass_kg_mu"]).sort_values("v")
    pk = pk.groupby(["trial_id", "plot"]).tail(1)
    pk["obs_biomass"] = pk["biomass_kg_mu"] * 15.0
    cmp = cmp.merge(pk[["trial_id", "plot", "obs_biomass"]], on=["trial_id", "plot"], how="left")

    la = obs_plant.dropna(subset=["lai_plot"]).sort_values("lai_plot")
    la = la.groupby(["trial_id", "plot"]).tail(1)
    cmp = cmp.merge(la[["trial_id", "plot", "lai_plot"]].rename(
        columns={"lai_plot": "obs_lai_max"}), on=["trial_id", "plot"], how="left")
    cmp["bias_yield"] = cmp["sim_yield"] - cmp["obs_yield"]
    cmp["bias_pct"] = cmp["bias_yield"] / cmp["obs_yield"] * 100
    return cmp


def gradient_table(cmp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tid, g in cmp.dropna(subset=["n_kg_ha"]).groupby("trial_id"):
        g = g.sort_values("n_kg_ha")
        if g["obs_yield"].notna().sum() >= 3:
            so = np.polyfit(g["n_kg_ha"], g["obs_yield"], 1)[0]
        else:
            so = np.nan
        ss = np.polyfit(g["n_kg_ha"], g["sim_yield"], 1)[0]
        rows.append(dict(
            trial_id=tid,
            n_levels=len(g),
            n_range_kg_ha=f"{g['n_kg_ha'].min():.0f}-{g['n_kg_ha'].max():.0f}",
            obs_yield_range=f"{g['obs_yield'].min():.0f}-{g['obs_yield'].max():.0f}",
            sim_yield_range=f"{g['sim_yield'].min():.0f}-{g['sim_yield'].max():.0f}",
            obs_dYield_dN_kg_per_kg=round(so, 2),
            sim_dYield_dN_kg_per_kg=round(ss, 2),
            obs_span_pct=round((g["obs_yield"].max() / g["obs_yield"].min() - 1) * 100, 1),
            sim_span_pct=round((g["sim_yield"].max() / g["sim_yield"].min() - 1) * 100, 1),
        ))
    return pd.DataFrame(rows)


def main() -> None:
    sim = load_sim()
    obs_plant = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    yld = pd.read_csv(DS / "yield.csv")
    fert = pd.read_csv(DS / "fert.csv")

    cmp = summary(sim, obs_plant, yld, fert)
    cmp.to_csv(SIM / "compare_obs_sim.csv", index=False)

    pd.set_option("display.width", 250)
    print("=== per-plot comparison (yield kg/ha) ===")
    print(cmp[["trial_id", "plot", "n_kg_ha", "obs_yield", "sim_yield", "bias_pct",
               "obs_biomass", "sim_biomass", "obs_lai_max", "sim_lai_max"]]
          .round(0).to_string(index=False))

    grad = gradient_table(cmp)
    grad.to_csv(SIM / "gradient_table.csv", index=False)
    print("\n=== N-gradients: observed vs simulated ===")
    print(grad.to_string(index=False))

    print("\n=== stress exposure ===")
    print(cmp[["trial_id", "plot", "n_kg_ha", "days_FN_below_09", "min_FN",
               "days_FW_below_09", "min_FW", "min_NO3_kgha"]]
          .round(3).to_string(index=False))

    print(f"\ntables -> {SIM / 'diagnose_plots.csv'}, {SIM / 'gradient_table.csv'}")


if __name__ == "__main__":
    main()
