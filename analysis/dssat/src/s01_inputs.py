"""Write the DSSAT input deck (weather, soil) and an input audit table.

The management file is generated per plot at run time by ``dssat_common``
because DSSATTools needs the soil object to fill the field/initial-condition
sections; a copy of the generated .MZX is kept by ``s02_run.py``.

Usage
-----
    python analysis/dssat/src/s01_inputs.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402


def _density(plot: str) -> float:
    """Stand count when it exists, otherwise the seed rate (SD_MY73)."""
    p = dc.load_plots().loc[plot]
    v = p["density_measured"]
    if pd.isna(v):
        v = p["density_planted"]
    return float(v) / dc.PLANTS_PER_MU


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    args = ap.parse_args()
    dc.set_trial(args.trial)
    print(f"trial {dc.TRIAL}  ({dc.CFG.site} {dc.CFG.variety})  -> {dc.OUT}")

    # --- weather ---------------------------------------------------------
    wth = dc.build_weather()
    wth.write(str(dc.WTH))
    met = pd.read_csv(dc.DATASET / "weather_daily.csv", parse_dates=["date"])
    met = met[met["zone"] == dc.CFG.zone]
    print(f"weather  {len(met)} days  {met.date.min():%Y-%m-%d}..{met.date.max():%Y-%m-%d}  "
          f"rain {met.rain.sum():.1f} mm -> {dc.WTH}")

    # --- soil ------------------------------------------------------------
    rows = []
    for plot in dc.ORDER:
        soil = dc.build_soil(plot)
        soil.write(str(dc.SOIL / f"{soil.id}.SOL"))
        top = soil.layers[min(soil.layers)]
        rows.append(dict(
            plot=plot, soil_id=soil.id, layers=soil.n_layers,
            depth_cm=soil.total_depth,
            slll=top["SLLL"], sdul=top["SDUL"], ssat=top["SSAT"],
            bd=top["SBDM"], sloc_top=top["SLOC"],
            awc_mm=sum((l["SDUL"] - l["SLLL"]) * (d - prev)
                       for (d, l), prev in zip(
                           sorted(soil.layers.items()),
                           [0] + sorted(soil.layers)[:-1])) * 10,
            no3_kg=dc.INIT_NO3_KG * dc.DEFAULT_NITRATE_SCALE,
            nh4_kg=dc.INIT_NH4_KG,
            density_plants_m2=_density(plot),
            n_rate=dc.N_RATE[plot],
        ))
    audit = pd.DataFrame(rows)
    audit.to_csv(dc.SIM / "inputs_audit.csv", index=False)
    print(audit.round(3).to_string(index=False))

    # --- management schedules (one row per event, for the report) --------
    fert, irrig = dc.load_schedule()
    sched = []
    for plot in dc.ORDER:
        for _, r in fert[fert["plot"] == plot].iterrows():
            sched.append(dict(plot=plot, date=r["date"].date(), op="urea",
                              amount=r["n_kg_ha"], unit="kg N/ha"))
        for _, r in irrig[irrig["plot"] == plot].iterrows():
            sched.append(dict(plot=plot, date=r["date"].date(), op="irrigation",
                              amount=r["mm"], unit="mm"))
    sch = pd.DataFrame(sched).sort_values(["plot", "date"])
    sch.to_csv(dc.SIM / "management_schedule.csv", index=False)
    print(f"management events: {len(sch)}")


if __name__ == "__main__":
    main()
