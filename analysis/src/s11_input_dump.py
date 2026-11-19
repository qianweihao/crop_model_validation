"""Step 11 - dump every model input exactly as APSIM receives it, for review.

Reads the generated .apsimx files (ground truth) plus the source tables, and
writes one flat table per input category to analysis/reports/inputs/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import APSIM_DIR, DS, REPORTS  # noqa: E402

OUT = REPORTS / "inputs"
OUT.mkdir(parents=True, exist_ok=True)


def field_of(sim):
    return next(c for c in sim["Children"] if c.get("Name") == "Field")


def rows_from_apsimx():
    rows = []
    for f in sorted(APSIM_DIR.glob("*.apsimx")):
        j = json.loads(f.read_text(encoding="utf-8-sig"))
        sim = j["Children"][0]
        fld = field_of(sim)
        clock = next(c for c in sim["Children"] if c.get("Name") == "Clock")
        wth = next(c for c in sim["Children"] if c.get("Name") == "Weather")
        sow = next(c for c in fld["Children"] if c.get("Name", "").startswith("Sow"))
        sch = next(c for c in fld["Children"] if c.get("Name", "").startswith("Nitrogen"))
        soil = next(c for c in fld["Children"] if c.get("Name") == "Soil")
        kids = {c["Name"]: c for c in soil["Children"]}
        p = {x["Key"]: x["Value"] for x in sow["Parameters"]}
        s = {x["Key"]: x["Value"] for x in sch["Parameters"]}
        ph, org, no3, nh4 = (kids["Physical"], kids["Organic"],
                             kids["NO3"], kids["NH4"])
        thick = np.array(ph["Thickness"], dtype=float)
        dul, ll = np.array(ph["DUL"], float), np.array(ph["LL15"], float)
        sat = np.array(ph["SAT"], float)
        rows.append(dict(
            sim=f.stem, clock_start=clock["Start"][:10], clock_end=clock["End"][:10],
            weather=Path(wth["FileName"]).name,
            sow_window=f"{p['StartDate']} .. {p['EndDate']}",
            population_m2=float(p["Population"]),
            population_per_mu=round(float(p["Population"]) * 666.667),
            cultivar=p["CultivarName"], sowing_depth_mm=float(p["SowingDepth"]),
            row_spacing_mm=float(p["RowSpacing"]),
            n_events=len([x for x in s["FertSchedule"].split(";") if x]),
            n_total_kg_ha=round(sum(float(x.split(":")[1])
                                    for x in s["FertSchedule"].split(";") if x), 1),
            irr_events=len([x for x in s["IrrigationSchedule"].split(";") if x]),
            irr_total_mm=round(sum(float(x.split(":")[1])
                                   for x in s["IrrigationSchedule"].split(";") if x), 1),
            profile_mm=float(thick.sum()),
            AWC_mm=round(float(((dul - ll) * thick).sum()), 0),
            drainable_mm=round(float(((sat - dul) * thick).sum()), 0),
            dul_top=round(float(dul[0]), 3), ll15_top=round(float(ll[0]), 3),
            sat_top=round(float(sat[0]), 3), bd_top=round(float(ph["BD"][0]), 3),
            carbon_top_pct=round(float(org["Carbon"][0]), 3),
            no3_kg_ha=round(float(sum(np.array(no3["InitialValues"], float)
                                      * np.array(ph["BD"], float) * thick * 0.01)), 1),
            nh4_kg_ha=round(float(sum(np.array(nh4["InitialValues"], float)
                                      * np.array(ph["BD"], float) * thick * 0.01)), 1),
        ))
    return pd.DataFrame(rows)


def main() -> None:
    inp = rows_from_apsimx()
    src = pd.read_csv(DS / "plots.csv")
    inp["trial_id"] = inp["sim"].str.rsplit("_", n=1).str[0]
    inp["plot"] = inp["sim"].str.rsplit("_", n=1).str[1]
    inp = inp.merge(src[["trial_id", "plot", "variety", "sow_date",
                         "density_planted", "density_measured"]],
                    on=["trial_id", "plot"], how="left")
    inp["density_source"] = np.where(inp["density_measured"].notna(),
                                     "measured", "planted (not measured)")
    inp.to_csv(OUT / "01_model_inputs.csv", index=False)

    # --- schedules in long form -------------------------------------------
    sched = []
    for f in sorted(APSIM_DIR.glob("*.apsimx")):
        j = json.loads(f.read_text(encoding="utf-8-sig"))
        fld = field_of(j["Children"][0])
        sch = next(c for c in fld["Children"] if c.get("Name", "").startswith("Nitrogen"))
        s = {x["Key"]: x["Value"] for x in sch["Parameters"]}
        for kind, key, unit in (("N", "FertSchedule", "kg N/ha"),
                                ("W", "IrrigationSchedule", "mm")):
            for item in [x for x in s[key].split(";") if x]:
                d, v = item.split(":")
                sched.append(dict(sim=f.stem, kind=kind, date=d,
                                  amount=float(v), unit=unit))
    sc = pd.DataFrame(sched)
    sc.to_csv(OUT / "02_management_schedule.csv", index=False)

    # --- weather ------------------------------------------------------------
    w = pd.read_csv(DS / "weather_daily.csv", parse_dates=["date"])
    wsum = w.groupby("zone").agg(start=("date", "min"), end=("date", "max"),
                                 days=("date", "size"), rain_mm=("rain", "sum"),
                                 radn_mean=("radn", "mean"),
                                 tmin=("tmin", "min"), tmax=("tmax", "max"))
    wsum.to_csv(OUT / "03_weather_summary.csv")
    w.to_csv(OUT / "04_weather_daily.csv", index=False)

    pd.set_option("display.width", 260)
    print("=" * 110)
    print("A. SOWING / CULTIVAR / DENSITY")
    print(inp[["sim", "sow_date", "cultivar", "population_per_mu", "density_planted",
               "density_measured", "density_source"]].to_string(index=False))
    print("\n" + "=" * 110)
    print("B. NITROGEN AND IRRIGATION TOTALS")
    print(inp[["sim", "n_events", "n_total_kg_ha", "irr_events", "irr_total_mm"]]
          .to_string(index=False))
    print("\n" + "=" * 110)
    print("C. SOIL (0-1800 mm profile)")
    print(inp[["sim", "profile_mm", "AWC_mm", "drainable_mm", "dul_top", "ll15_top",
               "sat_top", "bd_top", "carbon_top_pct", "no3_kg_ha", "nh4_kg_ha"]]
          .to_string(index=False))
    print("\n" + "=" * 110)
    print("D. WEATHER")
    print(wsum.round(1).to_string())
    print(f"\nfiles -> {OUT}")


if __name__ == "__main__":
    main()
