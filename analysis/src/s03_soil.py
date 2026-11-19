"""Step 03 - build the per-plot APSIM soil node from the field soil sampling.

Measured            : SAT, FC(=DUL), sand/silt/clay (0-20, 20-40 cm), OM, pH
Pedotransfer (S&R06): LL15, bulk density, KS
Assumed / calibrated: initial mineral N (no NO3-N / NH4-N was measured),
                      deep-layer values (only 0-40 cm was sampled)

Outputs one JSON per plot holding the fully-populated `Soil` node so step 04 can
substitute it into the template without re-deriving anything.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DS, SOIL  # noqa: E402

TEMPLATE = Path(__file__).resolve().parents[1] / "out" / "template" / "Maize.apsimx"

# layers (mm): 0-150, 150-300, 300-600, 600-900, 900-1200, 1200-1800
THICK = [150.0, 150.0, 300.0, 300.0, 300.0, 600.0]

# initial mineral N: total kg N/ha over 0-1800 mm and its vertical share
# Initial mineral N has never been measured at either site (only M3-P/K, OM and
# pH were analysed), so it is a calibration parameter.  A scan over the Shandong
# trial gives the minimum yield RRMSE at 110 kg N/ha:
#     40 -> RRMSE 14.2 %   A1 -45.7 %
#     70 -> RRMSE 11.5 %   A1 -30.8 %
#    110 -> RRMSE  8.4 %   A1  -0.8 %   <-- chosen
#    150 -> RRMSE 12.1 %   A1 +26.2 %
# Above ~110 kg N/ha the fertilised treatments stop changing (they are already
# N-saturated) and only the zero-N control keeps rising, so the optimum is set
# by A1 alone.
N_TOTAL_KG_HA = 110.0
N_SHARE = np.array([0.28, 0.22, 0.26, 0.14, 0.07, 0.03])
N_SHARE = N_SHARE / N_SHARE.sum()
NH4_FRACTION = 0.10

# pH by site (0-30 cm measured)
PH = {"SD_SF": 8.16, "XX_SF": 8.39, "XX_IRR": 8.39}


def find(nodes, name):
    for n in nodes:
        if n.get("Name") == name:
            return n
    raise KeyError(name)


def ppm_from_kgha(kg_ha: float, bd: float, thick_mm: float) -> float:
    """kg N/ha -> ppm (mg/kg) for one layer."""
    return kg_ha / (bd * thick_mm * 0.01)


def build(prof: pd.DataFrame, trial_id: str, plot: str) -> dict:
    g = prof[(prof["trial_id"] == trial_id) & (prof["plot"] == plot)].sort_values("layer")
    assert len(g) == len(THICK), f"{trial_id}/{plot}: expected {len(THICK)} layers, got {len(g)}"
    bd = g.bd.tolist()
    dul = g.dul.tolist()
    sat = g.sat.tolist()
    ll15 = [min(l, d - 0.01) for l, d in zip(g.ll15.tolist(), dul)]
    airdry = [round(l * 0.5, 4) for l in ll15]
    ks = g.ks.tolist()
    carbon = g.carbon_pct.tolist()

    tmpl = json.loads(TEMPLATE.read_text(encoding="utf-8-sig"))
    field = find(tmpl["Children"][0]["Children"], "Field")
    soil = json.loads(json.dumps(find(field["Children"], "Soil")))  # deep copy
    kids = {c["Name"]: c for c in soil["Children"]}

    # --- physical -------------------------------------------------------
    ph = kids["Physical"]
    ph.update(Thickness=THICK, BD=bd, AirDry=airdry, LL15=ll15, DUL=dul, SAT=sat, KS=ks)
    ph["ParticleSizeSand"] = g.sand.tolist()
    ph["ParticleSizeClay"] = g.clay.tolist()
    ph["ParticleSizeSilt"] = [100 - s - c for s, c in zip(g.sand, g.clay)]
    sc = find(ph["Children"], "MaizeSoil")
    sc.update(LL=ll15, KL=[0.10, 0.08, 0.08, 0.06, 0.04, 0.02],
              XF=[1.0] * 6, Thickness=THICK)

    # --- water ----------------------------------------------------------
    wat = kids["SoilWater"]
    wat["Thickness"] = THICK
    wat["SWCON"] = [0.55, 0.45, 0.35, 0.30, 0.25, 0.20]

    # --- organic --------------------------------------------------------
    org = kids["Organic"]
    org.update(Thickness=THICK, Carbon=carbon,
               SoilCNRatio=[12.0] * 6,
               FBiom=[0.04, 0.03, 0.02, 0.01, 0.01, 0.01],
               FInert=[0.40, 0.60, 0.80, 0.95, 1.00, 1.00],
               FOM=[350.0, 220.0, 150.0, 60.0, 30.0, 20.0])

    # --- chemical -------------------------------------------------------
    ch = kids["Chemical"]
    ch.update(Thickness=THICK, PH=[PH[trial_id]] * 6, PHUnits=0)

    # --- initial water: profile at field capacity at sowing ---------------
    # Uniform assumption requested for all plots: the profile starts at its
    # measured water content at -33 kPa (DUL), i.e. 100% of plant-available
    # water.  This is a standard pre-season assumption in crop modelling when
    # no sowing-date measurement exists (no such measurement for Shandong).
    init_theta = list(dul)
    w = kids["Water"]
    w.update(Thickness=THICK, InitialValues=[round(x, 4) for x in init_theta],
             InitialPAWmm=round(sum((t * (x - ll)) for t, x, ll in
                                    zip(THICK, init_theta, ll15)), 1),
             RelativeTo="LL15", FilledFromTop=False)

    # --- initial mineral N (calibration parameter) ----------------------
    n_layer = N_TOTAL_KG_HA * N_SHARE
    no3 = [round(ppm_from_kgha(v, b, t), 3)
           for v, b, t in zip(n_layer * (1 - NH4_FRACTION), bd, THICK)]
    nh4 = [round(ppm_from_kgha(v, b, t), 3)
           for v, b, t in zip(n_layer * NH4_FRACTION, bd, THICK)]
    kids["NO3"].update(Thickness=THICK, InitialValues=no3, InitialValuesUnits=0)
    kids["NH4"].update(Thickness=THICK, InitialValues=nh4, InitialValuesUnits=0)
    kids["Urea"].update(Thickness=THICK, InitialValues=[0.0] * 6)

    soil["Name"] = "Soil"
    soil["LocalName"] = f"{trial_id}_{plot}"
    soil["DataSource"] = "field sampling 2026 + Saxton-Rawls PTF"
    return soil


def main() -> None:
    prof = pd.read_csv(DS / "soil_profile.csv")
    keys = prof[["trial_id", "plot"]].drop_duplicates()
    for _, r in keys.iterrows():
        node = build(prof, r["trial_id"], r["plot"])
        (SOIL / f'{r["trial_id"]}_{r["plot"]}.json').write_text(
            json.dumps(node, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"wrote {len(keys)} soil nodes to {SOIL}")

    d = prof[prof.trial_id == "SD_SF"].groupby("layer").agg(
        ll15=("ll15", "mean"), dul=("dul", "mean"), sat=("sat", "mean"),
        bd=("bd", "mean"), ks=("ks", "mean"))
    print("\nSD_SF mean profile (volumetric):")
    print(d.round(3).to_string())
    print(f"\ninitial mineral N over 0-1800 mm = {N_TOTAL_KG_HA:.0f} kg N/ha "
          f"(NO3 {N_TOTAL_KG_HA*(1-NH4_FRACTION):.0f} + NH4 {N_TOTAL_KG_HA*NH4_FRACTION:.0f})")


if __name__ == "__main__":
    main()
