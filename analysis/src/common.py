"""Shared paths, trial registry, loaders and pedotransfer functions.

Project : zoomlion_crop_model - maize N x water stress simulation
Purpose : single source of truth for where data lives and how raw sheets map
          onto the canonical long tables used by every downstream step.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "产量预估输入数据"
OUT = ROOT / "analysis" / "out"
DS = OUT / "dataset"
MET = OUT / "met"
SOIL = OUT / "soil"
APSIM_DIR = OUT / "apsim"
SIM = OUT / "sim"
FIG = OUT / "fig"
REPORTS = ROOT / "analysis" / "reports"

for _d in (DS, MET, SOIL, APSIM_DIR, SIM, FIG, REPORTS):
    _d.mkdir(parents=True, exist_ok=True)

MGMT_XLSX = RAW / "玉米试验数据分类整理_施肥灌溉记录_产量.xlsx"
OBS_XLSX = RAW / "山东德州_河南新乡_单株叶面积生物量提取表.xlsx"
WTH_XLSX = RAW / "各区域各品种逐日气象信息_含播种前10天_截至2026-10-15.xlsx"
PHEN_XLSX = RAW / "品种生育期_汇总.xlsx"
NUT_XLSX = RAW / "种植前土壤养分测量.xlsx"
DSS_WTH = Path(r"E:/drought_smart_fertigation/drought_smart_fertigation/data/weather/standard_format.wth")

# --------------------------------------------------------------------------
# trial registry
#   trial_id -> everything needed to build a simulation unit
#   role: 'gradient'  -> N and/or water contrast, used for stress validation
#         'support'   -> usable for LAI / phenology only (no density, n=1)
# --------------------------------------------------------------------------
TRIALS = {
    "SD_SF": dict(
        site="山东德州", trial="水肥试验", variety="乐农87", weather="dezhou",
        obs_sheet="山东德州_乐农87", plot_col="实验地块", has_density=True,
        has_mgmt=True, role="gradient",
    ),
    "XX_SF": dict(
        site="河南新乡", trial="水肥试验", variety="MY73", weather="xinxiang",
        obs_sheet="河南新乡_MY73", plot_col="编号", has_density=True,
        has_mgmt=True, role="gradient",
    ),
    "XX_IRR": dict(
        site="河南新乡", trial="灌溉试验", variety="金华创826", weather="xinxiang",
        obs_sheet="河南新乡_金华创826", plot_col="编号", has_density=False,
        has_mgmt=True, role="support",
    ),
    "SD_MY73": dict(
        site="山东德州", trial="外围试验区", variety="MY73", weather="dezhou",
        obs_sheet="山东德州_MY73", plot_col=("实验小区", "MY73"), has_density=False,
        has_mgmt=False, role="support",
    ),
    "SD_DH997": dict(
        site="山东德州", trial="外围试验区", variety="登海997", weather="dezhou",
        obs_sheet="山东德州_登海997", plot_col=("实验小区", "DH997"), has_density=False,
        has_mgmt=False, role="support",
    ),
}

WEATHER_SHEETS = {
    "dezhou": "山东试验区_山东德州_乐农87",
    "xinxiang": "新乡灌溉_新乡_金华创826",
}
WEATHER_SITE = {
    "dezhou": dict(lat=37.4697, lon=116.9463, elev=22.0, name="Dezhou"),
    "xinxiang": dict(lat=35.1349, lon=113.7630, elev=75.0, name="Xinxiang"),
}

# canonical observation columns produced by s01
OBS_COLS = [
    "trial_id", "site", "trial", "variety", "plot", "date", "das",
    "plant_la_cm2", "leaf_dw", "stem_dw", "root_dw", "ear_dw", "grain_dw",
    "total_dw", "density_planted", "density_measured",
    "lai_plot", "biomass_kg_mu", "lai_std_nominal", "src_sheet", "src_row",
]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def norm_plot(x) -> str:
    """Normalise plot codes: 'A2-B1' -> 'A2B1'."""
    s = str(x).strip()
    if s.lower() == "nan":
        return ""
    return re.sub(r"[\s\-_]+", "", s).upper()


def read_header_sheet(path: Path, sheet: str, key_hint: str = "区域") -> pd.DataFrame:
    """Read an observation sheet whose header row is not row 0."""
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
    hr = next(
        i for i in range(8)
        if raw.iloc[i].astype(str).str.contains(key_hint).any()
        and raw.iloc[i].astype(str).str.contains("品种").any()
    )
    df = pd.read_excel(path, sheet_name=sheet, header=hr)
    df.columns = [str(c).strip() for c in df.columns]
    df = df[df["日期"].notna()].copy()
    df["日期"] = pd.to_datetime(df["日期"], errors="coerce")
    return df.dropna(subset=["日期"]).reset_index(drop=True)


def col(df: pd.DataFrame, prefix: str):
    """First column whose name starts with `prefix`."""
    for c in df.columns:
        if c.startswith(prefix):
            return c
    return None


def num(series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


# --------------------------------------------------------------------------
# pedotransfer: Saxton & Rawls (2006)
#   used only for LL15 (wilting point) and bulk density / KS, because the
#   field campaign measured SAT and FC but not wilting point.
# --------------------------------------------------------------------------
def saxton_rawls(sand_pct: float, clay_pct: float, om_pct: float) -> dict:
    A, C, OM = sand_pct / 100.0, clay_pct / 100.0, om_pct
    t1500 = (-0.024 * A + 0.487 * C + 0.006 * OM + 0.005 * A * OM
             - 0.013 * C * OM + 0.068 * A * C + 0.031)
    t33 = (-0.251 * A + 0.195 * C + 0.011 * OM + 0.006 * A * OM
           - 0.027 * C * OM + 0.452 * A * C + 0.299)
    ts33 = (0.278 * A + 0.034 * C + 0.022 * OM - 0.018 * A * OM
            - 0.027 * C * OM - 0.584 * A * C + 0.078)
    bd = (1.0 - ts33) * 2.65
    # saturated conductivity, Saxton & Rawls eq. (mm/h)
    lam = 1.0 / np.log(max(t33, 1e-6) / max(t1500, 1e-6))
    ks = 1930.0 * max(ts33 - t33, 1e-4) ** (3.0 - 1.0 / lam)
    return dict(theta_wp=t1500, theta_fc=t33, theta_sat=ts33, bulk_density=bd,
                ks_mm_h=ks)


def sr2006_ll15(sand_pct: float, clay_pct: float, om_pct: float) -> float:
    """Volumetric water content at -1500 kPa (wilting point), Saxton & Rawls (2006).

        theta_1500 = -0.024S + 0.487C + 0.006OM + 0.005(S*OM)
                     - 0.013(C*OM) + 0.068(S*C) + 0.031

    with S and C as fractions.  Validated against textbook values for sand
    (4.1 %), silt loam (11.4 %) and clay (32.4 %), and its predicted field
    capacity for this site (30.7 %) agrees with the measured 31.4-34.6 %.

    This replaces the previous ad-hoc LL15 = 0.46 x DUL, which had no local
    basis (the 0.46 was an artefact of dividing the driest Xinxiang sample by
    the mean field capacity across all plots).
    """
    S, C, OM = sand_pct / 100.0, clay_pct / 100.0, om_pct
    t1500 = (-0.024 * S + 0.487 * C + 0.006 * OM + 0.005 * S * OM
             - 0.013 * C * OM + 0.068 * S * C + 0.031)
    return float(np.clip(t1500, 0.04, 0.30))


def build_soil_profile(sand: float, clay: float, sat: float, fc: float,
                       om_pct: float, layers: tuple[int, ...] = (15, 30, 60, 90, 120, 180)) -> pd.DataFrame:
    """Depth-weighted APSIM soil layers for one plot.

    Measured SAT and FC are authoritative.  LL15 is set at 0.46 x DUL, which the
    2026 pre/post-irrigation moisture readings support: the lowest observed
    profile water (12.4-15.5 vol %) sits at 0.46-0.48 of the measured FC, so
    0.46 x DUL is a defensible wilting-point proxy for these silt loams.
    Bulk density follows from the measured porosity (capped at 1.60 g/cm3 because
    the reported SAT values are lower than textbook porosity for these textures).
    KS is a texture-class estimate, not a measurement.
    """
    bd = float(np.clip((1.0 - sat) * 2.65, 1.15, 1.60))
    # LL15 from texture (Saxton-Rawls 2006); measured SAT/FC stay authoritative
    ll15 = float(np.clip(sr2006_ll15(sand, clay, om_pct), 0.04, fc - 0.02))
    ks_by_layer = [40.0, 30.0, 20.0, 12.0, 8.0, 5.0]
    rows = []
    for i, bottom_cm in enumerate(layers):
        prev_cm = layers[i - 1] if i else 0
        thick = (bottom_cm - prev_cm) * 10.0      # layers are given in cm -> mm
        rows.append(dict(
            layer=i + 1, thickness_mm=thick,
            bottom_cm=bottom_cm,
            bd=round(bd, 3), ll15=round(ll15, 4), dul=round(fc, 4),
            sat=round(sat, 4), air_dry=round(ll15 * 0.5, 4),
            ks=ks_by_layer[min(i, len(ks_by_layer) - 1)],
            carbon_pct=round(om_pct / 1.724, 3),
            sand=sand, clay=clay,
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# APSIM constants
# --------------------------------------------------------------------------
MG_PER_MU_TO_KG_HA = 15.0      # 1 亩 = 666.7 m2 -> kg/亩 * 15 = kg/ha
# Farmer-reported grain yield in this project is at 14% moisture, while APSIM
# reports DRY grain weight.  Every observed-vs-simulated comparison must be done
# on one basis; we standardise on dry weight.
GRAIN_MOISTURE = 0.14
DRY_FACTOR = 1.0 - GRAIN_MOISTURE

# observation dates whose organ-level split is internally inconsistent
# (2026-08-14: stem dry weight falls 41->30 g while leaf jumps 27->55 g, which is
# physically impossible).  Total above-ground weight on that date is still usable,
# so only organ-level targets are suppressed.
ORGAN_SPLIT_SUSPECT = {("SD_SF", "2026-08-14")}
MU_TO_HA = 1.0 / 15.0
LAI_AREA_FACTOR = 6666667.0     # cm2/plant * plants/mu / 6666667 = m2/m2


def critical_n_pct(w_t_ha):
    """Critical shoot nitrogen concentration (%DM) for maize.

    Plenet & Lemaire (2000) dilution curve, valid for W >= 1 t/ha.
    """
    return 3.4 * np.power(np.asarray(w_t_ha, dtype=float), -0.37)


def nitrogen_nutrition_index(wt_g_m2, n_g_m2, min_w=100.0):
    """NNI = plant N% / critical N%.  ``min_w`` drops the divergent low-biomass tail.

    APSIM exposes no NNI variable and ``Photosynthesis.FN`` is hard-coded to 1.0,
    so this index is derived post-hoc from the standard dilution curve.
    """
    wt = np.asarray(wt_g_m2, dtype=float)
    n = np.asarray(n_g_m2, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        n_pct = np.where(wt >= min_w, n / wt * 100.0, np.nan)
        idx = n_pct / critical_n_pct(wt * 0.01)
    return idx
