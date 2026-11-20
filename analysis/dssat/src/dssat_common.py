"""DSSAT (CERES-Maize) validation harness for the SD_SF and XX_SF trials.

Wraps ``DSSATTools`` (https://github.com/daquinterop/Py_DSSATTools) around the
same standardised observation tables used by the APSIM pipeline in
``analysis/out/dataset/``.  Everything DSSAT specific lives under
``analysis/dssat/``.

Two trials are supported and kept completely separate on disk:

    SD_SF   山东德州 乐农87      out/SD_SF/...   fig/SD_SF/...
    XX_SF   河南新乡 MY73        out/XX_SF/...   fig/XX_SF/...

Every script selects a trial with ``--trial`` and calls :func:`set_trial` before
touching any path or parameter.

The module also works around a Windows quirk of DSSATTools 2.1.x: on import it
mirrors its ``static`` tree into ``%TEMP%/DSSAT048`` through ``os.symlink``,
which raises ``WinError 1314`` unless the process holds
``SeCreateSymbolicLinkPrivilege``.  We pre-clear that temp tree and degrade
``os.symlink`` to a real copy before importing the package.
"""
from __future__ import annotations

import os
import re
import shutil
from contextlib import contextmanager
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Windows import shim (must run before DSSATTools is imported)
# --------------------------------------------------------------------------
_DSSAT_TMP = os.path.join(tempfile.gettempdir(), "DSSAT048")
if os.path.isdir(_DSSAT_TMP):
    shutil.rmtree(_DSSAT_TMP, ignore_errors=True)

_REAL_SYMLINK = os.symlink


def _symlink_or_copy(src, dst, *args, **kwargs):
    try:
        _REAL_SYMLINK(src, dst, *args, **kwargs)
    except OSError:                       # WinError 1314 / 5
        if os.path.isdir(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dst)


os.symlink = _symlink_or_copy

from DSSATTools import (  # noqa: E402
    Crop, DSSAT, Management, SoilLayer, SoilProfile, Weather,
)
from DSSATTools.soil import DATA_FMT  # noqa: E402

# DSSATTools writes the @SITE record with F6.3 for LAT/LONG; Chinese longitudes
# (116.35) overflow that field and shift every following column, which DSSAT
# then reports as "More than 19 layers in the soil profile".  DSSAT's own soil
# files use F8.3.
DATA_FMT["location"] = 2 * ["A12"] + 2 * ["F8.3"] + ["A24"]

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[3]
DATASET = ROOT / "analysis" / "out" / "dataset"
BASE = ROOT / "analysis" / "dssat"
REPORT = BASE / "report"
PLANTS_PER_MU = 667.0

INIT_NO3_KG = 99.0      # reference total used by the nitrate_scale knob
INIT_NH4_KG = 11.0
N_LAYER_WEIGHTS = [0.28, 0.22, 0.26, 0.14, 0.07, 0.03]


# --------------------------------------------------------------------------
# trial registry
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class TrialCfg:
    name: str
    site: str
    variety: str
    zone: str
    plots: tuple
    n_rate: dict
    sow: datetime
    sim_start: datetime
    sta: dict
    cultivar: dict
    nitrate_scale: float
    sloc: tuple
    control: str
    obs_anthesis: pd.Timestamp | None
    obs_maturity: pd.Timestamp | None
    soil_prefix: str
    description: str
    vrname: str
    sloc_scale: float = 1.0
    notes: dict = field(default_factory=dict)
    # trials without their own soil / management records borrow documented
    # stand-ins instead of silently falling back to DSSAT defaults
    display: str = ""                         # 报告/图标题里的中文名
    soil_source: tuple | None = None          # (trial_id, plot) to copy
    assumed_n: tuple | None = None            # ((date, kg N/ha), ...)
    assumed_irrig: tuple | None = None        # ((date, mm), ...)


TRIALS: dict[str, TrialCfg] = {
    "SD_SF": TrialCfg(
        name="SD_SF", site="山东德州", variety="乐农87", zone="dezhou",
        display="德州水肥试验",
        plots=("A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"),
        n_rate=dict(zip(
            ("A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"),
            (0, 139, 139, 177, 177, 215, 215, 254, 254))),
        sow=datetime(2026, 6, 18), sim_start=datetime(2026, 6, 8),
        sta=dict(INSI="DZ01", LAT=37.45, LON=116.35, ELEV=25.0,
                 TAV=13.6, AMP=14.8, REFHT=2.0, WNDHT=2.0, CO2=420.0,
                 DESC="Shandong Dezhou 2026"),
        cultivar=dict(P1=235.0, P2=0.30, P5=880.0, G2=1300.0, G3=6.00, PHINT=60.0),
        nitrate_scale=0.24,
        sloc=(1.096, 0.85, 0.55, 0.35, 0.22, 0.12),
        control="A1",
        obs_anthesis=pd.Timestamp("2026-08-09"),
        obs_maturity=pd.Timestamp("2026-10-08"),
        soil_prefix="SD", description="Dezhou", vrname="LENONG87",
        sloc_scale=1.0,
        notes=dict(obs_maturity="实测成熟 2026-10-08"),
    ),
    "XX_SF": TrialCfg(
        name="XX_SF", site="河南新乡", variety="MY73", zone="xinxiang",
        display="新乡水肥试验",
        plots=("CK", "A1B1", "A1B2", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2"),
        n_rate=dict(zip(
            ("CK", "A1B1", "A1B2", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2"),
            (0, 144, 144, 190, 189, 228, 227, 276, 274))),
        sow=datetime(2026, 6, 18), sim_start=datetime(2026, 6, 8),
        sta=dict(INSI="XX01", LAT=35.30, LON=113.93, ELEV=75.0,
                 TAV=14.5, AMP=13.5, REFHT=2.0, WNDHT=2.0, CO2=420.0,
                 DESC="Henan Xinxiang 2026"),
        # Calibrated in s03_calibrate.py --trial XX_SF (see report).
        cultivar=dict(P1=235.0, P2=0.30, P5=900.0, G2=1300.0, G3=6.20, PHINT=50.0),
        nitrate_scale=0.55,
        sloc=(1.346, 1.05, 0.68, 0.43, 0.27, 0.15),
        sloc_scale=0.70,
        control="CK",
        obs_anthesis=pd.Timestamp("2026-08-08"),
        # 新乡 MY73 没有观测成熟期；用最后一次生物量观测（10-04 已开始下降）推定。
        obs_maturity=pd.Timestamp("2026-10-05"),
        soil_prefix="XX", description="Xinxiang", vrname="MY73",
        notes=dict(obs_maturity="实测无成熟期，按末次观测 10-04 推定为 10-05"),
    ),
    # ------------------------------------------------------------------
    # 山东潍坊 品种展示（MY73），只有生物量/产量/生育期，没有 LAI、土壤、
    # 施肥和灌溉记录。品种参数沿用 XX_SF 标定的 MY73，土壤借 SD_SF 的
    # A3B1 剖面，施肥与灌溉按下表假设（见 report/DSSAT_PARAMETERS_SD_MY73.md）。
    # ------------------------------------------------------------------
    "SD_MY73": TrialCfg(
        name="SD_MY73", site="山东潍坊", variety="MY73", zone="dezhou",
        display="潍坊品种展示",
        plots=("W1", "W2", "W3", "W4"),
        n_rate={p: 270 for p in ("W1", "W2", "W3", "W4")},
        sow=datetime(2026, 6, 16), sim_start=datetime(2026, 6, 8),
        sta=dict(INSI="WF01", LAT=36.71, LON=119.10, ELEV=30.0,
                 TAV=13.2, AMP=14.6, REFHT=2.0, WNDHT=2.0, CO2=420.0,
                 DESC="Shandong Weifang 2026"),
        # 品种参数沿用 XX_SF 的 MY73；只有 P5 按本站实测 R6（10-08）重定为 980
        # （新乡 900 对应 10-05）。
        cultivar=dict(P1=235.0, P2=0.30, P5=980.0, G2=1300.0, G3=6.20, PHINT=50.0),
        nitrate_scale=0.24,
        sloc=(1.096, 0.85, 0.55, 0.35, 0.22, 0.12),
        sloc_scale=1.0,
        control="W1",
        obs_anthesis=pd.Timestamp("2026-08-08"),
        obs_maturity=pd.Timestamp("2026-10-08"),
        soil_prefix="WF", description="Weifang", vrname="MY73",
        soil_source=("SD_SF", "A3B1"),
        assumed_n=(("2026-06-16", 108.0), ("2026-07-05", 81.0),
                   ("2026-07-25", 54.0), ("2026-08-05", 27.0)),
        assumed_irrig=(("2026-07-16", 25.0), ("2026-07-28", 25.0),
                       ("2026-08-05", 25.0), ("2026-08-21", 25.0)),
        notes=dict(management="无施肥/灌溉记录，按 270 kg N/ha + 100 mm 假设",
                   soil="无土壤剖面，借用 SD_SF A3B1",
                   lai="无 LAI 观测"),
    ),
}

# --------------------------------------------------------------------------
# active-trial globals (rebound by set_trial)
# --------------------------------------------------------------------------
TRIAL = "SD_SF"
ORDER: list[str] = []
N_RATE: dict = {}
SOW_DATE: datetime = datetime(2026, 6, 18)
SIM_START: datetime = datetime(2026, 6, 8)
STA: dict = {}
DEFAULT_CULTIVAR: dict = {}
SLOC_PROFILE: list = []
SLOC_SCALE: float = 1.0
DEFAULT_NITRATE_SCALE: float = 0.24
CONTROL_PLOT: str = "A1"
OBS_ANTHESIS: pd.Timestamp | None = None
OBS_MATURITY: pd.Timestamp | None = None
SOIL_PREFIX = "SD"
SITE_DESC = "Dezhou"
VRNAME = "LENONG87"
CFG: TrialCfg | None = None

OUT = BASE / "out"
SIM = OUT / "sim"
WTH = OUT / "weather"
SOIL = OUT / "soil"
MGMT = OUT / "mgmt"
FIG = BASE / "fig"

# --------------------------------------------------------------------------
# scenarios
# --------------------------------------------------------------------------
# "base"    measured plant stand  -> out/<trial>/sim,          figures D1..D5
# "nominal" planting density      -> out/<trial>/sim_nominal,  figures N1..N5
# The figure prefixes deliberately avoid C/L (APSIM, analysis/out/fig).
SCENARIOS = {
    "base": dict(subdir="sim", prefix="D", density="measured",
                 mgmt="mgmt", label="实测留苗密度", short="实测密度"),
    "nominal": dict(subdir="sim_nominal", prefix="N", density="nominal",
                    mgmt="mgmt_nominal", label="额定播种密度", short="额定密度"),
}
DISPLAY = ""              # 试验中文名（set_trial 时绑定）
SCENARIO = "base"
PREFIX = "D"
SCEN_LABEL = "实测留苗密度"
SCEN_SHORT = "实测密度"
DENSITY_MODE = "measured"


def set_trial(name: str) -> TrialCfg:
    """Activate a trial: rebinds paths and every parameter used downstream."""
    global TRIAL, ORDER, N_RATE, SOW_DATE, SIM_START, STA, DEFAULT_CULTIVAR
    global SLOC_PROFILE, SLOC_SCALE, DEFAULT_NITRATE_SCALE, CONTROL_PLOT, OBS_ANTHESIS
    global OBS_MATURITY, SOIL_PREFIX, SITE_DESC, VRNAME, CFG, DISPLAY
    global OUT, SIM, WTH, SOIL, MGMT, FIG

    if name not in TRIALS:
        raise SystemExit(f"unknown trial {name!r}; choose from {list(TRIALS)}")
    cfg = TRIALS[name]
    CFG = cfg
    TRIAL = cfg.name
    ORDER = list(cfg.plots)
    N_RATE = dict(cfg.n_rate)
    SOW_DATE = cfg.sow
    SIM_START = cfg.sim_start
    STA = dict(cfg.sta)
    DEFAULT_CULTIVAR = dict(cfg.cultivar)
    SLOC_PROFILE = list(cfg.sloc)
    SLOC_SCALE = cfg.sloc_scale
    DEFAULT_NITRATE_SCALE = cfg.nitrate_scale
    CONTROL_PLOT = cfg.control
    OBS_ANTHESIS = cfg.obs_anthesis
    OBS_MATURITY = cfg.obs_maturity
    SOIL_PREFIX = cfg.soil_prefix
    SITE_DESC = cfg.description
    VRNAME = cfg.vrname
    DISPLAY = cfg.display or cfg.name

    OUT = BASE / "out" / TRIAL
    SIM, WTH, SOIL, MGMT = OUT / "sim", OUT / "weather", OUT / "soil", OUT / "mgmt"
    FIG = BASE / "fig" / TRIAL
    for d in (SIM, WTH, SOIL, MGMT, FIG):
        d.mkdir(parents=True, exist_ok=True)
    set_scenario("base")
    return cfg


def set_scenario(name: str) -> dict:
    """Switch between the measured-stand and the nominal-density run."""
    global SCENARIO, SIM, MGMT, PREFIX, SCEN_LABEL, SCEN_SHORT, DENSITY_MODE
    if name not in SCENARIOS:
        raise SystemExit(f"unknown scenario {name!r}; choose from {list(SCENARIOS)}")
    scen = SCENARIOS[name]
    SCENARIO = name
    SIM = OUT / scen["subdir"]
    MGMT = OUT / scen["mgmt"]          # keep the .MZX audit trail per scenario
    PREFIX = scen["prefix"]
    SCEN_LABEL = scen["label"]
    SCEN_SHORT = scen["short"]
    DENSITY_MODE = scen["density"]
    for d in (SIM, MGMT):
        d.mkdir(parents=True, exist_ok=True)
    return scen


set_trial("SD_SF")


# --------------------------------------------------------------------------
# observations
# --------------------------------------------------------------------------
def load_plots() -> pd.DataFrame:
    p = pd.read_csv(DATASET / "plots.csv")
    p = p[p["trial_id"] == TRIAL].set_index("plot")
    return p.reindex(ORDER)


def load_obs() -> pd.DataFrame:
    o = pd.read_csv(DATASET / "obs_plant.csv", parse_dates=["date"])
    o = o[o["trial_id"] == TRIAL].copy()
    p = load_plots()
    # no stand count (SD_MY73) -> fall back to the seed rate, same rule as PPOP
    dens = p["density_measured"].fillna(p["density_planted"]).reindex(o["plot"]).values
    o["bio_kg_ha"] = o["total_dw"] * dens / PLANTS_PER_MU * 10.0
    return o


def load_yield() -> pd.DataFrame:
    y = pd.read_csv(DATASET / "yield.csv")
    return y[y["trial_id"] == TRIAL].set_index("plot").reindex(ORDER)


def load_schedule() -> tuple[pd.DataFrame, pd.DataFrame]:
    fert = pd.read_csv(DATASET / "fert.csv", parse_dates=["date"])
    irrig = pd.read_csv(DATASET / "irrig.csv", parse_dates=["date"])
    fert = fert[(fert["trial_id"] == TRIAL) & fert["n_kg_ha"].notna()
                & (fert["n_kg_ha"] > 0)]
    irrig = irrig[(irrig["trial_id"] == TRIAL) & (irrig["mm"] > 0)]
    if fert.empty and CFG is not None and CFG.assumed_n:
        fert = pd.DataFrame(
            [dict(trial_id=TRIAL, plot=p, date=pd.Timestamp(d), n_kg_ha=k)
             for p in ORDER for d, k in CFG.assumed_n])
    if irrig.empty and CFG is not None and CFG.assumed_irrig:
        irrig = pd.DataFrame(
            [dict(trial_id=TRIAL, plot=p, date=pd.Timestamp(d), mm=m)
             for p in ORDER for d, m in CFG.assumed_irrig])
    return fert, irrig


# --------------------------------------------------------------------------
# DSSAT input builders
# --------------------------------------------------------------------------
def build_weather() -> Weather:
    """Daily weather 2026-06-08..2026-10-15 from the standardised table."""
    w = pd.read_csv(DATASET / "weather_daily.csv", parse_dates=["date"])
    w = w[w["zone"] == CFG.zone].sort_values("date").set_index("date")
    df = pd.DataFrame({
        "tmax": w["tmax"].astype(float),
        "tmin": w["tmin"].astype(float),
        "srad": w["radn"].astype(float),
        "rain": w["rain"].astype(float),
        "rhum": w["rh"].astype(float),
        "wind": w["wind"].astype(float),
    })
    df.index.name = "date"
    pars = {"tmax": "TMAX", "tmin": "TMIN", "srad": "SRAD",
            "rain": "RAIN", "rhum": "RHUM", "wind": "WIND"}
    wth = Weather(df, pars, STA["LAT"], STA["LON"], STA["ELEV"],
                  tav=STA["TAV"], amp=STA["AMP"], refht=STA["REFHT"],
                  wndht=STA["WNDHT"], co2=STA["CO2"])
    wth.INSI = STA["INSI"]
    # _name is fixed in __init__ from the default INSI ("WSTA"); rebuild it so
    # the two trials do not both write WSTA2601.WTH.
    wth._name = f'{wth.INSI}{wth.data.index[0].year % 100:02d}'                 f'{wth.data.index[-1].year - wth.data.index[0].year + 1:02d}'
    wth.description = STA["DESC"]
    return wth


def build_soil(plot: str, sloc_scale: float | None = None) -> SoilProfile:
    """Six-layer profile: measured DUL/SAT, Saxton-Rawls LL15, measured texture."""
    tab = pd.read_csv(DATASET / "soil_profile.csv")
    src_trial, src_plot = (CFG.soil_source if CFG and CFG.soil_source
                           else (TRIAL, plot))
    tab = tab[(tab["trial_id"] == src_trial) & (tab["plot"] == src_plot)]
    tab = tab.sort_values("layer")
    assert len(tab) == 6, f"{plot}: expected 6 layers, got {len(tab)}"
    if sloc_scale is None:
        sloc_scale = SLOC_SCALE

    soil = SoilProfile(default_class="SCL")
    # Two DSSAT soil-reader quirks (both abort with a misleading message):
    #   * the profile ID must be exactly 10 characters, alphanumeric, no "_";
    #   * the description on the *SOILS line is capped at 18 characters.
    soil.id = f"{SOIL_PREFIX}{plot:0<8}"
    soil.description = f"{SITE_DESC} {plot}"[:18]
    soil.site = SITE_DESC
    soil.country = "CN"
    soil.lat, soil.lon = STA["LAT"], STA["LON"]
    soil.csc_family = "SCL"
    soil.layers = {}

    bottom = tab["bottom_cm"].values
    top = np.concatenate([[0], bottom[:-1]])
    total = np.cumsum(tab["thickness_mm"].values)[-1]
    srgf = np.exp(-4.0 * (top / total))

    for i, row in tab.reset_index(drop=True).iterrows():
        pars = dict(
            SLMH="-99", SLLL=float(row["ll15"]), SDUL=float(row["dul"]),
            SSAT=float(row["sat"]), SRGF=float(round(srgf[i], 4)),
            SSKS=float(row["ks"]) / 240.0,         # mm/d -> cm/h
            SBDM=float(row["bd"]), SLOC=round(SLOC_PROFILE[i] * sloc_scale, 3),
            SLCL=float(row["clay"]), SLSI=float(row["sand"]),
            SLCF=float(max(0.0, 100 - row["clay"] - row["sand"])),
            SLHW=8.16, SCEC=20.0,
        )
        soil.add_layer(SoilLayer(int(round(row["bottom_cm"])), pars))

    soil.SCOM = "DBR"
    soil.SALB = 0.12
    soil.SLU1 = 6.0        # stage-1 evaporation limit (mm)
    soil.SLDR = 0.40       # drainage fraction per day
    soil.SLRO = 75.0       # SCS runoff curve number
    soil.SLNF = 1.0        # no additional N fertility limitation
    soil.SLPF = 1.0        # no additional P fertility limitation
    return soil


def _initial_conditions(plot: str, soil: SoilProfile, ic_date: datetime):
    """Initial water at DUL and the calibrated nitrate/ammonium profile."""
    tab = pd.read_csv(DATASET / "soil_profile.csv")
    src_trial, src_plot = (CFG.soil_source if CFG and CFG.soil_source
                           else (TRIAL, plot))
    tab = tab[(tab["trial_id"] == src_trial) & (tab["plot"] == src_plot)]
    tab = tab.sort_values("layer")
    rows = []
    prev = 0
    for i, row in enumerate(tab.reset_index(drop=True).itertuples()):
        depth = float(row.bottom_cm)
        th = depth - prev
        mass_kg_ha = 1e5 * th * float(row.bd)         # kg soil per ha
        sno3 = INIT_NO3_KG * N_LAYER_WEIGHTS[i] * 1e6 / mass_kg_ha
        snh4 = INIT_NH4_KG * N_LAYER_WEIGHTS[i] * 1e6 / mass_kg_ha
        rows.append(dict(ICBL=depth, SH2O=float(row.dul),
                         SNH4=round(snh4, 4), SNO3=round(sno3, 4)))
        prev = depth
    return pd.DataFrame(rows)


def build_management(plot: str, soil: SoilProfile, *, cultivar: dict | None = None,
                     nitrate_scale: float | None = None, nitrogen: bool = True,
                     n_total_kg: float | None = None,
                     density: float | None = None) -> Management:
    """One treatment: plot density, reported irrigation, reported urea N."""
    from DSSATTools.base.sections import TabularSubsection

    if nitrate_scale is None:
        nitrate_scale = DEFAULT_NITRATE_SCALE
    if density is None:                       # plants/m2
        row = load_plots().loc[plot]
        col = ("density_planted" if DENSITY_MODE == "nominal"
               else "density_measured")
        if pd.isna(row[col]):                 # no stand count: use the seed rate
            col = "density_planted"
        dens = float(row[col]) / PLANTS_PER_MU
    else:
        dens = float(density)
    man = Management(planting_date=SOW_DATE, sim_start=SIM_START,
                     irrigation="R", fertilization="R", harvest="M",
                     organic_matter="G")

    # --- planting -------------------------------------------------------
    man.planting_details["PDATE"] = SOW_DATE.strftime("%y%j")
    man.planting_details["PPOP"] = round(dens, 2)
    man.planting_details["PPOE"] = round(dens, 2)
    man.planting_details["PLRS"] = 60
    man.planting_details["PLDP"] = 4
    man.planting_details["PLME"] = "S"
    man.planting_details["PLDS"] = "R"

    # --- initial conditions --------------------------------------------
    ic = man.sim_start.strftime("%y%j")
    man.initial_conditions["ICDAT"] = ic
    man.initial_conditions["ICRT"] = 4
    man.initial_conditions["ICND"] = len(soil.layers)
    man.initial_conditions["ICWD"] = 180
    man.initial_conditions["ICRES"] = 0
    man.initial_conditions["ICREN"] = 0
    man.initial_conditions["ICREP"] = 0
    man.initial_conditions["ICRIP"] = 0
    man.initial_conditions["ICNAME"] = f"{TRIAL}_{plot}"
    tab = _initial_conditions(plot, soil, man.sim_start)
    if n_total_kg is not None:
        # 90/10 NO3/NH4 split, same shape as the calibrated profile
        tab["SNO3"] = (tab["SNO3"] * 0.9 * n_total_kg / INIT_NO3_KG).round(4)
        tab["SNH4"] = (tab["SNH4"] * 0.1 * n_total_kg / INIT_NH4_KG).round(4)
    else:
        tab["SNO3"] = (tab["SNO3"] * nitrate_scale).round(4)
    man.initial_conditions["table"] = TabularSubsection(
        tab[["ICBL", "SH2O", "SNH4", "SNO3"]])

    # --- irrigation ------------------------------------------------------
    _, irrig = load_schedule()
    irr = irrig[irrig["plot"] == plot].sort_values("date")
    if len(irr):
        t = pd.DataFrame({
            "IDATE": irr["date"].dt.strftime("%y%j"),
            "IROP": "IR001",
            "IRVAL": irr["mm"].astype(float).round(1),
        })
    else:  # DSSAT needs at least one row
        t = pd.DataFrame({"IDATE": [ic], "IROP": ["IR001"], "IRVAL": [0]})
    man.irrigation["table"] = TabularSubsection(t)
    man.simulation_controls["IRRIG"] = "R"

    # --- nitrogen (urea) --------------------------------------------------
    fert, _ = load_schedule()
    fe = fert[fert["plot"] == plot].sort_values("date")
    # A 0 N record dated on the planting day makes DSSAT stop right after the
    # management events (only MgmtOps/MgmtEvent/OVERVIEW are written).  The
    # zero-N control therefore gets an unused placeholder one day later and
    # FERTI='N'.
    placeholder_date = ((SOW_DATE + timedelta(days=1)) if fe.empty
                        else SOW_DATE).strftime("%y%j")
    # Integers matter: DSSATTools derives the *TREATMENTS factor levels with
    # min(1, max(...)).  A float 0.0 is then written as "0.0", which is wider
    # than the 3-character I3 slot and corrupts the whole factor record (DSSAT
    # silently stops after the management events, writing no PlantGro.OUT).
    rows = [dict(FDATE=placeholder_date, FMCD="FE005", FACD="AP001", FDEP=5,
                 FAMN=0, FAMP=0, FAMK=0, FAMC=0, FAMO=0,
                 FOCD="-99", FERNAME="NONE")]
    for _, r in fe.iterrows():
        rows.append(dict(FDATE=r["date"].strftime("%y%j"), FMCD="FE005",
                         FACD="AP001", FDEP=5, FAMN=round(float(r["n_kg_ha"]), 2),
                         FAMP=0.0, FAMK=0.0, FAMC=0.0, FAMO=0.0,
                         FOCD="-99", FERNAME="UREA"))
    man.fertilizers["table"] = TabularSubsection(pd.DataFrame(rows))
    man.simulation_controls["FERTI"] = "R" if len(fe) else "N"

    # --- simulation controls ---------------------------------------------
    man.simulation_controls["WATER"] = "Y"
    man.simulation_controls["NITRO"] = "Y" if nitrogen else "N"
    man.simulation_controls["PHOTO"] = "C"
    man.simulation_controls["HARVS"] = "M"
    man.simulation_controls["GROUT"] = "Y"
    man.simulation_controls["CAOUT"] = "Y"
    man.simulation_controls["WAOUT"] = "Y"
    man.simulation_controls["NIOUT"] = "N"
    man.simulation_controls["FMOPT"] = "A"
    return man


CULTIVAR_KEYS = ["P1", "P2", "P5", "G2", "G3", "PHINT"]
# IXIM-Maize (MZIXM048) adds two cultivar-level leaf parameters that CERES
# does not have: AX (cm2 of the largest leaf) and ALL (leaf longevity, GDD).
IXIM_KEYS = CULTIVAR_KEYS + ["AX", "ALL"]
MODEL = "CERES"
DEFAULT_AX = 650.0        # cm2/leaf, largest leaf
DEFAULT_ALL = 800.0       # GDD, leaf longevity
# DSSATTools hard-codes the maize model and a 6-parameter cultivar record;
# both have to be widened for the IXIM cultivar file (P1..PHINT + AX + ALL).
_IXIM_HDR = "A5,1X,1X,A16,1X,A5,1X,A6,8(1X,A5)"
# AX is F5.1 in the stock file (max 999.9); MY73 needs ~1180 cm2 so the two
# leaf fields are widened to F6.1.  DSSAT reads the wider record fine.
_IXIM_ROW = ("A6,1X,A16,1X,A5,1X,A6,1X,F5.1,1X,F5.3,2(1X,F5.0),"
             "2(1X,F5.2),2(1X,F6.1)")


def set_model(name: str) -> str:
    """CERES (default) or IXIM."""
    global MODEL
    MODEL = name.upper()
    assert MODEL in ("CERES", "IXIM"), MODEL
    return MODEL


def build_crop(cultivar: dict | None = None,
               ecotype: dict | None = None) -> Crop:
    if MODEL == "IXIM":
        from DSSATTools.base import sections as _sec
        import DSSATTools.crop as _dcrop
        keep = (_sec.CULTIVAR_HEADER_FMT["Maize"],
                _sec.CULTIVAR_ROWS_FMT["Maize"],
                _dcrop.CROPS_MODULES["Maize"])
        _sec.CULTIVAR_HEADER_FMT["Maize"] = _IXIM_HDR
        _sec.CULTIVAR_ROWS_FMT["Maize"] = _IXIM_ROW
        _dcrop.CROPS_MODULES["Maize"] = "MZIXM"
        try:
            crop = _dcrop.Crop("maize", "990002")
        finally:
            (_sec.CULTIVAR_HEADER_FMT["Maize"],
             _sec.CULTIVAR_ROWS_FMT["Maize"],
             _dcrop.CROPS_MODULES["Maize"]) = keep
        # Section._row_writer looks the format up *at call time*, so it has to
        # be rebound or AX/ALL are silently dropped when the CUL is written.
        import fortranformat as _ff
        crop.cultivar._row_writer = (
            lambda x: _ff.FortranRecordWriter(_IXIM_ROW).write(x))
    else:
        crop = Crop("maize", "990002")
    if _SPECIES_PATCH:
        crop._spe_path = str(PATCHED_SPECIES)
    pars = dict(DEFAULT_CULTIVAR)
    if MODEL == "IXIM":
        pars.setdefault("AX", DEFAULT_AX)
        pars.setdefault("ALL", DEFAULT_ALL)
    if cultivar:
        keys = IXIM_KEYS if MODEL == "IXIM" else CULTIVAR_KEYS
        pars.update({k: v for k, v in cultivar.items() if k in keys})
    for k, v in pars.items():
        if k in crop.cultivar:
            crop.cultivar[k] = float(v)
    crop.cultivar["VRNAME.........."] = VRNAME
    if ecotype and hasattr(crop, "ecotype"):
        for k, v in ecotype.items():
            if k in crop.ecotype:
                crop.ecotype[k] = float(v)
    return crop


# --------------------------------------------------------------------------
# runner
# --------------------------------------------------------------------------
PHENO_KEYS = ["ADAT", "MDAT", "EDAT", "PDAT"]


def _read_summary(path: Path) -> dict:
    if not path.exists():
        return {}
    import re
    lines = path.read_text(encoding="cp437", errors="replace").splitlines()
    head = next((i for i, l in enumerate(lines) if l.strip().startswith("@")), None)
    if head is None:
        return {}
    header, row = lines[head], lines[head + 1]
    # Fields are right aligned, so slice each value by the *end* column of its
    # header name.  Splitting on whitespace misaligns on TNAM (names with
    # spaces) and shifts every later column by one.
    out, prev = {}, 0
    for m in re.finditer(r"\S+", header):
        name, end = m.group().strip(". "), m.end()
        token = row[prev:end].strip()
        prev = end
        try:
            out[name] = float(token)
        except ValueError:
            out[name] = token
    return out


def run_plot(plot: str, cultivar: dict | None = None, *,
             nitrate_scale: float | None = None, nitrogen: bool = True,
             n_total_kg: float | None = None, soil_n_output: bool = False,
             sloc_scale: float | None = None, density: float | None = None,
             ecotype: dict | None = None,
             save: bool = False, keep_dir: Path | None = None,
             verbose: bool = False) -> tuple[pd.DataFrame, dict]:
    """Run one plot and return (daily PlantGro, run summary).

    When ``save`` is set the DSSAT inputs/outputs are copied to
    ``out/<trial>/mgmt`` and ``out/<trial>/sim``.
    """
    soil = build_soil(plot, sloc_scale=sloc_scale)
    wth = build_weather()
    crop = build_crop(cultivar, ecotype)
    man = build_management(plot, soil, cultivar=cultivar,
                          nitrate_scale=nitrate_scale, nitrogen=nitrogen,
                          n_total_kg=n_total_kg, density=density)
    if soil_n_output:
        man.simulation_controls["NIOUT"] = "Y"

    dssat = DSSAT()
    dssat.run(soil, wth, crop, man, verbose=verbose)
    run_dir = Path(dssat._RUN_PATH)
    g = dssat.output["PlantGro"].copy()
    summary = _read_summary(run_dir / "Summary.OUT")
    if soil_n_output and "SoilNi" in dssat.output:
        sn = dssat.output["SoilNi"].copy()
        sn = sn.reset_index().rename(columns={"index": "date", "@YEAR": "year", "DOY": "doy"})
        sn["plot"] = plot
        sn["date"] = pd.to_datetime(sn["date"])
        summary = dict(summary, soil_ni=sn)

    if keep_dir is not None:
        # keep the *whole* run directory: .MZX + genotype files + SOIL.SOL +
        # Weather + DSSATPRO.V48, i.e. everything a third party needs to rerun
        dest = Path(keep_dir)
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(run_dir, dest)
    if save:
        keep_out = {"PlantGro.OUT", "Summary.OUT", "SoilWat.OUT",
                    "OVERVIEW.OUT", "WARNING.OUT"}
        for f in run_dir.iterdir():
            if f.suffix.upper() == ".MZX":
                shutil.copy2(f, MGMT / f"{TRIAL}_{plot}{f.suffix}")
            elif f.name in keep_out:
                shutil.copy2(f, SIM / f"{TRIAL}_{plot}_{f.name}")
    dssat.close()

    g = g.reset_index().rename(columns={"index": "date", "@YEAR": "year", "DOY": "doy"})
    g["plot"] = plot
    g["date"] = pd.to_datetime(g["date"])
    return g, summary


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def dssat_date(yyddd: float, default_year: int = 2026) -> pd.Timestamp | None:
    """DSSAT date (YYDDD, e.g. 26221) -> Timestamp."""
    try:
        v = int(round(yyddd))
    except (TypeError, ValueError):
        return None
    if v <= 0 or v in (-99, 999999):
        return None
    # Summary.OUT uses YYYYDDD, the management files use YYDDD.
    if v >= 1_000_000:
        year, ddd = divmod(v, 1000)
    else:
        yy, ddd = divmod(v, 1000)
        year = default_year // 100 * 100 + yy
    if not 1 <= ddd <= 366:
        return None
    try:
        return pd.Timestamp(datetime(year, 1, 1) + timedelta(days=ddd - 1))
    except ValueError:
        return None


def plateau_metrics(dates: pd.Series, lai: pd.Series, thresh: float = 0.95):
    """Peak LAI, the >=thresh*peak window, and the change across that window."""
    lai = pd.Series(np.asarray(lai, dtype=float), index=pd.to_datetime(dates))
    lai = lai.dropna()
    if lai.empty or lai.max() <= 0:
        return dict(peak=np.nan, peak_date=pd.NaT, days=0, change_pct=np.nan,
                    n_points=0, span_days=0, start=pd.NaT, end=pd.NaT)
    peak = float(lai.max())
    peak_date = lai.idxmax()
    win = lai[lai >= thresh * peak]
    start, end = win.index.min(), win.index.max()
    return dict(
        peak=peak, peak_date=peak_date,
        n_points=int(len(win)),
        span_days=int((end - start).days + 1),
        start=start, end=end,
        change_pct=(win.iloc[-1] / win.iloc[0] - 1) * 100,
    )


# --------------------------------------------------------------------------
# species-file patching (diagnostics only)
# --------------------------------------------------------------------------
SPECIES_FILE = Path(_DSSAT_TMP) / "Genotype" / "MZCER048.SPE"
# Crop.write() copies the *package* static file into the run directory, so a
# patched copy has to be handed to the Crop object through _spe_path.
PATCHED_SPECIES = Path(_DSSAT_TMP) / "MZCER048_patched.SPE"
_SPECIES_PATCH: dict = {}


def _species_text(overrides: dict) -> str:
    text = SPECIES_FILE.read_text(encoding="cp437")
    for key, val in overrides.items():
        pat = re.compile(rf"^(\s*{key}\s+)([-+0-9.]+)", re.M)
        if not pat.search(text):
            raise KeyError(f"species coefficient {key} not found")
        text = pat.sub(lambda m: f"{m.group(1)}{val}", text, count=1)
    return text


@contextmanager
def patch_species(**overrides):
    """Temporarily overwrite species coefficients for the next runs."""
    _SPECIES_PATCH.clear()
    _SPECIES_PATCH.update(overrides)
    PATCHED_SPECIES.write_text(_species_text(overrides), encoding="cp437")
    try:
        yield
    finally:
        _SPECIES_PATCH.clear()


# --------------------------------------------------------------------------
# experiment metrics
# --------------------------------------------------------------------------
def plot_metrics(g: pd.DataFrame, s: dict, density: float) -> dict:
    """Season metrics used by the diagnostic experiments (plants/m2 given)."""
    g = g.sort_values("date")
    adat = dssat_date(s.get("ADAT"))
    mat = g.iloc[-1]
    lai_peak = float(g["LAID"].max())
    bio_ant = np.nan
    if adat is not None:
        near = g[g["date"] <= adat]
        if len(near):
            bio_ant = float(near["CWAD"].iloc[-1])
    out = dict(
        density=density,
        lai_peak=lai_peak,
        la_per_plant=lai_peak / density,
        bio_anthesis=bio_ant,
        bio_final=float(mat["CWAD"]),
        yield_kg=float(mat["GWAD"]),
        hi=float(mat["HIAD"]),
        grains_m2=float(mat["G#AD"]),
        grain_w_mg=float(mat["GWGD"]),   # DSSAT GWGD is already mg/grain
        anthesis=adat,
        maturity=dssat_date(s.get("MDAT")),
    )
    out["grains_per_ear"] = (out["grains_m2"] / density
                             if density else np.nan)
    return out
