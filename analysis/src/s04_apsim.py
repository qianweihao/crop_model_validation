"""Step 04 - assemble one APSIM NG .apsimx per field plot from the Maize template.

Everything version-critical is done as plain JSON editing of the APSIM 2025.3
`Maize.apsimx` template, so the pipeline does not depend on apsimNGpy (which
rejects APSIM 2025.3 - it requires the 2026.2 build that ships APSIM.Core.dll).

Managers:
  Sowing    -> the template's own 'Sow using a variable rule', windowed to a
               single day with MinESW = MinRain = 0
  Nutrition -> a new manager that applies urea on the recorded dates (kg N/ha)
  Water     -> the same manager applies irrigation on the recorded dates (mm)
  Harvest   -> the template's own Harvest manager
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_env import resolve  # noqa: E402
from common import APSIM_DIR, DS, MET, SOIL  # noqa: E402

INSTALL = resolve()
MODELS_EXE = INSTALL.models_exe
TPL_DIR = Path(__file__).resolve().parents[1] / "out" / "template"
TPL_DIR.mkdir(parents=True, exist_ok=True)
# the template must come from the *running* build: its "Version" field gates
# whether Models.exe will read the file at all
TEMPLATE = TPL_DIR / f"Maize_{INSTALL.version}.apsimx"
if not TEMPLATE.exists():
    src = INSTALL.examples / "Maize.apsimx"
    if not src.exists():
        raise FileNotFoundError(f"no Maize.apsimx example under {INSTALL.examples}")
    shutil.copy2(src, TEMPLATE)

CLOCK_START = "2026-06-08T00:00:00"
CLOCK_END = "2026-10-15T00:00:00"
SOW_DEPTH = 40.0
ROW_SPACING = 600.0
MU_TO_M2 = 666.667

# variety -> (APSIM base cultivar, parameter commands applied on top)
CULTIVARS = {
    "乐农87": ("Dekalb_XL82", {}),
    "MY73": ("Dekalb_XL82", {}),
    "金华创826": ("Dekalb_XL82", {}),
}

REPORT_VARS = [
    "[Clock].Today",
    "[Maize].Phenology.CurrentStageName",
    "[Maize].Phenology.ThermalTime",
    "[Maize].Phenology.Stage",
    "[Maize].CoverGreen",
    "[Maize].Structure.LeafTipsAppeared",
    "[Maize].Leaf.ExpandedCohortNo",
    "[Maize].Leaf.LAI",
    "[Maize].Leaf.Wt",
    "[Maize].Stem.Wt",
    "[Maize].Root.Wt",
    "[Maize].AboveGround.Wt",
    "[Maize].AboveGround.N",
    "[Maize].Leaf.NConc",
    "[Maize].Grain.Total.Wt",
    "[Maize].Grain.Number",
    "[Maize].Grain.N",
    "[Maize].Leaf.Photosynthesis.FN",
    "[Maize].Leaf.Photosynthesis.FW",
    "[Maize].Leaf.Photosynthesis.RUE",
    "[Soil].SoilWater.PAW",
    "[Soil].SoilWater.ESW",
    "[Soil].SoilWater.SWmm",
    "[Soil].Nutrient.NO3.kgha",
    "[Irrigation].IrrigationApplied",
    "[Weather].Rain",
    "[Weather].Radn",
    "[Weather].MeanT",
]

SCHEDULE_CODE = [
    "using Models.Core;",
    "using Models.PMF;",
    "using Models.Soils;",
    "using System;",
    "",
    "namespace Models",
    "{",
    "    [Serializable]",
    "    public class Script : Model",
    "    {",
    "        [Link] private Clock Clock;",
    "        [Link] private Fertiliser Fertiliser;",
    "        [Link] private Irrigation Irrigation;",
    "",
    "        [Description(\"Fertiliser schedule 'yyyy-MM-dd:N_kg_ha;...'\")]",
    "        public string FertSchedule { get; set; }",
    "",
    "        [Description(\"Irrigation schedule 'yyyy-MM-dd:mm;...'\")]",
    "        public string IrrigationSchedule { get; set; }",
    "",
    "        private static DateTime D(string s)",
    "        {",
    "            return DateTime.ParseExact(s.Trim(), \"yyyy-MM-dd\", null);",
    "        }",
    "",
    "        [EventSubscribe(\"DoManagement\")]",
    "        private void OnDoManagement(object sender, EventArgs e)",
    "        {",
    "            if (!string.IsNullOrEmpty(IrrigationSchedule))",
    "            {",
    "                foreach (string item in IrrigationSchedule.Split(';'))",
    "                {",
    "                    if (item.Trim().Length == 0) continue;",
    "                    string[] kv = item.Split(':');",
    "                    if (D(kv[0]) == Clock.Today.Date)",
    "                        Irrigation.Apply(Convert.ToDouble(kv[1]));",
    "                }",
    "            }",
    "            if (!string.IsNullOrEmpty(FertSchedule))",
    "            {",
    "                foreach (string item in FertSchedule.Split(';'))",
    "                {",
    "                    if (item.Trim().Length == 0) continue;",
    "                    string[] kv = item.Split(':');",
    "                    if (D(kv[0]) == Clock.Today.Date)",
    "                        Fertiliser.Apply(Convert.ToDouble(kv[1]), Fertiliser.Types.UreaN);",
    "                }",
    "            }",
    "        }",
    "    }",
    "}",
]


def find(nodes, name):
    for n in nodes:
        if n.get("Name") == name:
            return n
    raise KeyError(name)


def make_manager(name: str, code: list[str], params: dict) -> dict:
    return {
        "$type": "Models.Manager, Models",
        "CodeArray": code,
        "Parameters": [{"Key": k, "Value": str(v)} for k, v in params.items()],
        "Name": name, "ResourceName": None, "Children": [],
        "Enabled": True, "ReadOnly": False,
    }


def build_sim(plot_row: pd.Series, fert: pd.DataFrame, irrig: pd.DataFrame) -> dict:
    tid, plot = plot_row["trial_id"], plot_row["plot"]

    tmpl = json.loads(TEMPLATE.read_text(encoding="utf-8-sig"))
    sim = tmpl["Children"][0]
    sim["Name"] = f"{tid}_{plot}"
    field = find(sim["Children"], "Field")

    # ---------- clock + weather -----------------------------------------
    find(sim["Children"], "Clock").update(Start=CLOCK_START, End=CLOCK_END)
    find(sim["Children"], "Weather")["FileName"] = str(MET / f"{plot_row['weather']}.met")

    # ---------- soil -----------------------------------------------------
    soil = json.loads((SOIL / f"{tid}_{plot}.json").read_text(encoding="utf-8"))
    idx = next(i for i, c in enumerate(field["Children"]) if c.get("Name") == "Soil")
    field["Children"][idx] = soil

    # ---------- irrigation node -----------------------------------------
    if not any(c.get("$type", "").startswith("Models.Irrigation") for c in field["Children"]):
        field["Children"].append({
            "$type": "Models.Irrigation, Models", "Name": "Irrigation",
            "ResourceName": None, "Children": [], "Enabled": True, "ReadOnly": False})

    # ---------- cultivar ------------------------------------------------
    variety = plot_row["variety"]
    base_cv, commands = CULTIVARS.get(variety, ("Dekalb_XL82", {}))
    cv_name = f"CN_{variety}" if commands else base_cv
    maize = find(field["Children"], "Maize")
    if commands:
        maize["Children"] = [{
            "$type": "Models.Core.Folder, Models", "Name": "Cultivars",
            "Children": [{
                "$type": "Models.PMF.Cultivar, Models",
                "Command": [f"{k} = {v}" for k, v in commands.items()],
                "Name": cv_name, "Children": [], "Enabled": True, "ReadOnly": False,
            }],
            "Enabled": True, "ReadOnly": False,
        }]
    else:
        maize["Children"] = []

    # ---------- management ----------------------------------------------
    pop = float(plot_row["density_measured"]) if pd.notna(plot_row["density_measured"]) \
        else float(plot_row["density_planted"])
    sow_date = pd.Timestamp(plot_row["sow_date"])
    dstr = f"{sow_date.day}-{sow_date.strftime('%b').lower()}"

    sow = find(field["Children"], "Sow using a variable rule")
    for p in sow["Parameters"]:
        if p["Key"] == "StartDate" or p["Key"] == "EndDate":
            p["Value"] = dstr
        elif p["Key"] == "MinESW":
            p["Value"] = "0"
        elif p["Key"] == "MinRain":
            p["Value"] = "-1"
        elif p["Key"] == "RainDays":
            p["Value"] = "1"
        elif p["Key"] == "CultivarName":
            p["Value"] = cv_name
        elif p["Key"] == "Population":
            p["Value"] = f"{pop / MU_TO_M2:.4f}"
        elif p["Key"] == "SowingDepth":
            p["Value"] = str(SOW_DEPTH)
        elif p["Key"] == "RowSpacing":
            p["Value"] = str(ROW_SPACING)

    # replace the template's at-sowing fertiliser rule with a schedule manager
    field["Children"] = [c for c in field["Children"]
                         if c.get("Name") != "Fertilise at sowing"]

    f = fert[(fert.trial_id == tid) & (fert["plot"] == plot) & (fert.n_kg_ha > 0)]
    w = irrig[(irrig.trial_id == tid) & (irrig["plot"] == plot) & (irrig.mm > 0)]
    fert_sched = ";".join(f"{pd.Timestamp(r['date']).strftime('%Y-%m-%d')}:{r['n_kg_ha']:.3f}"
                          for _, r in f.sort_values("date").iterrows())
    irr_sched = ";".join(f"{pd.Timestamp(r['date']).strftime('%Y-%m-%d')}:{r['mm']:.2f}"
                         for _, r in w.sort_values("date").iterrows())
    field["Children"].append(make_manager(
        "Nitrogen and irrigation schedule", SCHEDULE_CODE,
        {"FertSchedule": fert_sched, "IrrigationSchedule": irr_sched}))

    # ---------- remove the harvest rule ----------------------------------
    # The template's Harvest manager calls Crop.Harvest() + Crop.EndCrop() as
    # soon as the crop is ready, which wipes AboveGround.Wt and LAI to zero on
    # that single day.  That produces a vertical cliff at the end of every LAI
    # and biomass curve.  We keep the crop standing so the curves stay
    # continuous and the end-of-season values remain observable.
    field["Children"] = [c for c in field["Children"]
                         if not str(c.get("Name", "")).startswith("Harvest")]

    # ---------- report ---------------------------------------------------
    "".join(SCHEDULE_CODE)
    rep = find(field["Children"], "Report")
    rep["VariableNames"] = REPORT_VARS
    rep["EventNames"] = ["[Clock].EndOfDay"]
    rep["Enabled"] = True

    # ---------- drop plotting nodes -------------------------------------
    sim["Children"] = [c for c in sim["Children"]
                       if not c.get("$type", "").startswith("Models.Graph")]
    return sim


def main() -> None:
    plots = pd.read_csv(DS / "plots.csv", parse_dates=["sow_date"])
    fert = pd.read_csv(DS / "fert.csv", parse_dates=["date"])
    irrig = pd.read_csv(DS / "irrig.csv", parse_dates=["date"])

    simulable = plots[plots["trial_id"].isin(["SD_SF", "XX_SF", "XX_IRR"])]
    made = []
    for _, r in simulable.iterrows():
        sim = build_sim(r, fert, irrig)
        out = APSIM_DIR / f"{r['trial_id']}_{r['plot']}.apsimx"
        doc = json.loads(TEMPLATE.read_text(encoding="utf-8-sig"))
        doc["Children"] = [sim] + [c for c in doc["Children"]
                                   if not c.get("$type", "").startswith("Models.Core.Simulation")]
        out.write_text(json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")
        made.append(out)
    print(f"APSIM {INSTALL.version}  template={TEMPLATE.name}")
    print(f"wrote {len(made)} .apsimx to {APSIM_DIR}")
    print("simulable plots:", simulable.groupby("trial_id").size().to_dict())
    print("skipped (no management/soil records):",
          sorted(set(plots.trial_id) - set(simulable.trial_id)))


if __name__ == "__main__":
    main()
