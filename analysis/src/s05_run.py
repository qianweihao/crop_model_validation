"""Step 05 - run every .apsimx with Models.exe and collect the Report tables.

Writes a tidy long table (one row per simulation x day x variable) plus a wide
per-simulation table, so downstream steps never touch SQLite again.
"""
from __future__ import annotations

import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_env import resolve  # noqa: E402
from common import APSIM_DIR, SIM  # noqa: E402

INSTALL = resolve()
MODELS_EXE = INSTALL.models_exe


def run_all(force: bool = False) -> None:
    SIM.mkdir(parents=True, exist_ok=True)
    for f in sorted(APSIM_DIR.glob("*.apsimx")):
        dst = SIM / f.name
        if force or not dst.exists():
            shutil.copy2(f, dst)
    todo = [f for f in sorted(SIM.glob("*.apsimx"))
            if force or not f.with_suffix(".db").exists()]
    if not todo:
        print("all simulations already have a .db")
        return
    print(f"running {len(todo)} simulations with APSIM {INSTALL.version} ...")
    for f in todo:
        p = subprocess.run([str(MODELS_EXE), f.name], cwd=SIM,
                           capture_output=True, text=True)
        if "ERROR" in (p.stdout + p.stderr):
            print(f"  FAILED {f.stem}\n{(p.stdout + p.stderr)[:800]}")
    print("done")


def collect() -> pd.DataFrame:
    rows = []
    for db in sorted(SIM.glob("*.db")):
        con = sqlite3.connect(db)
        try:
            df = pd.read_sql("select * from Report", con)
        except Exception:
            con.close()
            continue
        con.close()
        df = df.drop(columns=[c for c in ("CheckpointID", "SimulationID", "Zone")
                              if c in df.columns])
        df.insert(0, "sim", db.stem)
        rows.append(df)
    if not rows:
        raise SystemExit("no Report tables found - run step 05 first")
    wide = pd.concat(rows, ignore_index=True)
    wide["Clock.Today"] = pd.to_datetime(wide["Clock.Today"])
    wide.to_csv(SIM / "sim_daily_wide.csv", index=False)

    id_cols = ["sim", "Clock.Today"]
    long = wide.melt(id_vars=id_cols, var_name="variable", value_name="value")
    long.to_csv(SIM / "sim_daily_long.csv", index=False)
    return wide


def main() -> None:
    run_all(force="--force" in sys.argv)
    w = collect()
    print(f"wide: {w.shape}  sims={w['sim'].nunique()}  days={w['Clock.Today'].nunique()}")
    print("columns:", [c for c in w.columns][:12], "...")

    # quick headline per simulation
    g = w.groupby("sim").agg(
        sim_yield_kg_ha=("Maize.Grain.Total.Wt", lambda s: s.max() * 10.0),
        sim_biomass_kg_ha=("Maize.AboveGround.Wt", lambda s: s.max() * 10.0),
        sim_lai_max=("Maize.Leaf.LAI", "max"),
        min_FN=("Maize.Leaf.Photosynthesis.FN", "min"),
        min_FW=("Maize.Leaf.Photosynthesis.FW", "min"),
        flowering=("Maize.Phenology.CurrentStageName",
                   lambda s: s.index[s.astype(str).str.contains("Flowering")].min()
                   if s.astype(str).str.contains("Flowering").any() else None),
    ).round(3)
    g.to_csv(SIM / "sim_summary.csv")
    print(g.to_string())


if __name__ == "__main__":
    main()
