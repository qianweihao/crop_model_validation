"""Programmatic parameter override + single-run helpers for APSIM NG.

Verified on APSIM 2025.3.7676.0 (see analysis/notes_apsim_overrides.md):

  * a custom cultivar is injected by putting a ``Models.Core.Folder`` named
    "Cultivars" under the plant, holding ``Models.PMF.Cultivar`` nodes whose
    ``Command`` list carries the overrides.  (``Models.PMF.CultivarFolder`` does
    not exist in 2025.3.)
  * function-valued parameters are written with a ``.FixedValue`` suffix,
    e.g. ``[Maize].Leaf.Photosynthesis.RUE.FixedValue = 2.4``.  Writing the bare
    property raises "Target is read only" / "Unable to convert ... to Constant".
"""
from __future__ import annotations

import json
import os
import sqlite3
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apsim_env import resolve  # noqa: E402

INSTALL = resolve()
MODELS_EXE = INSTALL.models_exe


def _max_workers() -> int:
    """Keep the machine responsive: APSIM_BIN_WORKERS or 1/3 of the CPU count, max 4."""
    try:
        return max(1, min(4, int(os.environ.get("APSIM_WORKERS", 0))
                           or max(1, (os.cpu_count() or 4) // 3)))
    except ValueError:
        return 2


def _run_low_priority(cmd: list[str], cwd: Path) -> tuple[str, str]:
    """Run a command with BELOW_NORMAL priority so the desktop stays responsive."""
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    if os.name == "nt":
        try:
            import psutil
            psutil.Process(proc.pid).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        except Exception:
            pass
    out, err = proc.communicate()
    return out, err


@dataclass
class SimCase:
    """One .apsimx document plus the cultivar overrides to apply."""
    doc: dict
    plant: str = "Maize"
    commands: dict[str, float] = field(default_factory=dict)

    @classmethod
    def load(cls, path: str | Path, plant: str = "Maize") -> "SimCase":
        doc = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        return cls(doc=doc, plant=plant)

    # -- structure helpers -------------------------------------------------
    def _field(self) -> dict:
        sim = next(c for c in self.doc["Children"]
                   if c.get("$type", "").startswith("Models.Core.Simulation"))
        return next(c for c in sim["Children"] if c.get("Name") == "Field")

    def _plant(self) -> dict:
        return next(c for c in self._field()["Children"] if c.get("Name") == self.plant)

    def set_commands(self, commands: dict[str, float]) -> "SimCase":
        self.commands = dict(commands)
        return self

    def write(self, path: str | Path) -> Path:
        doc = json.loads(json.dumps(self.doc))
        self.doc = doc
        if self.commands:
            plant = self._plant()
            plant["Children"] = [{
                "$type": "Models.Core.Folder, Models", "Name": "Cultivars",
                "Children": [{
                    "$type": "Models.PMF.Cultivar, Models",
                    "Command": [f"[{self.plant}].{k.lstrip('.')} = {v}"
                                if not k.startswith("[") else f"{k} = {v}"
                                for k, v in self.commands.items()],
                    "Name": "CAL", "Children": [], "Enabled": True, "ReadOnly": False,
                }],
                "Enabled": True, "ReadOnly": False,
            }]
            sow = next(c for c in self._field()["Children"]
                       if c.get("Name", "").startswith("Sow"))
            for p in sow["Parameters"]:
                if p["Key"] == "CultivarName":
                    p["Value"] = "CAL"
        path = Path(path)
        path.write_text(json.dumps(self.doc, ensure_ascii=False, indent=1), encoding="utf-8")
        return path


def run(case: SimCase, workdir: str | Path, tag: str = "run",
        keep_db: bool = False) -> pd.DataFrame:
    """Run one case and return its Report table."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    f = case.write(workdir / f"{tag}.apsimx")
    db = workdir / f"{tag}.db"
    if db.exists():
        db.unlink()
    out, err = _run_low_priority([str(MODELS_EXE), f.name], workdir)
    out += err
    if "ERROR" in out or not db.exists():
        raise RuntimeError(f"APSIM failed for {tag}:\n" + "\n".join(out.splitlines()[:12]))
    con = sqlite3.connect(db)
    try:
        df = pd.read_sql("select * from Report", con)
    finally:
        con.close()
    df["Clock.Today"] = pd.to_datetime(df["Clock.Today"])
    if not keep_db:
        db.unlink(missing_ok=True)
    return df


def run_many(jobs: list[tuple[SimCase, str]], workdir: str | Path,
             keep_db: bool = False) -> dict[str, pd.DataFrame]:
    """Run several cases concurrently, capped by _max_workers() and low priority.

    jobs = [(case, tag), ...]  ->  {tag: Report table}
    """
    from concurrent.futures import ThreadPoolExecutor
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out: dict[str, pd.DataFrame | None] = {}
    with ThreadPoolExecutor(max_workers=_max_workers()) as ex:
        futs = {ex.submit(run, c, workdir, t, keep_db): t for c, t in jobs}
        for fut, tag in futs.items():
            try:
                out[tag] = fut.result()
            except Exception:
                out[tag] = None      # caller decides how to score a failed run
    return out


def stage_dates(df: pd.DataFrame) -> dict:
    st = df["Maize.Phenology.CurrentStageName"].astype(str)
    out = {}
    for name in ("Emergence", "Flowering", "Maturity"):
        hit = df.loc[st.eq(name), "Clock.Today"]
        out[name] = hit.iloc[0] if len(hit) else None
    out["yield_kg_ha"] = float(df["Maize.Grain.Total.Wt"].max() * 10.0)
    out["biomass_kg_ha"] = float(df["Maize.AboveGround.Wt"].max() * 10.0)
    out["lai_max"] = float(df["Maize.Leaf.LAI"].max())
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="apply one parameter set and report")
    ap.add_argument("apsimx")
    ap.add_argument("--set", nargs="*", default=[], metavar="PATH=VALUE")
    a = ap.parse_args()
    case = SimCase.load(a.apsimx)
    cmds = {}
    for kv in a.set:
        k, v = kv.rsplit("=", 1)
        cmds[k] = float(v)
    with tempfile.TemporaryDirectory() as td:
        df = run(case.set_commands(cmds), td, "probe")
    print(a.apsimx, cmds)
    print(stage_dates(df))
