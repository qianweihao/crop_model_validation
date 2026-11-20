"""Run every DSSAT experiment file in this folder and collect the outputs.

    python run.py [--exe <dscsm048.exe>] [--keep]

Needs only the DSSAT 4.8 model executable; the support files are bundled in
_support/.  DSSAT uses fixed-length path buffers, so keep this folder at a
short path (e.g. C:\\dssat) - run.py warns if the path is too long.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAXPATH = 60          # DSSAT truncates longer paths silently


def find_exe(given: str | None) -> Path:
    if given:
        return Path(given).expanduser().resolve()
    env = os.environ.get("DSCSM048")
    if env:
        return Path(env)
    for c in (HERE / "_support" / "bin" / "dscsm048.exe",
              HERE / "_support" / "bin" / "dscsm048",
              Path(r"C:/DSSAT48/DSCSM048.EXE"),
              Path(r"C:/DSSAT48/bin/dscsm048.exe"),
              Path(r"C:/DSSAT48/Tools/dscsm048.exe")):
        if c.exists():
            return c
    found = shutil.which("dscsm048") or shutil.which("dscsm048.exe")
    if found:
        return Path(found)
    sys.exit("cannot find dscsm048; pass --exe <path to dscsm048.exe>")


def write_config() -> None:
    (HERE / "DSSATPRO.V48").write_text(
        f"WED    {HERE / 'Weather'}\n"
        f"MMZ    {HERE} dscsm048 MZCER048\n"
        f"CRD    {HERE}\n"
        f"PSD    {HERE / '_support' / 'Pest'}\n"
        f"SLD    {HERE / '_support' / 'Soil'}\n"
        f"STD    {HERE / '_support' / 'StandardData'}\n",
        encoding="ascii")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=None, help="path to dscsm048.exe")
    ap.add_argument("--keep", action="store_true",
                    help="keep the raw *.OUT files next to the .MZX")
    args = ap.parse_args()

    if len(str(HERE)) > MAXPATH:
        print(f"WARNING: {HERE} is {len(str(HERE))} characters long; "
              f"DSSAT silently truncates paths above ~{MAXPATH}. "
              f"Copy this folder to e.g. C:\\dssat and run it there.")
    exe = find_exe(args.exe)
    print(f"model : {exe}")
    print(f"folder: {HERE}")
    write_config()

    results = HERE / "results"
    results.mkdir(exist_ok=True)
    files = sorted(HERE.glob("*.MZX"))
    if not files:
        sys.exit("no .MZX files found")

    ok = 0
    for mzx in files:
        tag = mzx.stem
        for f in HERE.glob("*.OUT"):
            f.unlink()
        r = subprocess.run([str(exe), "C", mzx.name, "1"], cwd=str(HERE),
                           capture_output=True, text=True)
        summary = HERE / "Summary.OUT"
        if r.returncode == 0 and summary.exists():
            dest = results / tag
            dest.mkdir(exist_ok=True)
            for f in HERE.glob("*.OUT"):
                shutil.copy2(f, dest / f.name)
            if not args.keep:
                for f in HERE.glob("*.OUT"):
                    f.unlink()
            ok += 1
            print(f"  {tag:16s} ok   -> results/{tag}/")
        else:
            print(f"  {tag:16s} FAILED (rc={r.returncode})")
            print((HERE / "ERROR.OUT").read_text(errors="replace")[-400:]
                  if (HERE / "ERROR.OUT").exists() else r.stdout[-400:])
    print(f"\n{ok}/{len(files)} runs finished; results in {results}")


if __name__ == "__main__":
    main()
