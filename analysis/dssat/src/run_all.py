"""Run the whole DSSAT pipeline for both trials.

    python analysis/dssat/src/run_all.py            # s01,s02,s04,s05,s06
    python analysis/dssat/src/run_all.py --calibrate  # also re-run s03 (~10 min)

Each step is a separate process so a failure in one trial cannot leave the other
in a half-written state.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent
TRIALS = ["SD_SF", "XX_SF"]


def run(script: str, trial: str) -> None:
    cmd = [sys.executable, str(SRC / script), "--trial", trial]
    print(f"\n=== {script} --trial {trial} ===", flush=True)
    r = subprocess.run(cmd, cwd=SRC.parents[2])
    if r.returncode != 0:
        raise SystemExit(f"{script} failed for {trial}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true",
                    help="also re-run s03_calibrate (slow)")
    args = ap.parse_args()

    steps = ["s01_inputs.py", "s02_run.py"]
    if args.calibrate:
        steps.append("s03_calibrate.py")
    steps += ["s04_validate.py", "s05_figures.py", "s06_nstress.py"]

    for trial in TRIALS:
        for script in steps:
            run(script, trial)
    print("\nall trials done")


if __name__ == "__main__":
    main()
