"""Run the whole pipeline in order.

    python analysis/src/run_all.py            # data -> inputs -> sims -> diagnosis
    python analysis/src/run_all.py --calibrate  # ... plus cultivar calibration

Every step is idempotent; APSIM is auto-detected (see apsim_env.py), so after an
APSIM upgrade nothing needs editing - just rerun.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

SRC = Path(__file__).resolve().parent
STEPS = [
    ("s01_dataset.py", []),
    ("s02_weather.py", []),
    ("s03_soil.py", []),
    ("s04_apsim.py", []),
    ("s05_run.py", ["--force"]),
    ("s06_diagnose.py", []),
]


def call(script: str, args: list[str]) -> None:
    t0 = time.time()
    print(f"\n=== {script} {' '.join(args)}")
    p = subprocess.run([sys.executable, str(SRC / script), *args],
                       cwd=SRC.parents[1])
    if p.returncode != 0:
        raise SystemExit(f"{script} failed (exit {p.returncode})")
    print(f"--- {script} ok in {time.time() - t0:.1f}s")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--from", dest="start", default=None,
                    help="resume from this step, e.g. s05_run.py")
    a = ap.parse_args()

    steps = list(STEPS)
    if a.calibrate:
        steps.append(("s07_calibrate.py", ["--stage", "both"]))
    if a.start:
        idx = next((i for i, s in enumerate(steps) if s[0] == a.start), None)
        if idx is None:
            raise SystemExit(f"unknown step {a.start}")
        steps = steps[idx:]

    print(f"APSIM: {' '.join(str(x) for x in __import__('apsim_env').resolve().__dict__.values())[:80]}")
    for script, args in steps:
        call(script, args)
    print("\nall steps finished")


if __name__ == "__main__":
    main()
