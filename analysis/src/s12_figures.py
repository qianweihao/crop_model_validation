"""Step 12 - the L series, drawn with the same four functions as the C series.

Kept as a thin entry point so ``python analysis/src/s12_figures.py`` still works
and produces figures identical to the C set but under the L naming.

  L1  LAI              nine per-plot panels, observed vs simulated
  L2  地上部生物量       nine per-plot panels, observed vs simulated
  L3  产量              per-plot bars with bias, plus a 1:1 scatter
  L4  氮响应             yield / biomass / HI against nitrogen rate

Usage:  python analysis/src/s12_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s20_conservative_figures import load, make_all  # noqa: E402


def main() -> None:
    sim, _obs, val = load()
    print(f"sim rows={len(sim)}  plots={sim['sim'].nunique()}")
    for p in make_all("L"):
        print(" ", p)
    o, s = val["obs_yield_dry"], val["sim_yield"]
    print(f"\n产量 RRMSE {((o - s) ** 2).mean() ** 0.5 / o.mean() * 100:.2f}%")


if __name__ == "__main__":
    main()