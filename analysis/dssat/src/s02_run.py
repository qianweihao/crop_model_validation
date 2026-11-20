"""Run all nine SD_SF plots and store the daily and seasonal outputs.

Outputs
-------
out/sim/dssat_daily.csv    one row per plot per day (PlantGro columns)
out/sim/dssat_summary.csv  one row per plot (phenology, yield, biomass, N)
out/mgmt/*.MZX             copy of the generated management file per plot

Usage
-----
    python analysis/dssat/src/s02_run.py [--cultivar P1=235,...]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

KEEP = ["date", "plot", "DAP", "GSTD", "L#SD", "LAID", "LWAD", "SWAD", "GWAD",
        "CWAD", "VWAD", "G#AD", "GWGD", "HIAD", "NSTD", "P#AD", "RDPD", "SLAD"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    ap.add_argument("--cultivar", default=None,
                    help="JSON dict or k=v,k=v overrides, e.g. P1=235,P5=880")
    # default resolved after set_trial(): module defaults are bound at import,
    # i.e. before the second trial rebinds them.
    ap.add_argument("--scenario", default="base", choices=list(dc.SCENARIOS))
    ap.add_argument("--nitrate-scale", type=float, default=None)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    dc.set_trial(args.trial)
    dc.set_scenario(args.scenario)
    print(f"trial {dc.TRIAL}  ({dc.CFG.site} {dc.CFG.variety})  "
          f"scenario {dc.SCENARIO} ({dc.SCEN_LABEL})  -> {dc.SIM}")
    if args.nitrate_scale is None:
        args.nitrate_scale = dc.DEFAULT_NITRATE_SCALE
    cultivar = dict(dc.DEFAULT_CULTIVAR)
    if args.cultivar:
        if args.cultivar.strip().startswith("{"):
            cultivar.update(json.loads(args.cultivar))
        else:
            for kv in args.cultivar.split(","):
                k, v = kv.split("=")
                cultivar[k.strip()] = float(v)

    daily, summary = [], []
    for plot in dc.ORDER:
        g, s = dc.run_plot(plot, cultivar, nitrate_scale=args.nitrate_scale, save=True)
        daily.append(g[[c for c in KEEP if c in g.columns]])
        pl = dc.load_plots().loc[plot]
        dcol = ("density_planted" if dc.DENSITY_MODE == "nominal"
                else "density_measured")
        if pd.isna(pl[dcol]):                 # no stand count
            dcol = "density_planted"
        row = dict(plot=plot, n_rate=dc.N_RATE[plot],
                   density=float(pl[dcol]))
        row.update({k: s.get(k) for k in
                    ["PDAT", "EDAT", "ADAT", "MDAT", "HDAT", "HWAM", "CWAM",
                     "LAIX", "HIAM", "G#AM", "GWUM", "NUCM", "NICM", "NMINC",
                     "PRCM", "ETCM", "IRCM"]})
        row.update({f"par_{k}": v for k, v in cultivar.items()})
        summary.append(row)
        def num(k):
            v = s.get(k)
            return float("nan") if v is None else float(v)

        def day(k):
            d = dc.dssat_date(s.get(k))
            return "--" if d is None else f"{d:%m-%d}"

        print(f"  {plot}  LAI {num('LAIX'):.2f}  yield {num('HWAM'):.0f}  "
              f"anthesis {day('ADAT')}  maturity {day('MDAT')}")

    d = pd.concat(daily, ignore_index=True)
    d.to_csv(dc.SIM / "dssat_daily.csv", index=False)
    pd.DataFrame(summary).to_csv(dc.SIM / "dssat_summary.csv", index=False)
    (dc.SIM / "calibrated_cultivar.json").write_text(
        json.dumps({"cultivar": cultivar, "nitrate_scale": args.nitrate_scale,
                    "sloc_profile": dc.SLOC_PROFILE}, indent=2), encoding="utf-8")
    print(f"wrote {len(d)} daily rows -> {dc.SIM / 'dssat_daily.csv'}")


if __name__ == "__main__":
    main()
