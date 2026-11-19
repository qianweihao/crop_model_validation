"""Step 02 - write APSIM .met files for the two weather stations.

The five weather sheets reduce to two physical stations: all three Dezhou sheets
agree exactly on common dates, likewise the two Xinxiang sheets (verified in
review).  The peripheral Dezhou sheets simply start two days earlier.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DS, MET, WEATHER_SITE  # noqa: E402

# annual mean / amplitude for the soil-temperature model; taken from long-term
# station normals (the field record only covers June-October)
TAV = {"dezhou": 14.4, "xinxiang": 15.2}
AMP = {"dezhou": 14.6, "xinxiang": 13.9}

START = pd.Timestamp("2026-06-01")
END = pd.Timestamp("2026-11-15")


def write_met(zone: str, df: pd.DataFrame) -> Path:
    site = WEATHER_SITE[zone]
    out = MET / f"{zone}.met"
    lines = [
        "[weather.met.weather]",
        f"!station name = {site['name']} (2026 field campaign)",
        "!source = regional daily weather sheet, 2026-06-08 .. 2026-10-15",
        f"latitude = {site['lat']:.4f}",
        f"longitude = {site['lon']:.4f}",
        f"tav = {TAV[zone]:.1f} (oC)     ! annual average ambient temperature",
        f"amp = {AMP[zone]:.1f} (oC)     ! annual amplitude in mean monthly temperature",
        "",
        "Year   Day    radn    maxt    mint    rain",
        " ()     ()      ()      ()      ()      ()",
    ]
    for _, r in df.iterrows():
        lines.append(
            f"{r.date.year:4d} {r.date.dayofyear:5d} {r.radn:8.2f} {r.tmax:7.1f} "
            f"{r.tmin:7.1f} {r.rain:7.1f}")
    out.write_text("\n".join(lines) + "\n", encoding="ascii")
    return out


def main() -> None:
    w = pd.read_csv(DS / "weather_daily.csv", parse_dates=["date"])
    for zone, g in w.groupby("zone"):
        g = g[(g["date"] >= START) & (g["date"] <= END)].sort_values("date")
        # APSIM rejects maxt <= mint and negative radiation
        g = g.assign(maxt=g[["tmax", "tmin"]].max(axis=1),
                     mint=g[["tmax", "tmin"]].min(axis=1))
        g["radn"] = g["radn"].clip(lower=0.1)
        bad = g[["radn", "maxt", "mint", "rain"]].isna().any(axis=1).sum()
        p = write_met(zone, g)
        print(f"{zone:9s} {p.name}  {len(g)} days  "
              f"{g.date.min().date()}..{g.date.max().date()}  "
              f"rain {g.rain.sum():.0f} mm  radn {g.radn.mean():.1f} MJ/m2  "
              f"T {g.tavg.min():.1f}-{g.tavg.max():.1f}  missing={bad}")
        if (g["maxt"] <= g["mint"]).any():
            print("  WARNING: maxt <= mint after fix")


if __name__ == "__main__":
    main()
