"""D5 — nitrogen stress: when it starts, and what it costs.

DSSAT reports its own stress factor ``NSTD`` (0 = no stress, 1 = full stress) for
every day, so the figure is built around *timing* rather than around a binary
"limited / not limited" label:

  A  NSTD through the season, all nine plots coloured by N rate.  The onset
     moves later as N increases; only stress that starts before silking costs
     yield.
  B  stress intensity 1 - actual/potential, where "potential" is the same plot
     re-simulated with 400 kg N/ha of initial soil nitrogen.
  C  N uptake against N applied, actual vs potential: uptake saturates around
     170 kg N/ha while the potential ceiling is ~220.
  D  onset date (relative to silking) against the yield gap: the causal summary
     of panels A and B.

Outputs
-------
out/sim/nstress_potential_daily.csv   potential run, daily
out/sim/nstress_intensity.csv         per-plot ratios + NSTD onset
fig/D5_氮胁迫.png

Usage
-----
    python analysis/dssat/src/s06_nstress.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["SimHei", "Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 130, "savefig.dpi": 150,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    "axes.labelsize": 9.5, "axes.titlesize": 10, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
})

N_AMPLE = 400.0          # kg N/ha initial for the potential run
ONSET_THRESHOLD = 0.05   # NSTD above this counts as "stress has started"
ONSET_SUSTAIN = 7        # ... and must stay above it for this many days:
                         # XX_SF has a one-week transient around 07-29 in every
                         # plot that a single-day rule would mistake for onset.
A1COL = "#c0392b"
CMAP = "viridis"


def series(daily: pd.DataFrame, soil_ni: pd.DataFrame) -> pd.DataFrame:
    d = daily.set_index(daily["date"].dt.normalize())
    n = soil_ni.set_index(soil_ni["date"].dt.normalize())
    return pd.DataFrame({
        "lai": d["LAID"],
        "grainw": d["GWGD"].where(d["GWGD"] > 0),
        "biomass": d["CWAD"],
        "yield": d["GWAD"],
        "nstd": d["NSTD"],
        "uptake": n["NUCM"].reindex(d.index).ffill(),
        "no3": n["NITD"].reindex(d.index).ffill(),
    })


def _sustained_onset(nstd: pd.Series, thresh: float = ONSET_THRESHOLD,
                     n: int = ONSET_SUSTAIN) -> pd.Timestamp:
    """First day of the first >= n-day run with NSTD above ``thresh``."""
    above = (nstd > thresh).to_numpy()
    run = 0
    for i, flag in enumerate(above):
        run = run + 1 if flag else 0
        if run >= n:
            return nstd.index[i - n + 1]
    return pd.NaT


def collect() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    act_rows, pot_rows, ratios = [], [], []
    for plot in dc.ORDER:
        actual, sa = dc.run_plot(plot, soil_n_output=True)
        potential, sp = dc.run_plot(plot, n_total_kg=N_AMPLE, soil_n_output=True)
        sa_s = series(actual, sa.pop("soil_ni"))
        sp_s = series(potential, sp.pop("soil_ni"))
        act_rows.append(sa_s.assign(plot=plot))
        pot_rows.append(sp_s.assign(plot=plot))

        silking = dc.dssat_date(sa.get("ADAT"))
        onset = _sustained_onset(sa_s["nstd"])

        rec = dict(plot=plot, N=dc.N_RATE[plot],
                   silking=silking, onset=onset,
                   onset_vs_silking=((onset - silking).days
                                     if pd.notna(onset) else np.nan),
                   stress_days=int((sa_s["nstd"] > ONSET_THRESHOLD).sum()))
        for key in ("yield", "biomass", "grainw", "uptake", "lai"):
            rec[key] = float(np.nanmax(sa_s[key]))
            rec[key + "_pot"] = float(np.nanmax(sp_s[key]))
            rec[key + "_r"] = rec[key] / rec[key + "_pot"] if rec[key + "_pot"] else np.nan
        ratios.append(rec)
        print(f"  {plot}  N={rec['N']:3d}  胁迫起始 {onset:%m-%d} "
              f"({rec['onset_vs_silking']:+.0f} d vs 吐丝)  "
              f"产量差 {(1 - rec['yield_r']) * 100:4.1f}%")

    return (pd.concat(pot_rows).reset_index().rename(columns={"index": "date"}),
            pd.concat(act_rows).reset_index().rename(columns={"index": "date"}),
            pd.DataFrame(ratios))


def draw(act: pd.DataFrame, rt: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.4))
    nrates = sorted(rt["N"].unique())
    colours = {n: plt.get_cmap(CMAP)(i / max(1, len(nrates) - 1))
               for i, n in enumerate(nrates)}

    # ---- A: NSTD through the season ----
    ax = axes[0, 0]
    d = act.pivot_table(index="date", columns="plot", values="nstd").sort_index()
    for plot in dc.ORDER:
        n = dc.N_RATE[plot]
        ctrl = dc.CONTROL_PLOT
        ax.plot(d.index, d[plot], lw=2.4 if plot == ctrl else 1.3,
                color=A1COL if plot == ctrl else colours[n],
                alpha=1.0 if plot == ctrl else 0.75,
                zorder=5 if plot == ctrl else 3,
                label=f"{ctrl} (N=0)" if plot == ctrl else None)
    for n in nrates:
        ax.plot([], [], color=colours[n], lw=1.6, label=f"N={n}")
    silking = rt["silking"].iloc[0]
    ax.axvline(silking, color="k", ls=":", lw=1.2)
    ax.text(silking, 0.95, " 吐丝", transform=ax.get_xaxis_transform(),
            fontsize=8, va="top")
    ax.set_ylabel("氮胁迫因子 NSTD  (0 = 无胁迫)")
    ax.set_title("A  氮胁迫的时间进程：施氮越高，起始越晚")
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.tick_params(axis="x", rotation=45)
    ax.set_xlabel("日期 (2026)")

    # ---- B: stress intensity ----
    ax = axes[0, 1]
    inds = [("yield_r", "产量", "#c0392b"), ("biomass_r", "生物量", "#e67e22"),
            ("uptake_r", "吸氮量", "#2980b9"), ("lai_r", "LAI", "#16a085")]
    x = np.arange(len(rt))
    w = 0.2
    for i, (col, lab, colr) in enumerate(inds):
        ax.bar(x + (i - 1.5) * w, (1 - rt[col].to_numpy()) * 100, w,
               label=lab, color=colr, edgecolor="k", lw=0.3)
    ax.set_xticks(x, [f"{p}\nN={n}" for p, n in zip(rt["plot"], rt["N"])], fontsize=7)
    ax.set_ylabel("胁迫强度  1 - 实际/潜在 (%)")
    ax.set_title("B  氮胁迫强度（潜在 = 初始氮 400 kg N/ha）")
    ax.legend(frameon=False, ncol=2, fontsize=7.5)
    ax.set_ylim(0, max(35, float((1 - rt[[c for c, _, _ in inds]].min().min()) * 100) * 1.25))

    # ---- C: uptake vs applied N ----
    ax = axes[1, 0]
    ax.plot(rt["N"], rt["uptake_pot"], "o", ms=8, mfc="none", mec="#7f8c8d",
            mew=1.6, label="潜在吸氮（初始氮 400）")
    ax.plot(rt["N"], rt["uptake"], "o", ms=8, color=A1COL, label="实际吸氮")
    for r in rt.itertuples():
        ax.annotate(r.plot, (r.N, r.uptake), fontsize=7, xytext=(4, -10),
                    textcoords="offset points")
    ax.axhline(rt["uptake_pot"].max(), color="#7f8c8d", ls="--", lw=1)
    ax.text(rt["N"].min(), rt["uptake_pot"].max() + 4,
            f"上限 ≈ {rt['uptake_pot'].max():.0f} kg N/ha", fontsize=8, color="#555")
    ax.set_xlabel("施氮量 (kg N/ha)")
    ax.set_ylabel("季末地上部吸氮量 (kg N/ha)")
    ax.set_title("C  吸氮量对施氮量的响应：约 170 kg N/ha 后饱和")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")

    # ---- D: onset vs yield gap ----
    ax = axes[1, 1]
    gap = (1 - rt["yield_r"]) * 100
    ax.scatter(rt["onset_vs_silking"], gap, s=90,
               c=[colours[n] for n in rt["N"]], edgecolor="k", lw=0.5, zorder=4)
    for r, g in zip(rt.itertuples(), gap):
        ax.annotate(f"{r.plot}\nN={r.N}", (r.onset_vs_silking, g), fontsize=7,
                    xytext=(6, 4), textcoords="offset points")
    ax.axvline(0, color="k", ls=":", lw=1.2)
    ax.text(0, 0.95, " 吐丝", transform=ax.get_xaxis_transform(),
            fontsize=8, va="top")
    ax.set_xlabel(f"NSTD>{ONSET_THRESHOLD} 持续≥{ONSET_SUSTAIN} 天首次出现 减 吐丝日期 (d)")
    ax.set_ylabel("产量差  1 - 实际/潜在 (%)")
    ax.set_title("D  胁迫起始越早，产量损失越大")
    ax.set_ylim(-2, gap.max() * 1.22)

    fig.suptitle(f"{dc.PREFIX}5｜{dc.DISPLAY} {dc.CFG.variety}｜{dc.SCEN_LABEL}"
                 f"　氮胁迫：并非只有 {dc.CONTROL_PLOT} 受限，"
                 f"而是胁迫起始时间随施氮量后移",
                 fontsize=13.5, y=0.995)
    fig.tight_layout()
    p = dc.FIG / f"{dc.PREFIX}5_氮胁迫.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    ap.add_argument("--scenario", default="base", choices=list(dc.SCENARIOS))
    args = ap.parse_args()
    dc.set_trial(args.trial)
    dc.set_scenario(args.scenario)
    print(f"trial {dc.TRIAL}  ({dc.CFG.site} {dc.CFG.variety})  "
          f"scenario {dc.SCENARIO} ({dc.SCEN_LABEL})")

    pot, act, rt = collect()
    pot.to_csv(dc.SIM / "nstress_potential_daily.csv", index=False)
    rt.to_csv(dc.SIM / "nstress_intensity.csv", index=False)
    p = draw(act, rt)
    print(f"\n{p}  ({p.stat().st_size / 1024:.0f} kB)")
    print("\n胁迫起始与产量损失：")
    for r in rt.itertuples():
        print(f"  {r.plot:5s} N={int(r.N):3d}  起始 {r.onset:%m-%d} "
              f"({r.onset_vs_silking:+3.0f} d)  胁迫天数 {r.stress_days:3d}  "
              f"产量 {(1 - r.yield_r) * 100:5.1f}%  "
              f"吸氮 {(1 - r.uptake_r) * 100:5.1f}%")


if __name__ == "__main__":
    main()
