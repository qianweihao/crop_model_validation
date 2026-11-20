"""Comparison figures, produced twice under two names.

The same four drawings are emitted as the C series and as the L series, so either
naming can be cited from the reports without the content drifting apart:

  x1  LAI              nine per-plot panels, observed vs simulated
  x2  地上部生物量       nine per-plot panels, observed vs simulated
  x3  产量              per-plot bars with bias, plus a 1:1 scatter
  x4  氮响应             yield / biomass / HI against nitrogen rate

Usage
-----
    python analysis/src/s20_conservative_figures.py   # C series
    python analysis/src/s12_figures.py               # L series
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DS, FIG, SIM  # noqa: E402

plt.rcParams.update({
    "font.sans-serif": ["SimHei", "Microsoft YaHei", "DejaVu Sans"],
    "axes.unicode_minus": False, "figure.dpi": 130, "savefig.dpi": 150,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": "--",
    "axes.labelsize": 9.5, "axes.titlesize": 10, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
})

TRIAL = "SD_SF"
ORDER = ["A1", "A2B1", "A2B2", "A3B1", "A3B2", "A4B1", "A4B2", "A5B1", "A5B2"]
N_RATE = dict(zip(ORDER, [0, 139, 139, 177, 177, 215, 215, 254, 254]))
SEN = "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue"
SIMCOL = "#c0392b"
OBSCOL = "black"
PLANTS_PER_MU = 667.0        # density is plants per mu (1 mu = 667 m2)


# --------------------------------------------------------------------------
def load():
    sim = pd.read_csv(SIM / "sim_calibrated_daily.csv", parse_dates=["Clock.Today"])
    sim = sim[sim["trial_id"] == TRIAL]
    obs = pd.read_csv(DS / "obs_plant.csv", parse_dates=["date"])
    obs = obs[obs["trial_id"] == TRIAL].copy()
    obs["bio_kg_ha"] = obs["total_dw"] * obs["density_measured"] / PLANTS_PER_MU * 10.0
    val = pd.read_csv(SIM / "validate_plots.csv")
    return sim, obs, val[val["trial_id"] == TRIAL].set_index("plot")


def _series(sim, plot) -> pd.DataFrame:
    return sim[sim["sim"] == f"{TRIAL}_{plot}"].sort_values("Clock.Today")


# --------------------------------------------------------------------------
# nitrogen-stress reference run: every plot re-simulated with ample soil N, so
# that  stress = 1 - actual/potential  is a measurable, model-internal quantity
# --------------------------------------------------------------------------
N_SHARE = [0.28, 0.22, 0.26, 0.14, 0.07, 0.03]
N_AMPLE = 400.0                     # kg N/ha initial, far above any demand


def _params(sen: float | None = None) -> dict:
    import io
    import json
    cal = json.load(io.open(SIM / "phenology_calibrated.json", encoding="utf-8"))["乐农87"]["params"]
    gro = json.load(io.open(SIM / "calibrated_cultivars.json", encoding="utf-8"))["乐农87"]["growth"]
    p = {**cal, **gro}
    if sen is not None:
        p[SEN] = sen
    return p


def _set_initial_n(case, ntot: float):
    fld = case._field()
    soil = next(c for c in fld["Children"] if c.get("Name") == "Soil")
    k = {c["Name"]: c for c in soil["Children"]}
    bd, th = k["Physical"]["BD"], k["Physical"]["Thickness"]
    k["NO3"]["InitialValues"] = [round(ntot * s * 0.9 / (b * t * 0.01), 2)
                                  for s, b, t in zip(N_SHARE, bd, th)]
    k["NH4"]["InitialValues"] = [round(ntot * s * 0.1 / (b * t * 0.01), 2)
                                  for s, b, t in zip(N_SHARE, bd, th)]
    return case


def run_reference(sen: float | None = None) -> pd.DataFrame:
    """Simulate the nine plots with ample initial nitrogen (same management)."""
    import tempfile
    from apsim_override import SimCase, run_many

    params = _params(sen)
    jobs = []
    for plot in ORDER:
        case = SimCase.load(Path(SIM).parent / "apsim" / f"{TRIAL}_{plot}.apsimx")
        case = _set_initial_n(case.set_commands(params), N_AMPLE)
        jobs.append((case, f"pot__{plot}"))
    with tempfile.TemporaryDirectory() as td:
        tables = run_many(jobs, td)

    rows = []
    for plot in ORDER:
        df = tables.get(f"pot__{plot}")
        if df is None:
            continue
        df = df.copy()
        df["sim"] = f"{TRIAL}_{plot}"
        df["trial_id"], df["plot"] = TRIAL, plot
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def _ratio_table(sim, sim_pot) -> pd.DataFrame:
    """actual / potential for the five indicators that matter."""
    rows = []
    for plot in ORDER:
        a = _series(sim, plot).set_index("Clock.Today")
        b = _series(sim_pot, plot).set_index("Clock.Today")
        if a.empty or b.empty:
            continue
        ga = a["Maize.Grain.Number"].max()
        gb = b["Maize.Grain.Number"].max()
        rec = dict(plot=plot, N=N_RATE[plot])
        for key, col, f in (("yield", "Maize.Grain.Total.Wt", 10.0),
                            ("biomass", "Maize.AboveGround.Wt", 10.0),
                            ("uptake", "Maize.AboveGround.N", 10.0),
                            ("lai", "Maize.Leaf.LAI", 1.0)):
            rec[key] = a[col].max() * f
            rec[key + "_pot"] = b[col].max() * f
            rec[key + "_r"] = rec[key] / rec[key + "_pot"] if rec[key + "_pot"] else np.nan
        rec["grainw"] = a["Maize.Grain.Total.Wt"].max() / ga * 1000 if ga else np.nan
        rec["grainw_pot"] = b["Maize.Grain.Total.Wt"].max() / gb * 1000 if gb else np.nan
        rec["grainw_r"] = rec["grainw"] / rec["grainw_pot"] if rec["grainw_pot"] else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def _grid():
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 11.2), sharex=True, sharey=True)
    return fig, axes


def _finish(fig, axes, title, ylabel) -> None:
    for ax in axes[-1]:
        ax.set_xlabel("日期 (2026)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    for ax in axes.ravel():
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle(title, fontsize=13.5, y=0.995)


# --------------------------------------------------------------------------
def fig_lai(sim, obs, val, prefix: str, note: str = "") -> Path:
    fig, axes = _grid()
    for ax in axes.ravel():
        ax.set_ylim(0, 4.2)
    for ax, plot in zip(axes.ravel(), ORDER):
        o = obs[(obs["plot"] == plot)].dropna(subset=["lai_plot"])
        ax.plot(o["date"], o["lai_plot"], "o", ms=6.5, color=OBSCOL, zorder=5, label="实测")
        g = _series(sim, plot)
        ax.plot(g["Clock.Today"], g["Maize.Leaf.LAI"], lw=1.9, color=SIMCOL, label="模拟")
        ax.set_title(f"{plot}  N={N_RATE[plot]}  LAI RMSE={val.loc[plot, 'lai_rmse']:.2f}",
                     fontsize=9.5)
    axes[0, 0].legend(frameon=False, loc="upper left")
    _finish(fig, axes, f"{prefix}1  逐小区 LAI：实测 vs 模拟{note}", "LAI")
    fig.tight_layout()
    p = FIG / f"{prefix}1_LAI.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_biomass(sim, obs, prefix: str, note: str = "") -> Path:
    fig, axes = _grid()
    for ax in axes.ravel():
        ax.set_ylim(0, 20000)
    for ax, plot in zip(axes.ravel(), ORDER):
        o = obs[(obs["plot"] == plot) & obs["bio_kg_ha"].notna()]
        good = o[~o["organ_split_suspect"].fillna(False).astype(bool)]
        bad = o[o["organ_split_suspect"].fillna(False).astype(bool)]
        ax.plot(good["date"], good["bio_kg_ha"], "o", ms=6.5, color=OBSCOL,
                zorder=5, label="实测")
        if len(bad):
            ax.plot(bad["date"], bad["bio_kg_ha"], "^", ms=7, mfc="none",
                    mec="grey", mew=1.2, zorder=5, label="实测（器官分配可疑）")
        g = _series(sim, plot)
        ax.plot(g["Clock.Today"], g["Maize.AboveGround.Wt"] * 10, lw=1.9,
                color=SIMCOL, label="模拟")
        ax.set_title(f"{plot}  N={N_RATE[plot]}", fontsize=9.5)
    axes[0, 0].legend(frameon=False, loc="upper left", fontsize=7)
    _finish(fig, axes, f"{prefix}2  逐小区地上部生物量：实测 vs 模拟{note}", "生物量 (kg/ha)")
    fig.tight_layout()
    p = FIG / f"{prefix}2_生物量.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_yield(val, prefix: str, note: str = "") -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.0),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ob = val["obs_yield_dry"].reindex(ORDER).values
    sb = val["sim_yield"].reindex(ORDER).values

    ax = axes[0]
    x = np.arange(len(ORDER))
    ax.bar(x - 0.2, ob, 0.4, label="实测", color="#7f8c8d", edgecolor="k", lw=0.4)
    ax.bar(x + 0.2, sb, 0.4, label="模拟", color=SIMCOL, edgecolor="k", lw=0.4)
    for i, (a, b) in enumerate(zip(ob, sb)):
        ax.text(i, max(a, b) + 180, f"{(b / a - 1) * 100:+.0f}%", ha="center", fontsize=7.5)
    ax.set_xticks(x, [f"{p}\nN={N_RATE[p]}" for p in ORDER], fontsize=7.5)
    ax.set_ylabel("产量 (kg/ha, 干重)")
    ax.set_title("逐小区产量")
    ax.legend(frameon=False)
    ax.set_ylim(0, max(ob.max(), sb.max()) * 1.16)

    ax = axes[1]
    ax.scatter(ob, sb, s=85, c=SIMCOL, edgecolor="k", lw=0.5, zorder=3)
    lim = [min(ob.min(), sb.min()) * 0.9, max(ob.max(), sb.max()) * 1.06]
    ax.plot(lim, lim, "k--", lw=1.2)
    for p, a, b in zip(ORDER, ob, sb):
        ax.annotate(p, (a, b), fontsize=7.5, xytext=(6, 3), textcoords="offset points")
    rrmse = np.sqrt(np.mean((sb - ob) ** 2)) / ob.mean() * 100
    ax.text(0.04, 0.95, f"RRMSE = {rrmse:.2f}%\n偏差 = {(sb.sum() / ob.sum() - 1) * 100:+.2f}%",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(fc="white", ec="grey", alpha=0.9))
    ax.set_xlabel("实测产量 (kg/ha)")
    ax.set_ylabel("模拟产量 (kg/ha)")
    ax.set_title("1:1")

    fig.suptitle(f"{prefix}3  逐小区产量：实测 vs 模拟{note}", fontsize=13.5, y=1.01)
    fig.tight_layout()
    p = FIG / f"{prefix}3_产量.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_nitrogen_response(val, prefix: str, note: str = "") -> Path:
    """Yield, biomass and harvest index against nitrogen rate (the old L3 panel)."""
    n = np.array([N_RATE[p] for p in ORDER], dtype=float)
    panels = [
        ("产量 (kg/ha, 干重)", "obs_yield_dry", "sim_yield", "产量对氮的响应"),
        ("地上部生物量 (kg/ha)", "obs_bio_max", "sim_bio_max", "生物量对氮的响应"),
        ("收获指数 HI", "obs_HI", "sim_HI", "收获指数"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (ylabel, oc, sc, title) in zip(axes, panels):
        o = val[oc].reindex(ORDER).values
        s = val[sc].reindex(ORDER).values
        ax.plot(n, o, "o-", ms=8, lw=1.8, color=OBSCOL, label="实测", zorder=4)
        ax.plot(n, s, "s--", ms=7, lw=1.8, color=SIMCOL, label="模拟", zorder=4)
        ax.set_xlabel("施氮量 (kg N/ha)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False)
        lo = min(np.nanmin(o), np.nanmin(s))
        hi = max(np.nanmax(o), np.nanmax(s))
        pad = (hi - lo) * 0.18 or 1.0
        ax.set_ylim(lo - pad, hi + pad)
    d_o = (val["obs_yield_dry"].reindex(ORDER).iloc[-1] / val["obs_yield_dry"].reindex(ORDER).iloc[0] - 1) * 100
    d_s = (val["sim_yield"].reindex(ORDER).iloc[-1] / val["sim_yield"].reindex(ORDER).iloc[0] - 1) * 100
    axes[0].text(0.04, 0.94, f"0→254 kg N/ha\n实测 {d_o:+.0f}%  模拟 {d_s:+.0f}%",
                 transform=axes[0].transAxes, va="top", fontsize=8.5,
                 bbox=dict(fc="white", ec="grey", alpha=0.9))
    fig.suptitle(f"{prefix}4  氮素响应：模拟 vs 实测{note}", fontsize=13.5, y=1.02)
    fig.tight_layout()
    p = FIG / f"{prefix}4_氮响应.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_nstress(sim, sim_pot, prefix: str, note: str = "") -> Path:
    """Nitrogen stress: how much growth is lost to nitrogen shortage.

    Panel A quantifies the loss as 1 - actual/potential, where ``potential`` is the
    same plot re-simulated with 400 kg N/ha of initial soil nitrogen.  Panels B-D
    show when and through which process the gap opens, comparing A1 against the
    envelope of the eight fertilised plots.
    """
    rt = _ratio_table(sim, sim_pot)
    others = [p for p in ORDER if p != "A1"]

    def band(metric: str) -> pd.DataFrame:
        d = {}
        for plot in ORDER:
            g = _series(sim, plot).set_index("Clock.Today")
            if g.empty:
                continue
            if metric == "grainw":
                n = g["Maize.Grain.Number"]
                with np.errstate(divide="ignore", invalid="ignore"):
                    d[plot] = np.where(n > 0, g["Maize.Grain.Total.Wt"] / n * 1000, np.nan)
            elif metric == "uptake":
                d[plot] = g["Maize.AboveGround.N"] * 10.0
            else:  # topsoil nitrate, 0-30 cm
                cols = [c for c in g.columns if c.startswith("Soil.Nutrient.NO3.kgha")][:2]
                d[plot] = g[cols].sum(axis=1)
        return pd.DataFrame(d)

    fig, axes = plt.subplots(1, 4, figsize=(19.5, 4.6))

    # ---- A: stress intensity bars ----
    ax = axes[0]
    inds = [("yield_r", "产量", "#c0392b"), ("biomass_r", "生物量", "#e67e22"),
            ("grainw_r", "单粒重", "#8e44ad"), ("uptake_r", "吸氮量", "#2980b9")]
    x = np.arange(len(rt))
    w = 0.2
    for i, (col, lab, colr) in enumerate(inds):
        vals = (1.0 - rt[col].to_numpy()) * 100.0
        ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=colr, edgecolor="k", lw=0.3)
    ax.set_xticks(x, [f"{p}\nN={n}" for p, n in zip(rt["plot"], rt["N"])], fontsize=7)
    ax.set_ylabel("氮胁迫强度  1 − 实际/潜在 (%)")
    ax.set_title("A  氮胁迫强度（潜在 = 初始氮 400 kg N/ha）")
    ax.legend(frameon=False, ncol=2, fontsize=7.5)
    a1 = rt.index[rt["plot"] == "A1"][0]
    ax.annotate("仅 A1 受氮限制", xy=(a1, 20), xytext=(a1 + 1.4, 22),
                fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.9, color="grey"))
    ax.set_ylim(0, max(26, (1 - rt[[c for c, _, _ in inds]].min().min()) * 100 * 1.25))

    # ---- B-D: when the gap opens ----
    panels = [
        ("uptake", "地上部吸氮量 (kg N/ha)", "B  吸氮量：A1 vs 其余 8 小区范围",
         "约 7 月下旬开始分离"),
        ("grainw", "单粒重 (mg)", "C  单粒重：灌浆期分化", "9 月中旬开始分化"),
        ("no3", "表层 NO₃-N (kg/ha)", "D  表层土壤硝态氮", "A1 无追肥，供氮源仅为初始氮+矿化"),
    ]
    for ax, (metric, ylab, title, note_txt) in zip(axes[1:], panels):
        d = band(metric)
        if d.empty:
            continue
        o = d[others]
        ax.fill_between(d.index, o.min(axis=1), o.max(axis=1), color="#95a5a6",
                        alpha=0.30, lw=0, label="其余 8 小区范围")
        ax.plot(d.index, o.mean(axis=1), lw=1.4, color="#7f8c8d", label="其余 8 小区均值")
        ax.plot(d.index, d["A1"], lw=2.2, color="#c0392b", label="A1 (N=0)")
        ax.set_ylabel(ylab)
        ax.set_title(title)
        ax.set_xlabel("日期 (2026)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax.tick_params(axis="x", rotation=45)
        ax.legend(frameon=False, fontsize=7.5)
        ax.text(0.03, 0.06, note_txt, transform=ax.transAxes, fontsize=7.5,
                color="#555555")

    fig.suptitle(f"{prefix}5  氮胁迫：A1 与其余处理的差异{note}", fontsize=13.5, y=1.02)
    fig.tight_layout()
    p = FIG / f"{prefix}5_氮胁迫.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)

    rt.to_csv(SIM / f"nstress_intensity_{prefix}.csv", index=False)
    print(f"  氮胁迫强度（1-实际/潜在，%），{prefix} 系列：")
    for _, r in rt.iterrows():
        print(f"    {r['plot']:5s} N={int(r['N']):3d}  产量 {(1-r['yield_r'])*100:5.1f}  "
              f"生物量 {(1-r['biomass_r'])*100:5.1f}  单粒重 {(1-r['grainw_r'])*100:5.1f}  "
              f"吸氮 {(1-r['uptake_r'])*100:5.1f}  LAI {(1-r['lai_r'])*100:5.1f}")
    return p


def make_all(prefix: str, sim=None, obs=None, val=None, note: str = "",
             sim_pot=None, sen: float | None = None) -> list[Path]:
    """Draw the four figures.  Data can be injected so a second parameter set
    (the pre-optimisation baseline) can reuse exactly the same code path."""
    if sim is None or obs is None or val is None:
        sim, obs, val = load()
    if sim_pot is None:
        sim_pot = run_reference(sen)
    return [fig_lai(sim, obs, val, prefix, note), fig_biomass(sim, obs, prefix, note),
            fig_yield(val, prefix, note), fig_nitrogen_response(val, prefix, note),
            fig_nstress(sim, sim_pot, prefix, note)]


def main() -> None:
    sim, obs, val = load()
    print(f"sim rows={len(sim)}  plots={sim['sim'].nunique()}")
    for p in make_all("C", note="（sen=0.20，优化后）"):
        print(" ", p)
    o, s = val["obs_yield_dry"], val["sim_yield"]
    print(f"\n产量 RRMSE {np.sqrt(np.mean((s - o) ** 2)) / o.mean() * 100:.2f}%  "
          f"偏差 {(s.sum() / o.sum() - 1) * 100:+.2f}%  "
          f"LAI RMSE {val['lai_rmse'].mean():.3f}")


if __name__ == "__main__":
    main()