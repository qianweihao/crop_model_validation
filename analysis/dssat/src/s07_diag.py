"""Four diagnostics on top of the calibrated runs (no existing output touched).

Q1  measured stand vs the nominal planting density (5500/5640 plants per mu)
Q2  can the observed per-plant leaf area be reproduced?
Q3  does the biomass shortfall explain the high harvest index, and can it be fixed?
Q4  are the simulated grain number / grain weight in line with the 考种 data?

Outputs
-------
out/<trial>/diag/density.csv         Q1 per plot, measured vs nominal
out/<trial>/diag/density_stats.csv   Q1 accuracy under both densities
out/<trial>/diag/leafarea.csv        Q2 density x PHINT x P1 sweep + observed
out/<trial>/diag/biomass_hi.csv      Q3 scenarios
out/<trial>/diag/grain.csv           Q4 model vs 考种 (XX_SF only)
out/<trial>/diag/grain_refit.csv     Q4 G2/G3 grid scored on grains/ear + 100-grain wt
fig/<trial>/diag/E_诊断.png          2x2 summary

Usage
-----
    python analysis/dssat/src/s07_diag.py --trial XX_SF
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

KAOS = Path("data/产量预估输入数据/新乡玉米考种数据-1.xlsx")
NOMINAL_PLOT = None          # resolved after set_trial (control / mid / high N)
DIAG: Path = dc.OUT / "diag"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def measured_density(plot: str) -> float:
    return float(dc.load_plots().loc[plot, "density_measured"]) / dc.PLANTS_PER_MU


def nominal_density() -> float:
    return float(dc.load_plots()["density_planted"].iloc[0]) / dc.PLANTS_PER_MU


def obs_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    o = pd.read_csv(dc.DATASET / "obs_plant.csv", parse_dates=["date"])
    o = o[o.trial_id == dc.TRIAL]
    y = pd.read_csv(dc.DATASET / "yield.csv")
    y = y[y.trial_id == dc.TRIAL].set_index("plot")
    p = dc.load_plots()
    return o, y, p


def obs_lai_peak(o: pd.DataFrame, plots) -> dict:
    return {p: float(o.loc[o["plot"] == p, "lai_plot"].max()) for p in plots}


def obs_biomass(o: pd.DataFrame, p: pd.DataFrame, plots) -> dict:
    out = {}
    for pl in plots:
        dens = p.loc[pl, "density_measured"] / dc.PLANTS_PER_MU
        out[pl] = float((o.loc[o["plot"] == pl, "total_dw"] * dens * 10).max())
    return out


def load_kaos() -> pd.DataFrame:
    """考种 data: mean of the ten sampled ears per plot (XX_SF only)."""
    if dc.TRIAL != "XX_SF" or not KAOS.exists():
        return pd.DataFrame()
    rows = []
    for s in pd.ExcelFile(KAOS).sheet_names:
        d = pd.read_excel(KAOS, sheet_name=s, header=0).iloc[:10]
        d.columns = ["n", "rows_per_ear", "grain_w_ear", "hundred_gw",
                     "grains_per_ear", "grains_per_row", "theo_yield_mu",
                     "x", "theo_yield_mu2", "moist", "theo_yield_14"]
        rows.append(dict(plot=s, rows_per_ear=d.rows_per_ear.mean(),
                         grain_w_ear=d.grain_w_ear.mean(),
                         hundred_gw=d.hundred_gw.mean(),
                         grains_per_ear=d.grains_per_ear.mean(),
                         grains_per_row=d.grains_per_row.mean(),
                         theo_kg_ha=d.theo_yield_14.mean() * 15))
    return pd.DataFrame(rows).set_index("plot")


# --------------------------------------------------------------------------
# Q1  measured stand vs nominal planting density
# --------------------------------------------------------------------------
def q1_density() -> None:
    o, y, p = obs_frames()
    dn = nominal_density()
    rows = []
    for plot in dc.ORDER:
        dm = measured_density(plot)
        rec = dict(plot=plot, n_rate=dc.N_RATE[plot], dens_measured=dm,
                   dens_nominal=dn, obs_yield=y.loc[plot, "yield_kg_ha"] * 0.86)
        for tag, den in (("meas", dm), ("nom", dn)):
            g, s = dc.run_plot(plot, density=den)
            m = dc.plot_metrics(g, s, den)
            rec.update({f"{tag}_{k}": m[k] for k in
                        ("lai_peak", "la_per_plant", "bio_final", "yield_kg",
                         "hi", "grains_m2", "grain_w_mg", "grains_per_ear")})
        rec["obs_lai_peak"] = float(o.loc[o["plot"] == plot, "lai_plot"].max())
        rec["obs_bio"] = float((o.loc[o["plot"] == plot, "total_dw"]
                                * dm * 10).max())
        rows.append(rec)
    d = pd.DataFrame(rows)
    d.to_csv(DIAG / "density.csv", index=False)

    def stats(prefix: str) -> dict:
        ysim = d[f"{prefix}_yield_kg"].values
        yobs = d["obs_yield"].values
        e = ysim - yobs
        return dict(
            yield_rrmse=float(np.sqrt(np.mean(e ** 2)) / yobs.mean() * 100),
            yield_bias=float(e.sum() / yobs.sum() * 100),
            lai_rmse=float(np.sqrt(np.mean((d[f"{prefix}_lai_peak"]
                                            - d["obs_lai_peak"]) ** 2))),
            bio_bias=float((d[f"{prefix}_bio_final"].sum()
                            / d["obs_bio"].sum() - 1) * 100),
            hi_mean=float(d[f"{prefix}_hi"].mean()),
            grains_per_ear=float(d[f"{prefix}_grains_per_ear"].mean()),
        )

    st = pd.DataFrame({"measured": stats("meas"), "nominal": stats("nom")}).T
    st["obs_hi"] = float((d["obs_yield"] / d["obs_bio"]).mean())
    st.to_csv(DIAG / "density_stats.csv")
    print(st.round(3).to_string())


# --------------------------------------------------------------------------
# Q2  per-plant leaf area
# --------------------------------------------------------------------------
def q2_leafarea() -> None:
    o, y, p = obs_frames()
    # observed per-plant leaf area
    rows = []
    for plot in dc.ORDER:
        dens = measured_density(plot)
        lai = float(o.loc[o["plot"] == plot, "lai_plot"].max())
        rows.append(dict(plot=plot, density=dens, lai_peak=lai,
                         la_per_plant=lai / dens, source="observed"))
    # modelled: density x PHINT x P1
    ref = max(dc.ORDER, key=lambda q: float(
        o.loc[o["plot"] == q, "lai_plot"].max()))
    base = dict(dc.DEFAULT_CULTIVAR)
    for dens in (6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 11.0):
        for phint in (45.0, 50.0, 55.0):
            for p1 in (base["P1"], 275.0):
                cv = dict(base, PHINT=phint, P1=p1)
                g, s = dc.run_plot(ref, cv, density=dens)
                m = dc.plot_metrics(g, s, dens)
                rows.append(dict(plot=ref, density=dens, lai_peak=m["lai_peak"],
                                 la_per_plant=m["la_per_plant"],
                                 source=f"model P1={p1:.0f} PHINT={phint:.0f}",
                                 phint=phint, p1=p1, yield_kg=m["yield_kg"]))
    d = pd.DataFrame(rows)
    d.to_csv(DIAG / "leafarea.csv", index=False)
    mod = d[d.source != "observed"]
    print(f"  reference plot {ref}; model max per-plant leaf area "
          f"{mod.la_per_plant.max():.3f} m2, observed max "
          f"{d[d.source == 'observed'].la_per_plant.max():.3f} m2")
    print(f"  model LAI ceiling {mod.lai_peak.max():.2f} at density "
          f"{mod.loc[mod.lai_peak.idxmax(), 'density']:.1f}")


# --------------------------------------------------------------------------
# Q3  biomass and harvest index
# --------------------------------------------------------------------------
def q3_biomass() -> None:
    o, y, p = obs_frames()
    dn = nominal_density()
    scen = [
        ("base", {}, {}, None),
        ("nominal density", {}, {}, dn),
        ("PARSR 0.60", {}, dict(PARSR=0.60), None),
        ("PARSR 0.80", {}, dict(PARSR=0.80), None),
        ("initial N 400", {}, {}, None),
        ("PARSR 0.60 + N400", {}, dict(PARSR=0.60), None),
        ("PARSR 0.60, G3 -15%", None, dict(PARSR=0.60), None),
        ("G2 1100 G3 8.0", dict(G2=1100, G3=8.0), {}, None),
    ]
    base = dict(dc.DEFAULT_CULTIVAR)
    rows = []
    for tag, cv_ov, sp, den in scen:
        cv = dict(base)
        if cv_ov is None:                     # G3 -15 % relative to base
            cv["G3"] = base["G3"] * 0.85
        elif cv_ov:
            cv.update(cv_ov)
        n_tot = 400 if "N400" in tag or tag == "initial N 400" else None
        for plot in dc.ORDER:
            dens = den or measured_density(plot)
            if sp:
                with dc.patch_species(**sp):
                    g, s = dc.run_plot(plot, cv, density=dens,
                                       n_total_kg=n_tot)
            else:
                g, s = dc.run_plot(plot, cv, density=dens, n_total_kg=n_tot)
            m = dc.plot_metrics(g, s, dens)
            obs_bio = float((o.loc[o["plot"] == plot, "total_dw"]
                             * measured_density(plot) * 10).max())
            rows.append(dict(scenario=tag, plot=plot, n_rate=dc.N_RATE[plot],
                             density=dens, lai_peak=m["lai_peak"],
                             bio=m["bio_final"], yield_kg=m["yield_kg"],
                             hi=m["hi"], grain_w_mg=m["grain_w_mg"],
                             grains_per_ear=m["grains_per_ear"],
                             obs_yield=y.loc[plot, "yield_kg_ha"] * 0.86,
                             obs_bio=obs_bio))
    d = pd.DataFrame(rows)
    d.to_csv(DIAG / "biomass_hi.csv", index=False)
    agg = d.groupby("scenario").agg(
        bio=("bio", "mean"), yield_kg=("yield_kg", "mean"),
        hi=("hi", "mean"), grain_w_mg=("grain_w_mg", "mean"),
        grains_per_ear=("grains_per_ear", "mean"),
        obs_bio=("obs_bio", "mean"), obs_yield=("obs_yield", "mean"))
    agg["bio_bias_pct"] = (agg.bio / agg.obs_bio - 1) * 100
    agg["obs_hi"] = agg.obs_yield / agg.obs_bio
    print(agg.round(3).to_string())


# --------------------------------------------------------------------------
# Q4  grain number and grain weight
# --------------------------------------------------------------------------
def q4_grain() -> None:
    o, y, p = obs_frames()
    kaos = load_kaos()
    dn = nominal_density()
    base = dict(dc.DEFAULT_CULTIVAR)
    rows = []
    for plot in dc.ORDER:
        dm = measured_density(plot)
        rec = dict(plot=plot, n_rate=dc.N_RATE[plot], density=dm)
        for tag, den in (("meas", dm), ("nom", dn)):
            g, s = dc.run_plot(plot, density=den)
            m = dc.plot_metrics(g, s, den)
            rec[f"{tag}_grains_per_ear"] = m["grains_per_ear"]
            rec[f"{tag}_grain_w_mg"] = m["grain_w_mg"]
            rec[f"{tag}_grains_m2"] = m["grains_m2"]
        if not kaos.empty and plot in kaos.index:
            rec["obs_grains_per_ear"] = kaos.loc[plot, "grains_per_ear"]
            rec["obs_hundred_gw_mg"] = kaos.loc[plot, "hundred_gw"] * 10
            rec["obs_rows_per_ear"] = kaos.loc[plot, "rows_per_ear"]
            rec["obs_grains_per_row"] = kaos.loc[plot, "grains_per_row"]
            rec["obs_grain_w_ear"] = kaos.loc[plot, "grain_w_ear"]
        rows.append(rec)
    d = pd.DataFrame(rows)
    d.to_csv(DIAG / "grain.csv", index=False)
    if "obs_grains_per_ear" in d:
        print(d[["plot", "meas_grains_per_ear", "nom_grains_per_ear",
                 "obs_grains_per_ear", "meas_grain_w_mg", "obs_hundred_gw_mg"]]
              .round(1).to_string(index=False))
    else:
        print(d[["plot", "meas_grains_per_ear", "nom_grains_per_ear",
                 "meas_grain_w_mg"]].round(1).to_string(index=False))

    if kaos.empty:
        return
    # re-fit G2/G3 on the 考种 targets instead of the borrowed 240 mg anchor
    g_tgt = float(d["obs_grains_per_ear"].mean())
    w_tgt = float(d["obs_hundred_gw_mg"].mean())
    grid = []
    for g2 in (900.0, 1000.0, 1100.0, 1200.0, 1300.0):
        for g3 in (7.0, 8.0, 9.0, 10.0):
            cv = dict(base, G2=g2, G3=g3)
            rec = dict(G2=g2, G3=g3, grains_per_ear=[], grain_w_mg=[],
                       yield_kg=[], yield_obs=[])
            for plot in dc.ORDER:
                g, s = dc.run_plot(plot, cv)
                m = dc.plot_metrics(g, s, measured_density(plot))
                rec["grains_per_ear"].append(m["grains_per_ear"])
                rec["grain_w_mg"].append(m["grain_w_mg"])
                rec["yield_kg"].append(m["yield_kg"])
                rec["yield_obs"].append(y.loc[plot, "yield_kg_ha"] * 0.86)
            gp = np.array(rec["grains_per_ear"])
            gw = np.array(rec["grain_w_mg"])
            ys = np.array(rec["yield_kg"])
            yo = np.array(rec["yield_obs"])
            grid.append(dict(
                G2=g2, G3=g3,
                grains_per_ear=gp.mean(), grain_w_mg=gw.mean(),
                yield_rrmse=np.sqrt(np.mean((ys - yo) ** 2)) / yo.mean() * 100,
                score=(abs(gp.mean() - g_tgt) / g_tgt
                       + abs(gw.mean() - w_tgt) / w_tgt
                       + np.sqrt(np.mean((ys - yo) ** 2)) / yo.mean())))
    r = pd.DataFrame(grid).sort_values("score")
    r.to_csv(DIAG / "grain_refit.csv", index=False)
    print(f"\n  targets: grains/ear {g_tgt:.0f}, 100-grain weight {w_tgt:.0f} mg")
    print(r.head(5).round(3).to_string(index=False))


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------
def figure() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["SimHei"]
    plt.rcParams["axes.unicode_minus"] = False

    dens = pd.read_csv(DIAG / "density.csv")
    la = pd.read_csv(DIAG / "leafarea.csv")
    bio = pd.read_csv(DIAG / "biomass_hi.csv")
    gr = pd.read_csv(DIAG / "grain.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    x = np.arange(len(dens))
    ax.bar(x - 0.2, dens.meas_yield_kg, 0.4, label="实测密度", color="#4C72B0")
    ax.bar(x + 0.2, dens.nom_yield_kg, 0.4, label="额定密度", color="#DD8452")
    ax.scatter(x, dens.obs_yield, color="k", zorder=3, label="实测产量")
    ax.set_xticks(x)
    ax.set_xticklabels(dens["plot"], rotation=45)
    ax.set_ylabel("产量 kg/ha")
    ax.set_title("A  密度：实测留苗 vs 额定播种")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    obs = la[la.source == "observed"]
    ax.scatter(obs.density, obs.la_per_plant, color="k", zorder=3,
               label="实测")
    for tag, sub in la[la.source != "observed"].groupby("source"):
        ax.plot(sub.density, sub.la_per_plant, "o-", ms=3, lw=1, label=tag)
    ax.set_xlabel("密度 株/m2")
    ax.set_ylabel("峰值单株叶面积 m2/株")
    ax.set_title("B  单株叶面积：模型上限 vs 实测")
    ax.legend(fontsize=7)

    ax = axes[1, 0]
    agg = bio.groupby("scenario").agg(bio=("bio", "mean"),
                                      obs=("obs_bio", "mean"),
                                      hi=("hi", "mean"))
    order = agg.sort_values("bio").index
    ax.barh(np.arange(len(order)), agg.loc[order, "bio"], color="#4C72B0")
    ax.axvline(agg.obs.iloc[0], color="k", ls="--", label="实测生物量")
    for i, s in enumerate(order):
        ax.text(agg.loc[s, "bio"] + 200, i, f"HI {agg.loc[s,'hi']:.3f}",
                va="center", fontsize=7)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("成熟期地上部生物量 kg/ha")
    ax.set_title("C  生物量与收获指数：各情景")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    x = np.arange(len(gr))
    ax.bar(x - 0.2, gr.meas_grains_per_ear, 0.4, label="模拟 穗粒数",
           color="#4C72B0")
    if "obs_grains_per_ear" in gr:
        ax.bar(x + 0.2, gr.obs_grains_per_ear, 0.4, label="考种 穗粒数",
               color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels(gr["plot"], rotation=45)
    ax.set_ylabel("穗粒数 粒/穗")
    ax2 = ax.twinx()
    ax2.plot(x, gr.meas_grain_w_mg, "s-", color="#55A868", label="模拟百粒重")
    if "obs_hundred_gw_mg" in gr:
        ax2.plot(x, gr.obs_hundred_gw_mg, "o--", color="#C44E52",
                 label="考种百粒重")
    ax2.set_ylabel("百粒重 mg")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7)
    ax.set_title("D  穗粒数与百粒重")

    fig.suptitle(f"{dc.CFG.site} {dc.CFG.variety} — DSSAT 诊断", fontsize=14)
    fig.tight_layout()
    out = dc.FIG / "diag" / "E_诊断.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print(f"figure -> {out}")


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF")
    ap.add_argument("--section", default="all",
                    choices=["all", "density", "leafarea", "biomass", "grain",
                             "figure"])
    args = ap.parse_args()
    dc.set_trial(args.trial)
    global DIAG
    DIAG = dc.OUT / "diag"
    DIAG.mkdir(parents=True, exist_ok=True)
    print(f"trial {dc.TRIAL}  diag -> {DIAG}")

    if args.section in ("all", "density"):
        print("\n[Q1] density: measured vs nominal")
        q1_density()
    if args.section in ("all", "leafarea"):
        print("\n[Q2] per-plant leaf area")
        q2_leafarea()
    if args.section in ("all", "biomass"):
        print("\n[Q3] biomass and harvest index")
        q3_biomass()
    if args.section in ("all", "grain"):
        print("\n[Q4] grain number and grain weight")
        q4_grain()
    if args.section in ("all", "figure"):
        figure()


if __name__ == "__main__":
    main()
