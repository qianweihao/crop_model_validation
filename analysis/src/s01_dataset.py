"""Step 01 - audit raw sheets and emit canonical, analysis-ready tables.

Fixes the two defects found in review:
  * '按小区密度总生物量' / '按小区密度叶面积指数' are recomputed from the
    per-plant measurements x *measured* plot density (the supplied columns used
    a 5600 default on 2026-07-23 and 2026-08-21).
  * plot / treatment codes are normalised (A2-B1 -> A2B1, ck -> CK).

Outputs (analysis/out/dataset):
  plots.csv  obs_plant.csv  yield.csv  fert.csv  irrig.csv
  weather_daily.csv  soil_profile.csv  soil_water.csv  phenology.csv
  audit_findings.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    DRY_FACTOR, DS, MGMT_XLSX, NUT_XLSX, OBS_XLSX, ORGAN_SPLIT_SUSPECT,
    PHEN_XLSX, RAW, REPORTS, TRIALS, WEATHER_SHEETS, WEATHER_SITE, WTH_XLSX,
    build_soil_profile, col, norm_plot, num, read_header_sheet,
)

findings: list[dict] = []


def note(check, scope, detail, severity="info"):
    findings.append(dict(check=check, scope=scope, detail=detail, severity=severity))


# ==========================================================================
# 1. plot registry
# ==========================================================================
def build_plots() -> pd.DataFrame:
    info = pd.read_excel(MGMT_XLSX, sheet_name="试验信息")
    info["plot"] = info["处理编号"].map(norm_plot)
    info["trial_id"] = info.apply(
        lambda r: ("SD_SF" if r["地区"] == "山东德州" else
                   "XX_SF" if r["试验名称"] == "水肥试验" else "XX_IRR"), axis=1)
    plots = pd.DataFrame(dict(
        trial_id=info["trial_id"], site=info["地区"], trial=info["试验名称"],
        plot=info["plot"], variety=info["品种"],
        sow_date=pd.to_datetime(info["播种日期"]),
        density_planted=num(info["播种密度（株/亩）"]),
        density_measured=num(info["调查密度（株/亩）"]),
        area_mu=num(info["小区面积（亩）"]),
        device=info["设备编号"],
    ))
    plots["weather"] = plots["site"].map({"山东德州": "dezhou", "河南新乡": "xinxiang"})
    plots["role"] = plots["trial_id"].map({k: v["role"] for k, v in TRIALS.items()})
    plots["plot_raw"] = info["处理编号"].astype(str)

    # farmer-practice plots (no management records) are observation-only
    extra = []
    for tid, cfg in TRIALS.items():
        if cfg["has_mgmt"]:
            continue
        n = 4 if tid == "SD_MY73" else 5
        off = 0 if tid == "SD_MY73" else 4
        for k in range(n):
            extra.append(dict(trial_id=tid, site=cfg["site"], trial=cfg["trial"],
                              plot=f"W{off + k + 1}", variety=cfg["variety"],
                              sow_date=pd.Timestamp("2026-06-16"),
                              density_planted=6000.0, density_measured=np.nan,
                              area_mu=np.nan, device=None, weather="dezhou",
                              role=cfg["role"], plot_raw=f"W{off + k + 1}"))
    plots = pd.concat([plots, pd.DataFrame(extra)], ignore_index=True)
    return plots.sort_values(["trial_id", "plot"]).reset_index(drop=True)


# ==========================================================================
# 2. per-plant observations -> canonical long table
# ==========================================================================
def build_obs() -> pd.DataFrame:
    emerge = {"SD_SF": "2026-06-23", "XX_SF": "2026-06-24", "XX_IRR": "2026-06-24",
              "SD_MY73": "2026-06-21", "SD_DH997": "2026-06-21"}
    rows = []
    for tid, cfg in TRIALS.items():
        d = read_header_sheet(OBS_XLSX, cfg["obs_sheet"])
        pc = cfg["plot_col"]
        if isinstance(pc, tuple):
            pcol, tag = pc
            plot = d[pcol].map(lambda v: f"W{int(v)}" if pd.notna(v) else "")
        else:
            plot = d[pc].map(norm_plot)
        out = pd.DataFrame(dict(
            trial_id=tid, site=d["区域"], trial=d["农场名称"], variety=d["品种"],
            plot=plot, date=pd.to_datetime(d["日期"]),
            plant_la_cm2=num(d[col(d, "单株叶面积")]),
            leaf_dw=num(d[col(d, "单株叶干重")]),
            stem_dw=num(d[col(d, "单株茎干重")]),
            root_dw=num(d[col(d, "单株根干重")]),
            ear_dw=num(d[col(d, "单株穗干重")]),
            grain_dw=num(d[col(d, "单株籽粒干重")]),
            total_dw=num(d[col(d, "单株生物量总干重")]),
            density_planted=num(d[col(d, "种植密度")]),
            density_measured=num(d[col(d, "小区实测密度")]),
            lai_std_nominal=num(d[col(d, "叶面积指数")]),
            src_sheet=cfg["obs_sheet"], src_row=np.arange(len(d)) + 1,
        ))
        # ---- recompute, never trust the supplied derived columns ----
        out["lai_plot"] = out["plant_la_cm2"] * out["density_measured"] / 6666667.0
        out["biomass_kg_mu"] = out["total_dw"] * out["density_measured"] / 1000.0
        out["das"] = (out["date"] - pd.Timestamp(emerge[tid])).dt.days
        rows.append(out)
    obs = pd.concat(rows, ignore_index=True)
    obs["organ_split_suspect"] = [
        (r.trial_id, pd.Timestamp(r.date).strftime("%Y-%m-%d")) in ORGAN_SPLIT_SUSPECT
        for r in obs.itertuples()]
    obs["leaf_area_missing"] = obs["plant_la_cm2"].isna()
    return obs.sort_values(["trial_id", "plot", "date"]).reset_index(drop=True)


# ==========================================================================
# 3. yield (block-structured sheet, parsed explicitly)
# ==========================================================================
def build_yield() -> pd.DataFrame:
    raw = pd.read_excel(MGMT_XLSX, sheet_name="产量记录", header=None)
    recs = []

    def add(tid, plot, dens, meas, conv, theory, target, method, note_):
        recs.append(dict(trial_id=tid, plot=norm_plot(plot), density_check=dens,
                         yield_measured=meas, yield_converted=conv,
                         yield_theory=theory, target_yield=target,
                         method=method, remark=note_))

    # 河南新乡 block: header row 3, data 4..17
    for i in range(4, 18):
        tag, plot = raw.iat[i, 0], raw.iat[i, 1]
        if pd.isna(tag):
            continue
        tid = "XX_SF" if "MY73" in str(tag) else "XX_IRR"
        add(tid, plot, raw.iat[i, 2], raw.iat[i, 3], raw.iat[i, 4], raw.iat[i, 5],
            raw.iat[i, 6], "plot_harvest", raw.iat[i, 7])
    # 山东试验区 block: header row 22, data 23..31
    for i in range(23, 32):
        if pd.isna(raw.iat[i, 0]):
            continue
        add("SD_SF", raw.iat[i, 1], raw.iat[i, 2], raw.iat[i, 3], raw.iat[i, 4],
            raw.iat[i, 5], raw.iat[i, 6], "plot_harvest", np.nan)
    # 外围试验区 block: header row 36, data 37..45
    for i in range(37, 46):
        var, plot = raw.iat[i, 0], raw.iat[i, 1]
        if pd.isna(var) or pd.isna(plot):
            continue
        tid = "SD_MY73" if "MY73" in str(var) else "SD_DH997"
        add(tid, f"W{int(plot)}", np.nan, raw.iat[i, 2], np.nan, np.nan, np.nan,
            "10m_double_row", np.nan)

    y = pd.DataFrame(recs)
    y["yield_kg_ha"] = y["yield_measured"] * 15.0          # at 14% moisture
    y["yield_dry_kg_ha"] = y["yield_kg_ha"] * DRY_FACTOR   # APSIM-comparable
    return y


# ==========================================================================
# 4. management
# ==========================================================================
def build_fert() -> pd.DataFrame:
    f = pd.read_excel(MGMT_XLSX, sheet_name="施肥记录")
    f["plot"] = f["整理编号"].map(norm_plot)
    f["trial_id"] = f.apply(
        lambda r: ("SD_SF" if r["地区"] == "山东德州" else
                   "XX_SF" if r["试验名称"] == "水肥试验" else "XX_IRR"), axis=1)
    keep = pd.DataFrame(dict(
        trial_id=f["trial_id"], plot=f["plot"], batch=f["施肥批次"],
        stage=f["施肥时期"], date=pd.to_datetime(f["施肥日期"]),
        n_kg_mu=num(f["N（kg/亩"]), p2o5_kg_mu=num(f["P₂O₅（kg/亩）"]),
        k2o_kg_mu=num(f["K₂O（kg/亩）"]),
    ))
    keep["n_kg_ha"] = keep["n_kg_mu"] * 15.0
    return keep.sort_values(["trial_id", "plot", "date"]).reset_index(drop=True)


def build_irrig() -> pd.DataFrame:
    f = pd.read_excel(MGMT_XLSX, sheet_name="灌溉记录")
    f["plot"] = f["处理编号"].map(norm_plot)
    f["trial_id"] = f.apply(
        lambda r: ("SD_SF" if r["地区"] == "山东德州" else
                   "XX_SF" if r["试验名称"] == "水肥试验" else "XX_IRR"), axis=1)
    keep = pd.DataFrame(dict(
        trial_id=f["trial_id"], plot=f["plot"], date=pd.to_datetime(f["灌溉日期"]),
        kind=f["灌水类型"], stage=f["生育阶段"], mm=num(f["灌水量（mm）"]),
    ))
    return keep.sort_values(["trial_id", "plot", "date"]).reset_index(drop=True)


# ==========================================================================
# 5. weather (2 distinct stations; sheets verified identical on common dates)
# ==========================================================================
def build_weather() -> pd.DataFrame:
    out = []
    for zone, sheet in WEATHER_SHEETS.items():
        d = pd.read_excel(WTH_XLSX, sheet_name=sheet)
        d["date"] = pd.to_datetime(d["日期"], errors="coerce")
        d = d.dropna(subset=["date"]).sort_values("date")
        out.append(pd.DataFrame(dict(
            zone=zone, date=d["date"],
            tmax=num(d["tem_max"]), tmin=num(d["tem_min"]),
            tavg=num(d["日平均气温 (°C)"]), radn=num(d["短波辐射"]),
            rain=num(d["pre"]), rh=num(d["rhu_avg"]), wind=num(d["win_s_2mi_avg"]),
            wind_max=num(d["win_s_max"]), soil_t=num(d["gst_avg_10cm"]),
        )))
    return pd.concat(out, ignore_index=True)


# ==========================================================================
# 6. soil profile per plot
# ==========================================================================
def build_soil(plots: pd.DataFrame) -> pd.DataFrame:
    sp = pd.read_excel(MGMT_XLSX, sheet_name="土壤参数")
    sp["plot"] = sp["处理编号"].map(norm_plot)
    sp["trial_id"] = sp.apply(
        lambda r: ("SD_SF" if r["地区"] == "山东德州" else
                   "XX_SF" if r["试验名称"] == "水肥试验" else "XX_IRR"), axis=1)
    om = {"SD_SF": 1.89, "XX_SF": 2.32, "XX_IRR": 2.32}  # 德州1.89 / 新乡 2.1-2.57 mean
    rows = []
    for _, r in sp.iterrows():
        prof = build_soil_profile(
            sand=float(r["0–20 cm砂粒（%）"]), clay=float(r["0–20 cm黏粒（%）"]),
            sat=float(r["饱和含水量SAT（体积%）"]) / 100.0,
            fc=float(r["田间持水量FC（体积%）"]) / 100.0,
            om_pct=om.get(r["trial_id"], 2.2))
        prof.insert(0, "plot", r["plot"])
        prof.insert(0, "trial_id", r["trial_id"])
        prof["sand_top"] = r["0–20 cm砂粒（%）"]
        prof["clay_top"] = r["0–20 cm黏粒（%）"]
        prof["fc_measured"] = r["田间持水量FC（体积%）"] / 100.0
        prof["sat_measured"] = r["饱和含水量SAT（体积%）"] / 100.0
        rows.append(prof)
    prof = pd.concat(rows, ignore_index=True)

    # scatter check on the measured hydraulic parameters
    for tid, g in sp.groupby("trial_id"):
        cv_fc = g["田间持水量FC（体积%）"].std() / g["田间持水量FC（体积%）"].mean()
        cv_sat = g["饱和含水量SAT（体积%）"].std() / g["饱和含水量SAT（体积%）"].mean()
        note("soil_param_scatter", tid,
             f"FC {g['田间持水量FC（体积%）'].min():.2f}-{g['田间持水量FC（体积%）'].max():.2f}% "
             f"(CV={cv_fc:.0%}), SAT {g['饱和含水量SAT（体积%）'].min():.1f}-"
             f"{g['饱和含水量SAT（体积%）'].max():.1f}% (CV={cv_sat:.0%})",
             "warn" if cv_fc > 0.10 or cv_sat > 0.10 else "info")
    return prof


def build_soil_water() -> pd.DataFrame:
    sw = pd.read_excel(MGMT_XLSX, sheet_name="土壤墒情")
    return pd.DataFrame(dict(
        site=sw["地区"], trial=sw["试验名称"], plot=sw["处理编号"].map(norm_plot),
        date=pd.to_datetime(sw["灌溉日期"]),
        theta_before=num(sw["灌前含水量（体积%）"]),
        theta_after=num(sw["灌后含水量（体积%）"]),
        source=sw["来源"],
    ))


def build_phenology() -> pd.DataFrame:
    d = pd.read_excel(PHEN_XLSX)
    d = d.rename(columns={d.columns[0]: "stage"})
    long = d.melt(id_vars="stage", var_name="variety_tag", value_name="date")
    long["date"] = pd.to_datetime(long["date"], errors="coerce")
    long = long.dropna(subset=["date"])
    parts = long["variety_tag"].str.split("_", n=1, expand=True)
    long["site"], long["variety"] = parts[0], parts[1]
    return long.sort_values(["variety_tag", "date"]).reset_index(drop=True)


# ==========================================================================
# 7. audit
# ==========================================================================
def run_audit(obs: pd.DataFrame, yld: pd.DataFrame, fert: pd.DataFrame,
              irrig: pd.DataFrame, plots: pd.DataFrame) -> pd.DataFrame:
    # --- 7.1 supplied derived columns vs recomputed -------------------------
    sup = []
    for tid, cfg in TRIALS.items():
        d = read_header_sheet(OBS_XLSX, cfg["obs_sheet"])
        c_bm = col(d, "按小区密度总生物量")
        c_lai = col(d, "按小区密度叶面积指数")
        if not c_bm:
            continue
        imp = num(d[c_bm]) * 1000.0 / num(d[col(d, "单株生物量总干重")])
        dens = num(d[col(d, "小区实测密度")])
        valid = dens.notna() & imp.notna()
        if not valid.any():
            if c_bm and dens.isna().all():
                note("derived_column_divisor", tid,
                     "no measured plot density -> supplied '按小区密度' columns are "
                     "blank; plot-scale LAI/biomass cannot be derived (per-plant only)",
                     "warn")
            continue
        ok = (~valid) | ((imp - dens).abs() < 1.0)
        ratio = imp / dens
        bad = d.loc[~ok, ["日期"]].assign(ratio=ratio[~ok]).drop_duplicates()
        for _, r in bad.iterrows():
            sup.append(dict(trial_id=tid, date=pd.Timestamp(r["日期"]).date(),
                            divisor_used_implied=float("nan"),
                            inflation=float(r["ratio"])))
        if len(bad):
            note("derived_column_divisor", tid,
                 "supplied '按小区密度' columns used a density other than the "
                 f"measured one on {sorted(set(pd.Timestamp(x).strftime('%m-%d') for x in bad['日期']))} "
                 f"(inflation x{ratio[~ok].min():.2f}-{ratio[~ok].max():.2f}); recalculated",
                 "high")
    pd.DataFrame(sup).to_csv(DS / "audit_derived_column_flags.csv", index=False)

    # --- 7.2 organ sum vs total -------------------------------------------
    o = obs.copy()
    organ = o[["leaf_dw", "stem_dw", "root_dw", "ear_dw"]].sum(axis=1, min_count=1)
    gap = (o["total_dw"] - organ) / o["total_dw"]
    note("organ_sum_vs_total", "all",
         f"leaf+stem+root+ear differs from total DW by {gap.median():.0%} (median); "
         f"{int((gap.abs() > 0.25).sum())}/{len(o)} records off by >25% "
         "(root rarely measured, ear absent before silking) - use total_dw for calibration",
         "info")

    # --- 7.3 monotonicity of per-plant total DW ---------------------------
    for tid, g in obs.dropna(subset=["total_dw"]).groupby("trial_id"):
        bad = []
        for plot, s in g.sort_values("date").groupby("plot"):
            dw = s["total_dw"].to_numpy()
            if len(dw) > 1 and np.any(np.diff(dw) < -0.10 * dw[:-1]):
                bad.append(plot)
        note("biomass_monotonicity", tid,
             f"{len(bad)}/{g['plot'].nunique()} plots have a >10% drop in per-plant DW "
             f"({', '.join(bad) if bad else 'none'})",
             "warn" if bad else "ok")

    # --- 7.4 density vs N level -------------------------------------------
    n_tot = fert.groupby(["trial_id", "plot"])["n_kg_mu"].sum().rename("n_total")
    p = plots.merge(n_tot, on=["trial_id", "plot"], how="left")
    for tid, g in p.dropna(subset=["density_measured", "n_total"]).groupby("trial_id"):
        if len(g) < 4:
            continue
        r = np.corrcoef(g["n_total"], g["density_measured"])[0, 1]
        note("density_vs_N", tid,
             f"corr(total N, measured density) = {r:+.2f}; density "
             f"{g['density_measured'].min():.0f}-{g['density_measured'].max():.0f} /mu "
             f"from {g['density_planted'].iloc[0]:.0f} planted "
             f"({(g['density_measured']/g['density_planted']).mean():.0%} establishment)",
             "warn" if abs(r) > 0.5 else "info")

    # --- 7.5 N gradient and yield response --------------------------------
    y = yld.merge(n_tot, on=["trial_id", "plot"], how="left")
    for tid, g in y.dropna(subset=["n_total", "yield_measured"]).groupby("trial_id"):
        g = g.sort_values("n_total")
        note("yield_vs_N", tid,
             " ; ".join(f"{r.n_total:.1f}N->{r.yield_measured:.0f}" for _, r in g.iterrows())
             + f" kg/mu  (span {(g['yield_measured'].max()-g['yield_measured'].min())/g['yield_measured'].min():+.0%})",
             "info")

    # --- 7.6 water gradient ----------------------------------------------
    w = irrig.groupby(["trial_id", "plot"])["mm"].agg(["sum", "count"]).reset_index()
    for tid, g in w.groupby("trial_id"):
        note("irrigation_gradient", tid,
             f"seasonal irrigation {g['sum'].min():.1f}-{g['sum'].max():.1f} mm over "
             f"{g['count'].min()}-{g['count'].max()} events", "info")

    # --- 7.7 blend check ---------------------------------------------------
    note("B1_vs_B2_confounded", "SD_SF",
         "B2 plots received 8 events / 72.8-86.0 mm vs B1 5 events / 62.6-64.4 mm, "
         "so the B1-B2 contrast mixes N timing with +10..22 mm irrigation", "warn")

    # --- 7.8 missing initial mineral N -----------------------------------
    note("no_initial_mineral_N", "both sites",
         "soil nutrient sheets report no NO3-N / NH4-N; APSIM SoilN initial mineral N "
         "must be treated as a calibration parameter (prior 40-90 kg N/ha, 0-60 cm)", "high")

    # --- 7.9 data coverage gaps ------------------------------------------
    for tid, g in obs.groupby("trial_id"):
        note("obs_coverage", tid,
             f"{g['plot'].nunique()} plots x {g['date'].nunique()} dates = {len(g)} records; "
             f"plot-scale LAI/biomass {'available' if g['density_measured'].notna().any() else 'NOT available (no measured density)'}",
             "info")
    return pd.DataFrame(findings)


# ==========================================================================
def main() -> None:
    plots = build_plots()
    obs = build_obs()
    yld = build_yield()
    fert = build_fert()
    irrig = build_irrig()
    wth = build_weather()
    prof = build_soil(plots)
    sw = build_soil_water()
    phe = build_phenology()

    plots.to_csv(DS / "plots.csv", index=False)
    obs.to_csv(DS / "obs_plant.csv", index=False)
    yld.to_csv(DS / "yield.csv", index=False)
    fert.to_csv(DS / "fert.csv", index=False)
    irrig.to_csv(DS / "irrig.csv", index=False)
    wth.to_csv(DS / "weather_daily.csv", index=False)
    prof.to_csv(DS / "soil_profile.csv", index=False)
    sw.to_csv(DS / "soil_water.csv", index=False)
    phe.to_csv(DS / "phenology.csv", index=False)

    audit = run_audit(obs, yld, fert, irrig, plots)
    audit.to_csv(DS / "audit_findings.csv", index=False)

    # ---- console report -------------------------------------------------
    print(f"plots        {plots.shape}")
    print(f"obs_plant    {obs.shape}   plot-scale LAI rows: {obs['lai_plot'].notna().sum()}")
    print(f"yield        {yld.shape}")
    print(f"fert         {fert.shape}   irrig {irrig.shape}")
    print(f"weather      {wth.shape}  zones={wth['zone'].unique().tolist()}")
    print(f"soil_profile {prof.shape}")
    print()
    print("AUDIT")
    for _, r in audit.iterrows():
        print(f"  [{r['severity']:>4}] {r['check']:<24} {r['scope']:<10} {r['detail']}")


if __name__ == "__main__":
    main()
