"""Step 09 - build the observed leaf-appearance and stage series used to
constrain the phenology calibration.

Two data sources, different in structure:

  德州 (试验田_*.xlsx)  : daily records of 展开叶数量 / 可见叶数量 with a BBCH code
  新乡 (sheet per variety): BBCH-coded stage list with a date, plus the free-text
                            '生育期' field that often states the visible leaf count

The output is one row per (site, variety, date) with the constrained quantities:
  leaves_expanded, leaves_visible, stage_cn, bbch
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DS, RAW  # noqa: E402

DEZHOU = RAW / "原始记录" / "2026年德州玉米试验数据汇总.xlsx"
XINXIANG = RAW / "原始记录" / "2026新乡玉米试验数据汇总.xlsx"

# excel serial for 2026-06-19 etc. -> convert with origin 1899-12-30
EXCEL_ORIGIN = pd.Timestamp("1899-12-30")

CN_NUM = {"一": 1, "二": 2, "三": 3, "四": 4, "五": 5, "六": 6, "七": 7,
          "八": 8, "九": 9, "十": 10, "十一": 11, "十二": 12,
          "十三": 13, "十四": 14, "十五": 15, "十六": 16, "十七": 17,
          "十八": 18, "十九": 19, "二十": 20, "二十一": 21, "二十二": 22}


def parse_cn_leaves(text: str):
    """'6-7 叶一心, 八张可见叶' -> (expanded 6.5, visible 8)."""
    if not isinstance(text, str):
        return None, None
    visible = None
    m = re.search(r"([一二三四五六七八九十]+)\s*[张个]?可见", text)
    if m:
        visible = CN_NUM.get(m.group(1))
    expanded = None
    m = re.search(r"([一二三四五六七八九十\d]+)(?:\s*[-~到至]\s*([一二三四五六七八九十\d]+))?\s*叶一心", text)
    if m:
        def v(x):
            return CN_NUM.get(x, None) if x in CN_NUM else (int(x) if str(x).isdigit() else None)
        a, b = v(m.group(1)), v(m.group(2)) if m.group(2) else None
        expanded = a if b is None else (a + b) / 2.0 if a and b else a
    return expanded, visible


def excel_date(v):
    if pd.isna(v):
        return pd.NaT
    if isinstance(v, (pd.Timestamp,)):
        return v
    s = str(v).strip()
    if s.isdigit() and len(s) == 5:
        return EXCEL_ORIGIN + pd.Timedelta(days=int(s))
    m = re.match(r"(\d{1,2})月(\d{1,2})[日号]?", s)
    if m:
        return pd.Timestamp(2026, int(m.group(1)), int(m.group(2)))
    return pd.to_datetime(s, errors="coerce")


def read_dezhou() -> pd.DataFrame:
    rows = []
    xl = pd.ExcelFile(DEZHOU)
    for sheet in xl.sheet_names:
        if not sheet.startswith("试验田_"):
            continue
        variety = sheet.replace("试验田_", "")
        d = pd.read_excel(DEZHOU, sheet_name=sheet, header=0)
        d.columns = [str(c).strip() for c in d.columns]
        if "展开叶数量" not in d.columns:
            continue
        for _, r in d.iterrows():
            date = excel_date(r["日期"])
            if pd.isna(date):
                continue
            exp_ = pd.to_numeric(r.get("展开叶数量"), errors="coerce")
            vis = pd.to_numeric(r.get("可见叶数量"), errors="coerce")
            if pd.notna(exp_) or pd.notna(vis):
                rows.append(dict(site="山东德州", variety=variety, date=date,
                                 leaves_expanded=exp_, leaves_visible=vis,
                                 stage_cn=str(r.get("物候期阶段", "")).strip(),
                                 bbch=pd.to_numeric(r.get("BBCH"), errors="coerce")))
    return pd.DataFrame(rows)


def read_xinxiang() -> pd.DataFrame:
    rows = []
    xl = pd.ExcelFile(XINXIANG)
    for sheet in xl.sheet_names:
        d = pd.read_excel(XINXIANG, sheet_name=sheet, header=0)
        d.columns = [str(c).strip() for c in d.columns]
        # the 金华创826 sheet has two side-by-side blocks
        blocks = [(slice(0, 5), 0)]
        if d.shape[1] >= 12 and str(d.columns[7]).strip() == "品种":
            blocks.append((slice(7, 12), 7))
        for cols, off in blocks:
            sub = d.iloc[:, cols].copy()
            sub.columns = ["variety", "bbch", "stage", "desc", "date"][: sub.shape[1]]
            for _, r in sub.iterrows():
                var = str(r["variety"]).strip()
                if not var or var.lower() == "nan":
                    continue
                date = pd.to_datetime(r["date"], errors="coerce")
                if pd.isna(date):
                    continue
                exp_, vis = parse_cn_leaves(str(r.get("desc", "")))
                rows.append(dict(site="河南新乡", variety=var, date=date,
                                 leaves_expanded=exp_, leaves_visible=vis,
                                 stage_cn=str(r.get("stage", "")).strip(),
                                 bbch=pd.to_numeric(r["bbch"], errors="coerce")))
    return pd.DataFrame(rows)


def main() -> None:
    parts = [x for x in (read_dezhou(), read_xinxiang()) if len(x)]
    df = pd.concat(parts, ignore_index=True)
    df = df.sort_values(["site", "variety", "date"]).reset_index(drop=True)

    # keep only the varieties that are actually simulated
    keep = ["乐农87", "MY73", "金华创826"]
    df["keep"] = df["variety"].isin(keep)
    df.to_csv(DS / "leaf_series.csv", index=False)

    print(f"leaf_series {df.shape}  sites={df['site'].unique().tolist()}")
    for (site, var), g in df[df["keep"]].groupby(["site", "variety"]):
        g2 = g.dropna(subset=["leaves_visible"])
        print(f"\n{site} {var}: {len(g)} records, "
              f"{g2['leaves_visible'].notna().sum()} with visible leaf count")
        if len(g2):
            print("  visible leaves:",
                  ", ".join(f"{d:%m-%d}={v:.0f}" for d, v in
                            zip(g2["date"], g2["leaves_visible"])))
    bad = df[df["keep"] & df["leaves_visible"].isna() & df["leaves_expanded"].isna()]
    print(f"\nrecords without any leaf count: {len(bad)}/{int(df['keep'].sum())}")


if __name__ == "__main__":
    main()
