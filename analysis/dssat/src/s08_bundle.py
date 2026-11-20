"""Build a self-contained, runnable DSSAT bundle for one trial.

The bundle contains everything DSSAT needs except the model executable:

    README.md              how to run / how to import into DSSAT
    run.py, run.bat        one-click runner (writes DSSATPRO.V48, runs, collects)
    DSSATPRO.V48           template, rewritten by run.py with absolute paths
    MZCER048.CUL/.SPE/.ECO genotype files (CUL carries the calibrated cultivar)
    SOIL.SOL               all plot profiles in one file
    Weather/<INSI>2601.WTH
    <trial>_<plot>.MZX     one management file per plot
    _support/              StandardData, Soil, Pest (so only the exe is needed)

Usage
-----
    python analysis/dssat/src/s08_bundle.py --trial SD_SF
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dssat_common as dc  # noqa: E402

SUPPORT = ["StandardData", "Soil", "Pest"]


def _merge_soil(paths: list[Path], out: Path) -> int:
    """Concatenate single-profile SOIL.SOL files into one multi-profile file."""
    head, blocks = None, []
    for p in paths:
        text = p.read_text(encoding="cp437", errors="replace").splitlines()
        if head is None:
            head = text[0]                      # "*SOILS: File created with ..."
        start = next(i for i, l in enumerate(text) if l.startswith("*") and
                     not l.startswith("*SOILS"))
        blocks.append("\n".join(text[start:]).rstrip())
    out.write_text("\n\n".join([head, ""] + blocks) + "\n", encoding="cp437")
    return len(blocks)


OUT_README = """# {trial} — {site} {variety}（{nplot} 个小区）

## 想跑模拟 → 进 `{bundle}/`

```bash
cp -r {bundle} C:/dssat        # 一定要复制到短路径，见下
cd C:/dssat
python run.py --exe <dscsm048.exe>
```

`{bundle}/` 是自包含的（MZX + SOIL.SOL + WTH + 品种参数 + 支持文件），
里面有单独的 README。

⚠ DSSAT 的路径变量是定长的，超过约 60 字符会被静默截断，报
`File not found ... SIMCDE`。这个目录的路径本身就很长，**不要就地运行**。

## 其余目录都是"结果/证据"，不能直接跑

| 目录 | 内容 |
|---|---|
| `{bundle}/` | **可运行**的完整 DSSAT 文件包 |
| `sim/` | 实测留苗密度的 DSSAT 输出（`Summary.OUT`、`PlantGro.OUT`）+ 验证 CSV |
| `sim_nominal/` | 额定密度情景的输出 |
| `mgmt/` | 生成的 `.MZX`（文件名 >8 字符、且没有品种文件，单独跑会失败） |
| `mgmt_nominal/` | 额定密度情景的 `.MZX` |
| `soil/` | 单小区 `.SOL` |
| `weather/` | `.WTH` |
| `diag/` | 四个诊断问题的 CSV |

数据来源与标定过程见 `analysis/dssat/report/`。
"""

RUN_PY = '''"""Run every DSSAT experiment file in this folder and collect the outputs.

    python run.py [--exe <dscsm048.exe>] [--keep]

Needs only the DSSAT 4.8 model executable; the support files are bundled in
_support/.  DSSAT uses fixed-length path buffers, so keep this folder at a
short path (e.g. C:\\\\dssat) - run.py warns if the path is too long.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAXPATH = 60          # DSSAT truncates longer paths silently


def find_exe(given: str | None) -> Path:
    if given:
        return Path(given).expanduser().resolve()
    env = os.environ.get("DSCSM048")
    if env:
        return Path(env)
    for c in (HERE / "_support" / "bin" / "dscsm048.exe",
              HERE / "_support" / "bin" / "dscsm048",
              Path(r"C:/DSSAT48/DSCSM048.EXE"),
              Path(r"C:/DSSAT48/bin/dscsm048.exe"),
              Path(r"C:/DSSAT48/Tools/dscsm048.exe")):
        if c.exists():
            return c
    found = shutil.which("dscsm048") or shutil.which("dscsm048.exe")
    if found:
        return Path(found)
    sys.exit("cannot find dscsm048; pass --exe <path to dscsm048.exe>")


def write_config() -> None:
    (HERE / "DSSATPRO.V48").write_text(
        f"WED    {HERE / 'Weather'}\\n"
        f"MMZ    {HERE} dscsm048 MZCER048\\n"
        f"CRD    {HERE}\\n"
        f"PSD    {HERE / '_support' / 'Pest'}\\n"
        f"SLD    {HERE / '_support' / 'Soil'}\\n"
        f"STD    {HERE / '_support' / 'StandardData'}\\n",
        encoding="ascii")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=None, help="path to dscsm048.exe")
    ap.add_argument("--keep", action="store_true",
                    help="keep the raw *.OUT files next to the .MZX")
    args = ap.parse_args()

    if len(str(HERE)) > MAXPATH:
        print(f"WARNING: {HERE} is {len(str(HERE))} characters long; "
              f"DSSAT silently truncates paths above ~{MAXPATH}. "
              f"Copy this folder to e.g. C:\\\\dssat and run it there.")
    exe = find_exe(args.exe)
    print(f"model : {exe}")
    print(f"folder: {HERE}")
    write_config()

    results = HERE / "results"
    results.mkdir(exist_ok=True)
    files = sorted(HERE.glob("*.MZX"))
    if not files:
        sys.exit("no .MZX files found")

    ok = 0
    for mzx in files:
        tag = mzx.stem
        for f in HERE.glob("*.OUT"):
            f.unlink()
        r = subprocess.run([str(exe), "C", mzx.name, "1"], cwd=str(HERE),
                           capture_output=True, text=True)
        summary = HERE / "Summary.OUT"
        if r.returncode == 0 and summary.exists():
            dest = results / tag
            dest.mkdir(exist_ok=True)
            for f in HERE.glob("*.OUT"):
                shutil.copy2(f, dest / f.name)
            if not args.keep:
                for f in HERE.glob("*.OUT"):
                    f.unlink()
            ok += 1
            print(f"  {tag:16s} ok   -> results/{tag}/")
        else:
            print(f"  {tag:16s} FAILED (rc={r.returncode})")
            print((HERE / "ERROR.OUT").read_text(errors="replace")[-400:]
                  if (HERE / "ERROR.OUT").exists() else r.stdout[-400:])
    print(f"\\n{ok}/{len(files)} runs finished; results in {results}")


if __name__ == "__main__":
    main()
'''

README = """# DSSAT bundle — {site} {variety} ({trial})

DSSAT 4.8 / CERES-Maize (`MZCER048`).  {nplot} 个小区，{scenario}。

## 怎么跑

```bash
cd <这个文件夹>
python run.py                 # 自动写 DSSATPRO.V48 并跑全部 .MZX
```

或双击 `run.bat`。需要 `dscsm048.exe`（DSSAT 4.8 自带），默认在
`_support/bin/`、`C:\\DSSAT48\\` 下找，也可以 `python run.py --exe <路径>`。
用 Py_DSSATTools 的话，exe 在 `<python>\site-packages\DSSATTools\static\bin\dscsm048.exe`。

⚠ **路径要短**：DSSAT 的路径变量是定长的，超过约 60 字符会被静默截断，
表现为 `File not found ... SIMCDE` 错误。建议复制到 `C:\\dssat` 再跑，
`run.py` 会检查并提示。

## 文件

| 文件 | 说明 |
|---|---|
| `<plot>.MZX` | 每个小区一份 DSSAT 试验文件（FileX）；文件名 ≤8 字符是 DSSAT 的硬约束 |
| `SOIL.SOL` | {nplot} 个土壤剖面合并成一个文件 |
| `Weather/{wth}` | 气象文件（130 天） |
| `MZCER048.CUL` | 品种参数，含标定后的 `{vrname}` |
| `MZCER048.SPE` / `.ECO` | 物种 / 生态型参数 |
| `_support/` | StandardData / Soil / Pest，DSSAT 运行所需 |
| `results/<plot>/` | 跑完的 `Summary.OUT`、`PlantGro.OUT` 等 |

## 在 DSSAT 图形界面里导入

1. 把 `Weather/{wth}` 复制到 DSSAT 安装目录的 `Weather\\`；
2. 把 `SOIL.SOL` 里的剖面并入 `Soil\\SOIL.SOL`；
3. 把 `MZCER048.CUL` 里 `{vrname}` 那一行加到你的 `Genotype\\MZCER048.CUL`；
4. 用 XBuild 打开 `<plot>.MZX` 即可运行。

## 品种参数（{variety}）

| P1 | P2 | P5 | G2 | G3 | PHINT |
|---|---|---|---|---|---|
| {p1} | {p2} | {p5} | {g2} | {g3} | {phint} |

## 出处

由 `analysis/dssat/src/s08_bundle.py --trial {trial}` 生成，
原始结果见 `analysis/dssat/report/`。
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial", default="SD_SF", choices=list(dc.TRIALS))
    ap.add_argument("--scenario", default="base", choices=list(dc.SCENARIOS))
    ap.add_argument("--with-exe", action="store_true",
                    help="also copy dscsm048.exe into _support/bin (self-contained)")
    args = ap.parse_args()
    dc.set_trial(args.trial)
    dc.set_scenario(args.scenario)

    bundle = dc.OUT / "bundle" if dc.SCENARIO == "base" else dc.OUT / "bundle_nominal"
    if bundle.exists():
        shutil.rmtree(bundle)
    (bundle / "Weather").mkdir(parents=True)

    staging = dc.OUT / "_bundle_tmp"
    if staging.exists():
        shutil.rmtree(staging)

    print(f"trial {dc.TRIAL}  scenario {dc.SCENARIO}  -> {bundle}")
    for plot in dc.ORDER:
        dc.run_plot(plot, keep_dir=staging / plot)

    # --- genotype + config + weather, taken from any one run directory ----
    ref = staging / dc.ORDER[0]
    for f in ("MZCER048.CUL", "MZCER048.SPE", "MZCER048.ECO"):
        shutil.copy2(ref / f, bundle / f)
    wth = next((ref / "Weather").glob("*.WTH"))
    shutil.copy2(wth, bundle / "Weather" / wth.name)

    # --- one .MZX and one SOIL.SOL per plot -------------------------------
    soils = []
    for plot in dc.ORDER:
        d = staging / plot
        mzx = next(d.glob("*.MZX"))
        # DSSAT keeps only the first 8 characters of the FileX base name, so
        # "SD_SF_A5B2.MZX" would be truncated (Error key IPEXP).  Use the plot
        # code (<= 4 chars) instead.
        shutil.copy2(mzx, bundle / f"{plot}.MZX")
        soils.append(d / "SOIL.SOL")
    n = _merge_soil(soils, bundle / "SOIL.SOL")
    print(f"  merged {n} soil profiles into SOIL.SOL")

    # --- support files so only the exe is missing -------------------------
    static = Path(dc.__file__).resolve()
    import DSSATTools
    static = Path(DSSATTools.__file__).parent / "static"
    for name in SUPPORT:
        shutil.copytree(static / name, bundle / "_support" / name)
    print(f"  support: {', '.join(SUPPORT)}")
    if args.with_exe:
        import DSSATTools
        exe = Path(DSSATTools.__file__).parent / "static" / "bin" / "dscsm048.exe"
        (bundle / "_support" / "bin").mkdir(exist_ok=True)
        shutil.copy2(exe, bundle / "_support" / "bin" / exe.name)
        print(f"  exe    : {exe.name}")

    (bundle / "run.py").write_text(RUN_PY, encoding="utf-8")
    (bundle / "run.bat").write_text(
        "@echo off\r\npython \"%~dp0run.py\" %*\r\npause\r\n", encoding="utf-8")

    cv = dc.DEFAULT_CULTIVAR
    (bundle / "README.md").write_text(README.format(
        site=dc.CFG.site, variety=dc.CFG.variety, trial=dc.TRIAL,
        nplot=len(dc.ORDER), scenario=dc.SCEN_LABEL, wth=wth.name,
        vrname=dc.VRNAME, p1=cv["P1"], p2=cv["P2"], p5=cv["P5"],
        g2=cv["G2"], g3=cv["G3"], phint=cv["PHINT"]), encoding="utf-8")
    shutil.rmtree(staging)

    # --- static consistency checks (no exe needed) ------------------------
    soil_ids = {l.split()[0][1:] for l in
                (bundle / "SOIL.SOL").read_text(encoding="cp437").splitlines()
                if l.startswith("*") and not l.startswith("*SOILS")}
    cul = (bundle / "MZCER048.CUL").read_text(encoding="cp437")
    wsta = wth.stem[:6]
    rows = []
    for plot in dc.ORDER:
        txt = (bundle / f"{plot}.MZX").read_text(encoding="cp437")
        fld = [l.split() for l in txt.splitlines()
               if l.startswith(" 1 ") and "DR000" in l][0]
        sid = next(tok for tok in fld if len(tok) == 10 and tok.isalnum()
                   and tok.startswith(dc.SOIL_PREFIX))
        rows.append(dict(plot=plot, mzx=f"{plot}.MZX", soil=sid,
                         soil_ok=sid in soil_ids,
                         wsta_ok=wsta in txt,
                         cul_ok=dc.VRNAME in cul))
    chk = pd.DataFrame(rows)
    print(chk.to_string(index=False))
    if not (chk.soil_ok.all() and chk.wsta_ok.all() and chk.cul_ok.all()):
        raise SystemExit("bundle consistency check failed")
    print(f"  {len(soil_ids)} soil profiles, weather {wth.name}, "
          f"cultivar {dc.VRNAME} present")

    # --- index for the whole trial folder ---------------------------------
    (dc.OUT / "README.md").write_text(OUT_README.format(
        trial=dc.TRIAL, site=dc.CFG.site, variety=dc.CFG.variety,
        nplot=len(dc.ORDER), bundle=bundle.name), encoding="utf-8")

    print("")
    print("bundle ->", bundle)
    print("  run it with:  cd <bundle> && python run.py --exe <dscsm048.exe>")


if __name__ == "__main__":
    main()
