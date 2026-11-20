# DSSAT 校验（CERES-Maize，Py_DSSATTools）

三个试验的 DSSAT 复现，**输入、输出、图、报告全部按试验分开**：

试验目录名沿用数据集里的 `trial_id`（`SF` = 水肥试验，`SD_MY73` = 潍坊品种展示）：

| 试验 | 地点 / 品种 | 小区 | 输入输出 | 图 | 报告 |
|---|---|---|---|---|---|
| `SD_SF` | 山东德州 / 乐农87 | 9 | `out/SD_SF/` | `fig/SD_SF/` | [参数](report/DSSAT_PARAMETERS_SD_SF.md) · [验证](report/DSSAT_VALIDATION_SD_SF.md) |
| `XX_SF` | 河南新乡 / MY73 | 9 | `out/XX_SF/` | `fig/XX_SF/` | [参数](report/DSSAT_PARAMETERS_XX_SF.md) · [验证](report/DSSAT_VALIDATION_XX_SF.md) |
| `SD_MY73` | 山东潍坊 / MY73 | 4 | `out/SD_MY73/` | `fig/SD_MY73/` | [报告](report/DSSAT_SD_MY73.md) |

补充诊断：

- [四问诊断](report/DSSAT_DIAGNOSTICS.md)：密度、单株叶面积、生物量与 HI、穗粒数
- [MY73 两地三问](report/DSSAT_MY73_QUESTIONS.md)：潍坊密度敏感性、IXIM 叶面积参数、
  两地能否共用一套品种参数

产物在 `out/<trial>/diag/` 与 `fig/<trial>/diag/`，不影响上面的结果。

- 模型：DSSAT 4.8（`dscsm048`），CERES-Maize
- 驱动库：[Py_DSSATTools](https://github.com/daquinterop/Py_DSSATTools) 2.1.4
- 观测与土壤/气象数据复用 `analysis/out/dataset/` 的标准表

## 目录

```
analysis/dssat/
  src/
    dssat_common.py    试验注册表、导入垫片、输入构建、运行器
    s01_inputs.py      写 .WTH/.SOL，导出输入审计表与管理日程
    s02_run.py         9 小区批量运行 → out/<trial>/sim/dssat_daily.csv
    s03_calibrate.py   品种参数分阶段标定 → calibration_log.csv
    s04_validate.py    精度、物候、LAI 平台期 → validate_*.csv
    s05_figures.py     四张对比图 → fig/<trial>/D1..D4
    s06_nstress.py     氮胁迫 → fig/<trial>/D5_氮胁迫.png
    s07_diag.py        诊断：密度 / 单株叶面积 / 生物量·HI / 穗粒数
    s08_bundle.py      生成可直接运行的 DSSAT 文件包（out/<trial>/bundle）
    s09_my73_questions.py  MY73 三问：密度敏感性 / IXIM / 共用品种参数
    s10_leafarea_shared.py 单株叶面积定标成图 + 共用品种参数（产量/生物量优先）
    s11_ixim_eval.py   IXIM vs CERES 的产量/生物量精度对比
    s12_leafnumber.py  叶片数（L#SD）模拟 vs 实测 N 叶期
    s13_phint.py       PHINT × P1 扫描：调出叶时间的代价
    run_all.py         依次跑两个试验（s03 需 --calibrate）
  out/
    SD_SF/{weather,soil,mgmt,sim}
    XX_SF/{weather,soil,mgmt,sim}
    SD_MY73/{weather,soil,mgmt,sim}
  fig/
    SD_SF/D1..D5（实测留苗）  N1..N5（额定密度）  diag/E_诊断.png
    XX_SF/D1..D5（实测留苗）  N1..N5（额定密度）  diag/E_诊断.png
    SD_MY73/D1..D3（无 LAI 观测、单一氮水平，故无 D4/D5）
  report/
    DSSAT_PARAMETERS_SD_SF.md   DSSAT_VALIDATION_SD_SF.md
    DSSAT_PARAMETERS_XX_SF.md   DSSAT_VALIDATION_XX_SF.md
    DSSAT_SD_MY73.md            DSSAT_DIAGNOSTICS.md
    DSSAT_MY73_QUESTIONS.md
```

## 运行

```bash
# 单个试验
python analysis/dssat/src/s01_inputs.py --trial SD_SF
python analysis/dssat/src/s02_run.py    --trial SD_SF
python analysis/dssat/src/s03_calibrate.py --trial SD_SF   # 可选，约 4-5 分钟
python analysis/dssat/src/s04_validate.py  --trial SD_SF
python analysis/dssat/src/s05_figures.py   --trial SD_SF
python analysis/dssat/src/s06_nstress.py   --trial SD_SF

# 四个诊断（各 4-9 分钟，含 G2/G3 重标定网格）
python analysis/dssat/src/s07_diag.py --trial SD_SF

# 两个试验一起跑（s03 较慢，默认跳过）
python analysis/dssat/src/run_all.py
```

`--trial` 省略时默认 `SD_SF`。加 `--scenario nominal` 即用额定播种密度重跑，
结果写到 `out/<trial>/sim_nominal/`、图写到 `fig/<trial>/N1–N5`（D 系列不受影响）：

```bash
for s in s02_run s04_validate s05_figures s06_nstress; do
  python analysis/dssat/src/$s.py --trial SD_SF --scenario nominal
done
```
`s02_run.py` 支持临时覆盖品种参数与初始硝态氮：

```bash
python analysis/dssat/src/s02_run.py --trial XX_SF --cultivar P1=235,PHINT=50 --nitrate-scale 0.55
```

## 模拟文件包（给别人跑）

每个试验在 `out/<trial>/bundle/` 下有一套**自包含**的 DSSAT 文件，
包含品种参数，拷走即可运行：

```
out/SD_SF/bundle/
  A1.MZX … A5B2.MZX        9 份试验文件（FileX，文件名 ≤8 字符是 DSSAT 的硬约束）
  SOIL.SOL                 9 个土壤剖面合并
  Weather/DZ012601.WTH     气象（130 天）
  MZCER048.CUL             品种参数（含标定后的 LENONG87）
  MZCER048.SPE / .ECO      物种 / 生态型
  _support/                StandardData / Soil / Pest
  run.py / run.bat         一键运行
  README.md                说明
```

```bash
cd <bundle>                       # 建议放到短路径，如 C:\dssat
python run.py --exe <dscsm048.exe>
```

> DSSAT 的路径变量是定长的，**超过约 60 字符会被静默截断**（报
> `File not found ... SIMCDE`）。所以不要把 bundle 放在很深的目录里跑，
> `run.py` 会检查路径长度并提示。
>
> 生成/校验：`python analysis/dssat/src/s08_bundle.py --trial SD_SF`
> （检查每个 MZX 引用的 `ID_SOIL`、`WSTA`、品种名在包里都存在）。

## 图件系列

| 系列 | 位置 | 内容 |
|---|---|---|
| C / L | `analysis/out/fig/` | APSIM（既有，非本目录） |
| D1–D5 | `fig/<trial>/` | DSSAT，实测留苗密度（`out/<trial>/sim/`） |
| N1–N5 | `fig/<trial>/` | DSSAT，额定播种密度（`out/<trial>/sim_nominal/`） |
| E | `fig/<trial>/diag/` | 四个诊断问题汇总 |

每张图的标题都写成 `前缀｜试验｜密度　内容`，例如
`D1｜德州水肥试验 乐农87｜实测留苗密度　逐小区 LAI：实测 vs 模拟`，
单看图片也能分清系列和试验。`SD_MY73` 没有实测留苗，所以只有 D 系列（D1–D3）。

## 结果速览

| 指标 | SD_SF（乐农87） | XX_SF（MY73） | SD_MY73（MY73） |
|---|---|---|---|
| 产量 RRMSE | 5.16% | 6.29% | 15.7% |
| 产量偏差 | +1.52% | +5.07% | −15.7% |
| LAI RMSE | 0.697 | 0.539 | 无 LAI 观测 |
| 生物量 RRMSE | 9.28% | 20.9% | 39.4% |
| 平台期 模拟/实测 | 13.6 / 13.4 d | 19.0 / 10.1 d | 14.0 / — d |
| 吐丝 偏差 | +1 d | +2 d | −1 d |
| 成熟 偏差 | −5 d | +1 d（实测缺失，推定） | 0 d |
| 最终叶数 模拟/实测 | **20 / 20** | 20 / 21 | 20 / — |

> `SD_SF` 已于 2026-11 改用"方案 A"（`P1 250 / PHINT 51 / P5 880 / G3 6.8`）：
> 叶数从 17 补到 20、吐丝 08-10、生物量 −10.9% → −7.7%，
> 代价是产量 RRMSE 3.70% → 5.16%、LAI RMSE 0.44 → 0.70、成熟早 5 天。
>
> 两个 MY73 试验的物种参数 `PARSR` 已由 0.50 调到 **0.55**（仅这两个试验，
> 通过 `TrialCfg.species` 覆盖，SD_SF 不受影响）：新乡生物量 −22.1% → −19.2%、
> 产量 4.20% → 6.29%；潍坊产量 18.6% → 15.7%、生物量 −41.3% → −39.1%。
>
> `SD_MY73` 没有施肥/灌溉记录，管理是假设的，绝对误差不能当作模型精度，
> 详见其报告 §3。

> 顺带修正：D1/D2 原来用固定纵轴（LAI 0–4.4、生物量 0–20000），
> 会截掉 XX_SF 的实测峰值（LAI 5.63、生物量 24.8 t/ha）。现在纵轴按数据自动取值，
> D 系列因此重新渲染过一次，**数值与 CSV 完全未变**。

## 环境说明

DSSATTools 2.1.4 在 Windows 上有两个坑，`dssat_common.py` 已处理：

1. **导入**：库把 `static/` 镜像到 `%TEMP%/DSSAT048` 时用 `os.symlink`，无开发者模式会抛 `WinError 1314`。导入前清空该目录并把 `os.symlink` 降级为真实拷贝。
2. **输入格式**：土壤剖面 ID 必须 10 字符且不含下划线、`*SOILS` 描述限 18 字符、`@SITE` 需 `F8.3` 才装得下 113.93 的经度；`*TREATMENTS` 因子行必须是整数，若某处理无施肥（0 N 对照）库会写出 `0.0`，DSSAT 只输出 MgmtOps/OVERVIEW 而不报错。

另注意：模块级默认参数在 import 时绑定，切换试验后必须显式传参或在使用前解析
（`s02_run.py` 的 `--nitrate-scale` 已改为 `None` 再解析）。
