# DSSAT 校验（CERES-Maize，Py_DSSATTools）

三个试验的 DSSAT 复现，**输入、输出、图、结果全部按试验分开**：

| 试验 | 地点 / 品种 | 小区 | 输入输出 | 图 |
|---|---|---|---|---|
| `SD_SF` | 山东德州 / 乐农87 | 9 | `out/SD_SF/` | `fig/SD_SF/` |
| `XX_SF` | 河南新乡 / MY73 | 9 | `out/XX_SF/` | `fig/XX_SF/` |
| `SD_MY73` | 山东潍坊 / MY73 | 4 | `out/SD_MY73/` | `fig/SD_MY73/` |

> 目录名沿用数据集的 `trial_id`（`SF` = 水肥试验，`SD_MY73` = 潍坊品种展示）。

## 文档

| 文件 | 内容 |
|---|---|
| [`report/DSSAT_INPUTS_OUTPUTS.md`](report/DSSAT_INPUTS_OUTPUTS.md) | 模型输入（气象/土壤/初始条件/管理/品种参数）、输出文件、可移植文件包 |
| [`report/DSSAT_RESULTS.md`](report/DSSAT_RESULTS.md) | 精度、逐小区、物候、LAI、氮响应与胁迫、已知限制、图件、补充分析 |

- 模型：DSSAT 4.8（`dscsm048`），CERES-Maize
- 驱动库：[Py_DSSATTools](https://github.com/daquinterop/Py_DSSATTools) 2.1.4
- 观测与土壤/气象数据复用 `analysis/out/dataset/`

## 结果速览

| 指标 | SD_SF（乐农87） | XX_SF（MY73） | SD_MY73（MY73） |
|---|---|---|---|
| 产量 RRMSE | **5.80%** | 4.19% | 21.89% |
| 产量偏差 | **−0.72%** | +0.85% | −21.85% |
| LAI RMSE | **0.713** | 0.501 | 无观测 |
| 生物量 RRMSE | **7.29%** | 16.31% | 35.12% |
| 生物量偏差 | **−3.34%** | −13.61% | −34.74% |
| 吐丝 / 成熟 偏差 | +1 / −5 d | +2 / — | −1 / 0 d |
| 最终叶数 模拟/实测 | **20 / 20** | 20 / 21 | 20 / — |

> **2026-11 定稿**：
> 1. 施肥方式按记录判定（`FACD_MODE="auto"`）：与灌溉同日 → `AP005`（随水施肥），
>    播前基肥 → `AP008`（种肥同播），其余 → `AP001`。原来全按表施尿素算，会多算氨挥发。
>    灌溉水利用系数 **0.90**（记录的是毛水量，写入的是有效水量）。
> 2. 品种参数：SD_SF `P1 250 / PHINT 51 / G3 6.2`（20 片叶、吐丝 08-10）；
>    **两个 MY73 共用一套** `P1 235 / P5 980 / G2 1300 / G3 4.6 / PHINT 50`
>    （`P5` 由潍坊实测 R6 = 10-08 锚定，新乡成熟期因此为 10-12）。
> 3. 物种参数 `PARSR` 0.50 → 0.55（仅两个 MY73）；
>    德州 `CTCNP1` 1.52 → **2.00** + 初始 NO₃ 系数 0.24 → **0.65**，
>    把 N 饱和点从 110 提到 **168 kg N/ha**（田间实测 177），产量 RRMSE 同时从 7.40% 降到 **5.80%**。
> 4. 潍坊关闭氮限制（`NITRO=N`）：无施肥记录、实测 12.2 t/ha 说明氮充足，
>    原假设的 270 kg N/ha 会造出 `NSTD` 峰值 0.44 的假胁迫。
>
> `SD_MY73` 没有施肥/灌溉记录（管理是假设的），绝对误差不能当模型精度，见结果文档 §4。

## 目录

```
analysis/dssat/
  src/
    dssat_common.py        试验注册表、导入垫片、输入构建、运行器
    s01_inputs.py          写 .WTH/.SOL，导出输入审计表与管理日程
    s02_run.py             批量运行 → out/<trial>/sim/dssat_daily.csv
    s03_calibrate.py       品种参数分阶段标定 → calibration_log.csv
    s04_validate.py        精度、物候、LAI 平台期 → validate_*.csv
    s05_figures.py         D1–D4（+ N1–N4 额定密度）
    s06_nstress.py         D5 氮胁迫（+ N5）
    s07_diag.py            诊断：密度 / 单株叶面积 / 生物量·HI / 穗粒数
    s08_bundle.py          生成可直接运行的 DSSAT 文件包
    s09_my73_questions.py  MY73：密度敏感性 / IXIM / 共用品种参数
    s10_leafarea_shared.py 单株叶面积定标 + 共用品种参数（产量/生物量优先）
    s11_ixim_eval.py       IXIM vs CERES 精度对比
    s12_leafnumber.py      叶片数（L#SD）vs 实测可见叶
    s13_phint.py           PHINT × P1 扫描
    s14_sd_sf_phint55.py   SD_SF 出叶间隔优化
    s15_quadrat.py         德州 kg/亩 产量对比（小区实测 / 10 米双行 / 两种密度）
    s17_density_response.py 密度响应曲线（3 个氮水平 x 7 个密度）
    s18_nresponse_tune.py  氮响应重标定（按田间饱和点搜索 CTCNP1 x nitrate_scale）
    run_all.py             依次跑两个水肥试验
  out/<trial>/{weather,soil,mgmt,sim,bundle,diag}
  fig/<trial>/{D1..D5, N1..N5, diag/}
  report/
    DSSAT_INPUTS_OUTPUTS.md   DSSAT_RESULTS.md
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

# 两个水肥试验一起跑（s03 较慢，默认跳过）
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

## 图件系列

| 系列 | 位置 | 内容 |
|---|---|---|
| C / L | `analysis/out/fig/` | APSIM（既有，非本目录） |
| D1–D5 | `fig/<trial>/` | DSSAT，实测留苗密度（`out/<trial>/sim/`） |
| N1–N5 | `fig/<trial>/` | DSSAT，额定播种密度（`out/<trial>/sim_nominal/`） |
| E–L | `fig/<trial>/diag/` | 补充分析图：E 诊断、F 单株叶面积、G 共用品种参数、H 叶片数、I PHINT 优化、J 叶数 20 方案、**K 产量 kg/亩 对比**、**L 密度响应** |

每张图的标题写成 `前缀｜试验｜密度　内容`，例如
`D1｜德州水肥试验 乐农87｜实测留苗密度　逐小区 LAI：实测 vs 模拟`。
`SD_MY73` 没有实测留苗、只有单一氮水平，所以只有 D1–D3；
`K`/`L` 只对 `SD_SF` 生成（需要小区实测 + 多氮水平）。

## 环境说明

DSSATTools 2.1.4 在 Windows 上有几个坑，`dssat_common.py` 已处理：

1. **导入**：库把 `static/` 镜像到 `%TEMP%/DSSAT048` 时用 `os.symlink`，
   无开发者模式会抛 `WinError 1314`。导入前清空该目录并把 `os.symlink` 降级为真实拷贝。
2. **输入格式**：土壤剖面 ID 必须 10 字符且不含下划线、`*SOILS` 描述限 18 字符、
   `@SITE` 需 `F8.3` 才装得下 113.93 的经度；`*TREATMENTS` 因子行必须是整数，
   若某处理无施肥（0 N 对照）库会写出 `0.0`，DSSAT 只输出 MgmtOps/OVERVIEW 而不报错。
3. **IXIM**：DSSATTools 硬编码了 CERES 的 6 参数品种记录，用 IXIM 时需要
   把 `CULTIVAR_HEADER_FMT/ROWS_FMT` 扩到 8 列（多 `AX`/`ALL`），
   并重绑 `_row_writer`，否则参数写不进 CUL。

另注意：模块级默认参数在 import 时绑定，切换试验后必须显式传参或在使用前解析
（`s02_run.py` 的 `--nitrate-scale` 已改为 `None` 再解析）。
