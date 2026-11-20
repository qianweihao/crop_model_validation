# DSSAT 校验（CERES-Maize，Py_DSSATTools）

两个试验的 DSSAT 复现，**输入、输出、图、报告全部按试验分开**：

| 试验 | 地点 / 品种 | 小区 | 输入输出 | 图 | 报告 |
|---|---|---|---|---|---|
| `SD_SF` | 山东德州 / 乐农87 | 9 | `out/SD_SF/` | `fig/SD_SF/` | [参数](report/DSSAT_PARAMETERS_SD_SF.md) · [验证](report/DSSAT_VALIDATION_SD_SF.md) |
| `XX_SF` | 河南新乡 / MY73 | 9 | `out/XX_SF/` | `fig/XX_SF/` | [参数](report/DSSAT_PARAMETERS_XX_SF.md) · [验证](report/DSSAT_VALIDATION_XX_SF.md) |

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
  out/
    SD_SF/{weather,soil,mgmt,sim}
    XX_SF/{weather,soil,mgmt,sim}
  fig/
    SD_SF/D1..D5   XX_SF/D1..D5
  report/
    DSSAT_PARAMETERS_SD_SF.md   DSSAT_VALIDATION_SD_SF.md
    DSSAT_PARAMETERS_XX_SF.md   DSSAT_VALIDATION_XX_SF.md
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

# 两个试验一起跑（s03 较慢，默认跳过）
python analysis/dssat/src/run_all.py
```

`--trial` 省略时默认 `SD_SF`。`s02_run.py` 支持临时覆盖品种参数与初始硝态氮：

```bash
python analysis/dssat/src/s02_run.py --trial XX_SF --cultivar P1=235,PHINT=50 --nitrate-scale 0.55
```

## 结果速览

| 指标 | SD_SF（乐农87） | XX_SF（MY73） |
|---|---|---|
| 产量 RRMSE | 3.70% | 4.20% |
| 产量偏差 | +0.11% | +0.09% |
| LAI RMSE | 0.436 | 0.504 |
| 生物量 RRMSE | 12.3% | 23.5% |
| 平台期 模拟/实测 | 14.1 / 13.4 d | 19.0 / 10.1 d |
| 吐丝 偏差 | +3 d | +2 d |
| 成熟 偏差 | 0 d | +1 d（实测缺失，推定） |

## 环境说明

DSSATTools 2.1.4 在 Windows 上有两个坑，`dssat_common.py` 已处理：

1. **导入**：库把 `static/` 镜像到 `%TEMP%/DSSAT048` 时用 `os.symlink`，无开发者模式会抛 `WinError 1314`。导入前清空该目录并把 `os.symlink` 降级为真实拷贝。
2. **输入格式**：土壤剖面 ID 必须 10 字符且不含下划线、`*SOILS` 描述限 18 字符、`@SITE` 需 `F8.3` 才装得下 113.93 的经度；`*TREATMENTS` 因子行必须是整数，若某处理无施肥（0 N 对照）库会写出 `0.0`，DSSAT 只输出 MgmtOps/OVERVIEW 而不报错。

另注意：模块级默认参数在 import 时绑定，切换试验后必须显式传参或在使用前解析
（`s02_run.py` 的 `--nitrate-scale` 已改为 `None` 再解析）。
