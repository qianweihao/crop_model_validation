# DSSAT bundle — 山东德州 乐农87 (SD_SF)

DSSAT 4.8 / CERES-Maize (`MZCER048`).  9 个小区，实测留苗密度。

## 怎么跑

```bash
cd <这个文件夹>
python run.py                 # 自动写 DSSATPRO.V48 并跑全部 .MZX
```

或双击 `run.bat`。需要 `dscsm048.exe`（DSSAT 4.8 自带），默认在
`_support/bin/`、`C:\DSSAT48\` 下找，也可以 `python run.py --exe <路径>`。
用 Py_DSSATTools 的话，exe 在 `<python>\site-packages\DSSATTools\staticin\dscsm048.exe`。

⚠ **路径要短**：DSSAT 的路径变量是定长的，超过约 60 字符会被静默截断，
表现为 `File not found ... SIMCDE` 错误。建议复制到 `C:\dssat` 再跑，
`run.py` 会检查并提示。

## 文件

| 文件 | 说明 |
|---|---|
| `<plot>.MZX` | 每个小区一份 DSSAT 试验文件（FileX）；文件名 ≤8 字符是 DSSAT 的硬约束 |
| `SOIL.SOL` | 9 个土壤剖面合并成一个文件 |
| `Weather/DZ012601.WTH` | 气象文件（130 天） |
| `MZCER048.CUL` | 品种参数，含标定后的 `LENONG87` |
| `MZCER048.SPE` / `.ECO` | 物种 / 生态型参数 |
| `_support/` | StandardData / Soil / Pest，DSSAT 运行所需 |
| `results/<plot>/` | 跑完的 `Summary.OUT`、`PlantGro.OUT` 等 |

## 在 DSSAT 图形界面里导入

1. 把 `Weather/DZ012601.WTH` 复制到 DSSAT 安装目录的 `Weather\`；
2. 把 `SOIL.SOL` 里的剖面并入 `Soil\SOIL.SOL`；
3. 把 `MZCER048.CUL` 里 `LENONG87` 那一行加到你的 `Genotype\MZCER048.CUL`；
4. 用 XBuild 打开 `<plot>.MZX` 即可运行。

## 品种参数（乐农87）

| P1 | P2 | P5 | G2 | G3 | PHINT |
|---|---|---|---|---|---|
| 250.0 | 0.3 | 880.0 | 1300.0 | 6.2 | 51.0 |

## 出处

由 `analysis/dssat/src/s08_bundle.py --trial SD_SF` 生成，
原始结果见 `analysis/dssat/report/`。
