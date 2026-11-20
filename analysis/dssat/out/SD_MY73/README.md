# SD_MY73 — 山东潍坊 MY73（4 个小区）

## 想跑模拟 → 进 `bundle/`

```bash
cp -r bundle C:/dssat        # 一定要复制到短路径，见下
cd C:/dssat
python run.py --exe <dscsm048.exe>
```

`bundle/` 是自包含的（MZX + SOIL.SOL + WTH + 品种参数 + 支持文件），
里面有单独的 README。

⚠ DSSAT 的路径变量是定长的，超过约 60 字符会被静默截断，报
`File not found ... SIMCDE`。这个目录的路径本身就很长，**不要就地运行**。

## 其余目录都是"结果/证据"，不能直接跑

| 目录 | 内容 |
|---|---|
| `bundle/` | **可运行**的完整 DSSAT 文件包 |
| `sim/` | 实测留苗密度的 DSSAT 输出（`Summary.OUT`、`PlantGro.OUT`）+ 验证 CSV |
| `sim_nominal/` | 额定密度情景的输出 |
| `mgmt/` | 生成的 `.MZX`（文件名 >8 字符、且没有品种文件，单独跑会失败） |
| `mgmt_nominal/` | 额定密度情景的 `.MZX` |
| `soil/` | 单小区 `.SOL` |
| `weather/` | `.WTH` |
| `diag/` | 四个诊断问题的 CSV |

数据来源与标定过程见 `analysis/dssat/report/`。
