# 玉米水氮胁迫模拟（zoomlion_crop_model）

用 APSIM Next Generation 复现 2026 年山东德州水肥试验（品种乐农87，9 个小区，
4 个氮水平 × 2 种施肥方式 + 零氮对照）的 LAI、生物量和产量。

**技术参考（完整输入与输出）→ [`reports/MODEL_INPUTS_OUTPUTS.md`](reports/MODEL_INPUTS_OUTPUTS.md)**

---

## 环境

| 项 | 状态 |
|---|---|
| APSIM | **2026.2.7980.0**（`C:/Program Files/APSIM2026.2.7980.0`） |
| apsimNGpy | 可用（升级后 `APSIM.Core.dll` 就位） |
| 版本选择 | `apsim_env.py` 自动选最高版本，可用 `APSIM_BIN` 覆盖 |
| 方案 | JSON 直编 `.apsimx` + `Models.exe` 批跑 + sqlite3 读结果 |
| 并发 | `APSIM_WORKERS=3`，`Models.exe` 低于正常优先级，不卡桌面 |

## 结果（山东乐农87，9 小区）

| 指标 | 值 |
|---|---|
| 产量 RRMSE | **8.47%** |
| 产量偏差 | **+2.81%** |
| LAI RMSE | **0.342** |
| 生物量 RRMSE | 21.92%（偏差 −19%） |
| HI 模拟/实测 | 0.619 / 0.486 |

河南两个试验本轮未纳入（`XX_IRR` 全部为 0），原因见技术参考 §3.6。

## 运行

```bash
python analysis/src/run_all.py                    # 全流程
python analysis/src/s10_phenology.py --variety 乐农87      # 物候标定
python analysis/src/s07_calibrate.py --variety 乐农87 --stage growth
python analysis/src/s08_validate.py               # 验证（只出表，不出图）
python analysis/src/s20_conservative_figures.py   # C 系列图（sen=0.20，优化后）
python analysis/src/s12_figures.py                # L 系列图（sen=0.33 重跑，优化前）
```

一键：`python analysis/src/run_all.py [--calibrate]`

## 图件

`analysis/out/fig/`，共 8 张。C 与 L 图式完全一致，仅参数版本不同，可并排对照。

| C 系列（**优化后** sen=0.20） | L 系列（**优化前** sen=0.33） | 内容 |
|---|---|---|
| `C1_LAI.png` | `L1_LAI.png` | 逐小区 LAI |
| `C2_生物量.png` | `L2_生物量.png` | 逐小区地上部生物量 |
| `C3_产量.png` | `L3_产量.png` | 逐小区产量 + 1:1 |
| `C4_氮响应.png` | `L4_氮响应.png` | 产量/生物量/HI 对施氮量 |

C 与 L 只差一个参数 `SenescenceAcceleration`（0.20 vs 0.33）：

| 指标 | L 优化前 | C 优化后 |
|---|---|---|
| 产量 RRMSE | 8.98% | **8.47%** |
| LAI RMSE | 0.375 | **0.342** |
| LAI 平台期 | −8.06% | **−5.24%**（实测 −3.2%） |
| 生物量 RRMSE | 22.53% | 21.92% |
| A1 季末衰减 | −45% | −29% |

## 标定流程

两阶段顺序标定，顺序不可颠倒（先物候后生长，否则生长参数会吸收物候误差）。

1. **物候**：用 28 个叶片出现观测 + 抽雄/开花/成熟三个日期联合标定。
   仅用后两个日期会不可辨识（曾标出 30 和 207，差 7 倍）。
   结果：抽雄 −1 d、开花 0 d、成熟 0 d、叶片数 RMSE 1.36。
2. **生长**：RUE、籽粒库容、冠层形态、衰老速率。
3. **土壤氮**：以零氮对照 A1 产量标定初始无机氮 → 110 kg N/ha。

## 已知限制（详见技术参考 §3）

1. **氮响应过弱**：只有 A1 受氮限制，139 kg N/ha 及以上完全饱和（加氮到 400 kg N/ha 产量变化 +0）。
2. **生物量偏低 19%**：缺口随施氮量放大（A1 −5% → A5B2 −32%），是氮响应斜率问题而非尺度问题。
3. **HI 偏高 0.13**：产量准但生物量低，结构性偏差。
4. **平台期与季末陡降不可兼得**：APSIM 玉米衰老只有一个恒定加速度参数。
5. **FN（光合氮胁迫因子）恒为 1.0**，是硬编码常量，不能用于判断氮胁迫；
   氮限制实际作用于叶面积扩张与灌浆充实。

## 目录

```
analysis/src/         s01-s23 管线脚本
analysis/reports/     MODEL_INPUTS_OUTPUTS.md（技术参考）
analysis/out/dataset/ 标准化观测表
analysis/out/apsim/   23 个 .apsimx
analysis/out/sim/     模拟结果与验证表
analysis/out/fig/     图件（8 张）
analysis/out/_archive/ 历史结果快照
```
