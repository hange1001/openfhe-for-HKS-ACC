# P1 消除无效塔清零与装载报告 — 2026-09-04

P1 目标：不增加并行度，消除 alpha 以外的无用整塔处理。
预测节省：两 digit 主体少约 12288 次写零 + 1536 次装载迭代 ≈ 13824 周期（约 9.9%）。

## 一、结论摘要

- RTL perf 周期：暖态 139734 → **125910**，实际节省 **13824 周期（9.89%）**，
  与模型预测偏差 **0%**（4608 + 9216，逐 digit 精确命中）。
- 资源：BRAM 688 / DSP 1160 / URAM 96 不变；FF 78323→76984（-1339）、
  LUT 147680→147916（+236）；Fmax 190.48 MHz（5.250ns）不变，蝶形 II=1 不变，
  结构审计 PASS（整机仍 20 MultMod、无 AUTO、AXI 256）。
- INIT 295063 不变。
- 功能：native Top 18/18、HKS 22 用例、BConv 独立 9/9（新增毒化用例 TC-7/8）、
  RTL smoke 35/35、perf 40960 residues 精确一致、OpenFHE Release/ASan
  472 checks / 1523712 residues、ASan 零泄漏。

## 二、实现改动

| 文件 | 改动 |
|---|---|
| `include/bconv_systolic.h` | `Compute_BConv_Systolic` 增加 `int active_q` 参数 |
| `src/bconv_systolic.cpp` | Feed_X 对 `q>=active_q` 显式注入 0（不读未装载的 local_in_x 旧缓存）；Load_X 只装载 `l<active_q` 行（tripcount 1..LIMB_Q）；`active_q` 走 s_axilite |
| `src/top.cpp` | `Prepare_HKS_BConv_Input` 循环界改为 `q<alpha`，去掉无效行写零分支（tripcount 1..LIMB_Q）；HKS 调用传 `alpha`，独立 OP_BCONV 传 `LIMB_Q` |
| `testbench/bconv_systolic_tb.cpp` | golden 模型支持 active_q；新增 TC-7（active_q=2 行2毒化）与 TC-8（active_q=1 行1-2毒化） |

无效行语义：Prepare 不再整塔清零（省写），Load_X 不再装载无效行（省读），
BConv Feed_X 在消费者入口显式注入 0（保证数学等价，不依赖"旧值乘零"）。

## 三、周期结果（真实 OpenFHE fixture，RTL co-sim）

| 调用 | P0 基线 | P1 | 节省 |
|---|---:|---:|---:|
| INIT（冷启动） | 295063 | 295063 | 0 |
| digit alpha=2 start=0 | 70131 | 65523 | 4608 |
| digit alpha=1 start=2 | 69603 | 60387 | 9216 |
| 两个 digit 暖态合计 | 139734 | 125910 | **13824（9.89%）** |

模型核对：每 digit 节省 = (3-alpha) 行 × (4096 写零 + 512 装载) = (3-alpha)×4608。
alpha=2：4608 ✓；alpha=1：9216 ✓；合计 13824 ✓，偏差 0%。

## 四、资源与时序

| 指标 | P0 基线 | P1 | 判定 |
|---|---:|---:|---|
| BRAM_18K / DSP / URAM | 688 / 1160 / 96 | 688 / 1160 / 96 | 不变 ✓ |
| FF | 78323 | 76984 | -1339（去掉无效行守卫逻辑） |
| LUT | 147680 | 147916 | +236（Feed_X 注入 0 的选择逻辑） |
| 估算时钟 | 5.250 ns | 5.250 ns | II 不退化 ✓ |
| 结构审计 | PASS | PASS | 20 MultMod、无 AUTO、256-bit AXI、蝶形 II=1 ✓ |

Prepare_HKS_BConv_Input 的 HLS 报告延迟仍显示 12312 为 tripcount max=3 的保守估计；
实际运行按 alpha 缩短（见第三节实测）。

## 五、验证矩阵

| 层级 | 结果 |
|---|---|
| native Top / HKS | 18/18 通过；22 合法用例 0 失败（含 alpha 交替、非法描述符与哨兵） |
| BConv 独立模块 | 9/9 通过（新增毒化行 TC-7/8 验证 Feed_X 注入 0 而非读旧值） |
| OpenFHE Release | 472 checks / 1523712 residues 精确一致 |
| OpenFHE ASan | 472 checks / 1523712 residues，泄漏检查 0 错误 |
| RTL smoke | 35/35 调用、4 种 digit、0 失败 |
| RTL perf | 40960 residues 与 CPU oracle 精确一致 |
| 结构/调度 | bank 数/类型、实例数、II、AXI 位宽与 P0 契约一致 |

## 六、证据文件

- `top-csynth.rpt` / `rtl-audit.json`：mem_p1_r1 综合与结构审计
- `perf-cosim.rpt` / `perf-vitis-output.txt`、`smoke-cosim.rpt` / `smoke-vitis-output.txt`
- `native-test.txt`、`bconv-module-csim.txt`、`openfhe-*`、`asan-*`
- `source-diff-stat.txt`：P1 源码改动统计（相对 P0 基线提交 4b7eabf）
