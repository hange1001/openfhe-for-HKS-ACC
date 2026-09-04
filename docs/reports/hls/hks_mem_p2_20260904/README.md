# P2 BConv 直接读写工作区报告 — 2026-09-04

P2 目标：去掉 `poly_buffer_1→local_in_x→local_out_x→poly_buffer_1` 的两次整块复制，
systolic 核心直接从 Q 区 Feed_X、直接向补基区 Collect。

## 一、结论摘要

- RTL perf 周期：暖态 125910 → **120822**，P2 节省 **5088 周期**（模型约 5128，
  偏差 0.8%，小于 10% 验收线）；INIT 295063 不变。
- HLS BRAM_18K：688 → **440（-248）**。其中 local_in_x(48)+local_out_x(80)=128 为
  直接消除对象；另外 120 来自去掉 8 路 Load_X/Store_X 访问视图后，Vitis 将
  poly_buffer_1 按"每 limb 1R+1W"重新映射（每 limb 15→? BRAM，少复制 120）。
- Systolic_Loop II=1；BConv 实例 BRAM=0、延迟 8230→4145 周期。
- DSP 1160 不变、Fmax 190.48 MHz（5.250ns）不变、结构审计 PASS（20 MultMod、
  无 AUTO、256-bit AXI、蝶形 II=1）。
- LUT 147916→154852（+6936）：直接寻址的 bank 选择/写使能 MUX 代价，如实报告；
  FF 76984→76770（-214）。
- 功能：native Top 18/18、HKS 22、BConv 9/9（含毒化）、smoke 35/35、
  perf 40960 residues 精确一致、OpenFHE Release/ASan 472/1523712、ASan 零泄漏。

## 二、实现改动

| 文件 | 改动 |
|---|---|
| `src/bconv_systolic.cpp` | `bconv_systolic_core` 接口改为直接收 `in_x[MAX_LIMBS][SQRT][SQRT]` 工作区：Feed_X 读 `in_x[q][(t-q)>>LOG_SQRT][(t-q)&(SQRT-1)]`（负地址用 safe_idx 规避、值由 valid 选择），Collect 写 `in_x[LIMB_Q+p][(t-3-p)>>LOG_SQRT][(t-3-p)&(SQRT-1)]`；删除 local_in_x/local_out_x 与 Load_X/Store_X；保留 local_w/mod/S/m_barrett 与阵列流水寄存器 |
| `check_systolic_banks.py`（新增） | bank 访问检查器：遍历全部 4103 拍 × active_q(1..3) × sizeP(1..5)，验证每 limb 每周期至多 1 读+1 写、读写不同 limb 实例、地址范围合法——15 组合全部无冲突 |

## 三、周期结果（真实 OpenFHE fixture，RTL co-sim）

| 调用 | P1 | P2 | P2 节省 |
|---|---:|---:|---:|
| INIT（冷启动） | 295063 | 295063 | 0 |
| digit alpha=2 start=0 | 65523 | 62979 | 2544 |
| digit alpha=1 start=2 | 60387 | 57843 | 2544 |
| 两个 digit 暖态合计 | 125910 | 120822 | **5088** |

模型：每 digit 节省 = Load_X(alpha×512+2) + Store_X(sizeP×512+2)。
alpha=2：1024+1536+4 ≈ 2564；alpha=1：512+2048+4 ≈ 2564；合计约 5128。
实测 5088，偏差 0.8%。历史基线把"两次 BConv 加载/写回"归属为 6664 周期的
独立上界项；P1 已先删掉其中无效行装载（不计入 P2），P2 实测 5088 不与 P1 的
13824 重复相加。

## 四、资源与时序

| 指标 | P1 | P2 | 判定 |
|---|---:|---:|---|
| BRAM_18K | 688 | **440** | 严格下降 ✓（-248，超出 128 目标） |
| DSP | 1160 | 1160 | 不增加 ✓ |
| URAM | 96 | 96 | 不变 |
| FF | 76984 | 76770 | -214 |
| LUT | 147916 | 154852 | +6936（直接寻址 MUX，如实报告） |
| 估算时钟 | 5.250 ns | 5.250 ns | 频率未退化 ✓ |
| Systolic_Loop II | 1 | 1 | ✓ |
| BConv 实例 | 128 BRAM / 8230 周期 | 0 BRAM / 4145 周期 | ✓ |

LUT 增长是"存储减少换选择逻辑"的代价；本阶段验收标准（II=1、DSP≤基线、
BRAM 严格下降、BConv 周期下降、频率不降）全部满足，不因 LUT 上升视为失败，
但如实记录并留给 P5 物理实现观察布线影响。

## 五、验证矩阵

| 层级 | 结果 |
|---|---|
| 独立数学参考 | BConv golden（128-bit 精确）9/9 一致；HKS 22 用例（scalar CRT/NTT 参考）0 失败 |
| native/C-sim | Top 18/18、BConv 9/9（含 TC-7/8 无效行毒化）、hks-csim 0 错误 |
| bank 检查器 | 4103 拍 × 15 组合（active_q×sizeP）逐周期无冲突 |
| OpenFHE Release/ASan | 472 checks / 1523712 residues 精确一致；ASan 零泄漏 |
| RTL smoke | 35/35 调用、4 种 digit、0 失败 |
| RTL perf | 40960 residues 与 CPU oracle 精确一致 |
| 结构/调度 | 审计 PASS：单套 transform 链、4 共享模乘、整机 20 MultMod、AXI 256、蝶形 II=1 |

## 六、证据文件

- `top-csynth.rpt` / `rtl-audit.json`：mem_p2_r1 综合与结构审计
- `systolic-bank-check.txt`：bank 检查器全组合输出
- `perf-cosim.rpt` / `perf-vitis-output.txt`、`smoke-cosim.rpt` / `smoke-vitis-output.txt`
- `native-test.txt`、`bconv-module-native.txt`、`openfhe-*`、`asan-*`
- `source-diff-stat.txt`：P1+P2 累计源码改动统计（相对基线提交 4b7eabf）
