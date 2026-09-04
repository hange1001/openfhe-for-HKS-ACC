# P0 基线冻结与同源码重跑报告 — 2026-09-04

P0 目标：建立可复现、同源码的前后比较。对比锚点为历史基线
`hks_no_auto_20260904`（139734 周期），当前源码重跑数字用于解释与历史基线的差异。

## 一、结论摘要

- 当前源码（AUTO 删除 + 地址位选择改写）重跑 native/C-sim/HLS 综合/结构审计/OpenFHE
  集成的结果与历史基线**逐项一致**，位操作改写对功能、资源、时序估计零影响。
- 基线源码已按用户要求提交 Git：`4b7eabf`（防改错检查点）。
- P0 期间发现并修正计划文档一处事实错误：BConv 是**单实例** wrapper
  （128 BRAM_18K），不是"两套各 128"。
- RTL perf co-sim 周期结果见第五节（仿真完成后回填）。

## 二、版本指纹与执行环境

| 项目 | 值 |
|---|---|
| 基线提交（本轮新建） | `4b7eabf3568d34c704db8b2a660470612dba4ad4` |
| 删除 AUTO 前检查点 | `480dc91c81ec9a4b401c739316e8ff8a7d1d5bde` |
| OpenFHE fixture SHA-256 | `f8979c6688fb1ae78ab2d1f99119b32e05d16552cab26c3af33ffe2c0357e348` |
| 实验标签 | `mem_p0_r1`（Vitis HKS_PROJECT_TAG） |
| Vitis HLS | 2023.2 (Build 4023990, Oct 11 2023) |
| 目标器件 | U55C `xcu55c-fsvh2892-2L-e` |
| 平台 | `xilinx_u55c_gen3x16_xdma_3_202210_1` |
| 时钟 | 6ns 目标 + 0.75ns setup uncertainty（物理 7ns 情景另算） |
| g++ | 11.4.0 (Ubuntu 22.04) |

源码状态复原记录：未提交改动此前被 `git stash -u`（stash@{0}
checkpoint-before-hks-storage-rework-20260904）连同计划文档与报告目录一并打包，
导致工作树一度回到含 AUTO 的 `480dc91`。经用户确认后 `git stash pop` 恢复源码改动，
并从 stash 第三提交恢复 `hks_no_auto_20260904/` 报告与计划文档；随后经用户授权提交
`4b7eabf`。stash 未删除，保留为额外保险。

## 三、重跑结果 vs 历史基线

| 验证层 | 历史基线 | mem_p0_r1 重跑 | 一致 |
|---|---|---:|---|
| native Top | 18/18 通过 | 18/18 通过 | ✓ |
| native HKS | 22 用例 0 失败 | 22 用例 0 失败 | ✓ |
| C-sim | 0 错误 | 0 错误 | ✓ |
| HLS BRAM_18K | 688 | 688 | ✓ |
| HLS DSP | 1160 | 1160 | ✓ |
| HLS FF | 78323 | 78323 | ✓ |
| HLS LUT | 147680 | 147680 | ✓ |
| HLS URAM | 96 | 96 | ✓ |
| 估算时钟 | 5.250 ns | 5.250 ns（Fmax 190.48 MHz） | ✓ |
| 结构审计（rtl-audit.json） | PASS | PASS，JSON 逐字段一致 | ✓ |
| RTL smoke | 35 调用 4 种 digit | 35/35，HKS 4 用例 0 失败 | ✓ |
| OpenFHE Release | 472 checks / 1523712 residues | 472 / 1523712 精确一致 | ✓ |
| OpenFHE ASan | 472 / 1523712，零泄漏 | 472 / 1523712，零泄漏 | ✓ |
| XRT 语法编译 | 0 输出 0 退出 | 0 输出 0 退出 | ✓ |
| RTL perf 周期 | 139734（70131+69603） | 139734（70131+69603） | ✓ |

结构审计（`check_shared_transform.py --axi-width 256 --lanes 4
--total-multipliers 20 --no-auto`）关键项：单套 `Execute_Transform_Operation →
Execute_Transform → CG_Transform_Banks` 链各 1 实例；无 AUTO 可达硬件；变换内 4 个
共享模乘（提升到 Top 的物理池）；整机 20 个 MultMod；GMEM0/1/2 均 256 位；正反向蝶形
II=1。与本目录复制的 `top-csynth.rpt` 及历史 `rtl-audit.json` 逐字段一致。

## 四、数组/实例清单（生产者-消费者）

| 数组 | 布局/RAM | 生产者 | 消费者 |
|---|---|---|---|
| `poly_buffer_1[8][64][64]` | cyclic=4、ram_2p BRAM（8×15=120） | Prepare_HKS_BConv_Input（compact 结果）、BConv Store_X | TRANSFORM_LOAD（complement 源）、BConv in_x、ADD/SUB/MULT 操作数1 |
| `poly_buffer_2[8][64][64]` | cyclic=4、ram_2p BRAM（4×29=116） | Load(AXI)、Load_Transform_Tower、TRANSFORM_STORE(store_coeff) | TRANSFORM_LOAD（非 complement）、Prepare 源、ADD/SUB/MULT 操作数2 |
| `result_buffer[8][64][64]` | cyclic=4、ram_2p BRAM（4×29=116） | EVAL 旁路（输入加载时）、TRANSFORM_STORE(!store_coeff)、ADD/SUB/MULT 结果 | Store → mem_out |
| `NTT/INTTTwiddleFactor[8][12][2048]` | complete dim2、cyclic=4、URAM（各 32） | INIT | 变换读 |
| `bank_a/bank_b[RING_DIM]`（Execute_Transform 局部） | cyclic=8、ram_t2p BRAM（32） | TRANSFORM_LOAD | CG_Transform_Banks 12 级 ping-pong → TRANSFORM_STORE |
| `local_in_x[3][4096]`（BConv 局部） | complete dim1、cyclic=8、ram_1wnr BRAM（48） | Load_X | systolic Feed_X |
| `local_out_x[5][4096]`（BConv 局部） | complete dim1、cyclic=8、ram_s2p BRAM（80） | systolic Collect | Store_X |

BConv 修正：`Compute_BConv_Systolic` 在 Top 中为**单实例**（`_32_33` 调用点），
BRAM_18K=128（local_in_x 24 组×2=48 + local_out_x 40 组×2=80），0 DSP；
OP_BCONV 与 HKS_DIGIT 两个调用点共用。`Compute_BConv`/`Compute_BConv_Naive` 不可达。
Load_X 延迟 1538 周期、Store_X 2562 周期，为 P2 的直接消除对象。

## 五、RTL perf 周期（真实 OpenFHE fixture）

perf co-sim 完成：3 次调用（INIT + 两个 digit），40,960 个模数元素与 CPU oracle
精确一致。周期与历史基线逐项相同：

| 调用 | 历史基线 | mem_p0_r1 |
|---|---:|---:|
| INIT（冷启动） | 295063 | 295063 |
| digit alpha=2 start=0 | 70131 | 70131 |
| digit alpha=1 start=2 | 69603 | 69603 |
| 两个 digit 暖态合计 | 139734 | 139734 |
| 总执行（含 INIT） | 434797 | 434797 |

**差异为 0**：地址位选择改写不影响 RTL 周期，无需差异解释，P0 验收达成，可进入 P1。

## 六、证据文件

- `top-csynth.rpt`：mem_p0_r1 综合原始报告（自 csynth 工程复制）
- `rtl-audit.json`：结构审计输出（与历史 rtl-audit.json 逐字段一致）
- `native-test.txt`：native Top 18/18 + HKS 22 用例输出
- `perf-cosim.rpt` / `perf-vitis-output.txt`：真实 fixture RTL co-sim 周期与完整日志
- `smoke-cosim.rpt` / `smoke-vitis-output.txt`：35 调用扩展 smoke
- `openfhe-*` / `asan-*` / `xrt-compile-output.txt`：OpenFHE Release/ASan 集成与 XRT 语法检查
- 对照基线证据在 `../hks_no_auto_20260904/`（提交 4b7eabf 已纳入版本管理）。
