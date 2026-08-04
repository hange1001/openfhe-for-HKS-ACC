# 项目状态总览（OpenFHE for HKS Accelerator）

> 笔记入口。每次新对话/新阶段从这里开始，避免重复读 13 份过期笔记。
> 上次更新：2026-06-29。

---

## 1. 项目背景

基于 [OpenFHE](https://github.com/openfheorg/openfhe-development) 的 RNS-CKKS 实现，将密钥切换（KeySwitch）链路中计算最重的算子卸载到 FPGA（Xilinx Alveo U55C，备用 U250）；同时在 Host 侧设计 **Hybrid KeySwitching（HKS）三种数据调度策略**（DC / MP / OC），权衡延迟、峰值缓冲与可扩展性。

**两个研究轴**：
1. **FPGA 微架构**（[src/fpga_backend/](../../src/fpga_backend/)）— BConv 脉动阵列 + CG-NTT 常几何蝶形 + Barrett 模乘 + 512-bit AXI 突发。
2. **Host 调度**（[src/pke/lib/keyswitch/keyswitch-hybrid.cpp](../../src/pke/lib/keyswitch/keyswitch-hybrid.cpp) + [src/pke/include/keyswitch/hks_strategy.h](../../src/pke/include/keyswitch/hks_strategy.h)）— DC/MP/OC 三策略 + `MemoryTracker` + `HKSStats`。

## 2. 当前进度

### FPGA 侧

| 模块 | 状态 | 关键指标 / 备注 |
|------|------|---------------|
| CG-NTT 替换标准 NTT | ✅ 完成 | 单 limb 延迟 112,671 → 15,701 cycles（7.17× 加速） |
| AXI 地址 / BRAM 索引解耦 | ✅ 完成 | 引入 `axi_l = l - mod_idx_offset`，修复 OOB |
| CG-NTT bit-reversed 输出兼容 | ✅ Host 软件补丁 | `NttForwardOffload` 出口 + `NttInverseOffload` 入口做 BitReverse |
| BConv Systolic vs Naive | ✅ 6.5× 加速实测 | Naive 53,613 → Systolic 8,232 cycles |
| Barrett 替换 128-bit `%` | ✅ 完成 | 消除 30 个 urem，LUT 从 198K 释放 |
| **L1**：512-bit AXI 总线拓宽 | ✅ 完成 | I/O 8×、单 Limb 3.0×、绝对时间 1.733→0.586 ms |
| **L2**：DATAFLOW Limb 流水 | ⏳ 待实现 | 预期 8 Limbs 117K → 87K cycles |
| **L3**：URAM 预加载 TW + BUTTERFLY II=1 | ⏳ 待实现 | 预期 87K → ~46K cycles，累计 7.5× |
| Top FSM 子函数拆解 | ✅ 完成 | `Load_Init_Params` / `Execute_NTT` / `Execute_INTT` 隔离 |
| 时序收敛 @ 200 MHz | ⚠️ HLS Slack -0.33ns（估计偏差） | 待 Vivado P&R 验证 |

### Host 侧

| 模块 | 状态 | 备注 |
|------|------|------|
| `FpgaManager`（XRT 封装） | ✅ 完成 | `NttForwardOffload` / `BConvOffload` 等接口 |
| `BuildCGTwiddle` Host 预计算 | ✅ 完成 | `MathUtils::BitReverse` + perfect-shuffle 追踪 |
| DC 策略（Digit-Centric） | ✅ 完成 | [keyswitch-hybrid.cpp:413](../../src/pke/lib/keyswitch/keyswitch-hybrid.cpp) |
| MP 策略（Max-Parallel） | ✅ 完成 | [keyswitch-hybrid.cpp:396](../../src/pke/lib/keyswitch/keyswitch-hybrid.cpp) |
| OC 策略（Output-Centric, sizeP=1） | ✅ Host 完成 | [keyswitch-hybrid.cpp:468](../../src/pke/lib/keyswitch/keyswitch-hybrid.cpp)，但当前 `BConvOffload(sizeP=1)` 走 CPU 后丢余 tower |
| `HKSStats` 计时插桩 | ✅ 完成 | INTT/BConv/NTT/ModMul 分项计时 |
| `MemoryTracker` 内存轨迹 | ✅ 已实施 | `dump_csv` 输出 step/op/digit/tower/delta/watermark |
| `hks-benchmark` 入口 | ✅ 完成 | [hks-benchmark.cpp](../../src/pke/examples/hks-benchmark.cpp)；`--strategy DC/MP/OC --repeat N` |

### 三策略实测（N=4096, sizeQl=3, sizeP=2, 2 digits）

| 指标 | DC | MP | OC |
|------|-----|-----|-----|
| Precompute 总延迟 (ms) | **4.47** | 5.17 | 5.99 |
| BConv 调用次数 | 2 | 2 | **4** |
| 峰值缓冲 (KB) | 64 | 128 | **32** |
| 传输开销占比 | 32.7% | 34.3% | 32.8% |

> 单 kernel 架构下 DC ≈ MP；OC 用调用次数换最小峰值缓冲；DMA 固定延迟（~20–30 μs/call）是当前瓶颈。

## 3. 笔记导引（保留的 5 份）

| 文件 | 用途 |
|------|------|
| [CG-NTT-Migration-Report.md](CG-NTT-Migration-Report.md) | **当前 FPGA 主架构权威文档**。3 阶段迁移（替换内核 / AXI 解耦 / bit-reverse 修复）+ 数据流图 |
| [cg_ntt_optimization.md](cg_ntt_optimization.md) | **三级递进优化路线**。L1 已落地，L2/L3 是下一步工作的入口文档 |
| [OC_strategy_gap_analysis.md](OC_strategy_gap_analysis.md) | **OC 策略 vs CiFlow 论文 5 项 gap 对比**。解释为什么实测 OC > DC，给出 #2/#3/#4/#5 修复路径与可行性论证 |
| [MemoryTracker_Plan.md](MemoryTracker_Plan.md) | 三策略峰值 SRAM 对比的探针设计（已实施记录） |
| [ntt_path.md](ntt_path.md) | `is_ntt` 模板特化技巧（编译期分支消除），L3 优化的子项 |
| [role.md](role.md) | 协作角色定义（数字 IC 架构师 / 苏格拉底导师），供其他 AI 工具加载 |

[archive/](archive/) 下是 8 份过期笔记，参考 [archive/README.md](archive/README.md) 了解归档原因。

## 4. 论文产出（[../papers/](../papers/)）

| 文件 | 用途 |
|------|------|
| [content.md](../papers/content.md) | 毕业论文 6 章完整大纲（§1 绪论 → §6 总结） |
| [实验规划.md](../papers/实验规划.md) | 实验 1a/1b/1c + 实验 2 + 实验 3 + 三章实测数据 + 图表规划 |
| [fsm_vs_cpu_scheduling.md](../papers/fsm_vs_cpu_scheduling.md) | "CPU 调度 vs FPGA FSM" 决策框架（开销比判据 + 当前 opcode-RPC 模型为何最优） |
| `fig5-2_memory_trace.{py,pdf,png}` | 动态内存轨迹图（DC/MP/OC 三策略） |
| `fig5-3_*.{py,pdf,png}` | 端到端延迟分解柱状图 |
| `rns_fig.{tex,pdf}` | RNS 表示 TikZ 配图 |

## 5. 下一步待办

**FPGA 优化**
- [ ] L2：`#pragma HLS DATAFLOW` 拆 `load_task / compute_task / store_task`，让 Limb 间 Load/Compute/Store 重叠
- [ ] L3：8 Limbs 旋转因子全量预加载到 URAM；NTT/INTT BUTTERFLY 拆分独立循环，消除条件分支重新达 II=1
- [ ] Vivado P&R 实测时序，确认 HLS 估计的 -0.33ns 是否会在后端被吸收

**OC 策略对齐 CiFlow**（详见 [OC_strategy_gap_analysis.md](OC_strategy_gap_analysis.md)）
- [ ] Gap #3 + #4：Section1 bypass + 循环范围改 `[0, sizeQl + sizeP)`（纯 Host 改动，最先验证）
- [ ] Gap #2：OpenFHE 加 `ApproxSwitchCRTBasisSingleTower(target_idx)` 重载，CPU 路径与 FPGA path 都直接受益（FPGA hook 已支持 sizeP=1，HLS 不用改）
- [ ] Gap #5：OC 内层加 tower 累加器，partial sum 在 Host 内存累加再回写
- [ ] 修复后回归 `hks-benchmark --strategy OC`，更新 [实验规划.md 表 5-4](../papers/实验规划.md)，同步 [content.md §5.3](../papers/content.md) 与答辩 PPT Slide 19 论述

**论文实验补充**（见 [../papers/实验规划.md](../papers/实验规划.md) §"缺失实验 / 待补充数据"）
- [ ] 补 4-4：与 CraterLake / BTS / FAB / SHARP 横向对比表（高优先级，纯文献调研）
- [ ] 补 5-1 真实 FPGA 跑通 OC 策略：把 `BConvOffload(sizeP=1)` 走真硬件而非 CPU
- [ ] 补 3-1：综合一版用硬件除法器的 BConv，作为 Barrett LUT 节省的 baseline 对照
