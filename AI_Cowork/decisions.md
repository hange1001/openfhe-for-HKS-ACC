# 架构决策记录 (ADR)

> 记录项目中已做的关键设计决策，避免 AI 或新贡献者重复讨论相同问题。
> 格式：日期 + 决策内容 + 原因 + 影响范围。
> **最新决策在最前面**，编号从大到小逆序排列。

---

## ADR-010: OC 策略采用 `Section1 bypass → single-tower BConv → on-chip 累加` 的渐进式修复顺序

- **日期**: 2026-06-29
- **背景**: 当前 OC 实测 5.99 ms > DC 4.47 ms，并未达到 CiFlow 论文宣称的"延迟最优"。根因是当前实现仍调用全 sizeP 的 `ApproxSwitchCRTBasis()` 然后只挑 1 个 tower，等价于 sizeP 倍冗余 BConv 计算。`docs/notes/OC_strategy_gap_analysis.md` 把差距拆成 5 条 gap（#1 INTT 复用已做、#2-#5 未做）。
- **决策**: 按以下顺序逐条修复：
  1. **Gap #3 + #4（一起做）**：纯 Host 改动。给 OC 加 Section1 bypass（算 mod Q 输出时第 `p/alpha` 个 digit 跳过 BConv 直接复用原 EVAL 系数）+ 把外层循环从 `for p in sizeP` 改为 `for p in [0, sizeQl + sizeP)`
  2. **Gap #2**：给 OpenFHE 加 `ApproxSwitchCRTBasisSingleTower(target_idx)` 重载，让 BConv 只算 1 个 output tower
  3. **Gap #5**：ModUp P5 的 partial sum 在 Host 内存累加，最后一次性回写
- **原因**:
  - **先 #3+#4 是因为零库改动**，最快验证逻辑正确性。Section1 bypass 的收益是省 `1/dnum` 的 BConv+NTT
  - **#2 才是最大收益**（BConv 计算量 5×、D2H 流量 5× 缩减），但需要新增 OpenFHE API，改动面更大
  - **#2 的 FPGA 侧零成本**：`FpgaManager.h:614` 的 hook 判断 `sizeP <= KERNEL_MAX_OUT_COLS=5` 已包含 sizeP=1，HLS 不用改
  - **#5 是次要优化**，先做不动主体逻辑
- **影响**: `src/pke/lib/keyswitch/keyswitch-hybrid.cpp` OC 分支；`src/core/include/lattice/hal/default/dcrtpoly{.h,-impl.h}`；后续 `docs/papers/实验规划.md` 表 5-4 与 `docs/papers/content.md` §5.3 论述同步
- **替代方案**:
  - **一次性把 #2-#5 全做**：风险大，回归时难二分定位
  - **直接做 #2 跳过 #3/#4**：BConv 计算节省了，但 Section1 bypass 的 1/dnum 收益丢失
- **参考**: [../docs/notes/OC_strategy_gap_analysis.md](../docs/notes/OC_strategy_gap_analysis.md)（5 条 gap 完整定义）、[task.yaml](task.yaml) `oc-align-with-ciflow`

---

## ADR-009: 项目笔记结构化为 PROJECT_STATUS.md 入口 + archive/ 归档

- **日期**: 2026-06-29
- **决策**: 在 `docs/notes/` 根建立 `PROJECT_STATUS.md` 作为唯一笔记入口；过时的 8 份旧 NTT 笔记 `git mv` 到 `docs/notes/archive/`，保留可追溯但不再维护；新 AI 会话从 PROJECT_STATUS.md 开始而非从 git log 重建上下文
- **原因**:
  - `docs/notes/` 累积了 13 份笔记，5 份是关于已废弃的标准 NTT 实现（NTT_COMPREHENSIVE_REPORT 等），AI 接手时浪费上下文窗口
  - 直接删除会丢失早期问题分析的可追溯性（如 HW_XCLBIN_ANALYSIS 的 13 个 P0/P1 综合失败诊断）
  - git mv 保留历史，archive/README.md 标注每份的归档原因和被谁替代
- **影响**: `docs/notes/` 目录结构；后续 AI 协作的"先读什么"约定
- **替代方案**:
  - **git rm 删除**：彻底但破坏可追溯性
  - **保持原样**：上下文窗口压力持续，未来更难收拾

---

## ADR-008: opcode-RPC 调度模型，不在 FPGA 内部建顶层 FSM

- **日期**: 2026-05-20
- **决策**: 当前 `top.cpp` 是 opcode-driven RPC dispatcher（CPU 逐个发 opcode），不演化为 FPGA 内部完整流水 FSM。下一阶段如果有需要，最多做"复合 opcode"（如 `OP_NTT_MULT_INTT_CHAIN`）
- **原因**:
  - 决策判据：开销比 = CPU 调度开销 / 单次计算周期。NTT (~300K cyc) / BConv (~500K cyc) 等粗粒度算子的开销比 < 1%，CPU 调度无感
  - FPGA FSM 灵活性差：改顺序要重综合（小时级），CPU 改一行源码秒级
  - 复合 opcode 是阶段 2 的性价比最高路径——不动 top.cpp 顶层结构，只在 `switch(opcode)` 中加新 case 串联已有子函数
- **影响**: `src/fpga_backend/src/top.cpp` 整体架构；Host 端 `FpgaManager` API 的设计风格；未来 Batch HKS 也按这个原则走
- **替代方案**:
  - **完整 FPGA FSM**：极致性能，但每次改流程都要重综合，灵活性丧失
  - **完全 CPU 控制 + 细粒度 op**：通信开销 > 计算开销，开销比爆炸
- **参考**: [../docs/papers/fsm_vs_cpu_scheduling.md](../docs/papers/fsm_vs_cpu_scheduling.md)

---

## ADR-007: 拆解顶层 FSM 为子函数 + `#pragma HLS INLINE off`

- **日期**: 2026-04-26
- **决策**: 把 `OP_INIT` / `OP_BCONV` 中直接裸写的 AXI for 循环全部提取到独立子函数（`Load_Init_Params` / `Load_BConv_Params` / `Execute_NTT` / `Execute_INTT`），并打 `#pragma HLS INLINE off` 强制不内联
- **原因**:
  - 综合后顶层 `Top` 模块 Slack = −0.33 ns，根因是 HLS 为多个 case 共享的 AXI master 生成巨型 MUX 导致控制信号扇出爆炸
  - INLINE off 让 HLS 把 AXI 握手封装在子模块内部，顶层只剩纯函数调用 dispatcher，时序大幅缓解
  - 同时降低了 CG-NTT 蝶形并行度（从 16 PE 降到 8 PE）以匹配可承受的布线密度
- **影响**: `src/fpga_backend/src/top.cpp`；`src/fpga_backend/src/cg_ntt.cpp` 并行度
- **替代方案**:
  - **降时钟到 5ns (200MHz)**：放弃 4ns (250MHz) 目标，已采用作为兜底
  - **不拆函数靠 Vivado Implementation 自己 retiming**：风险大，HLS 估计的 -0.33ns 在后端不一定收敛
- **参考**: [../docs/notes/archive/fix.md](../docs/notes/archive/fix.md)
- **后续修订（2026-08-19）**: 时钟最终降到的是 **6ns (166 MHz)**，不是本条写的 5ns —— 见
  `csynth.tcl` 的 `create_clock -period 6ns`。本条记录保留原样（ADR 是历史快照），
  但引用时以 tcl 为准。同日重跑 `make csynth MODULE=Top` 的基线：
  Slack **−0.33 ns** / Fmax 179.19 MHz / LUT 14% / FF 2% / BRAM 15% / DSP 14% / URAM 10%，
  顶层 FSM 最大扇出 **5496**（本条说的"扇出爆炸"仍在）。数字见
  [../docs/reports/summary.csv](../docs/reports/summary.csv)。

---

## ADR-006: Barrett 模约减替换 BConv 中的 128-bit 硬件除法器

- **日期**: 2026-04-26
- **决策**: BConv 脉动阵列 PE 中的 `(x * w) % mod_p` 全部用 MultMod（Barrett 算法 II=1 流水线）替换，删除所有 `%` 运算
- **原因**:
  - 原 `bconv.cpp` 用 `ap_uint<128> % mod_p` 让 HLS 综合出 30 个 `urem_128ns_64ns` 硬件除法器
  - 每个除法器 ~8651 FF + 6607 LUT → 单此一处占 LUT 15% / FF 10%，是时序违例的另一重要来源
  - Barrett 算法的 II=1 流水线已在 NTT 蝶形单元中验证可用，复用同一模块
- **影响**: `src/fpga_backend/src/bconv.cpp`；`src/fpga_backend/src/arithmetic.cpp` 的 MultMod 模块；BConv 接口（需要从 Host 端传入 Barrett 常数 `k_half` 和 `m`）
- **替代方案**:
  - **保留 `%` 等 Vivado 后端优化**：尝试过，LUT 仍爆炸
  - **Montgomery 约减**：精度更高但需要更多预计算，性价比不如 Barrett

---

## ADR-005: CG-NTT bit-reverse 输出由 Host 软件补丁修复，不重综合 HLS

- **日期**: 2026-04-26
- **决策**: CG-NTT 经 STAGE 次 perfect shuffle 后输出 bit-reversed 顺序，与 OpenFHE 期望的标准顺序不兼容。修复方案是在 Host 端 `FpgaManager::NttForwardOffload`（出口）和 `NttInverseOffload`（入口）做 bit-reverse 重排，不修改 HLS 内核
- **原因**:
  - Bit-reversal 是自反操作（执行两次=恒等），Host 拦截即可
  - 单次 N=4096 的 bit-reversal ~10 μs，远小于 PCIe 传输延迟（>40 μs），开销可忽略
  - 改 HLS 内核要重综合（小时级），影响其他正在进行的优化工作
- **影响**: `src/core/include/FpgaManager.h` 的 `NttForwardOffload` / `NttInverseOffload`；CKKS 端到端测试（密文加密、EvalMult、EvalRotate 全部依赖）
- **替代方案**:
  - **硬件内嵌 reorder**：可在 STORE_OUT 阶段以零额外延迟完成，但需重综合 CG-NTT 内核；保留为未来 L3 优化的子项
- **参考**: [../docs/notes/CG-NTT-Migration-Report.md](../docs/notes/CG-NTT-Migration-Report.md) §5

---

## ADR-004: 512-bit AXI 总线拓宽（L1 优化）

- **日期**: 2026-04-23
- **决策**: 把 `Compute_CG_NTT` 顶层接口的三个大数组参数从 `uint64_t` 多维数组改为 `ap_uint<512>*` 一维指针，强制 HLS 推断出 512-bit burst
- **原因**:
  - 综合报告中 m_axi 全部退化为 64-bit（"Could not widen since type i64 size is greater than or equal to alignment"），实际只用了 U55C HBM 12.5% 的带宽
  - 改前 LOAD_IN 4097 cycles → 改后 513 cycles（8× 缩减），单 Limb 总延迟从 43,316 降到 14,644（3× 加速）
- **影响**: `src/fpga_backend/src/cg_ntt.cpp` 接口；Host 端 buffer 准备代码；BRAM 略增（+78 块，FIFO 缓冲）
- **替代方案**:
  - **保持多维数组靠 HLS 自动推断**：已证明失败，HLS 拒绝拓宽
  - **手动写 burst pragma**：等价但代码可读性差
- **参考**: [../docs/notes/cg_ntt_optimization.md](../docs/notes/cg_ntt_optimization.md) §3

---

## ADR-003: CG-NTT（恒定几何 NTT）取代标准变几何 NTT

- **日期**: 2026-04-15
- **决策**: 引入 `cg_ntt.cpp` 新内核取代原 `ntt_kernel.cpp` 的标准变几何 NTT 作为主力 NTT 实现
- **原因**:
  - 标准变几何 NTT 每 stage 跨度变化，HLS 生成大量 MUX 和动态地址生成，布线复杂度高
  - 标准 NTT 综合报告显示 BRAM 超标 206%（旋转因子 [MAX_LIMBS][PE_PARALLEL][RING_DIM] = 2.1 MB×2）
  - CG-NTT 通过 perfect shuffle 让每 stage 的读写地址固定（i, i+N/2 → 2i, 2i+1），消除 MUX，地址生成硬连线
  - 单 limb 实测延迟：变几何 112,671 cycles → CG-NTT 15,701 cycles（7.17× 加速），DSP -75% / FF -89% / LUT -57%
- **影响**: `src/fpga_backend/src/cg_ntt.cpp`（新）；`src/fpga_backend/src/top.cpp` 调度路径替换；Host 端 `BuildCGTwiddle` 取代 `GenerateTwiddleIndices`；旋转因子布局从 `[MAX_LIMBS][PE_PARALLEL][RING_DIM]` 变为 `[MAX_LIMBS][STAGE][CG_HALF_N]`
- **替代方案**:
  - **优化标准 NTT 内 MUX**：试过删 BU_NUM 冗余维度等，能解 BRAM 但 MUX 复杂度根因没动
  - **混合基 NTT (mixed radix)**：实现复杂度高，无明确收益
- **参考**: [../docs/notes/CG-NTT-Migration-Report.md](../docs/notes/CG-NTT-Migration-Report.md)

---

## ADR-002: PE_PARALLEL 全局常量统一所有 partition factor + TwiddleFactor PE 副本 URAM

- **日期**: 2026-04-14
- **决策**: 在 `define.h` 新增 `static const int PE_PARALLEL = 8`，让所有 `ARRAY_PARTITION cyclic factor` / `UNROLL factor` / Twiddle 副本数都引用此常量；旋转因子建立 PE_PARALLEL 个物理副本，绑定 URAM
- **原因**:
  - 早期 partition factor 散落各文件，调整并行度要改 10+ 处
  - 旋转因子降维为 1D 后 8 个 PE 读相同 TwiddleIndex 触发 bank 冲突，II 被拖到 4
  - 建立 PE_PARALLEL 个独立副本（每个 PE 读自己的副本）即使地址相同也零冲突；URAM 容量充裕（4MB 占 U55C 12%）
- **影响**: `src/fpga_backend/include/define.h`；`top.cpp` / `ntt_kernel.cpp` / `interleave.cpp` / `mod_*_kernel.cpp` 所有 partition pragma
- **替代方案**:
  - **保持 cyclic factor 散落**：维护成本高
  - **TwiddleFactor 不复制单副本**：bank 冲突无解，II > 1
- **参考**: [../docs/notes/archive/综合后修改方案.md](../docs/notes/archive/综合后修改方案.md) FIX-B

---

## ADR-001: HKS 三策略（DC / MP / OC）作为本课题的研究框架

- **日期**: 2026-03-28
- **决策**: 引入 `HKSStrategy::{DC, MP, OC}` 枚举 + `SetHKSStrategy()` 全局切换；`keyswitch-hybrid.cpp::EvalKeySwitchPrecomputeCore` 内分三个分支实现；`HKSStats` + `MemoryTracker` 提供操作计数、延迟、峰值缓冲采集；DC 为默认（与未修改 OpenFHE 行为一致）
- **原因**:
  - 直接对标 CiFlow (Neda et al., ISPASS 2024) 论文的三种 dataflow 分类，便于在论文 §5 章作差异化对比
  - DC 是 OpenFHE 原生路径（per-digit INTT→BConv→NTT），改动最少作为 baseline
  - 三策略走同一 FPGA kernel（不重综合），仅 Host 端调度顺序不同——符合 ADR-008 的 opcode-RPC 模型
- **影响**: `src/pke/include/keyswitch/hks_strategy.h`（新）；`src/pke/lib/keyswitch/keyswitch-hybrid.cpp` 大改动；`src/pke/examples/hks-benchmark.cpp`（新）作为命令行入口；本课题的论文第 5 章主线
- **替代方案**:
  - **只实现 OC**：失去对比基准
  - **每个策略一份 FPGA bitstream**：违反 ADR-008，灵活性丧失
- **参考**: [../docs/papers/content.md](../docs/papers/content.md) §5；[../docs/papers/实验规划.md](../docs/papers/实验规划.md) §补 5-1

---

<!-- ADR 模板（新增决策时复制此段到本文件最顶部，编号 +1）
## ADR-NNN: 决策标题

- **日期**: YYYY-MM-DD
- **决策**: 简述做了什么决定
- **原因**:
  - 原因1
  - 原因2
- **影响**: 影响哪些模块/文件
- **替代方案**: 曾考虑但未采用的其他方案（如有）
- **参考**: 相关笔记 / 论文章节 / commit
-->
