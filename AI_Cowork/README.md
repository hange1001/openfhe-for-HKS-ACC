# openfhe-for-HKS-ACC 项目背景（AI 上下文）

> 本文件是给 AI Agent 看的项目入口。详细进度索引另见 [../docs/notes/PROJECT_STATUS.md](../docs/notes/PROJECT_STATUS.md)。

## 项目概述

- **项目**：openfhe-for-HKS-ACC — 基于 OpenFHE 的全同态加密 (CKKS) Hybrid Key-switching 软硬协同加速器
- **语言**：C++17（OpenFHE 本体）+ HLS C++（FPGA backend）
- **工具链**：CMake / Vitis HLS 2024.1 / XRT / Xilinx Alveo U55C（备用 U250）
- **License**：BSD-2-Clause（继承自 OpenFHE）
- **当前分支**：main

## 两个研究轴

1. **FPGA 微架构**（[src/fpga_backend/](../src/fpga_backend/)）— BConv 脉动阵列 + CG-NTT 常几何蝶形 + Barrett 模乘 + 512-bit AXI 突发
2. **Host 调度**（[src/pke/lib/keyswitch/keyswitch-hybrid.cpp](../src/pke/lib/keyswitch/keyswitch-hybrid.cpp) + [src/pke/include/keyswitch/hks_strategy.h](../src/pke/include/keyswitch/hks_strategy.h)）— 三种 HKS 数据流策略（DC / MP / OC）+ MemoryTracker + HKSStats

## 代码目录结构

```
openfhe-for-HKS-ACC/
├── src/
│   ├── core/                  # OpenFHE 核心（含 FpgaManager.h, dcrtpoly-impl.h 等被改的位置）
│   ├── pke/                   # Public Key Encryption（含 keyswitch-hybrid.cpp）
│   │   ├── examples/hks-benchmark.cpp   # HKS 三策略实测入口
│   │   └── include/keyswitch/hks_strategy.h
│   ├── binfhe/                # （未改动）
│   └── fpga_backend/          # ⭐ HLS kernels + Makefile + .cfg
│       ├── src/               # cg_ntt.cpp / bconv.cpp / top.cpp 等
│       ├── include/           # define.h（PE_PARALLEL, RING_DIM, LIMB_Q 等常量）
│       └── testbench/         # C-sim testbench
├── docs/
│   ├── notes/                 # ⭐ 项目笔记（PROJECT_STATUS.md 是入口）
│   │   └── archive/           # 过时笔记
│   └── papers/                # ⭐ 论文产出（content.md 大纲 / 实验规划.md / 图表）
├── AI_Cowork/                 # ⭐ 本文件夹：AI 协作基础设施
└── benchmark/, test/, third-party/   # （较少改动）
```

## 关键架构决策（详见 [decisions.md](decisions.md)）

- **CG-NTT 取代标准 NTT**：固定几何消除变 stride MUX，单 limb 延迟 7.17× 加速
- **Barrett 模约减取代硬件除法器**：消除 BConv 30 个 urem，释放 ~198K LUT
- **opcode-RPC 调度而非 FPGA 内部 FSM**：粗粒度 FHE 算子计算量远大于 CPU 调度开销（详见 [../docs/papers/fsm_vs_cpu_scheduling.md](../docs/papers/fsm_vs_cpu_scheduling.md)）
- **512-bit AXI 总线**：L1 优化已完成，I/O 8× 缩减

## 关键约定

- **FPGA 参数**（[src/fpga_backend/include/define.h](../src/fpga_backend/include/define.h)）：
  - `RING_DIM = 1 << 12` (4096)、`SQRT = 1 << 6` (64)、`STAGE = 12`、`PE_PARALLEL = 8`
  - `LIMB_Q = 3`、`LIMB_P = 2`、`MAX_LIMBS = 8`
- **opcode**：`OP_INIT=0 / OP_ADD=1 / OP_SUB=2 / OP_MULT=3 / OP_NTT=4 / OP_INTT=5 / OP_BCONV=6 / OP_AUTO=7`
- **Host 端 HKS 策略**：`HKSStrategy::{DC, MP, OC}` 枚举 + `SetHKSStrategy()`
- **不做的事**：不修改 OpenFHE 第三方依赖、不在 Top 顶层裸写 AXI for 循环（必须封装到 `INLINE off` 子函数）
- **协作角色**：详见 [../docs/notes/role.md](../docs/notes/role.md)（数字 IC 架构师 / 苏格拉底导师）

## 当前进度速览（详见 [../docs/notes/PROJECT_STATUS.md](../docs/notes/PROJECT_STATUS.md)）

| 模块 | 状态 |
|------|------|
| CG-NTT 迁移 + bit-reverse 修复 | ✅ |
| BConv Systolic（6.5× 加速） | ✅ |
| 512-bit AXI L1 优化 | ✅ |
| HKS 三策略 Host 端 | ✅ 实测：DC 4.47ms / MP 5.17ms / OC 5.99ms |
| OC 策略对齐 CiFlow 论文 | ⚠️ 当前实现是"DC+sizeP 冗余"，详见 [../docs/notes/OC_strategy_gap_analysis.md](../docs/notes/OC_strategy_gap_analysis.md) |
| L1 算子划分 + L2 三性能模型 | ✅ 推导 v1 完成：[../docs/notes/L1L2_推导v1.md](../docs/notes/L1L2_推导v1.md)（算子账 / 三模型实例化 / F2 口径对账），方法论见[指南](../docs/notes/L1L2_算子划分与性能模型指南.md) |
| L2 DATAFLOW / L3 URAM 预加载 | ⏳ |
| 论文撰写 | ⏳ 大纲完整，部分实验数据待补 |

## 构建与运行

```bash
# OpenFHE + Host 侧（无 FPGA）
cd build && cmake .. -DWITH_FPGA=OFF && make -j hks-benchmark
./bin/examples/pke/hks-benchmark --strategy DC --repeat 100

# 启用 FPGA 链路
cd build && cmake .. -DWITH_FPGA=ON && make -j hks-benchmark

# FPGA HLS 综合
cd src/fpga_backend && vitis_hls -f csynth.tcl

# FPGA xclbin 链接（U55C）
cd src/fpga_backend && make u55c
```

## 给 AI Agent 的提示

- **新对话开始时**：先读 [../docs/notes/PROJECT_STATUS.md](../docs/notes/PROJECT_STATUS.md) 而不是从 git log 重建上下文
- **接到具体任务**：在 [task.yaml](task.yaml) 描述清楚再开始，做完后在 [log.md](log.md) 追加会话记录
- **遇到关键设计抉择**：写一条 ADR 到 [decisions.md](decisions.md)
- **避免做的事**：不要再写一遍 NTT_COMPREHENSIVE_REPORT 类的早期探索文档（已归档在 `docs/notes/archive/`），现在的 NTT 主线是 CG-NTT
