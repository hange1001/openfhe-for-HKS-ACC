# fpga_backend — HKS 加速器的 HLS Kernel

本目录是 [openfhe-for-HKS-ACC](../../README.md) 的 FPGA 侧全部内容：Vitis HLS C++ 写的
kernel、C 仿真 testbench、以及驱动 HLS/Vitis 两条流程的 Makefile 与 tcl 脚本。

- **目标板卡**：Xilinx Alveo U55C（`xcu55c-fsvh2892-2L-e`）
- **目标时钟**：6 ns（约 166.7 MHz），`set_clock_uncertainty 0.75ns`；是否达成须看对应版本 P&R
- **工具链**：Vitis HLS **2023.2** + XRT 2.16（WSL `/tools/Xilinx/Vitis_HLS/2023.2`）
- **HLS 顶层函数**：`Top`（[src/top.cpp](src/top.cpp)）

> 项目整体进度、实测加速比、待办见 [docs/notes/PROJECT_STATUS.md](../../docs/notes/PROJECT_STATUS.md)；
> 架构决策的来龙去脉见 [AI_Cowork/decisions.md](../../AI_Cowork/decisions.md)。

---

## 1. 数据流：Host 怎么调到 Kernel

软硬件边界是一个**极薄的 opcode-RPC**：Host 传一个 opcode + 两块输入缓冲，
FPGA 侧 `Top` 用一个 `switch` 分发；函数内仍有 HLS 生成的控制状态机。
当前以 HKS 为目标，已移除 FPGA AUTO，自同构保留在 OpenFHE CPU 路径。

融合入口为 `EvalKeySwitchPrecomputeCore` → `TryHKSDigitOffload` → C-model 或
`HksDigitTransfer` → `Top(OP_HKS_DIGIT)`。目前融合的是单 digit ModUp，不是完整 HKS。

| 保留指令 | 用途 |
|---|---|
| OP_INIT=0 | 初始化模数、Barrett 参数与双方向 twiddle |
| OP_ADD=1 / OP_SUB=2 / OP_MULT=3 | HKS 相关模运算基础模块与回归入口 |
| OP_NTT=4 / OP_INTT=5 | 共用一套双向 CG_Transform_Banks 引擎 |
| OP_BCONV=6 | 脉动阵列基转换 |
| 编号 7 | 退役保留，不映射到其他指令；Top 不读写外部内存，主机 Execute 抛出异常 |
| OP_HKS_DIGIT=8 | INTT、预缩放、BConv、NTT 融合与原 EVAL digit 旁路 |

删除前的完整版本已保存为 `480dc91`。源码没有 FPGA AUTO 的实现、元数据加载或卸载 API；
历史报告不改写，最新验证见 [移除 AUTO 报告](../../docs/reports/hls/hks_no_auto_20260904/README.md)。

**AXI 端口**（[top.cpp](src/top.cpp)）：`mem_in1`→gmem0、`mem_in2`→gmem1、`mem_out`→gmem2，数据端口均为 256 位；
标量参数走 s_axilite。约定：**顶层不得裸写 AXI 循环**，所有 AXI 访问必须封装进
`#pragma HLS INLINE off` 的子函数，否则 HLS 会在顶层生成巨型 MUX 导致时序崩溃
（这是 [decisions.md](../../AI_Cowork/decisions.md) ADR-007 的结论）。

---

## 2. 文件角色

⚠️ **只有约一半的 `.cpp` 在 `Top` 的调用路径上。** 其余是论文的对照基线，
故意保留但不参与综合 —— 读代码时先看这张表，别从 `bconv.cpp` 入手。

| 文件 | 角色 | 说明 |
|---|---|---|
| [src/top.cpp](src/top.cpp) | 🟢 **产线** | HLS 顶层，opcode 分发 + 片上 buffer / 模数 / 旋转因子的静态存储 |
| [src/cg_ntt.cpp](src/cg_ntt.cpp) | 🟢 **产线** | `CG_Transform_Banks` — 一套运行时双向 NTT/INTT（恒定几何蝶形） |
| [src/bconv_systolic.cpp](src/bconv_systolic.cpp) | 🟢 **产线** | `Compute_BConv_Systolic` — 部署的 BConv（3×5 脉动阵列） |
| [src/arithmetic.cpp](src/arithmetic.cpp) | 🟢 **产线** | `MultMod`（Barrett）/ `AddMod` / `Karatsuba`，被所有 kernel 复用 |
| [src/load.cpp](src/load.cpp) | 🟢 **产线** | `Load` / `Store` — DDR ↔ 片上 `[MAX_LIMBS][SQRT][SQRT]` 搬运 |
| [src/mod_{add,sub,mult}_kernel.cpp](src/) | 🟢 **产线** | 逐元素模加/减/乘 |
| [src/ntt_kernel.cpp](src/ntt_kernel.cpp) | ⚪ **基线** | 标准 NTT。CG-NTT 7.17× 加速比的对照组，**不在 `Top` 路径** |
| [src/bconv_naive.cpp](src/bconv_naive.cpp) | ⚪ **基线** | 朴素 BConv。Systolic 6.5× 加速比的对照组 |
| [src/bconv.cpp](src/bconv.cpp) | ⚪ **中间版本** | 向量内积阵列（**不是**脉动阵列），介于 naive 与 systolic 之间 |
| `Compute_CG_NTT`（[cg_ntt.cpp](src/cg_ntt.cpp) 内） | ⚪ **独立顶层** | 512-bit AXI 包装器，仅供单模块 cosim；`Top` 不调用它 |
| `cg_ntt_reorder`（同上） | ⚪ 仅 testbench | bit-reversal 重排的参考实现 |
| [src/interleave.cpp](src/interleave.cpp) | 🔴 **死代码** | `InterLeave` 无任何调用者，可删 |
| [include/memory.h](include/memory.h) | 🔴 **死代码** | 全文注释掉，仍被 `top.h` include |

**已知的重复定义**（改动时三处都要同步，目前无单一真值源）：

- opcode 在 [include/define.h:57](include/define.h#L57)、[include/opcode.h:11](include/opcode.h#L11)、
  [FpgaManager.h:41](../core/include/FpgaManager.h#L41) 各定义一遍（注意 `define.h` 叫 `OP_MULT`，`opcode.h` 叫 `OP_MUL`，同值不同名）
- `LIMB_Q` / `MAX_OUT_COLS` 在 [define.h:37](include/define.h#L37)、[FpgaManager.h:611](../core/include/FpgaManager.h#L611)、
  [dcrtpoly-impl.h:1171](../core/include/lattice/hal/default/dcrtpoly-impl.h#L1171) 各写一遍

---

## 3. 参数

全部在 [include/define.h](include/define.h) —— **改参数只改那里，本文不复制数值**，
避免又一处会过期的副本。当前口径（`RING_DIM` / `LIMB_Q` / `LIMB_P` / `PE_PARALLEL` 等）以该文件为准。

派生关系值得单独记一笔：

| 常量 | 定义处 | 由什么决定 |
|---|---|---|
| `MAX_LIMBS` | define.h | `LIMB_Q + MAX_OUT_COLS`，**不是**独立可调参数 |
| `MAX_OUT_COLS` | define.h | `LIMB_Q + LIMB_P`，脉动阵列的列数 |
| `CG_HALF_N` | [cg_ntt.h](include/cg_ntt.h) | `RING_DIM / 2`，蝶形跨度 |
| `CG_PE_NUM` | cg_ntt.h | `= PE_PARALLEL` |
| `CG_BUF_PARTITION` | cg_ntt.h | `= 2 * PE_PARALLEL`，由**写端口**决定（每 PE 写 2 个连续地址），推导见该处注释 |

> 注意 [bconv_systolic.cpp](src/bconv_systolic.cpp) 的 `LOAD_PAR = 8` 是**独立于 `PE_PARALLEL`** 的量：
> 它表示局部缓存搬运的 8×64-bit 并行量，不代表 Top 外部 AXI 位宽，不要跟着 `PE_PARALLEL` 一起改。

---

## 4. 构建与验证

两条独立流程。日常迭代用 HLS 流程（快），出 xclbin 用 Vitis 流程（慢）。

### HLS 流程：单模块 csim / csynth / cosim

`MODULE` 取值 = HLS 顶层函数名，注册表在 [Makefile](Makefile) 的「模块注册表」一节。

```bash
make csim   MODULE=Compute_CG_NTT          # C 仿真，验证算法正确性
make csynth MODULE=Top                     # C 综合，看资源与时序
make cosim  MODULE=Compute_BConv_Systolic  # C/RTL 协同仿真 + 波形
```

已注册的 `MODULE`：`Top`、`Compute_CG_NTT`、`NTT_Kernel`、`Compute_NTT`、
`Compute_BConv`、`Compute_BConv_Naive`、`Compute_BConv_Systolic`。
拼错会被 `check_module` 挡下并提示。

工程目录（两条流程**不共用**，`make report` 已按此适配）：

- `make csynth` → `Solution/$(MODULE)/solution1/`
- `make cosim` → `Solution/cosim_$(MODULE)/solution1/`

把报告提拔进 git（`docs/reports/hls/` 是有意跟踪的）：

```bash
make report MODULE=Compute_BConv_Systolic
```

`make report` 会经 [trim_report.sh](trim_report.sh) 处理，产出两样东西：

| 产物 | 内容 |
|---|---|
| `docs/reports/hls/<MODULE>_csynth.rpt` | 裁剪版。丢掉 `SW I/O Information` 与 `Bind Op Report`（后者单独占 ~56%，是逐算子 DSP 绑定明细，diff 时全是噪声），保留资源/时序、AXI 接口、Storage、Pragma 四段 |
| `docs/reports/summary.csv` | **一个模块一行**的顶层数字（slack / BRAM / DSP / FF / LUT / URAM）。跨 commit 看 delta 只 diff 这个文件 |

`summary.csv` 的 `clk_ns` 由 `latency_ns / latency_cycles` 反推——用来抓时钟口径漂移
（历史报告有 5ns 综合的，当前工程是 6ns）。要调整保留哪几段：`HLS_RPT_DROP="段名1|段名2" make report ...`

### Vitis 流程：出 xclbin

```bash
make sw_emu    # C model 仿真
make hw_emu    # RTL 仿真
make hw        # 综合出比特流
```

### 切换板卡

器件型号由 Makefile 顶部的 `PLATFORM` + `HLS_PART` **单点驱动**，`HLS_PART` 通过环境变量
传给三个 tcl 的 `set_part`。切板卡时两个变量一起改即可，不要再去改 tcl。

---

## 5. Testbench

在 [testbench/](testbench/) 下，均为 HLS C 仿真（`make csim`），不依赖真实硬件。

| 文件 | 覆盖 |
|---|---|
| [cg_ntt_tb.cpp](testbench/cg_ntt_tb.cpp) | CG-NTT 蝶形、单 limb 往返、`Compute_CG_NTT` 多 limb、重排 |
| [ntt_kernel_tb.cpp](testbench/ntt_kernel_tb.cpp) | 基线 NTT：`Configurable_PE`、`NTT_Kernel`、`Compute_NTT` |
| [bconv_tb.cpp](testbench/bconv_tb.cpp) | `Compute_BConv` / `Compute_BConv_Naive`（用 `-DBCONV_USE_NAIVE` 切换，共用用例） |
| [bconv_systolic_tb.cpp](testbench/bconv_systolic_tb.cpp) | 脉动阵列 BConv |
| [top_tb.cpp](testbench/top_tb.cpp) | 顶层 opcode 分发全路径 |
| [auto_test.cpp](testbench/auto_test.cpp) | 自同态 |

⚠️ 例外：[host_overhead_bench.cpp](testbench/host_overhead_bench.cpp) **不是 HLS testbench**，不参与 `make csim`。
它是实验 D（task.yaml step 1.2）的离线分量——量 `bo.write/read` 的 memcpy 与 BConv hook 的
host 侧 gather/scatter 代价，不需要板卡也不需要 XRT：

```bash
g++ -O2 -std=c++17 -o /tmp/hob testbench/host_overhead_bench.cpp && /tmp/hob
```

---

## 6. 修改本目录时的约定

1. **顶层不裸写 AXI**。新增 opcode 分支时，把 AXI 读写封进 `#pragma HLS INLINE off` 子函数
   （照抄 `Load_Init_Params` / `Load_BConv_Params` 的写法）。
2. **partition / unroll factor 引用具名常量**，不写字面量。需要非 `PE_PARALLEL` 的倍数时，
   像 `CG_BUF_PARTITION` 那样先定义一个有推导注释的派生常量。
3. **注释里写死数值前先想清楚**。本目录曾出现整段注释按 `PE_PARALLEL=8` 描述而代码已是 4 的情况，
   比没注释更误导。优先写「为什么是这个值」而不是「这个值是多少」。
4. **动了基线文件（⚪）要同步论文**：那几个加速比数字直接引自它们的综合/仿真报告。
