# HKS dataflow 软件仿真器

实施方案：[`doc/HKS_dataflow_simulation_plan.md`](../../doc/HKS_dataflow_simulation_plan.md)

在**完全相同的硬件资源约束**下比较 DC / MP / OC 三种 dataflow。不修改任何 RTL。

标定基线 **P4 = `e948a69`**，证据目录 `docs/reports/hls/hks_mem_p4_20260905/`。

---

## 快速开始

```bash
cd tools/hks_flow_sim

# P4 复现（校准点）：DC 必须精确等于实测 91242 cycles
python3 hks_sim.py --config configs/p4_reproduction.yaml

# 完整 HKS 边界 + output tile 全宽度对比（M4 Pareto）
python3 hks_sim.py --config configs/p4_reproduction.yaml \
                   --config configs/keymult_nominal.yaml \
                   --boundary full_hks --invocation per_hks --oc-tile-sweep

# M4 完整分析：Pareto + 容量边界曲线 + 敏感性
python3 analyze_m4.py --output results/m4

# 测试
python3 -m unittest discover -s tests -t .
```

`--config` 可重复，按顺序深合并（base 配置 + KeyMult 档位配置）。

---

## 标定与验证

### P4 零残差

`transform=6481`、`scale=1047`、`bconv=4145` 三个数直接来自 HLS 报告，不是拟合。
只留两个自由参数，用两个 digit 事务解出：

| 参数 | 值 | 来源 |
|---|---:|---|
| `axi_tower_cycles` | 1053 | = 4096·64/256 拍 + 29，**calibrated** |
| `control_per_digit` | 656 | 每次 kernel 调用固定开销，**calibrated** |

```
alpha=2: 5×6481 + 2×1047 + 4145 + 7×1053 + 656 = 46671   实测 46671
alpha=1: 5×6481 + 1×1047 + 4145 + 6×1053 + 656 = 44571   实测 44571
```

### P3 hold-out

P3 的预乘是**独立硬件、1 系数/拍**（P4 才改成复用 4 路 lane），结构与 P4 不同。
把 P4 标定出的 `axi`/`control` 原样搬过去，只解预乘一项：两个事务分别给出
**4101** 和 **4119**，互差 0.4%。取中值反代，两 digit 合计预测 100431 vs
实测 100422，误差 **+0.009%**（M3 门槛 5%）。顺带印证：4110 ÷ 1047 = 3.93 ≈ lane 数 4。

断言在 `tests/test_cost_model.py`。

---

## 相对方案的四处改动

方案有四个地方会**静默改变结论**，且其中三个**标定点分辨不出来**。

### 1. tie-break 必须固定并计入 hash

方案 §5.2 没规定多个 ready 事件争同一引擎时选谁。这个自由度会系统性偏袒某种
dataflow（LIFO 偏袒 DC 的深度优先，FIFO 偏袒 MP 的广度优先）。默认
`tie_break="trace_order"`（严格按 trace 发射序，唯一不引入额外自由度的选择），
并计入 `hardware_config_hash`。

*当前结果对它不敏感*（三种 tie-break 排序一致），但这是**测出来的**，不是假设。

### 2. `allow_engine_overlap` —— 会翻转结论

P4 的 Top 是 opcode-driven RPC dispatcher（ADR-008），BConv 与变换引擎由顶层
**顺序调用**，没有 DATAFLOW。方案 §5.2 把它们列为两个独立容量池，隐含了可并发。

DC 的单次 BConv 必然先于它的全部 NTT，永远不会有两个引擎同时 ready，所以
**P4 的 91242 在开关两种取值下都能对上**；OC 的 single-tower 粒度会让 BConv/NTT 重叠。

| | 排序 |
|---|---|
| `False`（默认，P4 忠实） | dc < mp < oc-w5 < oc-w3 < oc-w4 < oc-w2 < oc-w1 |
| `True`（假设已做 DATAFLOW） | **oc-w1 < oc-w2 < oc-w3 < oc-w4 < oc-w5 < dc < mp** |

完全翻转。主结论固定用 `False`，`True` 只作敏感性结果。

### 3. AXI 是单一共享池，容量 = `axi_channels`（默认 1）

**最初的建模是错的**：按用途拆成 `dma_h2d` / `dma_d2h` / `dma_spill` 三个各自
容量 1 的池，等于白送 3 倍并发带宽，而且按「用途」划分通道本身就是臆造的
——真实绑定是按 gmem 端口。

实测能证伪它：P4 单 digit 事务 46671 只有在 `alpha` 个 load 塔加 `L+K` 个
store 塔**全部串行**（7×1053）时才对得上；任意两笔重叠，实测值都会更小。

修正前后的差别不小，两个此前报告过的"现象"都是这个 bug 的产物：

- DC 靠增量回写藏掉 8424 拍传输而 MP 不能 → 修正后差值归零，DC = MP = 145566；
- OC 的 latency 在 w=4 处非单调 → 修正后单调不增；
- tie-break 敏感（fifo/lifo 会换排序）→ 修正后三种 tie-break 排序一致。

（收尾报告里 INIT 的两张 twiddle 表确实在 gmem0/gmem1 上重叠，但那是另一条
代码路径，且 295063 是直接实测的常数，整体计入，不由本池建模。）

### 4. 调用边界是串行的

XRT 路径完全阻塞（`bo.sync → run.wait → bo.sync`，无 async、无双缓冲，
task.yaml q8），下一次 kernel 调用不能在上一次返回前开始。不加这道 barrier，
两个 digit 会算成 89533 而不是实测的 91242。

`invocation` 粒度（`per_hks` / `per_digit` / `per_output_tower`）可配置：
方案 §12 点名的「OC 每个 target tower 一次调用会被固定开销吃掉收益」
在这里是**可证伪的假设**，不是写死的结论。

---

## 一条必须先看清的恒等式

```
sum_d beta_d = D(L+K) - L = L(D-1) + KD
```

左边是 DC 逐 digit 的 BConv 输出塔总数，右边是方案 §6.3 给 OC 的闭式——**它们相等**。
INTT 也都是 L 次，NTT 也都是 L(D-1)+KD 次。

**三种 dataflow 的算术量完全相同**，差异 100% 来自 lifetime / 存储 / 搬运 / 调度。
与推导v1 §2.3「每算子一次 offload 下三策略 AI 恒等」是同一件事。

`tests/test_trace_counts.py` 断言三策略的 INTT/SCALE/NTT/KeyMult/Accumulate
计数**逐项相等**，BConv 输出塔总数也相等。允许不同的只有 BConv 的**调用次数**。

---

## M4：output-tiled OC

固定 3×5 阵列下，真正缺的不是把阵列变窄，而是把 OC 从 single-tower 扩成
**output-tiled**：同一套硬件，只改调度粒度。

```
OC-w1  一次处理 1 个 output tower（= 原始真 OC 基线）
...
OC-w5  一次最多利用 5 个 BConv columns
```

一个 tile 内同时分配 `2w` 个 accumulator；对每个 digit，把该 tile 中需要 BConv 的
target towers **合并成一次调用**，成本按实际 `len(non_native_targets)`；
native Q towers 继续 bypass；KeyMult/Accumulate 仍逐 tower 执行。

`oc_output_tile_width` 是**调度**参数：不进 `hardware_config_hash`
（否则 OC-w1..w5 之间会显示成"配置不同"而失去可比性），单独进
`strategy_config_hash` 与结果元数据。

**结果一律写 `OC-w1`～`OC-w5`，不会静默挑最好的那个写成"OC"。**

### P4 参数（N=4096, L=3, K=2, D=2；KeyMult nominal）

| 策略 | cycles | BConv 调用 | peak | acc 塔 | Pareto |
|---|---:|---:|---:|---:|---|
| dc | **145566** | 2 | 576 KB | 10 | front |
| mp | 145566 | 2 | 1600 KB | 10 | dominated |
| oc-w1 | 166291 | 7 | **384 KB** | 2 | front |
| oc-w2 | 158001 | 5 | 480 KB | 4 | front |
| oc-w3 | 153856 | 4 | 544 KB | 6 | front |
| oc-w4 | 153856 | 4 | 640 KB | 8 | dominated |
| oc-w5 | 145566 | 2 | 736 KB | 10 | dominated |

Pareto front = **{dc, oc-w1, oc-w2, oc-w3}**。

- **DC 与 MP 的算子数量和理论 arithmetic workload 相同，当前建模下周期也相同；
  但 MP 需要 2.8 倍 peak memory**（barrier 把中间塔生命周期拉长到跨阶段）。
  周期相等只在串行 AXI 下成立，`axi_channels=3` 时两者立刻分开——
  `tests/test_mp_dc_delta.py` 同时锁住这两条。
- w=5 时 BConv 调用退化到与 DC 相同（2 次），compute 也相等；代价是 accumulator
  从 2 塔涨到 10 塔，于是被 DC 支配。
- w3 == w4：两者 tile 数相同（`[0,1,2][3,4]` vs `[0,1,2,3][4]`），BConv 调用都是 4。

### 工业参数（N=16384, L=12, K=4, D=4）

| 策略 | cycles | BConv 调用 | peak | Pareto |
|---|---:|---:|---:|---|
| dc | **3623460** | 4 | 6528 KB | front |
| mp | 3623460 | 4 | **36864 KB** | dominated |
| oc-w1 | 4281172 | 52 | **3840 KB** | front |
| oc-w2 | 3886780 | 28 | 4224 KB | front |
| oc-w3 | 3755316 | 20 | 4608 KB | front |
| oc-w4 | 3689584 | 16 | 4992 KB | front |
| oc-w5 | 3689584 | 16 | 5376 KB | dominated |

一个结构性上限：`L+K=16 > bconv_cols=5`，所以 OC **永远达不到** DC 的 4 次
BConv 调用，最好是 16 次（w=4 或 5），仍是 DC 的 4 倍。**当 L+K 超过阵列列数时，
output tiling 只能压缩 OC 的冗余，不能消除它。**

MP 的 peak 达 36 MB，对照 U55C 全片约 43 MiB 原始位容量（可用数据容量更低），
基本不可行。

### 容量边界：交叉点在这里

以一个 tower 为步长扫描（P4 tower=32 KB / 工业 tower=128 KB）：

| 容量 | P4：dc | P4：最好的 OC | 工业：dc | 工业：最好的 OC |
|---|---:|---:|---:|---:|
| 无约束 | **145566** | 153856 (w3) | **3623460** | 3689584 (w4) |
| 最紧可行 | 228097 | **187351** (w1) | 5524429 | **4010410** (w3) |

**容量宽裕时 DC 赢；容量收紧后 OC 反超**——P4 下 OC-w1 快 17.9%，
工业参数下 OC-w3 快 **27.4%**。这正是方案 §12 待验证判断里
「当 (L,K,D) 增大并使 DC/MP 超过片上容量时，OC 优势可能扩大」的定量版本。

### 敏感性（非主结论）

- **tie-break**：三种排序完全一致，结论对调度优先级**不敏感**。
- **`allow_engine_overlap=True`**：排序完全翻转，OC 全面领先 DC。
  这是「做没做 DATAFLOW 拆分」的问题，不是 dataflow 的问题。

⚠️ 全部为 RTL 周期下界换算，**不含 PCIe / 驱动 / buffer 分配 / Host 调度**，无板卡。

---

## PCIe 参数怎么处理

方案 §4.2 把 `pcie_fixed_us` / `pcie_bandwidth_GBps` 标为"待板卡实测"。但没有板卡，
且 task.yaml q8 已确认现有的 182.8 μs/call **是下界不是墙钟**（漏计 `xrt::bo`
分配、`bo.write/read` 的 host memcpy 与析构）。

所以它们是**扫描维，不是标定量**，默认 0（不建模）：

```bash
--sweep hardware.pcie_fixed_us=0,50,182.8,400
```

---

## 目录结构

```
tools/hks_flow_sim/
  hks_sim.py            CLI
  analyze_m4.py         M4 分析：Pareto + 容量曲线 + 敏感性
  hks_sim/
    config.py           schema / 校验 / 三种 config hash
    workload.py         digit 布局推导（alpha_d, beta_d, D）
    op.py               Op / Buffer / TraceBuilder
    resources.py        资源池与绑定（AXI 为单一共享池）
    cost_model.py       周期模型（measured / calibrated / projected）
    memory.py           lifetime / peak / spill 容量 pre-pass
    scheduler.py        事件驱动调度器 + makespan 四分解
    trace_common.py     三种 trace 共用的发射器
    trace_dc.py         Digit-Centric
    trace_mp.py         Maximum-Parallel
    trace_oc.py         Output-Centric（含 output tiling）
    engine.py           编排（方案 §8 没有这一层，从 report.py 分出来的）
    report.py           CSV / JSON / 事件 trace / 文本摘要
  configs/              p4_reproduction / industrial_projection / keymult 三档
  tests/                90 个断言
```

---

## 里程碑状态

| | 状态 | 说明 |
|---|---|---|
| M0 配置与事件规格 | ✅ | (alpha_d, beta_d, D) 推导 + 非法参数拒绝 + config hash |
| M1 静态 trace 与存储统计 | ✅ | 不依赖周期参数即可出 operation-count / memory 报告 |
| M2 资源约束调度 | ✅ | ready queue / 引擎互斥 / DMA overlap / makespan 四分解 |
| M3 P4 校准与 KeyMult 区间 | ✅ | P4 零残差、P3 hold-out 0.009%、三档区间 |
| M4 Tiling 与 spill | ✅ | OC-w1..w5、容量边界曲线、Pareto、两档敏感性 |
| M5 Functional reference | ✅ | 真 OC 的 test-only C++ reference，8 组形状 524288 residue 精确一致 |
| M6 Industrial projection | ⬜ | 扫参框架就绪，尚未系统跑 |

M4 尚未做的：BConv 的**行**方向分块（`ceil(alpha/rows)`，当前所有参数点都是 1），
以及 `allow_tiling` 目前只做校验、未改变 BConv 成本模型。

### M5 已完成：真 OC 的功能正确性

`src/pke/unittest/UnitTestHKSTrueOC.cpp`（**只供测试使用**，不参与 production dispatch）。
它按真 OC 的顺序重算 `cTilda0/cTilda1`：

```
for output tower t in [0, L+K):
    for digit j:
        native Q tower -> bypass（直接复用原 EVAL 系数）
        non-native     -> 调用已有 full BConv 后只取目标塔 -> NTT
        立即 KeyMult
        立即 Accumulate
```

non-native 分支调用全宽 BConv **只为提供正确 residue**，其 CPU 内存与运行时间
不进入任何性能模型——本文件验证的是真 OC 的功能正确性，不是它的效率。

6 个测试全过：

| 测试 | 覆盖 |
|---|---|
| `MatchesDCReference` | 与 OpenFHE 的 DC 路径逐 residue 一致（49152 个） |
| `MatchesAcrossShapesAndLevels` | 8 组形状、L=4..7、alpha=1/2/3、D=2/3/4/6，共 **524288 个 residue**；含末位 digit 不完整与每 digit 单塔 |
| `AllSchedulesAgree` | DC / MP / 伪OC / 真OC 四者结果全等 |
| `OperatorCountsMatchSimulatorClosedForm` | bypass=L、逐目标塔 BConv=L(D−1)+KD、KeyMult=D(L+K)，在真实 OpenFHE 参数集上验证仿真器的闭式 |
| `ScalarSingleTargetBConvMatchesFullBConv` | 按 BConv 定义逐系数重算单列（32768 个 residue），证明「只算一列」= 全宽结果的那一列 |
| `EndToEndDecryptionIsCorrect` | 三种 production 策略的 EvalMult 解密结果正确 |

**运行方式**（注意不是 `pke_tests`）：

```bash
cd build && make hks_true_oc_tests -j8 && ./unittest/hks_true_oc_tests
```

`pke_tests` 目前在**静态初始化阶段**就会 abort：
`unittest/utckksrns/UnitTestInteractiveBootstrap.cpp` 的
`static ... = getTestData(__FILE__);` 在 main 之前读 CSV 并 `std::stoul`，
抛出的 `std::invalid_argument` 逃逸出静态初始化，连 `--gtest_list_tests` 都跑不起来
（既存问题，见 MAP.md）。该问题与 HKS 无关，本轮没有去修它，只是把 HKS 的功能
验证隔离到一个能跑的独立 target 里。等它修好，`hks_true_oc_tests` 这个 target
即可删除——同一个 `.cpp` 已经在 `pke_tests` 的 glob 里。
