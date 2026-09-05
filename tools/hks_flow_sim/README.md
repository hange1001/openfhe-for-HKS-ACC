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

# 完整 HKS 边界，KeyMult 三档必须一起跑
python3 hks_sim.py --config configs/p4_reproduction.yaml \
                   --config configs/keymult_nominal.yaml \
                   --boundary full_hks --invocation per_hks

# 扫参（每个点都完整跑三种策略）
python3 hks_sim.py --config configs/industrial_projection.yaml \
                   --sweep workload.ring_dimension=4096,8192,16384 \
                   --output results/n_sweep --emit-trace

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

两参数解两方程，没有剩余自由度，所以必须留 hold-out。

### P3 hold-out

P3 的预乘是**独立硬件、1 系数/拍**（P4 才改成复用 4 路 lane），结构与 P4 不同。
把 P4 标定出的 `axi`/`control` 原样搬过去，只解预乘一项：两个事务分别给出
**4101** 和 **4119**，互差 0.4%。取中值反代：

| | 预测 | 实测 | 误差 |
|---|---:|---:|---:|
| P3 alpha=2 | 52797 | 52779 | +0.034% |
| P3 alpha=1 | 47634 | 47643 | −0.019% |
| P3 两 digit | 100431 | 100422 | **+0.009%** |

M3 门槛是 5%。顺带独立印证了一件已知的事：解出的 4110 ÷ P4 的 1047 = 3.93 ≈ lane 数 4。

断言写在 `tests/test_cost_model.py`。

---

## 相对方案的三处改动

方案 §5.2 / §4.2 有三个地方会静默改变结论，实现时补上了：

### 1. tie-break 必须固定并计入 hash

方案只规定「依赖完成且资源可用即可启动」，没规定多个 ready 事件争同一引擎时选谁。
这个自由度会**系统性偏袒某种 dataflow**：LIFO 偏袒 DC 的深度优先，FIFO 偏袒 MP
的广度优先。默认 `tie_break="trace_order"`（严格按 trace 发射序，是唯一不引入额外
自由度的选择），并计入 `hardware_config_hash`。

### 2. `allow_engine_overlap` —— 标定点分辨不出、却能翻转结论的一条

P4 的 Top 是 opcode-driven RPC dispatcher（ADR-008），BConv 与变换引擎由顶层
**顺序调用**，没有 DATAFLOW。但方案 §5.2 把它们列为两个独立容量池，隐含了可并发。

DC 的单次 BConv 必然先于它的全部 NTT，永远不会有两个引擎同时 ready，所以
**P4 的 91242 在开关的两种取值下都能对上**；而 OC 的 single-tower 粒度会让
BConv/NTT 重叠。结果：

| | DC | MP | OC | |
|---|---:|---:|---:|---|
| `allow_engine_overlap=False`（默认，P4 忠实） | 91242 | 91242 | 109861 | OC 慢 20.4% |
| `allow_engine_overlap=True`（假设已做 DATAFLOW） | 91242 | 91242 | 84991 | OC 快 6.9% |

默认取 `False`。置 `True` 等于假设已完成 phase4 step 4.1 的 DATAFLOW 拆分，
那是**另一个硬件配置点**，hash 会变。

### 3. 调用边界是串行的

XRT 路径完全阻塞（`bo.sync → run.wait → bo.sync`，无 async、无双缓冲，
见 task.yaml q8），下一次 kernel 调用不能在上一次返回前开始。不加这道 barrier，
调度器会让下一个 digit 的 H2D 与上一个 digit 的 D2H 重叠，两个 digit 会算成
89533 而不是实测的 91242。

`invocation` 粒度（`per_hks` / `per_digit` / `per_output_tower`）是可配置的：
方案 §12 点名的「OC 每个 target tower 一次调用会被固定开销吃掉收益」，
在这里是一个**可证伪的假设**，不是写死的结论。

---

## 一条必须先看清的恒等式

```
sum_d beta_d = D(L+K) - L = L(D-1) + KD
```

左边是 DC 逐 digit 的 BConv 输出塔总数，右边是方案 §6.3 给 OC 的闭式——**它们相等**。
INTT 也都是 L 次，NTT 也都是 7 次（P4 参数下）。

**三种 dataflow 的算术量完全相同**，差异 100% 来自 lifetime / 存储 / 搬运 / 调度。
这与推导v1 §2.3「每算子一次 offload 下三策略 AI 恒等」是同一件事。

所以：
- M1 的 operation-count 报告**没有区分度**，真正的区分在 peak memory 和 stall breakdown；
- `tests/test_trace_counts.py` 断言三策略的 INTT/SCALE/NTT/KeyMult/Accumulate
  计数**逐项相等**，BConv 输出塔总数也相等。允许不同的只有 BConv 的**调用次数**
  （OC 是 single-tower，调用更多但每次产出更少）。任何一方算子数不同，就是 trace 写错了。

### 为什么 OC 的 single-tower BConv 不省时间

成本模型是 `T = N·⌈alpha/rows⌉·⌈beta/cols⌉ + 49`。固定 3×5 阵列**并行产出全部 5 列**，
所以 beta=1 和 beta=5 同价。OC 把 BConv 调用次数从 2 涨到 7，每次仍是满价 4145，
净多付 20725 拍。

这正是 ADR-010 / `OC_strategy_gap_analysis` 说的「sizeP 倍冗余」的周期形式。
**要让 OC 真正省下 BConv 时间，必须缩窄阵列或做列分块**（`allow_tiling`）。

---

## 当前结果（P4 参数，N=4096/L=3/K=2/D=2）

`boundary=full_hks`，`invocation=per_hks`，引擎顺序调用：

| KeyMult 档 | DC | MP | OC |
|---|---:|---:|---:|
| optimistic | **134132** | 142556 | 154857 |
| nominal | **137142** | 145566 | 157867 |
| pessimistic | **165682** | 174106 | 186407 |
| peak memory | 576 KB | 1600 KB | **384 KB** |

**排序在三档下不变** —— 结论不依赖未实现模块的参数，这是个稳健性结果。

对照方案 §12 的待验证判断：

- ✅ 单 engine + 单 BConv array 下，MP 不减少 compute cycles，反而把 peak 抬到 2.8×
  （barrier 拉长生命周期）；DC 与 MP 除 peak 外每项指标相同。
- ✅ DC 在当前小参数下 latency 最好。
- ✅ OC 的优势是 peak memory（0.67× DC），不是 compute latency。
- ⏳ 「(L,K,D) 增大使 DC/MP 超过片上容量时 OC 优势扩大」—— 需要 M4 的容量边界扫描。
- ⏳ 「OC digit input 无法驻留时收益被重读抵消」—— 收紧 `sram_capacity_bytes`
  即可触发 SpillLoad，已实现但尚未系统扫描。

⚠️ 全部为 RTL 周期下界换算，**不含 PCIe / 驱动 / buffer 分配 / Host 调度**，无板卡。

---

## PCIe 参数怎么处理

方案 §4.2 把 `pcie_fixed_us` / `pcie_bandwidth_GBps` 标为「待板卡实测」。但项目没有
板卡，而且 task.yaml q8 已经确认现有的 182.8 μs/call **是下界不是墙钟**（漏计
`xrt::bo` 分配、`bo.write/read` 的 host memcpy 与析构）。

所以这两个参数在这里是**扫描维，不是标定量**，默认 0（不建模）。要回答方案 §1 的
问题 2（瓶颈在计算、片上存储、DDR 还是 PCIe），用：

```bash
--sweep hardware.pcie_fixed_us=0,50,182.8,400
```

输出敏感性曲线，全部标 `projected`。

---

## 目录结构

```
tools/hks_flow_sim/
  hks_sim.py            CLI
  hks_sim/
    config.py           schema / 校验 / config hash
    workload.py         digit 布局推导（alpha_d, beta_d, D）
    op.py               Op / Buffer / TraceBuilder
    resources.py        资源池与绑定
    cost_model.py       周期模型（measured / calibrated / projected）
    memory.py           lifetime / peak / spill 容量 pre-pass
    scheduler.py        事件驱动调度器 + makespan 四分解
    trace_common.py     三种 trace 共用的发射器
    trace_dc.py         Digit-Centric
    trace_mp.py         Maximum-Parallel
    trace_oc.py         Output-Centric
    engine.py           编排（方案 §8 没有这一层，从 report.py 分出来的）
    report.py           CSV / JSON / 事件 trace / 文本摘要
  configs/              p4_reproduction / industrial_projection / keymult 三档
  tests/                59 个断言
```

`engine.py` 是相对方案 §8 唯一新增的文件：把编排从 `report.py` 分出来，
让 `report.py` 只负责输出格式，测试可以直接调 `run_all` 而不经过 CLI。

---

## 里程碑状态

| | 状态 | 说明 |
|---|---|---|
| M0 配置与事件规格 | ✅ | (alpha_d, beta_d, D) 推导 + 非法参数拒绝 + config hash |
| M1 静态 trace 与存储统计 | ✅ | 不依赖周期参数即可出 operation-count / memory 报告 |
| M2 资源约束调度 | ✅ | ready queue / 引擎互斥 / DMA overlap 开关 / makespan 四分解 |
| M3 P4 校准与 KeyMult 区间 | ✅ | P4 零残差、P3 hold-out 0.009%、三档区间 |
| M4 Tiling 与 spill | 🟡 | 容量 pre-pass 与 SpillStore/Load 已实现；容量边界曲线未系统扫描；BConv 列分块未做 |
| M5 Functional reference | ⬜ | **待定：路线需拍板，见下** |
| M6 Industrial projection | ⬜ | 扫参框架就绪，尚未系统跑 |

### M5 需要拍板的事

项目里已经有 DC/MP/OC 的 C++ 实现（`hks_strategy.h` + `keyswitch-hybrid.cpp`），
但按 ADR-010 / `OC_strategy_gap_analysis`，**现有 OC 是「伪 OC」**（DC + sizeP 冗余），
跟本仿真器 §6.3 建模的真 OC（single-tower BConv + 立即 KeyMult/Accumulate）不是一回事。

两条路：

- **(a)** 复用现有 C++ 三策略做 reference —— 便宜，但 OC 那一路对不上要验的东西；
- **(b)** 另写一份能验真 OC 的 functional model —— 贵，但 M5 的完成标准才真正成立。

在拍板之前 M5 不动。
