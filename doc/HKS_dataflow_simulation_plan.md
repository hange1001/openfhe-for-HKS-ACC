# HKS 三种数据流同构硬件软件仿真实施方案

## 1. 目标

在不修改当前 RTL 的前提下，建立面向 Hybrid Key Switching（HKS）的软件仿真器，在**完全相同的硬件资源约束**下比较 Digit-Centric（DC）、Maximum-Parallel（MP）和 Output-Centric（OC）三种 dataflow。

仿真器需要回答：

1. 三种 dataflow 的 latency、peak memory、data movement 和资源利用率分别是多少；
2. 参数从当前小规模扩展到更大 CKKS/HKS workload 后，瓶颈位于计算、片上存储、DDR 还是 PCIe；
3. OC 获得收益需要哪些条件，例如 input residency、single-tower BConv、KeyMult/Accumulate 融合与 double buffering。

本文档只规定实施和验证方法。最终哪种 dataflow 更优，应由数据决定。

## 2. 仿真边界：方案 B

仿真覆盖：

```text
P1 INTT -> P2 BConv -> P3 NTT -> P4 KeyMult -> P5 Accumulate/Sum
```

ModDown 第一版作为三种 dataflow 共同的后处理，可计入总时间，但不参与第一阶段的调度差异分析。

当前已实现的 INTT、SCALE、BConv、NTT 使用 HLS/RTL 周期校准。尚未实现的 KeyMult/Accumulate 使用统一资源模型和延迟模型，三种 dataflow 必须共享同一组参数。

默认不增加额外 DSP：KeyMult 复用现有 4-lane modular multiplication 资源，并与 NTT/SCALE 互斥。后续可测试独立 KeyMult engine，但每个硬件配置点必须同时运行 DC、MP、OC。

## 3. 公平性规则

固定硬件配置 (H)，只改变调度 trace：

```text
T(strategy) = Simulate(Trace(strategy), HardwareConfig)
strategy in {DC, MP, OC}
```

以下条件在同一对比组内保持一致：

- clock period；
- Transform PE 数和 pipeline 参数；
- BConv array 的 row/column 数；
- KeyMult/Accumulate 的资源绑定和延迟；
- SRAM/URAM 容量、bank 和端口数；
- DDR/PCIe 带宽及固定启动开销；
- DMA 是否允许与计算重叠；
- kernel invocation boundary 和输入输出格式。

若 live data 超过片上容量，仿真器必须插入 `SpillStore/SpillLoad`，不能为某种策略隐式增加存储。每份结果输出 `hardware_config_hash`，证明三种策略使用同一配置。

实验分成两类：

1. **P4 reproduction**：严格使用当前 `N=4096, L=3, K=2, PE=4` 和现有 buffer/twiddle 限制，校准已实现算子。
2. **Industrial projection**：允许更大的 (N,L,K,D)，显式模拟 tiling、spill 和带宽限制。结果必须标为 projected，不能当作当前 bitstream 实测。

## 4. 可配置工业参数

### 4.1 Workload 参数

| 参数 | 含义 | 建议扫描值 |
|---|---|---|
| `ring_dimension / N` | polynomial ring dimension | 4096, 8192, 16384, 32768, 65536 |
| `q_towers / L` | 当前 level 的 Q tower 数 | 3, 6, 8, 12, 16, 24, 30 |
| `p_towers / K` | special-prime tower 数 | 1, 2, 3, 4, 5, 6, 8 |
| `num_part_q` | Q 的 digit partition 数 | 1, 2, 3, 4, 5, 6, 8 |
| `alpha` | 每个 digit 的目标 Q tower 数 | 默认由 `ceil(L/num_part_q)` 推导 |
| `limb_bits` | 单 modulus 位宽 | 50, 54, 59, 60, 64 |
| `batch_count` | 连续 HKS 数量 | 1, 2, 4, 8, 16 |

最后一个 digit 可能不完整，必须按以下公式推导：

```text
alpha_d = min(alpha, L - d*alpha)
beta_d = L + K - alpha_d
```

建议增加 OpenFHE parameter exporter，从实际 `CryptoContext` 导出：

```text
N, sizeQl, sizeP, numPartQ, modulus_bits[], level, digit_layout[]
```

手写大参数只表示 workload shape，不代表已经通过特定 CKKS 安全性或精度验证。正式 industrial case 优先使用 OpenFHE 真实上下文。

### 4.2 Hardware 参数

| 参数 | P4 默认值/来源 | 含义 |
|---|---|---|
| `clock_period_ns` | 7.0 ns | cycle 转换时间 |
| `transform_lanes` | 4 | INTT/NTT/SCALE 吞吐 |
| `bconv_rows` | 3 | BConv reduction tile |
| `bconv_cols` | 5 | BConv output tile |
| `keymult_binding` | `shared_transform_mul` | 默认不新增 DSP |
| `modmul_ii` | HLS/RTL 校准 | modular multiply 吞吐 |
| `sram_capacity_bytes` | HLS buffer map | live data 容量 |
| `sram_banks` | array partition/map | bank conflict |
| `ddr_bandwidth_GBps` | 待校准/扫描 | spill 和 evk 访问 |
| `pcie_fixed_us` | 待板卡实测 | invocation 固定开销 |
| `pcie_bandwidth_GBps` | 待板卡实测 | H2D/D2H |
| `dma_compute_overlap` | false | 是否通信计算重叠 |
| `allow_tiling` | projection 开启 | 处理大参数 |

硬件参数可以做敏感性扫描，但每个硬件配置点都必须完整运行三种 dataflow。

推荐分三步扫描：

1. **One-factor-at-a-time**：分别改变 (N,L,K,	ext{numPartQ})；
2. **Real-context sweep**：导入 3–5 组真实 OpenFHE `CryptoContext`；
3. **Capacity-boundary sweep**：在 peak memory 刚好超过 SRAM 附近密集取点。

## 5. 两层仿真架构

### 5.1 Functional reference

复用 OpenFHE modular arithmetic，根据 dataflow 改变执行顺序和对象生命周期。验证：

- DC、MP、OC 的最终 `cTilda0/cTilda1` 与 reference 逐 tower、逐 coefficient 相同；
- 正确处理不完整 digit；
- 正确处理 OC 的 native Q tower bypass；
- single-tower BConv 与 full BConv 对应 tower 一致。

### 5.2 Trace-driven performance simulator

不使用 CPU wall-clock 模拟 FPGA。生成硬件事件 DAG，再按资源、依赖和存储条件计算 makespan。

核心事件：

```text
Op {
    id
    kind                  # H2D/INTT/SCALE/BConv/NTT/KeyMult/Accumulate/Spill/D2H
    deps[]
    resource_requirements # engine, lane, memory port
    work                  # N, alpha, beta, tower count
    reads[] / writes[]
    allocs[] / frees[]
}
```

调度规则：

1. 依赖完成后事件才 ready；
2. engine 和 memory port 同时可用后才能启动；
3. allocation 超容量时插入 spill；
4. 仅在配置允许时重叠 DMA 和计算；
5. 记录每个事件的 start/end cycle 和 stall reason。

核心实现建议使用 Python 标准库 `dataclasses + heapq`。CSV/JSON 为必选输出，matplotlib/pandas 仅作可选报告依赖。

默认资源：

```text
transform_engine      capacity = 1
transform_mul_lanes   capacity = 4
bconv_array           capacity = 1, shape = 3 x 5
keymult_engine        binding  = transform_mul_lanes
accumulator_ports     fixed read/write ports
dma_h2d               capacity = 1
dma_d2h               capacity = 1
external_memory       bandwidth + access latency
on_chip_memory        capacity + banks + ports
```

KeyMult 初始模型：

```text
T_KM(d,t) =
  2*ceil(N/PE_mul)*II_mul
  + T_pipe + T_acc + T_mem
```

因子 2 对应 `c_j*b_j` 和 `c_j*a_j`。必须同时运行 optimistic、nominal、pessimistic 三组参数，避免对未实现模块给出伪精确单点结果。

## 6. 三种 trace

### 6.1 DC

```text
for each digit d:
    INTT/SCALE native towers
    BConv all complement towers
    NTT required towers
    KeyMult all output towers
    Accumulate into cTilda0/cTilda1
    release digit temporaries
```

DC 可以释放当前 digit 的 ModUp 临时数据，但要跨 digit 保存 output accumulators；放不下时必须模拟写回和重读。

### 6.2 MP

```text
all digits INTT/SCALE
barrier
all digits BConv
barrier
all digits NTT
barrier
all digits KeyMult
Accumulate
```

MP 和 DC 的基本算术量应一致，差异主要来自 barrier、并行机会和 live working set。单 transform engine、单 BConv array 下不得通过虚假并行缩短 MP 时间。

### 6.3 OC

```text
retain/reload digit inputs under the common SRAM limit

for each target output tower:
    allocate one cTilda0/cTilda1 tower accumulator
    for each required digit:
        bypass native Q tower when valid
        otherwise BConv only the target tower
        NTT if required
        KeyMult immediately
        Accumulate immediately
    write completed output tower
    release accumulator
```

OC 的 single-output BConv 数量：

```text
N_BConv_OC = L*(D - 1) + K*D
```

当前 (L=3,K=2,D=2) 时为 7。若 digit input 不能驻留，每个 target tower 的重新读取成本必须计入。

## 7. 存储模型

至少跟踪：

- digit native/coefficient-domain towers；
- BConv/NTT temporaries；
- full `partsCtExt`；
- evaluation-key A/B towers；
- `cTilda0/cTilda1` accumulators；
- twiddle/constants/metadata；
- DMA input/output buffers。

64-bit limb 下：

```text
B_tower = 8*N bytes
```

```text
B_OC_acc = 2*8*N bytes
```

```text
B_full_acc = 2*(L + K)*8*N bytes
```

Peak memory 必须由 trace 的真实 lifetime 计算，不能使用预填固定结果。

## 8. 建议代码结构

```text
tools/hks_flow_sim/
  README.md
  hks_sim.py
  hks_sim/
    config.py
    workload.py
    op.py
    resources.py
    cost_model.py
    memory.py
    scheduler.py
    trace_dc.py
    trace_mp.py
    trace_oc.py
    report.py
  configs/
    p4_reproduction.yaml
    industrial_projection.yaml
    keymult_optimistic.yaml
    keymult_nominal.yaml
    keymult_pessimistic.yaml
  tests/
    test_workload.py
    test_trace_counts.py
    test_fairness.py
    test_memory.py
    test_scheduler.py
```

CLI：

```text
python tools/hks_flow_sim/hks_sim.py \
  --workload configs/workload.yaml \
  --hardware configs/hardware.yaml \
  --strategies dc,mp,oc \
  --output results/run_name
```

参数扫描示例：

```text
--sweep N=4096,8192,16384,32768,65536
--sweep L=3,6,8,12,16,24,30
--sweep sram_mb=1,2,4,8,16,32
--strict-p4
--enable-tiling
--emit-trace
--emit-gantt
```

## 9. 输出指标

每次运行至少输出：

```text
strategy
hardware_config_hash
workload_config_hash
total_cycles / latency_us
compute_cycles
memory_stall_cycles
dependency_stall_cycles
h2d_bytes / d2h_bytes
ddr_read_bytes / ddr_write_bytes
spill_read_bytes / spill_write_bytes
peak_live_bytes
transform_utilization
bconv_utilization
keymult_utilization
```

事件 trace 输出：

```text
event_id, strategy, kind, digit, tower, resource,
start_cycle, end_cycle, bytes_read, bytes_written, stall_reason
```

最终报告生成 latency、peak memory、traffic、utilization、SRAM-capacity curve、参数热力图和典型 workload 的 Gantt trace。

## 10. 实施里程碑

### M0：配置和事件规格

- 定义 workload/hardware schema；
- 定义 `Op`、`Buffer`、`Resource`；
- 固化 P4 reproduction；
- 增加 config hash 和公平性检查。

完成标准：正确推导 (alpha_d, beta_d, D)，拒绝非法参数。

### M1：静态 trace 和存储统计

- 生成 DC/MP/OC trace；
- 统计 INTT/BConv/NTT/KeyMult 工作量；
- 跟踪 buffer lifetime 和 peak memory；
- 检查当前参数下 OC BConv 数为 7。

完成标准：不依赖周期参数即可生成 operation-count/memory 报告。

### M2：资源约束调度

- 实现 ready queue；
- 实现 engine 互斥和端口冲突；
- 实现 DMA overlap 开关；
- 输出 makespan、utilization、stall breakdown。

完成标准：单 engine 配置下，MP 不会因虚假并行获得额外吞吐。

### M3：P4 校准和 KeyMult 区间

- 填入已有算子的 RTL co-sim 周期；
- 用 fused digit kernel 总周期校验子模型；
- 区分 compute 和 PCIe/XRT 固定开销；
- 建立三组 KeyMult 参数。

完成标准：P1–P3 路径预测与可信 RTL 数据的目标偏差不超过 5%；数据不足时输出区间。

### M4：Tiling 和 spill

- fixed-capacity allocator；
- BConv row/output tiling；
- digit input residency/reload；
- accumulator spill；
- capacity-boundary curve。

完成标准：超容量 workload 仍可运行，所有额外搬运可追溯。

### M5：Functional reference 和参数导入

- 导出真实 `CryptoContext`；
- 实现/连接三种 functional reference；
- 逐 tower、逐 coefficient 对比；
- 覆盖不完整 digit 和多个 level。

完成标准：三种调度与 OpenFHE reference 数学结果一致。

### M6：Industrial projection

- 执行 (N,L,K,D) 扫描；
- 执行 SRAM、bandwidth、KeyMult 敏感性扫描；
- 导入真实 workload；
- 输出 CSV/JSON、图表和 Gantt；
- 标记 measured、calibrated、projected。

完成标准：得到带条件的结论，而不是脱离参数的简单排名。

## 11. 下午优先收集的数据

P0，已实现算子：

- 单 limb INTT/NTT cycle；
- SCALE cycle；
- `alpha=2,beta=3` 和 `alpha=1,beta=4` 的 BConv cycle；
- single-output BConv cycle 或 pipeline 参数；
- buffer 的 BRAM/URAM mapping、bank 和 port；
- fused `OP_HKS_DIGIT` 的阶段周期；若只有总周期，先用于整体校准。

P1，主机和存储：

- kernel launch 固定开销；
- 1/2/4/8/16 limb 的 H2D/D2H 时间；
- DDR 有效带宽；
- XRT 是否完全同步；
- DMA/compute overlap 上限。

P2，未实现 KeyMult：

- modular multiplier latency 和 II；
- `a/b` 两路乘法顺序或并行关系；
- accumulate 能否与 write-back 融合；
- evk 从片上、DDR 或 host 提供时的访问成本。

## 12. 验证、风险与候选判断

必须建立：

- **Trace count test**：检查算子数量和 bypass；
- **Fairness test**：三种策略 hardware hash 相同；
- **Dependency test**：KeyMult 不得早于对应 tower；
- **Capacity test**：容量不足时 spill 或 infeasible；
- **Monotonicity test**：降低带宽/SRAM 不应无故变快；
- **Calibration hold-out**：预留数据点验证误差；
- **Functional test**：三种调度与 OpenFHE reference 一致。

风险控制：

1. KeyMult 未实现时报告敏感性区间，不报伪精确单点；
2. 所有数据标记 `measured/calibrated/projected`；
3. 不只统计算术量，必须跟踪 lifetime 和 bytes；
4. 正确处理最后一个不完整 digit；
5. 先做单变量和容量边界扫描，再做组合扫描。

待验证的候选判断：

- 单 transform engine 和单 BConv array 下，MP 可能不减少 compute cycles，反而增加 peak memory；
- DC 在当前小参数下可能具有最好或接近最好的 latency；
- OC 更稳定的优势可能是 peak memory 和 external traffic，不一定是当前参数下的 compute latency；
- 当 (L,K,D) 增大并使 DC/MP 超过片上容量时，OC 优势可能扩大；
- 若 OC 的 digit input 无法驻留，或每个 target tower 都产生新的 host/kernel 调用，其收益可能被重读和固定开销抵消。

最终报告应给出条件化结论，例如“当 SRAM 小于 X MB 且 DDR 带宽低于 Y GB/s 时 OC 更优”，而不是宣称某种 dataflow 始终最优。

