# CG-NTT 性能优化技术文档

## 1. 背景与动机

### 1.1 设计目标

本项目为 HKS 全同态加密（FHE）硬件加速器实现高性能 NTT（Number Theoretic Transform）模块。目标器件为 Xilinx Alveo U55C（xcu55c-fsvh2892-2L-e），时钟频率 200 MHz（5 ns 周期）。

CG-NTT（Constant Geometry NTT）采用恒定几何蝶形拓扑，消除了传统 NTT 中的 MUX 交叉开关，使 PE 阵列的物理连线固定，适合 FPGA 流水线实现。

### 1.2 基线性能分析

综合报告（Vitis HLS 2023.2）显示，`CG_NTT_Kernel` 核心计算仅需 10,499 个周期，但单次 Limb 迭代被拉长到 43,316 个周期，8 个 Limb 总延迟高达 346,529 个周期（1.733 ms @ 200 MHz）。

**基线延迟分解（单 Limb = 43,316 cycles）：**

| 阶段 | 延迟 (cycles) | 占比 |
|------|--------------|------|
| LOAD_IN（加载多项式数据） | 4,097 | 9.5% |
| LOAD_TW（加载旋转因子） | 24,577 | 56.7% |
| CG_NTT_Kernel（蝶形计算） | 10,499 | 24.2% |
| STORE_OUT（写回结果） | 4,097 | 9.5% |

I/O 搬运占总延迟的 75.8%，核心计算仅占 24.2%。这是典型的"计算快、搬运慢"的 HLS 系统级瓶颈。

### 1.3 瓶颈根因

综合报告中 m_axi 端口全部退化为 64-bit 位宽：

```
m_axi_gmem0: 64 -> 64 bits
m_axi_gmem1: 64 -> 64 bits
m_axi_gmem2: 64 -> 64 bits
```

HLS 报出多个 Widen Fail 警告：

> Could not widen since type i64 size is greater than or equal to alignment 1(bytes)

原因：顶层接口使用多维数组指针（`uint64_t in_data[MAX_LIMBS][RING_DIM]`），HLS 无法证明多维索引计算后的地址对齐性，拒绝自动拓宽。U55C 的 HBM/DDR 物理 AXI 端口支持 512-bit，但实际只用了 64-bit，带宽利用率仅 12.5%。

---

## 2. 优化方案总览

针对上述瓶颈，提出三级递进优化方案：

| 优化级别 | 方案 | 目标 | 预期加速 |
|---------|------|------|---------|
| L1 | 512-bit AXI 总线拓宽 | 消除 I/O 带宽瓶颈 | 单 Limb 3.0x |
| L2 | DATAFLOW 任务级流水 | Limb 间 Load/Compute/Store 重叠 | 吞吐量质变 |
| L3 | 旋转因子 URAM 预加载 + BUTTERFLY II=1 修复 | 消除循环内 TW 加载 + 计算提速 | 进一步 2-3x |

---

## 3. L1 优化：512-bit AXI 总线拓宽（已实现）

### 3.1 实现方法

将 `Compute_CG_NTT` 顶层接口的三个大数组参数从 `uint64_t` 多维数组改为 `ap_uint<512>*` 一维指针：

```cpp
// 改前
void Compute_CG_NTT(
    uint64_t in_data[MAX_LIMBS][RING_DIM],
    const uint64_t cg_ntt_twiddle[MAX_LIMBS][STAGE][CG_HALF_N],
    const uint64_t cg_intt_twiddle[MAX_LIMBS][STAGE][CG_HALF_N],
    ...
);

// 改后
void Compute_CG_NTT(
    ap_uint<512> *in_data,
    const ap_uint<512> *cg_ntt_twiddle,
    const ap_uint<512> *cg_intt_twiddle,
    ...
);
```

片上局部缓冲保持 `uint64_t`，在 Load/Store 阶段手动 pack/unpack：

```cpp
// 512-bit burst 读入，每拍解包 8 个 uint64_t
LOAD_IN:
for (int i = 0; i < PACKED_RING_DIM; i++) {  // 512 次迭代（原 4096 次）
    #pragma HLS PIPELINE II=1
    ap_uint<512> pack = in_data[l * PACKED_RING_DIM + i];
    for (int j = 0; j < 8; j++) {
        #pragma HLS UNROLL
        local_in_data[i * 8 + j] = pack.range((j+1)*64-1, j*64);
    }
}
```

小参数表（`modulus[8]`、`K_HALF[8]`、`M_barrett[8]`）保持 `uint64_t*`，仅 8 个元素，拓宽无意义。

### 3.2 综合结果对比

**M_AXI 端口位宽变化：**

| 端口 | 改前 (SW→HW) | 改后 (SW→HW) | 用途 |
|------|-------------|-------------|------|
| gmem0 | 64 → 64 | 512 → 512 | 多项式数据 (in_data) |
| gmem1 | 64 → 64 | 512 → 512 | NTT 旋转因子 |
| gmem2 | 64 → 64 | 512 → 512 | INTT 旋转因子 |
| gmem3 | 64 → 64 | 64 → 64 | 模数参数（不变） |

**Burst 推断结果：**

| 端口 | 方向 | 改前 (长度×宽度) | 改后 (长度×宽度) | 数据量 |
|------|------|-----------------|-----------------|--------|
| gmem0 | read | 4,096 × 64-bit | 512 × 512-bit | 32 KB |
| gmem0 | write | 4,096 × 64-bit | 512 × 512-bit | 32 KB |
| gmem1 | read | 24,576 × 64-bit | 3,072 × 512-bit | 192 KB |
| gmem2 | read | 24,576 × 64-bit | 3,072 × 512-bit | 192 KB |

**延迟对比：**

| 阶段 | 改前 (cycles) | 改后 (cycles) | 加速比 |
|------|--------------|--------------|--------|
| LOAD_IN | 4,097 | 513 | 8.0x |
| LOAD_TW (NTT/INTT) | 24,577 | 3,073 | 8.0x |
| CG_NTT_Kernel | 10,499 | 10,499 | 1.0x |
| STORE_OUT | 4,097 | 513 | 8.0x |
| **单 Limb 迭代** | **43,316** | **14,644** | **3.0x** |
| **8 Limbs 总计** | **346,529** | **117,153** | **3.0x** |

**绝对时间：** 1.733 ms → 0.586 ms @ 200 MHz

**资源开销变化：**

| 资源 | 改前 | 改后 | 变化 | U55C 占比 |
|------|------|------|------|----------|
| BRAM | 96 (2%) | 174 (4%) | +78 | 4% |
| DSP | 348 (3%) | 348 (3%) | 0 | 3% |
| FF | 30,827 (1%) | 40,713 (1%) | +9,886 | 1% |
| LUT | 51,412 (3%) | 56,498 (4%) | +5,086 | 4% |
| URAM | 8 (~0%) | 8 (~0%) | 0 | ~0% |

BRAM 增长主要来自 512-bit m_axi FIFO 缓冲（每个 gmem 端口从 4 BRAM 增至 30 BRAM）。总资源占比仍然很低。

---

## 4. L2 优化：DATAFLOW 任务级流水（待实现）

### 4.1 问题分析

当前 `LIMB_LOOP` 内部严格串行执行：

```
Limb 0: [LOAD_IN] → [LOAD_TW] → [CG_NTT_Kernel] → [STORE_OUT]
Limb 1:                                               [LOAD_IN] → [LOAD_TW] → ...
```

前一个 Limb 完全处理完毕后，下一个 Limb 才开始加载数据。Limb 间没有任何重叠。

### 4.2 优化方案

使用 `#pragma HLS DATAFLOW` 将 Load、Compute、Store 拆分为独立子函数，使 HLS 自动插入 Ping-Pong（PIPO）缓冲，实现 Limb 间流水重叠：

```
Limb 0: [LOAD] → [COMPUTE] → [STORE]
Limb 1:          [LOAD]     → [COMPUTE] → [STORE]
Limb 2:                       [LOAD]     → [COMPUTE] → [STORE]
```

```cpp
void load_task(const ap_uint<512> *in_data, const ap_uint<512> *twiddle,
               uint64_t local_data[RING_DIM], uint64_t local_tw[STAGE][CG_HALF_N],
               int l, bool is_ntt);
void compute_task(uint64_t local_data[RING_DIM], const uint64_t local_tw[STAGE][CG_HALF_N],
                  uint64_t mod, uint64_t k, uint64_t m, bool is_ntt);
void store_task(ap_uint<512> *in_data, const uint64_t local_data[RING_DIM], int l);

LIMB_LOOP:
for (int l = mod_idx_offset; l < mod_idx_offset + num_active_limbs; l++) {
    #pragma HLS DATAFLOW
    load_task(...);
    compute_task(...);
    store_task(...);
}
```

### 4.3 注意事项

- DATAFLOW 要求子函数间通过局部数组传递数据，HLS 会将其综合为 PIPO 缓冲
- PIPO 会使 `local_in_data` 和 `local_twiddle` 的 BRAM/URAM 用量翻倍（当前 16 BRAM + 8 URAM → 32 BRAM + 16 URAM），在 U55C 上完全可接受
- `in_data` 在 gmem0 上既读又写，但 AXI 读写通道独立，Limb N 的 STORE 和 Limb N+1 的 LOAD 可以在同一端口上重叠
- `is_ntt` 条件分支选择不同的 twiddle 指针（gmem1 vs gmem2），需要在 DATAFLOW 区域外预先选好指针，避免条件分支破坏 DATAFLOW 推断

### 4.4 预期收益

DATAFLOW 后，稳态吞吐量由最慢的子任务决定。L1 优化后各阶段延迟：

| 子任务 | 延迟 (cycles) |
|--------|--------------|
| load_task (LOAD_IN + LOAD_TW) | ~3,586 |
| compute_task (CG_NTT_Kernel) | 10,499 |
| store_task (STORE_OUT) | 513 |

瓶颈为 compute_task（10,499 cycles）。稳态下每个 Limb 的吞吐间隔 ≈ 10,499 cycles，8 Limbs 总延迟 ≈ 10,499 × 8 + 流水线填充/排空 ≈ **87,000 cycles**（较 L1 的 117,153 再降 ~26%）。

---

## 5. L3 优化：旋转因子 URAM 预加载 + BUTTERFLY II 修复（待实现）

### 5.1 旋转因子 URAM 全量预加载

**问题：** 当前每个 Limb 迭代都从 DDR 加载完整的旋转因子表（3,073 cycles），不同 Limb 使用不同模数对应不同的旋转因子，无法简单提到循环外。

**方案：** 在 `LIMB_LOOP` 之前，一次性将所有 Limb 的旋转因子全部预加载到片上 URAM：

```cpp
// 预加载所有 Limb 的旋转因子到 URAM
uint64_t all_twiddle[MAX_LIMBS][STAGE][CG_HALF_N];
#pragma HLS BIND_STORAGE variable=all_twiddle type=ram_1p impl=uram

PRELOAD_ALL_TW:
for (int i = 0; i < MAX_LIMBS * PACKED_TW_SIZE; i++) {
    #pragma HLS PIPELINE II=1
    ap_uint<512> pack = twiddle_ptr[i];
    // ... unpack to all_twiddle
}
```

**资源开销：** 8 limbs × 12 stages × 2048 × 8 bytes = 1.5 MB。U55C 有 320 个 URAM（每个 288 Kbit），总容量约 11.5 MB。1.5 MB 仅占 ~13%。

**收益：** LIMB_LOOP 内部彻底消除 LOAD_TW 阶段。结合 DATAFLOW，load_task 只剩 LOAD_IN（513 cycles），compute_task（10,499 cycles）仍为瓶颈，但流水线填充更快。

### 5.2 BUTTERFLY_LOOP II=3 → II=1 修复

**问题：** 综合报告显示 BUTTERFLY_LOOP 的 Initiation Interval 为 3，未达到目标 II=1：

```
BUTTERFLY_LOOP | II | 785 cycles | iteration latency=21 | II=3 | trip count=256
```

256 次迭代 × II=3 = 768 + 流水线深度 = 785 cycles/stage。若能达到 II=1，每 stage 仅需 ~277 cycles，12 个 stage 的 STAGE_LOOP 从 9,468 降至 ~3,324 cycles，CG_NTT_Kernel 从 10,499 降至 ~4,350 cycles。

**根因分析：** 8 个 PE 每拍需要 2 读 + 2 写 = 32 次 BRAM 访问。`buf_A`/`buf_B` 使用 `cyclic factor=16` + `ram_2p` 理论上提供 32 端口，但 `is_ntt` 条件分支导致 NTT/INTT 两种读写模式的端口映射不同，HLS 调度器无法在单拍内同时满足两种模式。

**方案：** 将 NTT 和 INTT 的蝶形循环拆分为两个独立的 BUTTERFLY_LOOP，消除条件分支，让 HLS 对每个分支独立做端口调度：

```cpp
STAGE_LOOP:
for (int stage = 0; stage < STAGE; stage++) {
    int actual_stage = is_ntt ? stage : (STAGE - 1 - stage);
    if (is_ntt) {
        NTT_BUTTERFLY:
        for (int i = 0; i < CG_HALF_N / CG_PE_NUM; i++) {
            #pragma HLS PIPELINE II=1
            // 仅 NTT 读写模式，无条件分支
        }
    } else {
        INTT_BUTTERFLY:
        for (int i = 0; i < CG_HALF_N / CG_PE_NUM; i++) {
            #pragma HLS PIPELINE II=1
            // 仅 INTT 读写模式，无条件分支
        }
    }
}
```

### 5.3 L3 综合预期

| 场景 | 单 Limb (cycles) | 8 Limbs (cycles) | 时间 @ 200MHz |
|------|-----------------|------------------|---------------|
| 基线 | 43,316 | 346,529 | 1.733 ms |
| +L1 (512-bit) | 14,644 | 117,153 | 0.586 ms |
| +L2 (DATAFLOW) | 吞吐 ~10,500 | ~87,000 | 0.435 ms |
| +L3 (URAM + II=1) | 吞吐 ~5,400 | ~46,000 | 0.230 ms |

从基线到 L3 全部落地：**346,529 → ~46,000 cycles，约 7.5x 系统级加速**。

---

## 6. 总结

CG-NTT 核心计算（PE 阵列 + 乒乓缓冲 + 恒定几何寻址）的设计已经非常高效，瓶颈完全在系统级 I/O 和任务调度上。三级优化方案按风险和收益递进：

- **L1（已完成）：** 512-bit 总线拓宽，I/O 延迟 8x 缩减，总延迟 3.0x 加速，零功能风险
- **L2（待实现）：** DATAFLOW 流水，Limb 间重叠执行，吞吐量由计算瓶颈决定
- **L3（待实现）：** URAM 预加载消除循环内 TW 搬运 + 蝶形 II 修复，计算本身再提速 2.4x

三级优化互相正交，可独立实施和验证。
