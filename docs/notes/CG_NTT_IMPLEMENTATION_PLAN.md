# CG-NTT（恒定几何 NTT）实现计划

## Context

当前 FPGA NTT 使用 2D [64×64] 布局 + 复杂的位运算索引生成（`generate_input_index` / `generate_output_index`），导致大量 MUX 和路由资源消耗。CG-NTT 通过"每层固定读写模式 + Perfect Shuffle 重排"，将地址生成简化为 `i`, `i+N/2`, `2i`, `2i+1`，彻底消灭交叉开关，释放 LUT 资源。

本计划创建 3 个新文件，**不修改任何现有文件**。

---

## 文件清单

| 文件 | 用途 |
|------|------|
| `src/fpga_backend/include/cg_ntt.h` | CG-NTT 函数声明 |
| `src/fpga_backend/src/cg_ntt.cpp` | CG-NTT 核心实现 |
| `src/fpga_backend/testbench/cg_ntt_tb.cpp` | 独立测试台（NTT→INTT 往返验证） |

---

## 1. cg_ntt.h — 函数声明

```cpp
// 核心：单 limb CG-NTT
void CG_NTT_Kernel(
    uint64_t in_data[RING_DIM],          // 输入/输出（1D 平铺）
    const uint64_t modulus,
    const uint64_t K_HALF,
    const uint64_t M_barrett,
    const uint64_t cg_twiddle[STAGE][RING_DIM/2],  // 预计算的 CG 旋转因子
    bool is_ntt
);

// 多 limb 包装器
void Compute_CG_NTT(
    uint64_t in_data[MAX_LIMBS][RING_DIM],
    const uint64_t cg_ntt_twiddle[MAX_LIMBS][STAGE][RING_DIM/2],
    const uint64_t cg_intt_twiddle[MAX_LIMBS][STAGE][RING_DIM/2],
    const uint64_t modulus[MAX_LIMBS],
    const uint64_t K_HALF[MAX_LIMBS],
    const uint64_t M_barrett[MAX_LIMBS],
    bool is_ntt,
    int num_active_limbs,
    int mod_idx_offset
);

// 2D ↔ 1D 转换辅助（与现有 [SQRT][SQRT] 布局互转）
void flatten_2d_to_1d(const uint64_t src[SQRT][SQRT], uint64_t dst[RING_DIM]);
void reshape_1d_to_2d(const uint64_t src[RING_DIM], uint64_t dst[SQRT][SQRT]);

// CG-NTT 输出重排（bit-reversal + unshuffle）
void cg_ntt_reorder(uint64_t data[RING_DIM]);
```

---

## 2. cg_ntt.cpp — 核心实现

### 2.1 存储架构

```cpp
static const int HALF_N = RING_DIM / 2;  // 2048
static const int CG_PE_NUM = 8;          // 并行 PE 数

uint64_t buf_A[RING_DIM];
uint64_t buf_B[RING_DIM];

// 16 bank cyclic 分区：8 PE × 读 2 个数 = 16 端口无冲突
#pragma HLS ARRAY_PARTITION variable=buf_A cyclic factor=16 dim=1
#pragma HLS ARRAY_PARTITION variable=buf_B cyclic factor=16 dim=1
```

**为什么 factor=16？**
- 写端：`buf[2*global_i]` 和 `buf[2*global_i+1]` 是相邻地址，落在相邻 bank
- 读端：`buf[global_i]` 和 `buf[global_i + 2048]`，当 factor=16 时，两者的 `%16` 值相差 `2048%16=0`，看似冲突
- **但** global_i 的 8 个并行实例 `i*8+0, i*8+1, ..., i*8+7` 分别落在 bank 0~7 和 bank 0~7，读操作本身不冲突
- 实际上读端 `global_i` 和 `global_i+2048` 落在同一 bank（因为 2048%16=0），这意味着需要双端口 RAM（2P）。每个 bank 只被同一个 PE 读两次（一次读 u，一次读 v），2P BRAM 可以在一个周期完成
- 综合策略：`#pragma HLS BIND_STORAGE variable=buf_A type=ram_2p impl=bram`

### 2.2 主循环 — CG_NTT_Kernel

```
函数 CG_NTT_Kernel(in_data, modulus, K_HALF, M, cg_twiddle, is_ntt):
    
    // 初始化：in_data → buf_A
    INIT: for i in 0..RING_DIM-1: buf_A[i] = in_data[i]
    
    // 主循环
    STAGE_LOOP: for stage = 0..STAGE-1:
        // NTT 正序，INTT 逆序
        actual_stage = is_ntt ? stage : (STAGE - 1 - stage)
        
        BUTTERFLY_LOOP: for i = 0..HALF_N/CG_PE_NUM-1:
            #pragma HLS PIPELINE II=1
            
            PE_UNROLL: for p = 0..CG_PE_NUM-1:
                #pragma HLS UNROLL
                
                global_i = i * CG_PE_NUM + p    // 0 ~ 2047
                
                // ① 固定跨度读取：永远读 i 和 i + N/2
                u = buf_read[global_i]
                v = buf_read[global_i + HALF_N]
                
                // ② 顺序读取旋转因子
                tf = cg_twiddle[actual_stage][global_i]
                
                // ③ 蝶形运算（复用现有 Configurable_PE）
                Configurable_PE(u, v, tf, out_u, out_v, modulus, K_HALF, M, is_ntt)
                
                // ④ 完美洗牌写回：永远写到 2i 和 2i+1
                buf_write[2 * global_i]     = out_u
                buf_write[2 * global_i + 1] = out_v
        
        // 乒乓切换：偶数 stage 读 A 写 B，奇数读 B 写 A
    
    // 回写
    WRITEBACK: for i in 0..RING_DIM-1: in_data[i] = buf_final[i]
```

**关键点：**
- 每个 stage 256 个周期（2048/8），12 个 stage 共 3072 周期
- 加上初始化和回写各 4096/8=512 周期 → **总计约 4096 周期**（vs 现有的 12×64=768 周期，但现有方案 MUX 资源消耗极大）

### 2.3 旋转因子预计算（Host 端软件，Testbench 中实现）

CG-NTT 的旋转因子不能直接复用标准 NTT 的 twiddle 表。必须在 Host 端模拟 CG-NTT 的数据流动，记录每层每个 PE 实际需要的旋转因子：

```
函数 build_cg_twiddle_table(tf_out[STAGE][HALF_N], mod, root_2N):
    
    // 标准 NTT 的旋转因子：stage s, butterfly group g, position within group
    // 在 CG-NTT 中，数据经过 perfect shuffle 后位置被打乱
    // 我们需要追踪每个数据元素的"逻辑位置"
    
    perm[RING_DIM]  // 追踪当前物理位置 → 逻辑位置的映射
    初始化: perm[i] = i
    
    for stage = 0..STAGE-1:
        // 在当前排列下，物理位置 global_i 和 global_i+HALF_N 
        // 对应的逻辑位置是 perm[global_i] 和 perm[global_i+HALF_N]
        
        // 标准 NTT stage s 中，蝶形对 (a, a+stride) 的旋转因子为：
        //   stride = N >> (s+1)
        //   group = a / stride
        //   tf = root_2N ^ (group * stride)  (这个取决于具体的 NTT 变体)
        
        // 简化方法：直接模拟标准 NTT，记录每个蝶形对使用的旋转因子，
        // 然后按 CG-NTT 的物理位置重新排列
        
        for i = 0..HALF_N-1:
            logical_a = perm[i]
            logical_b = perm[i + HALF_N]
            tf_out[stage][i] = 计算该蝶形对的旋转因子(logical_a, logical_b, stage)
        
        // 模拟 perfect shuffle：更新 perm
        new_perm[RING_DIM]
        for i = 0..HALF_N-1:
            new_perm[2*i]     = perm[i]         // out_u 的逻辑位置
            new_perm[2*i + 1] = perm[i + HALF_N]  // out_v 的逻辑位置
        perm = new_perm
```

### 2.4 输出重排（cg_ntt_reorder）

CG-NTT 跑完 12 层后，数据处于特定乱序状态。需要一个 bit-reversal + unshuffle 模块把数据导回正常顺序：

```
函数 cg_ntt_reorder(data[RING_DIM]):
    temp[RING_DIM]
    
    // 模拟 12 次 perfect shuffle 后的最终排列
    perm[RING_DIM]
    初始化 perm[i] = i
    for s = 0..STAGE-1:
        new_perm[RING_DIM]
        for i = 0..HALF_N-1:
            new_perm[2*i]     = perm[i]
            new_perm[2*i + 1] = perm[i + HALF_N]
        perm = new_perm
    
    // perm[physical] = original_logical
    // 逆映射：temp[perm[i]] = data[i]
    for i = 0..RING_DIM-1:
        temp[perm[i]] = data[i]
    
    拷贝 temp → data
```

这个重排在 HLS 中可以用一个简单的查找表实现（编译时确定的常量排列）。

---

## 3. cg_ntt_tb.cpp — 测试台

### 测试项：

| # | 测试 | 说明 |
|---|------|------|
| 1 | `test_cg_shuffle_permutation` | 验证 perfect shuffle 排列的数学性质 |
| 2 | `test_cg_twiddle_build` | 验证旋转因子表的正确性（与暴力计算对比） |
| 3 | `test_cg_ntt_roundtrip` | 单 limb NTT→INTT 往返验证 |
| 4 | `test_compute_cg_ntt_roundtrip` | 多 limb NTT→INTT 往返验证 |
| 5 | `test_cg_vs_standard_ntt` | CG-NTT 结果（重排后）与标准 NTT 结果一致 |

### 编译方式：
```bash
g++ -std=c++14 -O2 -DFPGA_STANDALONE_TEST \
    -I../include \
    cg_ntt_tb.cpp \
    ../src/cg_ntt.cpp \
    ../src/arithmetic.cpp \
    -o cg_ntt_tb && ./cg_ntt_tb
```

### 关键：旋转因子生成

测试台中需要实现完整的软件参考 CG-NTT：
1. 使用 NTT-friendly 质数 `p = 3221225473 (3×2^30+1)`，原根 `g = 5`
2. 计算 2N 次本原根 `root_2N = g^((p-1)/(2N))`
3. 用上述 `build_cg_twiddle_table` 算法生成旋转因子
4. 对 INTT，使用逆根 `inv_root_2N = root_2N^(p-2)` 生成逆旋转因子表

---

## 4. HLS 优化要点

| 优化项 | 策略 |
|--------|------|
| buf_A/buf_B 分区 | `cyclic factor=16` + `ram_2p bram` |
| PE 并行 | `#pragma HLS UNROLL` (factor=CG_PE_NUM) |
| 流水线 | `#pragma HLS PIPELINE II=1` 在 BUTTERFLY_LOOP |
| 依赖声明 | `#pragma HLS DEPENDENCE variable=buf_A/B inter false` |
| 旋转因子存储 | `#pragma HLS BIND_STORAGE type=rom_2p impl=bram` 或 URAM |
| 旋转因子分区 | 按 CG_PE_NUM cyclic 分区，保证 8 个 PE 同时读无冲突 |
| 重排查找表 | 编译时常量数组，`complete` 分区或流水线化 |

---

## 5. 验证步骤

1. **独立编译测试台**：`g++ -DFPGA_STANDALONE_TEST ...` → 5 项测试全部 PASS
2. **Vitis HLS C Simulation**：将 `cg_ntt.cpp` 加入 `csim.tcl`，验证功能
3. **关键指标**：NTT→INTT 往返恢复原始数据、CG-NTT 输出重排后与标准 NTT 一致
