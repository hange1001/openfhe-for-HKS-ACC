# CG-NTT 迁移技术报告

## 1. 背景与动机

本项目将 FPGA 加速器中的标准 Cooley-Tukey NTT 替换为 CG-NTT（Constant Geometry NTT，恒定几何 NTT）。CG-NTT 的核心优势在于每一层蝶形单元的读写地址模式完全固定（读 `[i, i+N/2]`，写 `[2i, 2i+1]`），消除了传统 NTT 中随 stage 变化的交叉开关（MUX），使地址生成变为硬连线，显著降低了 FPGA 布线复杂度。

迁移过程中暴露并修复了两个架构级缺陷：
- **AXI 地址与 BRAM 索引耦合**：导致 `mod_index > 0` 时越界崩溃
- **CG-NTT 输出顺序不兼容**：导致 OpenFHE 端到端 CKKS 测试全部失败

---

## 2. 变更文件清单

| 文件 | 变更类型 | 涉及阶段 |
|------|----------|----------|
| `src/fpga_backend/src/top.cpp` | 重写 | Phase 1 + 2 |
| `src/fpga_backend/include/top.h` | 修改 | Phase 1 |
| `src/fpga_backend/src/cg_ntt.cpp` | 修改 | Phase 2 |
| `src/fpga_backend/src/load.cpp` | 已验证无需改动 | Phase 2 |
| `src/fpga_backend/testbench/top_tb.cpp` | 新增 | Phase 1 |
| `src/core/include/FpgaManager.h` | 重写 | Phase 1 + 2 + 3 |

---

## 3. Phase 1：NTT → CG-NTT 内核替换

### 3.1 删除的组件

| 组件 | 说明 |
|------|------|
| `NTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM]` | 片上 URAM 旋转因子存储（NTT） |
| `INTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM]` | 片上 URAM 旋转因子存储（INTT） |
| `InterLeave()` 调用 | 标准 NTT 的前/后置重排 |
| `Compute_NTT()` 调用 | 标准多 limb NTT 包装器 |
| `#include "ntt_kernel.h"` / `#include "interleave.h"` | top.cpp 中的旧头文件 |

### 3.2 新增的数据流

**旧流程（标准 NTT）：**
```
OP_NTT:  Load → InterLeave(true) → Compute_NTT → Store
OP_INTT: Load → Compute_NTT(false) → InterLeave(false) → Store
```

**新流程（CG-NTT）：**
```
OP_NTT:  burst_copy(mem_in1 → mem_out) → Compute_CG_NTT(mem_out, mem_in2, ..., true)
OP_INTT: burst_copy(mem_in1 → mem_out) → Compute_CG_NTT(mem_out, mem_in2, ..., false)
```

关键变化：
- `Compute_CG_NTT` 自行管理 DDR ↔ 片上 BRAM 的搬运，不再依赖外部 `Load`/`Store`
- 旋转因子不再存片上 URAM，改为每次调用时通过 `mem_in2`（AXI gmem1）从 DDR 传入
- `OP_INIT` 不再加载旋转因子，仅保留标量参数（MODULUS、K_HALF、M）

### 3.3 旋转因子架构变化

| 属性 | 标准 NTT | CG-NTT |
|------|----------|--------|
| 存储位置 | 片上 URAM（OP_INIT 时加载） | DDR（每次调用传入） |
| 布局 | `[MAX_LIMBS][PE_PARALLEL][RING_DIM]` | `[MAX_LIMBS][STAGE][CG_HALF_N]` |
| 生成方式 | `psi^i` 幂次 + `GenerateTwiddleIndices` 排列 | `BuildCGTwiddle`：追踪 bit-reverse + perfect shuffle 排列 |
| 大小/limb | `PE_PARALLEL × RING_DIM = 8 × 4096 = 32K` uint64 | `STAGE × CG_HALF_N = 12 × 2048 = 24K` uint64 |

### 3.4 FpgaManager.h 主机侧变更

- 新增 `MathUtils::BitReverse()` 和 `MathUtils::BuildCGTwiddle()` 静态方法
- `InitModuli()` 中旋转因子生成改为调用 `BuildCGTwiddle`，结果存入 `m_ntt_twiddle` / `m_intt_twiddle` 成员变量
- 新增 `ExecuteNTT()` 方法，与通用 `Execute()` 分离，独立管理旋转因子 XRT Buffer
- `NttForwardOffload` / `NttInverseOffload` 改为调用 `ExecuteNTT`

---

## 4. Phase 2：AXI 地址与 BRAM 索引解耦

### 4.1 问题根因

`mod_index` 在旧代码中身兼两职：

1. **BRAM 参数索引**：`MODULUS[mod_index]`、旋转因子表 `twiddle[mod_index][...]` → 正确
2. **AXI 物理偏移**：`mem_in[mod_index * PACKED_RING_DIM + i]` → **错误**

主机通过 XRT 分配的 Buffer 大小为 `num_active_limbs × RING_DIM`，数据从偏移 0 开始。当 `mod_index = 3` 时，FPGA 试图从 `3 × 512 = 1536` 号 512-bit 字开始读取，远超 Buffer 边界，触发 OOB 段错误。

### 4.2 修复原则

**AXI 读写强行归零，BRAM 查表保留索引。**

引入局部索引 `axi_l = l - mod_idx_offset`：

| 访问目标 | 索引 | 原因 |
|----------|------|------|
| `in_data[...]` / `mem_out[...]` | `axi_l`（局部） | Host Buffer 从 0 开始，仅含 active limbs |
| `cg_ntt_twiddle[...]` / `cg_intt_twiddle[...]` | `l`（全局） | Host 上传完整旋转因子表（所有 limbs） |
| `modulus[l]` / `K_HALF[l]` / `M_barrett[l]` | `l`（全局） | 片上参数表在 OP_INIT 时按全局索引填充 |

### 4.3 具体代码变更

**`top.cpp` OP_NTT/OP_INTT：**
```cpp
// Before:
int base = mod_index * (RING_DIM / PW);  // BUG
mem_out[base + i] = mem_in1[base + i];

// After:
mem_out[i] = mem_in1[i];  // AXI 从 0 开始
```

**`cg_ntt.cpp` Compute_CG_NTT：**
```cpp
// Before:
ap_uint<512> pack = in_data[l * PACKED_RING_DIM + i];       // BUG: l 是全局索引

// After:
int axi_l = l - mod_idx_offset;                              // 局部索引
ap_uint<512> pack = in_data[axi_l * PACKED_RING_DIM + i];   // 修复
// 旋转因子仍用全局 l：cg_ntt_twiddle[l * PACKED_TW_SIZE + i]
// 片上参数仍用全局 l：modulus[l], K_HALF[l], M_barrett[l]
```

### 4.4 已验证无需改动的模块

- **`Load()` / `Store()`**：已使用 `(l - mod_index)` 作为 AXI 地址，`l` 作为 BRAM 索引
- **OP_BCONV**：`Load(mem_in1, ..., LIMB_Q, 0)` 固定 mod_index=0，无越界风险
- **OP_AUTO**：通过 Load/Store 间接访问，已正确解耦

---

## 5. Phase 3：CG-NTT 输出顺序兼容性修复

### 5.1 问题根因

CG-NTT 经过 STAGE 次 perfect shuffle 后，输出数据处于 **bit-reversed 顺序**，而 OpenFHE 期望 NTT 输出为**标准顺序**。

在 CKKS 加密流程中：
```
Encrypt: c₀ = A · s + e + m
```
- `A` 由 OpenFHE 在主机侧生成，处于标准顺序 NTT 域
- `s` 经 FPGA NTT 后返回 bit-reversed 顺序
- 主机执行 `A[i] * s[rev(i)]`，数学上完全错误
- 密文从生成起就已损坏，后续所有 CKKS 操作（EvalAdd、EvalMult、EvalRotate）必然失败

### 5.2 为什么单元测试未发现

`top_tb.cpp` 的 `test_op_ntt_vs_sw()` 在验证时隐式调用了 `cg_ntt_reorder(poly_hw)`：
```cpp
Top(pin, g_ntt_tw, pout, OP_NTT, 1, mod_idx);
unpack_512_to_u64(pout, poly_hw, N);
cg_ntt_reorder(poly_hw);  // ← 测试平台做了重排，但硬件内核没有！
```
测试通过了，但实际部署路径（`FpgaManager::NttForwardOffload`）直接返回了未重排的裸数据。

### 5.3 修复方案：主机侧软件补丁

无需重编 HLS。Bit-reversal 是自反操作（执行两次等于恒等），在主机侧 C++ 中拦截即可：

**NttForwardOffload（NTT 后重排）：**
```cpp
ExecuteNTT(OP_NTT, in, m_ntt_twiddle.data(), out, 1, mod_idx);

// CG-NTT 输出为 bit-reversed，还原为标准顺序
int bits = (int)std::log2(n);
std::vector<uint64_t> temp(n);
for (size_t i = 0; i < n; i++) {
    int rev = MathUtils::BitReverse((int)i, bits);
    temp[rev] = out[i];
}
std::memcpy(out, temp.data(), n * sizeof(uint64_t));
```

**NttInverseOffload（INTT 前预打乱）：**
```cpp
// OpenFHE 给的是标准顺序，CG-INTT 期望 bit-reversed 输入
int bits = (int)std::log2(n);
std::vector<uint64_t> scrambled_in(n);
for (size_t i = 0; i < n; i++) {
    int rev = MathUtils::BitReverse((int)i, bits);
    scrambled_in[rev] = in[i];
}
ExecuteNTT(OP_INTT, scrambled_in.data(), m_intt_twiddle.data(), out, 1, mod_idx);
```

### 5.4 对 KeySwitching 链路的影响

EvalRotate 依赖 KeySwitching，其底层数据流为：
```
INTT → BConv → NTT → Mult
```
NTT/INTT 的顺序错误会在 KeySwitching 中产生天文数字级噪声，导致 `Decode()` 抛出近似误差。修复 NTT 端到端时序后，EvalRotate 自然恢复正常。

---

## 6. 架构总览

```
┌─────────────────────────────────────────────────────────┐
│                    Host (FpgaManager)                    │
│                                                         │
│  NttForwardOffload:                                     │
│    ExecuteNTT(OP_NTT) → BitReverse reorder (SW)         │
│                                                         │
│  NttInverseOffload:                                     │
│    BitReverse scramble (SW) → ExecuteNTT(OP_INTT)       │
│                                                         │
│  InitModuli:                                            │
│    BuildCGTwiddle → m_ntt_twiddle / m_intt_twiddle      │
│    OP_INIT: 仅标量参数 (MODULUS, K_HALF, M)             │
└──────────────────────┬──────────────────────────────────┘
                       │ PCIe / AXI
┌──────────────────────▼──────────────────────────────────┐
│                    FPGA (Top Kernel)                     │
│                                                         │
│  OP_NTT/OP_INTT:                                        │
│    burst_copy(mem_in1→mem_out, base=0)                  │
│    Compute_CG_NTT(                                      │
│      in_data:  AXI 局部索引 axi_l                       │
│      twiddle:  AXI 全局索引 l                            │
│      modulus:  片上参数 全局索引 l                        │
│    )                                                    │
│                                                         │
│  OP_ADD/SUB/MULT:                                       │
│    Load(AXI局部) → Compute → Store(AXI局部)             │
│    BRAM buffer 使用全局索引                               │
└─────────────────────────────────────────────────────────┘
```

---

## 7. 遗留事项

| 项目 | 状态 | 说明 |
|------|------|------|
| 主机侧 BitReverse 开销 | 可接受 | N=4096 的 bit-reversal 约 ~10μs，远小于 PCIe 传输延迟 |
| 硬件内嵌 reorder | 待评估 | 可在 `Compute_CG_NTT` 的 STORE_OUT 阶段以零额外延迟完成 bit-reversal 写回，消除主机侧补丁 |
| Makefile 中 `interleave.cpp` | 保留 | 虽然 top.cpp 不再使用，但未确认其他模块是否依赖，暂不移除 |
| `SRC_Top` 中 `cg_ntt.cpp` 重复 | 已知 | `KERNEL_SRC` 和 `SRC_Top` 均包含 `cg_ntt.cpp`，可能导致 HLS 链接警告 |
