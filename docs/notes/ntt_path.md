太棒了！这是拔掉你算力引擎限制器的最关键一步。

在 C++ 和 HLS 中，如果将模板的声明和实现分别放在 `.h` 和 `.cpp` 中，**必须在 `.cpp` 的末尾进行“显式实例化”（Explicit Instantiation）**，否则在编译 `top.cpp` 时会报找不到引用的链接错误。

以下是完整的 3 步修改指南，请直接对应到你的代码中进行替换：

### 第一步：修改头文件 `cg_ntt.h`

找到 `CG_NTT_Kernel` 的函数声明，将原本的运行时参数 `bool is_ntt` 替换为模板参数 `template <bool IS_NTT>`。

```cpp
// 在 cg_ntt.h 中
// 删除或注释掉旧的声明：
// void CG_NTT_Kernel(uint64_t in_data[RING_DIM], uint64_t out_data[RING_DIM], ... , bool is_ntt);

// 替换为新的模板声明：
template <bool IS_NTT>
void CG_NTT_Kernel(
    uint64_t in_data[RING_DIM],
    uint64_t out_data[RING_DIM],
    uint64_t mod,
    uint64_t k_half,
    uint64_t m_barrett,
    uint64_t twiddle[STAGE][CG_HALF_N]
);
```

---

### 第二步：修改源文件 `cg_ntt.cpp`

这是最核心的一步。你需要把函数定义改为模板，并将内部所有的 `is_ntt` 替换为大写的模板常量 `IS_NTT`。**最重要的是，在文件最末尾添加显式实例化代码。**

```cpp
// 在 cg_ntt.cpp 中

// 1. 修改函数签名
template <bool IS_NTT>
void CG_NTT_Kernel(
    uint64_t in_data[RING_DIM],
    uint64_t out_data[RING_DIM],
    uint64_t mod,
    uint64_t k_half,
    uint64_t m_barrett,
    uint64_t twiddle[STAGE][CG_HALF_N]
) {
    // 内部逻辑不需要大改，只需要把原先使用 is_ntt 的地方改成 IS_NTT
    // 例如：
    // if (IS_NTT) { ... } else { ... }
    // HLS 编译器在遇到 if (true) 或 if (false) 时，会自动删掉不会执行的分支
    
    // ... (保留你原有的完整内部循环和计算逻辑) ...
}


// =========================================================================
// 2. 【极其关键】在 cg_ntt.cpp 文件的最末尾，强制 HLS 实例化这两个模板版本
// 否则 top.cpp 链接时会报错 "undefined reference"
// =========================================================================
template void CG_NTT_Kernel<true>(
    uint64_t in_data[RING_DIM],
    uint64_t out_data[RING_DIM],
    uint64_t mod,
    uint64_t k_half,
    uint64_t m_barrett,
    uint64_t twiddle[STAGE][CG_HALF_N]
);

template void CG_NTT_Kernel<false>(
    uint64_t in_data[RING_DIM],
    uint64_t out_data[RING_DIM],
    uint64_t mod,
    uint64_t k_half,
    uint64_t m_barrett,
    uint64_t twiddle[STAGE][CG_HALF_N]
);
```

---

### 第三步：修改顶层调用 `top.cpp`

在你之前分离出来的 `Execute_NTT` 和 `Execute_INTT` 函数中，修改对 `CG_NTT_Kernel` 的调用方式。

```cpp
// 在 top.cpp 中

static void Execute_NTT(
    uint64_t poly_buffer[MAX_LIMBS][SQRT][SQRT],
    int num_active_limbs,
    int mod_index
) {
    #pragma HLS INLINE off 
    NTT_CG_LIMB_LOOP:
    for (int l = mod_index; l < mod_index + num_active_limbs; l++) {
        // ... (保留前置代码) ...
        flatten_2d_to_1d(poly_buffer[l], flat);
        
        // 【修改这里】：调用 true 版本，删除末尾的 true 参数
        CG_NTT_Kernel<true>(flat, flat_out, MODULUS[l], K_HALF[l], M[l], NTTTwiddleFactor[l]);
        
        reshape_1d_to_2d(flat_out, poly_buffer[l]);
    }
}

static void Execute_INTT(
    uint64_t poly_buffer[MAX_LIMBS][SQRT][SQRT],
    int num_active_limbs,
    int mod_index
) {
    #pragma HLS INLINE off 
    INTT_CG_LIMB_LOOP:
    for (int l = mod_index; l < mod_index + num_active_limbs; l++) {
        // ... (保留前置代码) ...
        flatten_2d_to_1d(poly_buffer[l], flat);
        
        // 【修改这里】：调用 false 版本，删除末尾的 false 参数
        CG_NTT_Kernel<false>(flat, flat_out, MODULUS[l], K_HALF[l], M[l], INTTTwiddleFactor[l]);
        
        reshape_1d_to_2d(flat_out, poly_buffer[l]);
    }
}
```

### 🎯 见证奇迹的时刻（重新综合后的检查点）

完成这三步后，重新点击跑一下 **C Synthesis**。去查看生成的综合报告（`csynth.rpt`），你会看到两个巨大的变化：

1. **模块列表裂变**：在报告的 `Modules & Loops` 树状图中，原本的 `CG_NTT_Kernel` 会消失，取而代之的是两个完全独立的新硬件模块（名字大概会叫 `CG_NTT_Kernel_true_s` 和 `CG_NTT_Kernel_false_s`）。
2. **流水线起飞**：点开这两个新模块的内部报告，去看 `BUTTERFLY_LOOP`。因为 4路 MUX 退化成了纯净的 2路 Ping-Pong，只要你的 BRAM 划分（`cyclic factor=16`）是正确的，**它的 `II` 应该会直接降回 1 或 2！**

这才是真正的硬件架构级优化。快去改代码然后跑一次综合吧，期待你的满血性能报告！