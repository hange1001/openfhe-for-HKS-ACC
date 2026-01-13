#ifndef BCONV_H
#define BCONV_H

#include <ap_int.h>
#include <hls_stream.h>
#include "define.h"
#include "arithmetic.h"


// ------------------------------------------------
// [CONFIG] 维度定义 (Weight Stationary)
// ------------------------------------------------
// RING_DIM (M): 序列长度 (Time Steps / Batch Size) -> 16
// SIZE_Q   (K): 输入特征数 (Input Channels)        -> 4
// SIZE_P   (N): 输出特征数 (Output Channels)       -> 8
//static const int RING_DIM = 65536;
static const int SIZE_Q   = 4;
static const int SIZE_P   = 4;

// 顶层函数声明 (保持 4 个参数接口以匹配 Testbench)
extern "C" void bconv_systolic(
    const ap_uint<64>* x_linear,    // [RING_DIM * SIZE_Q]
    const ap_uint<64>* w_linear,    // [SIZE_Q * SIZE_P]
    const ap_uint<64> mod_linear,  // [SIZE_P]
    ap_uint<64>* out_linear         // [RING_DIM * SIZE_P]
);
void Compute_BConv(
     const uint64_t in_x[MAX_LIMBS][SQRT][SQRT], // 输入 X (片上 BRAM)
    const uint64_t in_w[MAX_LIMBS][SQRT][SQRT], // 输入 W (片上 BRAM)
    uint64_t       out[MAX_LIMBS][SQRT][SQRT],  // 输出 Out (片上 BRAM)
    const uint64_t MODULUS[MAX_LIMBS],
    int num_active_limbs,                // 当前有效层数 (比如 44)
    int mod_idx_offset               // 如果有模数偏移 (通常是 0)
);

#endif