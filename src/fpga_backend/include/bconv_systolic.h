#ifndef BCONV_SYSTOLIC_H
#define BCONV_SYSTOLIC_H

#include <ap_int.h>
#include <hls_stream.h>
#include "define.h"

// BConv 脉动阵列实现：LIMB_Q × MAX_OUT_COLS PE 网格
// 用硬件整数除法（ap_uint % mod）取代 Barrett 约减，便于功能验证

extern "C" {
    // in_x  : [0..LIMB_Q-1] 为输入 limb，[LIMB_Q..LIMB_Q+sizeP-1] 为输出槽
    // in_w  : 权重矩阵 [LIMB_Q][MAX_OUT_COLS]
    // out_mod: 各输出列模数，无效列填非零占位值（避免除零）
    // sizeP : 有效输出列数（1 ~ MAX_OUT_COLS）
    void Compute_BConv_Systolic(
        uint64_t in_x[MAX_LIMBS][SQRT][SQRT],
        const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],
        const uint64_t out_mod[MAX_OUT_COLS],
        int sizeP
    );
}

#endif
