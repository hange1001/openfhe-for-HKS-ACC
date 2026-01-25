#ifndef BCONV_H
#define BCONV_H

#include <ap_int.h>
#include <hls_stream.h>
#include "define.h"
#include "arithmetic.h"

// 顶层函数声明 - 固定维度 LIMB_Q x LIMB_P
extern "C" {
    // Host wrapper: uses uint64_t for C++ compatibility
    // in_w is 2D: weights[q][p] for each (input limb, output limb) pair
    void Compute_BConv(
        uint64_t in_x[LIMB_Q + LIMB_P][SQRT][SQRT], 
        const uint64_t in_w[LIMB_Q][LIMB_P], 
        const uint64_t MODULUS[LIMB_Q + LIMB_P],
        int num_active_limbs,
        int mod_idx_offset = 0
    );
}

#endif