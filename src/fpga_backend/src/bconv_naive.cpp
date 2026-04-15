#include "../include/bconv_naive.h"

Compute_BConv_Naive(
    uint64_t in_x[MAX_LIMBS][SQRT][SQRT],
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],
    const uint64_t out_mod[MAX_OUT_COLS],
    int sizeP
){
    #pragma HLS INLINE
    
    // 朴素BConv实现：直接在2D数组上进行三重循环计算
    // 对于每个输出模数p，计算所有系数位置(r,c)的线性组合
    
    Output_Loop: for (int p = 0; p < sizeP; ++p) {
        Row_Loop: for (int r = 0; r < SQRT; ++r) {
            Col_Loop: for (int c = 0; c < SQRT; ++c) {
                #pragma HLS PIPELINE II=1
                
                uint64_t sum = 0;
                
                // 对每个输入模数q进行乘积累加
                Limb_Loop: for (int q = 0; q < LIMB_Q; ++q) {
                    // 使用128位乘法避免溢出
                    uint128_t prod = (uint128_t)in_x[q][r][c] * in_w[q][p];
                    
                    // 模运算：先对prod取模，然后累加
                    uint64_t prod_mod = prod % out_mod[p];
                    sum = (sum + prod_mod) % out_mod[p];
                }
                
                // 将结果写回输出位置
                in_x[LIMB_Q + p][r][c] = sum;
            }
        }
    }
}