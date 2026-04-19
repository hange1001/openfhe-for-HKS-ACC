#include "../include/arithmetic.h"

void AddMod(
    uint64_t &a,
    const uint64_t &b,
    const uint64_t &mod,
    const bool &is_add
){
    #pragma HLS INLINE
    unsigned __int128 temp_res;
    if (is_add) { 
        temp_res = (unsigned __int128)a + b;
        if (temp_res >= mod) {
            temp_res -= mod;
        }
        a = (uint64_t)temp_res;
    } 
    else {
        if (a >= b) {
            a = a - b;
        } else {
            temp_res = (unsigned __int128)a + mod - b;
            a = (uint64_t)temp_res;
        }
    }
}

void Karatsuba(
    const uint64_t &a,      
    const uint64_t &b,      
    uint128_t &result      
){
    const uint64_t MASK_32 = 0xFFFFFFFF;
    
    // 1. 拆分高低位
    uint64_t a_low  = a & MASK_32;
    uint64_t a_high = a >> 32;
    uint64_t b_low  = b & MASK_32;
    uint64_t b_high = b >> 32;

    // 2. 计算三个中间乘积
    uint64_t z0 = a_low * b_low;
    uint64_t z2 = a_high * b_high;
    uint128_t z1 = (uint128_t)(a_low + a_high) * (b_low + b_high);

    // 3. 计算中间项 mid = z1 - z2 - z0
    uint128_t mid = z1 - z2 - z0;

    // 4. 移位拼接
    result = ((uint128_t)z2 << 64) + (mid << 32) + z0;
}


// =================================================
// Barrett Modular Multiplication (全精度宽位宽版)
// =================================================
void MultMod(
    const uint64_t &a,
    const uint64_t &b,
    const uint64_t &mod,
    const uint64_t &m,
    const uint64_t &S,      // 总移位量 S = bitwidth(mod) + 62
    uint64_t &res_mod
){
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1

    if (mod == 0) {
        res_mod = 0;
        return;
    }

    // 1. 真实乘积（最高约 120-bit）
    ap_uint<128> res_mult = (ap_uint<128>)a * b;
    #pragma HLS BIND_OP variable=res_mult op=mul impl=dsp latency=4

    ap_uint<64> res_mult_H = res_mult.range(127,64);
    ap_uint<64> res_mult_L = res_mult.range(63,0);

    // 2. 全精度：不截断任何低位，直接 128-bit × 64-bit → 192-bit
    //ap_uint<192> res_mult_shift = (ap_uint<192>)res_mult * m;

    ap_uint<128> p_high = (ap_uint<128>)res_mult_H * m;
    ap_uint<128> p_low  = (ap_uint<128>)res_mult_L * m;


    #pragma HLS BIND_OP variable=p_high op=mul impl=dsp latency=4
    #pragma HLS BIND_OP variable=p_low op=mul impl=dsp latency=4

    // 3. 右移 S 位得到精准商 q (提前在外部声明作用域)
    ap_uint<128> q = 0;
    
    if (S > 64) {
        // 高位部分右移 S-64，低位部分右移 S 
        q = (p_high >> (S - 64)) + (p_low >> S);
    } 
    else if (S < 64) {
        // S较小时，高位部分需要补偿左移，低位部分右移 S
        q = (p_high << (64 - S)) + (p_low >> S);
    } 
    else { 
        // S == 64 的边界情况
        q = p_high + (p_low >> 64);
    }

    // 4. 余数 r = z - q * mod
    ap_uint<128> q_times_mod = q * (ap_uint<128>)mod;
    ap_uint<128> r_full = res_mult - q_times_mod;
    uint64_t r = (uint64_t)r_full;

    // 5. 校正（全精度下误差极小，最多 3 次）
    if (r >= mod) { r -= mod; }
    if (r >= mod) { r -= mod; }
    if (r >= mod) { r -= mod; }

    res_mod = r;
}