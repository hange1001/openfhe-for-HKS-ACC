#include "../include/arithmetic.h"

void AddMod(uint64_t &a, const uint64_t &b, const uint64_t &mod, const bool &is_add) {
    #pragma HLS inline
    if(is_add) {
        uint64_t sum = a + b;
        #pragma HLS bind_op variable=sum op=add impl=fabric latency=1
        uint64_t sub = sum - mod;
        a = (sum >= mod) ? sub : sum;
    } else {
        uint64_t diff = a - b;
        #pragma HLS bind_op variable=diff op=sub impl=fabric latency=1
        uint64_t add_mod = a + mod;
        uint64_t diff_mod = add_mod - b;
        a = (a >= b) ? diff : diff_mod;
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
    #pragma HLS LATENCY min=2

    if (mod == 0) {
        res_mod = 0;
        return;
    }

#ifdef FPGA_STANDALONE_TEST
    unsigned __int128 res_mult = (unsigned __int128)a * b;
    uint64_t res_mult_H = (uint64_t)(res_mult >> 64);
    uint64_t res_mult_L = (uint64_t)(res_mult);

    unsigned __int128 p_high = (unsigned __int128)res_mult_H * m;
    unsigned __int128 p_low  = (unsigned __int128)res_mult_L * m;

    unsigned __int128 q = 0;
    if (S > 64) {
        q = (p_high >> (S - 64)) + (p_low >> S);
    } else if (S < 64) {
        q = (p_high << (64 - S)) + (p_low >> S);
    } else {
        q = p_high + (p_low >> 64);
    }

    unsigned __int128 q_times_mod = q * (unsigned __int128)mod;
    unsigned __int128 r_full = res_mult - q_times_mod;
    uint64_t r = (uint64_t)r_full;

    // 5. 校正（全精度下误差极小，最多 3 次）
    if (r >= mod) { r -= mod; }
    if (r >= mod) { r -= mod; }

    res_mod = r;
#else
    // 1. 计算乘积
    ap_uint<128> res_mult = (ap_uint<128>)a * (ap_uint<128>)b;
    #pragma HLS bind_op variable=res_mult op=mul impl=dsp latency=4

    // 2. 拆分高低 64 位
    ap_uint<64> res_mult_H = (ap_uint<64>)(res_mult >> 64);
    ap_uint<64> res_mult_L = (ap_uint<64>)(res_mult);

    // 3. 计算 p_high 和 p_low（使用 Barrett 预计算因子 m）
    ap_uint<128> p_high = (ap_uint<128>)res_mult_H * (ap_uint<128>)m;
    #pragma HLS bind_op variable=p_high op=mul impl=dsp latency=4

    ap_uint<128> p_low = (ap_uint<128>)res_mult_L * (ap_uint<128>)m;
    #pragma HLS bind_op variable=p_low op=mul impl=dsp latency=4

    // 4. 计算 q（S = bitwidth(mod) + 62）
    //    S 按 case 离散化：每个 case 内移位量为编译期常数，HLS 综合为硬连线，
    //    彻底消除 128-bit 动态桶形移位器。default 分支保留通用兜底。
    //    OpenFHE 常用模数位宽：60-bit → S=122，61-bit → S=123，59-bit → S=121
    ap_uint<128> shift_high = 0;
    ap_uint<128> shift_low  = 0;
    switch (S) {
        case 122: shift_high = p_high >> 58;  shift_low = p_low >> 122; break; // 60-bit mod
        case 123: shift_high = p_high >> 59;  shift_low = p_low >> 123; break; // 61-bit mod
        case 121: shift_high = p_high >> 57;  shift_low = p_low >> 121; break; // 59-bit mod
        case 120: shift_high = p_high >> 56;  shift_low = p_low >> 120; break; // 58-bit mod
        case 112: shift_high = p_high >> 48;  shift_low = p_low >> 112; break; // 50-bit mod
        default:
            if      (S > 64) { shift_high = p_high >> (S - 64); shift_low = p_low >> S; }
            else if (S < 64) { shift_high = p_high << (64 - S); shift_low = p_low >> S; }
            else             { shift_high = p_high;              shift_low = p_low >> 64; }
            break;
    }

    ap_uint<128> q = shift_high + shift_low;

    // 5. 计算 q * mod
    ap_uint<128> q_times_mod = q * (ap_uint<128>)mod;
    #pragma HLS bind_op variable=q_times_mod op=mul impl=dsp latency=4

    // 6. 强行给 128-bit 大减法打一拍解决时序
    ap_uint<128> r_full = res_mult - q_times_mod;
    #pragma HLS bind_op variable=r_full op=sub impl=fabric latency=2

    uint64_t r_stg1 = (uint64_t)r_full;

    // 7. 第一次模校正，强制打拍
    uint64_t sub1 = r_stg1 - mod;
    #pragma HLS bind_op variable=sub1 op=sub impl=fabric latency=2
    uint64_t r_stg2 = (r_stg1 >= mod) ? sub1 : r_stg1;

    // 8. 第二次模校正，强制打拍
    uint64_t sub2 = r_stg2 - mod;
    #pragma HLS bind_op variable=sub2 op=sub impl=fabric latency=2
    uint64_t r_stg3 = (r_stg2 >= mod) ? sub2 : r_stg2;

    res_mod = r_stg3;
#endif
}