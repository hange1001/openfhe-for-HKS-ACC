#include "../include/bconv_systolic.h"
#include "../include/arithmetic.h"
#include <hls_stream.h>

static const int TOTAL_CYCLES = LIMB_Q + RING_DIM + MAX_OUT_COLS - 1;

// ==================================================================
// P2：无局部中转。脉动核心直接读写工作区 in_x（top 的 poly_buffer_1）：
//   Feed_X  读  in_x[q][(t-q)>>LOG_SQRT][(t-q)&(SQRT-1)]           q<active_q
//   Collect 写  in_x[LIMB_Q+p][(t-3-p)>>LOG_SQRT][(t-3-p)&(SQRT-1)]  p<sizeP
// 每 limb 是独立 ram_2p 实例，每周期每 limb 至多 1 读+1 写；
// 逐周期 bank 无冲突由 check_systolic_banks.py 验证，并用调度报告复核 II=1。
// ==================================================================
static_assert((SQRT & (SQRT - 1)) == 0, "SQRT must be power of 2");
static constexpr int SQRT_MASK = SQRT - 1;  // SQRT=64 时取模用

// =============================================================
// bconv_systolic_core: 脉动阵列核心计算（直接读写工作区）
// =============================================================
static void bconv_systolic_core(
    uint64_t in_x[MAX_LIMBS][SQRT][SQRT],
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],
    const uint64_t out_mod[MAX_OUT_COLS],
    const uint64_t out_S[MAX_OUT_COLS],
    const uint64_t out_m_barrett[MAX_OUT_COLS],
    int sizeP,
    int active_q
) {
    #pragma HLS INLINE

    #pragma HLS ARRAY_PARTITION variable=in_w          complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_mod       complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_S         complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_m_barrett complete dim=0

    ap_uint<64> x_reg  [LIMB_Q]      [MAX_OUT_COLS + 1];
    uint64_t    sum_reg[MAX_OUT_COLS] [LIMB_Q + 1];
    #pragma HLS ARRAY_PARTITION variable=x_reg   complete dim=0
    #pragma HLS ARRAY_PARTITION variable=sum_reg complete dim=0

    Init_X_Reg: for (int q = 0; q < LIMB_Q; ++q) {
    #pragma HLS UNROLL
        for (int p = 0; p <= MAX_OUT_COLS; ++p) {
    #pragma HLS UNROLL
            x_reg[q][p] = 0;
        }
    }

    Init_Sum_Reg: for (int p = 0; p < MAX_OUT_COLS; ++p) {
    #pragma HLS UNROLL
        for (int q = 0; q <= LIMB_Q; ++q) {
    #pragma HLS UNROLL
            sum_reg[p][q] = 0;
        }
    }

    Systolic_Loop: for (int t = 0; t < TOTAL_CYCLES; ++t) {
    #pragma HLS PIPELINE II=1

        ap_uint<64> x_curr  [LIMB_Q]      [MAX_OUT_COLS + 1];
        uint64_t    sum_curr[MAX_OUT_COLS] [LIMB_Q + 1];
        #pragma HLS ARRAY_PARTITION variable=x_curr   complete dim=0
        #pragma HLS ARRAY_PARTITION variable=sum_curr complete dim=0

        Save_X: for (int q = 0; q < LIMB_Q; ++q) {
        #pragma HLS UNROLL
            for (int p = 0; p <= MAX_OUT_COLS; ++p) {
        #pragma HLS UNROLL
                x_curr[q][p] = x_reg[q][p];
            }
        }

        Save_Sum: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS UNROLL
            for (int q = 0; q <= LIMB_Q; ++q) {
        #pragma HLS UNROLL
                sum_curr[p][q] = sum_reg[p][q];
            }
        }

        Feed_X: for (int q = 0; q < LIMB_Q; ++q) {
        #pragma HLS UNROLL
            int data_idx = t - q;
            // 无效行显式注入 0：不靠旧值乘零，也不读未装载的 local_in_x 旧缓存。
            // 地址负值用 safe_idx 规避，值由 valid 选择，保证寻址范围合法。
            bool x_valid = (q < active_q) && (data_idx >= 0) && (data_idx < RING_DIM);
            int safe_idx = x_valid ? data_idx : 0;
            x_reg[q][0] = x_valid
                          ? (ap_uint<64>)in_x[q][safe_idx >> LOG_SQRT][safe_idx & SQRT_MASK]
                          : (ap_uint<64>)0;
        }

        Init_Sum: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS UNROLL
            sum_reg[p][0] = 0;
        }

        PE_Row: for (int q = 0; q < LIMB_Q; ++q) {
        #pragma HLS UNROLL
            PE_Col: for (int p = 0; p < MAX_OUT_COLS; ++p) {
            #pragma HLS UNROLL
                uint64_t x_in    = (uint64_t)x_curr[q][p];
                uint64_t sum_in  = sum_curr[p][q];
                uint64_t mod_p   = out_mod[p];

                uint64_t prod = 0;
                MultMod(x_in, in_w[q][p], mod_p, out_m_barrett[p], out_S[p], prod);

                // prod < mod_p，sum_in < mod_p → 和至多 2*mod_p，一次条件减即可
                uint64_t sum_out = sum_in + prod;
                if (sum_out >= mod_p) sum_out -= mod_p;

                x_reg[q][p + 1]   = x_in;
                sum_reg[p][q + 1] = sum_out;
            }
        }

        // 直接从 t 计算输出下标，消除 valid_count 的跨迭代依赖；写回工作区补基塔。
        Collect: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS UNROLL
            int out_idx = t - (LIMB_Q + p);
            if (p < sizeP && out_idx >= 0 && out_idx < RING_DIM) {
                in_x[LIMB_Q + p][out_idx >> LOG_SQRT][out_idx & SQRT_MASK] = sum_reg[p][LIMB_Q];
            }
        }
    }
}

// =================================================
// Compute_BConv_Systolic: 顶层接口
// =================================================
void Compute_BConv_Systolic(
    uint64_t in_x[MAX_LIMBS][SQRT][SQRT],
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],
    const uint64_t out_mod[MAX_OUT_COLS],
    const uint64_t out_S[MAX_OUT_COLS],
    const uint64_t out_m_barrett[MAX_OUT_COLS],
    int sizeP,
    int active_q
) {
    // in_x 是 top.cpp 里的片上 poly_buffer（BRAM），不是 DDR；
    // dim=3 按 8 路交错，Feed_X/Collect 的直接寻址逐周期无 bank 冲突。
    #pragma HLS ARRAY_PARTITION variable=in_x type=cyclic factor=8 dim=3

    // 标量/小数组参数走 AXI-Lite 控制（仅 cosim 顶层需要；top.cpp 调用时会被忽略）
    #pragma HLS INTERFACE s_axilite port=out_mod bundle=control
    #pragma HLS INTERFACE s_axilite port=sizeP   bundle=control
    #pragma HLS INTERFACE s_axilite port=active_q bundle=control
    #pragma HLS INTERFACE s_axilite port=return  bundle=control

    // 只保留小型 weights/modulus 寄存器与阵列流水寄存器；
    // local_in_x/local_out_x 及 Load_X/Store_X 全量中转已删除。
    uint64_t local_w  [LIMB_Q][MAX_OUT_COLS];
    uint64_t local_mod[MAX_OUT_COLS];
    uint64_t local_S  [MAX_OUT_COLS];
    uint64_t local_m_barrett[MAX_OUT_COLS];
    #pragma HLS ARRAY_PARTITION variable=local_w          complete dim=0
    #pragma HLS ARRAY_PARTITION variable=local_mod        complete dim=0
    #pragma HLS ARRAY_PARTITION variable=local_S          complete dim=0
    #pragma HLS ARRAY_PARTITION variable=local_m_barrett  complete dim=0

    Load_W: for (int q = 0; q < LIMB_Q; ++q) {
        for (int p = 0; p < MAX_OUT_COLS; ++p) {
            #pragma HLS PIPELINE II=1
            local_w[q][p] = in_w[q][p];
        }
    }

    Load_Mod: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS PIPELINE II=1
        local_mod[p]        = out_mod[p];
        local_S[p]          = out_S[p];
        local_m_barrett[p]  = out_m_barrett[p];
    }

    bconv_systolic_core(in_x, local_w, local_mod, local_S, local_m_barrett, sizeP, active_q);
}

