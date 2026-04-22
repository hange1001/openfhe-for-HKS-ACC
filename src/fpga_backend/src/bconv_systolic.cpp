#include "../include/bconv_systolic.h"
#include <hls_stream.h>

static const int TOTAL_CYCLES = LIMB_Q + RING_DIM + MAX_OUT_COLS - 1;

// ==================================================================
// I/O 带宽参数：匹配 Alveo 512-bit AXI（8 × 64-bit/拍）
// ==================================================================
constexpr int LOAD_PAR = 8;
static_assert((SQRT & (SQRT - 1)) == 0, "SQRT must be power of 2");
static_assert((RING_DIM % LOAD_PAR) == 0, "RING_DIM must be multiple of LOAD_PAR");
static constexpr int SQRT_MASK = SQRT - 1;  // SQRT=64 时取模用

// =============================================================
// bconv_systolic_core: 脉动阵列核心计算
// =============================================================
static void bconv_systolic_core(
    const uint64_t in_x[LIMB_Q][RING_DIM],
    uint64_t out_x[MAX_OUT_COLS][RING_DIM],
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],
    const uint64_t out_mod[MAX_OUT_COLS],
    int sizeP
) {
    #pragma HLS INLINE

    #pragma HLS ARRAY_PARTITION variable=in_w    complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_mod complete dim=0

    ap_uint<64>  x_reg  [LIMB_Q]      [MAX_OUT_COLS + 1];
    ap_uint<128> sum_reg[MAX_OUT_COLS] [LIMB_Q + 1];
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

        ap_uint<64>  x_curr  [LIMB_Q]      [MAX_OUT_COLS + 1];
        ap_uint<128> sum_curr[MAX_OUT_COLS] [LIMB_Q + 1];
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
            x_reg[q][0] = (data_idx >= 0 && data_idx < RING_DIM)
                          ? (ap_uint<64>)in_x[q][data_idx]
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
                ap_uint<64>  x_in   = x_curr[q][p];
                ap_uint<128> sum_in = sum_curr[p][q];
                ap_uint<64>  mod_p  = out_mod[p];

                // 守卫零模数（非活跃列 out_mod 可能为 0）
                ap_uint<128> prod    = (mod_p > 0)
                    ? (ap_uint<128>)(((ap_uint<128>)x_in * (ap_uint<128>)in_w[q][p]) % mod_p)
                    : (ap_uint<128>)0;
                ap_uint<128> sum_out = (mod_p > 0)
                    ? (ap_uint<128>)((sum_in + prod) % mod_p)
                    : (ap_uint<128>)0;

                x_reg[q][p + 1]    = x_in;
                sum_reg[p][q + 1]  = sum_out;
            }
        }

        // 直接从 t 计算输出下标，消除 valid_count 的跨迭代依赖
        Collect: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS UNROLL
            int out_idx = t - (LIMB_Q + p);
            if (p < sizeP && out_idx >= 0 && out_idx < RING_DIM) {
                out_x[p][out_idx] = (uint64_t)sum_reg[p][LIMB_Q];
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
    int sizeP
) {
    #pragma HLS INTERFACE m_axi port=in_x   bundle=gmem0 \
        max_read_burst_length=64 max_write_burst_length=64 \
        num_read_outstanding=8   num_write_outstanding=8
    #pragma HLS INTERFACE m_axi port=in_w   bundle=gmem1 \
        max_read_burst_length=64 num_read_outstanding=8
    #pragma HLS INTERFACE s_axilite port=out_mod bundle=control
    #pragma HLS INTERFACE s_axilite port=sizeP   bundle=control
    #pragma HLS INTERFACE s_axilite port=return  bundle=control

    // 环维度按 LOAD_PAR 交错，使 Load/Store 单拍可并行 8 路访问
    uint64_t local_in_x[LIMB_Q][RING_DIM];
    #pragma HLS BIND_STORAGE variable=local_in_x type=ram_1wnr impl=bram
    #pragma HLS ARRAY_PARTITION variable=local_in_x type=complete dim=1
    #pragma HLS ARRAY_PARTITION variable=local_in_x type=cyclic factor=LOAD_PAR dim=2

    uint64_t local_out_x[MAX_OUT_COLS][RING_DIM];
    #pragma HLS BIND_STORAGE variable=local_out_x type=ram_s2p impl=bram
    #pragma HLS ARRAY_PARTITION variable=local_out_x type=complete dim=1
    #pragma HLS ARRAY_PARTITION variable=local_out_x type=cyclic factor=LOAD_PAR dim=2

    uint64_t local_w  [LIMB_Q][MAX_OUT_COLS];
    uint64_t local_mod[MAX_OUT_COLS];
    #pragma HLS ARRAY_PARTITION variable=local_w   complete dim=0
    #pragma HLS ARRAY_PARTITION variable=local_mod complete dim=0

    Load_W: for (int q = 0; q < LIMB_Q; ++q) {
        for (int p = 0; p < MAX_OUT_COLS; ++p) {
            #pragma HLS PIPELINE II=1
            local_w[q][p] = in_w[q][p];
        }
    }

    Load_Mod: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS PIPELINE II=1
        local_mod[p] = out_mod[p];
    }

    // Load_X: 单拍突发读 LOAD_PAR 个连续系数（512-bit AXI burst）
    Load_X: for (int l = 0; l < LIMB_Q; ++l) {
        for (int i = 0; i < RING_DIM; i += LOAD_PAR) {
            #pragma HLS PIPELINE II=1
            Load_Burst: for (int pp = 0; pp < LOAD_PAR; ++pp) {
                #pragma HLS UNROLL
                int idx = i + pp;
                local_in_x[l][idx] = in_x[l][idx >> LOG_SQRT][idx & SQRT_MASK];
            }
        }
    }

    bconv_systolic_core(local_in_x, local_out_x, local_w, local_mod, sizeP);

    // Store_X: 单拍突发写回 LOAD_PAR 个连续系数
    Store_X: for (int p = 0; p < sizeP; ++p) {
        #pragma HLS LOOP_TRIPCOUNT min=1 max=MAX_OUT_COLS avg=3
        for (int i = 0; i < RING_DIM; i += LOAD_PAR) {
            #pragma HLS PIPELINE II=1
            Store_Burst: for (int pp = 0; pp < LOAD_PAR; ++pp) {
                #pragma HLS UNROLL
                int idx = i + pp;
                in_x[LIMB_Q + p][idx >> LOG_SQRT][idx & SQRT_MASK] = local_out_x[p][idx];
            }
        }
    }
}
