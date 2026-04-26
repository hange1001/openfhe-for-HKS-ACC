#include "../include/top.h"
#include "../include/load.h"
#include "../include/arithmetic.h"
#include "../include/ntt_kernel.h"
#include "../include/cg_ntt.h"
#include "../include/interleave.h"
#include "../include/mod_mult_kernel.h"
#include "../include/mod_add_kernel.h"
#include "../include/mod_sub_kernel.h"
#include "../include/bconv.h"
#include "../include/bconv_systolic.h"
#include "../include/auto.h"


// -------------------------
// Store the Memory
// -------------------------
static uint64_t poly_buffer_1[MAX_LIMBS][SQRT][SQRT];
static uint64_t poly_buffer_2[MAX_LIMBS][SQRT][SQRT];
static uint64_t result_buffer[MAX_LIMBS][SQRT][SQRT];


// -------------------------
// Store the Modulus
// -------------------------
static uint64_t MODULUS[MAX_LIMBS];
static uint64_t K_HALF[MAX_LIMBS];
static uint64_t M[MAX_LIMBS];

// ------------------------
// Store the CG-NTT TwiddleFactor
// ------------------------
static uint64_t NTTTwiddleFactor[MAX_LIMBS][STAGE][CG_HALF_N];
static uint64_t INTTTwiddleFactor[MAX_LIMBS][STAGE][CG_HALF_N];


// =========================================================================
// 将 OP_INIT 的所有 AXI 裸读操作提取成独立子函数，隔离顶层状态机。
// INLINE off 强制 HLS 生成独立 RTL 模块，AXI arbiter MUX 下沉至子模块内部，
// 顶层只保留 ap_start/ap_done 握手，消除 gmem0/gmem1 的顶层巨型 MUX。
// =========================================================================
static void Load_Init_Params(const uint64_t *mem_in1, const uint64_t *mem_in2) {
    #pragma HLS INLINE off

    static const int CG_TF_SIZE   = STAGE * CG_HALF_N;
    static const int NTT_TF_BASE  = LIMB_Q * 3;
    static const int INTT_TF_BASE = LIMB_P * 3;

    init_Q_MOD:
    for (int i = 0; i < LIMB_Q; i++) {
        #pragma HLS PIPELINE II=1
        MODULUS[i] = mem_in1[i];
    }
    init_Q_KHALF:
    for (int i = 0; i < LIMB_Q; i++) {
        #pragma HLS PIPELINE II=1
        K_HALF[i] = mem_in1[LIMB_Q + i];
    }
    init_Q_M:
    for (int i = 0; i < LIMB_Q; i++) {
        #pragma HLS PIPELINE II=1
        M[i] = mem_in1[LIMB_Q * 2 + i];
    }
    init_P_MOD:
    for (int j = 0; j < LIMB_P; j++) {
        #pragma HLS PIPELINE II=1
        MODULUS[LIMB_Q + j] = mem_in2[j];
    }
    init_P_KHALF:
    for (int j = 0; j < LIMB_P; j++) {
        #pragma HLS PIPELINE II=1
        K_HALF[LIMB_Q + j] = mem_in2[LIMB_P + j];
    }
    init_P_M:
    for (int j = 0; j < LIMB_P; j++) {
        #pragma HLS PIPELINE II=1
        M[LIMB_Q + j] = mem_in2[LIMB_P * 2 + j];

        #ifndef __SYNTHESIS__
        int idx = LIMB_Q + j;
        std::cout << "[FPGA Init] P[" << j << "] (idx=" << idx << "): MOD=" << MODULUS[idx]
                  << ", K=" << K_HALF[idx] << ", M=" << M[idx] << std::endl;
        #endif
    }

    init_NTTTwiddle_Loop:
    for (int l = 0; l < MAX_LIMBS; l++) {
        for (int s = 0; s < STAGE; s++) {
            for (int t = 0; t < CG_HALF_N; t++) {
                #pragma HLS PIPELINE II=1
                NTTTwiddleFactor[l][s][t] =
                    mem_in1[NTT_TF_BASE + l * CG_TF_SIZE + s * CG_HALF_N + t];
            }
        }
    }
    init_INTTTwiddle_Loop:
    for (int l = 0; l < MAX_LIMBS; l++) {
        for (int s = 0; s < STAGE; s++) {
            for (int t = 0; t < CG_HALF_N; t++) {
                #pragma HLS PIPELINE II=1
                INTTTwiddleFactor[l][s][t] =
                    mem_in2[INTT_TF_BASE + l * CG_TF_SIZE + s * CG_HALF_N + t];
            }
        }
    }
}

// =========================================================================
// 将 OP_BCONV 的权重和参数加载提取成独立子函数，隔离顶层状态机。
// =========================================================================
static void Load_BConv_Params(
    const uint64_t *mem_in2,
    uint64_t in_w[LIMB_Q][MAX_OUT_COLS],
    uint64_t out_mod[MAX_OUT_COLS],
    uint64_t out_S[MAX_OUT_COLS],
    uint64_t out_m_barrett[MAX_OUT_COLS]
) {
    #pragma HLS INLINE off

    load_bconv_w:
    for (int q = 0; q < LIMB_Q; q++) {
        for (int p = 0; p < MAX_OUT_COLS; p++) {
            #pragma HLS PIPELINE II=1
            in_w[q][p] = mem_in2[q * MAX_OUT_COLS + p];
        }
    }

    int mod_offset   = LIMB_Q * MAX_OUT_COLS;
    int khalf_offset = mod_offset   + MAX_OUT_COLS;
    int m_offset     = khalf_offset + MAX_OUT_COLS;

    load_bconv_meta:
    for (int p = 0; p < MAX_OUT_COLS; p++) {
        #pragma HLS PIPELINE II=1
        out_mod[p]       = mem_in2[mod_offset   + p];
        out_S[p]         = mem_in2[khalf_offset + p];
        out_m_barrett[p] = mem_in2[m_offset     + p];
    }
}

// =========================================================================
// 将 OP_AUTO 的 meta 读取提取成独立子函数，消除 gmem1 的残余直接访问点。
// =========================================================================
static void Load_Auto_Meta(const uint64_t *mem_in2, uint32_t &k, uint32_t &kinv) {
    #pragma HLS INLINE off
    k    = (uint32_t)mem_in2[0];
    kinv = (uint32_t)mem_in2[1];
}

// =========================================================================
// 提取 NTT/INTT 核心循环（含 flatten/CG_NTT_Kernel/reshape），
// 消除 Top FSM 对 poly_buffer_1 BRAM 的直接访问，降低顶层 BRAM Port MUX。
// =========================================================================
static void Execute_NTT(
    uint64_t poly_buffer[MAX_LIMBS][SQRT][SQRT],
    int num_active_limbs,
    int mod_index
) {
    #pragma HLS INLINE off

    NTT_CG_LIMB_LOOP:
    for (int l = mod_index; l < mod_index + num_active_limbs; l++) {
        #pragma HLS LOOP_TRIPCOUNT min=1 max=5 avg=3
        uint64_t flat[RING_DIM];
        uint64_t flat_out[RING_DIM];
        #pragma HLS ARRAY_PARTITION variable=flat     cyclic factor=PE_PARALLEL dim=1
        #pragma HLS ARRAY_PARTITION variable=flat_out cyclic factor=PE_PARALLEL dim=1
        flatten_2d_to_1d(poly_buffer[l], flat);
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
        #pragma HLS LOOP_TRIPCOUNT min=1 max=5 avg=3
        uint64_t flat[RING_DIM];
        uint64_t flat_out[RING_DIM];
        #pragma HLS ARRAY_PARTITION variable=flat     cyclic factor=PE_PARALLEL dim=1
        #pragma HLS ARRAY_PARTITION variable=flat_out cyclic factor=PE_PARALLEL dim=1
        flatten_2d_to_1d(poly_buffer[l], flat);
        CG_NTT_Kernel<false>(flat, flat_out, MODULUS[l], K_HALF[l], M[l], INTTTwiddleFactor[l]);
        reshape_1d_to_2d(flat_out, poly_buffer[l]);
    }
}


// =========================================================================
// 顶层主函数：只做状态机分发，不直接访问任何 AXI 端口
// =========================================================================
void Top(
    const uint64_t *mem_in1,
    const uint64_t *mem_in2,
    uint64_t *mem_out,
    const uint8_t opcode,
    const int num_active_limbs,
    const int mod_index
){
    #pragma HLS INTERFACE m_axi port=mem_in1  offset=slave bundle=gmem0 depth=196617
    #pragma HLS INTERFACE m_axi port=mem_in2  offset=slave bundle=gmem1 depth=196614
    #pragma HLS INTERFACE m_axi port=mem_out  offset=slave bundle=gmem2 depth=32768

    #pragma HLS INTERFACE s_axilite port=mem_in1          bundle=control
    #pragma HLS INTERFACE s_axilite port=mem_in2          bundle=control
    #pragma HLS INTERFACE s_axilite port=mem_out          bundle=control
    #pragma HLS INTERFACE s_axilite port=opcode           bundle=control
    #pragma HLS INTERFACE s_axilite port=num_active_limbs bundle=control
    #pragma HLS INTERFACE s_axilite port=mod_index        bundle=control
    #pragma HLS INTERFACE s_axilite port=return           bundle=control

    #pragma HLS ARRAY_PARTITION variable=poly_buffer_1 cyclic factor=PE_PARALLEL dim=3
    #pragma HLS ARRAY_PARTITION variable=poly_buffer_2 cyclic factor=PE_PARALLEL dim=3
    #pragma HLS ARRAY_PARTITION variable=result_buffer cyclic factor=PE_PARALLEL dim=3

    #pragma HLS BIND_STORAGE variable=poly_buffer_1 type=ram_2p impl=bram
    #pragma HLS BIND_STORAGE variable=poly_buffer_2 type=ram_2p impl=bram
    #pragma HLS BIND_STORAGE variable=result_buffer type=ram_2p impl=bram

    #pragma HLS ARRAY_PARTITION variable=NTTTwiddleFactor complete dim=2
    #pragma HLS ARRAY_PARTITION variable=NTTTwiddleFactor cyclic factor=PE_PARALLEL dim=3
    #pragma HLS ARRAY_PARTITION variable=INTTTwiddleFactor complete dim=2
    #pragma HLS ARRAY_PARTITION variable=INTTTwiddleFactor cyclic factor=PE_PARALLEL dim=3
    #pragma HLS BIND_STORAGE variable=NTTTwiddleFactor type=ram_2p impl=uram
    #pragma HLS BIND_STORAGE variable=INTTTwiddleFactor type=ram_2p impl=uram

    switch(opcode) {
        case OP_INIT: {
            std::cout << "[FPGA] Initializing Modulus Parameters..." << std::endl;
            Load_Init_Params(mem_in1, mem_in2);
            break;
        }

        case OP_ADD:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Load(mem_in2, poly_buffer_2, num_active_limbs, mod_index);
            Compute_Add(poly_buffer_1, poly_buffer_2, result_buffer, MODULUS, num_active_limbs, mod_index);
            Store(result_buffer, mem_out, num_active_limbs, mod_index);
            break;

        case OP_SUB:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Load(mem_in2, poly_buffer_2, num_active_limbs, mod_index);
            Compute_Sub(poly_buffer_1, poly_buffer_2, result_buffer, MODULUS, num_active_limbs, mod_index);
            Store(result_buffer, mem_out, num_active_limbs, mod_index);
            break;

        case OP_MULT:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Load(mem_in2, poly_buffer_2, num_active_limbs, mod_index);
            Compute_Mult(poly_buffer_1, poly_buffer_2, result_buffer, MODULUS, K_HALF, M, num_active_limbs, mod_index);
            Store(result_buffer, mem_out, num_active_limbs, mod_index);
            break;

        case OP_NTT: {
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Execute_NTT(poly_buffer_1, num_active_limbs, mod_index);
            Store(poly_buffer_1, mem_out, num_active_limbs, mod_index);
            break;
        }

        case OP_INTT: {
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Execute_INTT(poly_buffer_1, num_active_limbs, mod_index);
            Store(poly_buffer_1, mem_out, num_active_limbs, mod_index);
            break;
        }

        case OP_BCONV: {
            int sizeP = num_active_limbs;
            Load(mem_in1, poly_buffer_1, LIMB_Q, 0);

            static uint64_t in_w[LIMB_Q][MAX_OUT_COLS];
            static uint64_t out_mod[MAX_OUT_COLS];
            static uint64_t out_S[MAX_OUT_COLS];
            static uint64_t out_m_barrett[MAX_OUT_COLS];

            Load_BConv_Params(mem_in2, in_w, out_mod, out_S, out_m_barrett);

            #ifndef __SYNTHESIS__
            std::cout << "[BCONV] sizeP=" << sizeP << std::endl;
            for (int p = 0; p < sizeP; p++) {
                std::cout << "  out_mod[" << p << "] = " << out_mod[p]
                          << ", S=" << out_S[p] << ", m_barrett=" << out_m_barrett[p] << std::endl;
            }
            #endif

            Compute_BConv_Systolic(poly_buffer_1, in_w, out_mod, out_S, out_m_barrett, sizeP);
            Store(poly_buffer_1, mem_out, sizeP, LIMB_Q);
            break;
        }

        case OP_AUTO: {
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            uint32_t k, kinv;
            Load_Auto_Meta(mem_in2, k, kinv);
            Compute_Auto(poly_buffer_1, k, kinv, result_buffer, MODULUS, num_active_limbs, mod_index);
            Store(result_buffer, mem_out, num_active_limbs, mod_index);
            break;
        }

        default:
            std::cout << "[FPGA] Unknown opcode: " << opcode << std::endl;
            break;
    }
}
