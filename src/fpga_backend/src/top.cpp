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
#include "../include/hks_digit.h"


// -------------------------
// Store the Memory
// -------------------------
static uint64_t poly_buffer_1[MAX_LIMBS][SQRT][SQRT];
static uint64_t poly_buffer_2[MAX_LIMBS][SQRT][SQRT];
static uint64_t result_buffer[MAX_LIMBS][SQRT][SQRT];

// Flattened coefficient indices are decoded as wiring, not integer division.
static_assert(SQRT == (1 << LOG_SQRT),
              "SQRT must match the power-of-two address width");


// -------------------------
// Store the Modulus
// -------------------------
static uint64_t MODULUS[MAX_LIMBS];
static uint64_t S[MAX_LIMBS];
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
    init_Q_S:
    for (int i = 0; i < LIMB_Q; i++) {
        #pragma HLS PIPELINE II=1
        S[i] = mem_in1[LIMB_Q + i];
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
    init_P_S:
    for (int j = 0; j < LIMB_P; j++) {
        #pragma HLS PIPELINE II=1
        S[LIMB_Q + j] = mem_in2[LIMB_P + j];
    }
    init_P_M:
    for (int j = 0; j < LIMB_P; j++) {
        #pragma HLS PIPELINE II=1
        M[LIMB_Q + j] = mem_in2[LIMB_P * 2 + j];

        #ifndef __SYNTHESIS__
        int idx = LIMB_Q + j;
        std::cout << "[FPGA Init] P[" << j << "] (idx=" << idx << "): MOD=" << MODULUS[idx]
                  << ", S=" << S[idx] << ", M=" << M[idx] << std::endl;
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
    int S_offset = mod_offset   + MAX_OUT_COLS;
    int m_offset     = S_offset + MAX_OUT_COLS;

    load_bconv_meta:
    for (int p = 0; p < MAX_OUT_COLS; p++) {
        #pragma HLS PIPELINE II=1
        out_mod[p]       = mem_in2[mod_offset   + p];
        out_S[p]         = mem_in2[S_offset + p];
        out_m_barrett[p] = mem_in2[m_offset     + p];
    }
}

// =========================================================================
// One runtime call site, shared work A bank and exactly one scratch tower.
// Physical slot and modulus index are distinct (compact digit vs global basis).
// =========================================================================
static void Execute_Transform(int limb, bool is_ntt, int source_limb,
                              bool scale_only, uint64_t scale_factor) {
    #pragma HLS INLINE off
    uint64_t scratch[RING_DIM];
    #pragma HLS ARRAY_PARTITION variable=scratch cyclic factor=CG_BUF_PARTITION dim=1
    #pragma HLS BIND_STORAGE variable=scratch type=ram_t2p impl=bram
    CG_Transform_Work(poly_buffer_1, source_limb, scratch, MODULUS[limb], S[limb], M[limb],
                     NTTTwiddleFactor[limb], INTTTwiddleFactor[limb], is_ntt,
                     scale_only, scale_factor);
}


// Preserve EVAL bypass values during the original AXI load, not in another pass.
static void Load_Transform_Tower(const uint64_t* input, int src, int global_limb,
                                 bool preserve_eval) {
    #pragma HLS INLINE off
    LOAD_TRANSFORM_WORD: for (int k = 0; k < RING_DIM; ++k) {
        #pragma HLS PIPELINE II=1
        #pragma HLS UNROLL factor=PE_PARALLEL
        const uint64_t value = input[src * RING_DIM + k];
        const int row = k >> LOG_SQRT;
        const int col = k & (SQRT - 1);
        poly_buffer_1[src][row][col] = value;
        if (preserve_eval) result_buffer[global_limb][row][col] = value;
    }
}

static void Load_Transform_Input(const uint64_t* input, int count, int start,
                                  bool preserve_eval) {
    #pragma HLS INLINE off
    LOAD_TRANSFORM_LIMB: for (int l = 0; l < count; ++l) {
        Load_Transform_Tower(input, l, start + l, preserve_eval);
    }
}

// Scatter directly from physical work slots to the requested output basis.
// Both input towers and bypass copies are already resident before any AXI write,
// preserving the existing input/output alias contract.
static void Store_Transform_Tower(uint64_t* output, int dst, int slot, bool bypass) {
    #pragma HLS INLINE off
    STORE_TRANSFORM_WORD: for (int k = 0; k < RING_DIM; ++k) {
        #pragma HLS PIPELINE II=1
        #pragma HLS UNROLL factor=PE_PARALLEL
        const int row = k >> LOG_SQRT;
        const int col = k & (SQRT - 1);
        output[dst * RING_DIM + k] = bypass ? result_buffer[dst][row][col]
                                            : poly_buffer_1[slot][row][col];
    }
}

static void Store_Transform_Output(uint64_t* output, bool fused, int alpha, int start) {
    #pragma HLS INLINE off
    const int count = fused ? MAX_OUT_COLS : alpha;
    STORE_TRANSFORM_LIMB: for (int l = 0; l < count; ++l) {
        const bool bypass = fused && l >= start && l < start + alpha;
        const int p = l < start ? l : l - alpha;
        const int slot = fused ? (bypass ? 0 : LIMB_Q + p) : l;
        Store_Transform_Tower(output, l, slot, bypass);
    }
}

// Keep metadata AXI traffic inside a non-inlined helper, like OP_BCONV.
static void Load_HKS_Params(
    const uint64_t* meta, int alpha, int digit_start,
    uint64_t weights[LIMB_Q][MAX_OUT_COLS], uint64_t inv[LIMB_Q],
    uint64_t out_mod[MAX_OUT_COLS], uint64_t out_s[MAX_OUT_COLS],
    uint64_t out_m[MAX_OUT_COLS]
) {
    #pragma HLS INLINE off
    const int count = MAX_OUT_COLS - alpha;
    HKS_LOAD_WEIGHTS: for (int q = 0; q < LIMB_Q; ++q) {
        for (int p = 0; p < MAX_OUT_COLS; ++p) {
            #pragma HLS PIPELINE II=1
            weights[q][p] = (q < alpha && p < count)
                ? meta[q * MAX_OUT_COLS + p] : 0;
        }
    }
    HKS_LOAD_INV: for (int q = 0; q < LIMB_Q; ++q) {
        #pragma HLS PIPELINE II=1
        inv[q] = q < alpha ? meta[HKS_INV_OFFSET + q] : 0;
    }
    HKS_OUTPUT_PARAMS: for (int p = 0; p < MAX_OUT_COLS; ++p) {
        #pragma HLS PIPELINE II=1
        if (p < count) {
            const int global_limb = p < digit_start ? p : p + alpha;
            out_mod[p] = MODULUS[global_limb];
            out_s[p] = S[global_limb];
            out_m[p] = M[global_limb];
        } else {
            out_mod[p] = 0;
            out_s[p] = 0;
            out_m[p] = 0;
        }
    }
}

static void Execute_Transform_Operation(
    const uint64_t* mem_in1, const uint64_t* mem_in2, uint64_t* mem_out,
    uint8_t opcode, int alpha, int digit_start
) {
    #pragma HLS INLINE off
#ifdef HKS_DISABLE_FUSED_OPCODE
    const bool fused = false;
#else
    const bool fused = opcode == OP_HKS_DIGIT;
#endif
    // Reject before any AXI access or state mutation; avoid signed overflow.
    if (fused) {
        if (alpha < 1 || alpha > LIMB_Q || digit_start < 0 ||
            digit_start > LIMB_Q - alpha) return;
        HKS_VALIDATE_MODULI: for (int l = 0; l < MAX_OUT_COLS; ++l) {
            if (MODULUS[l] <= 1 || !(MODULUS[l] & 1) ||
                MODULUS[l] >= (uint64_t(1) << 62)) return;
        }
    } else if (alpha < 1 || alpha > MAX_LIMBS || digit_start < 0 ||
               digit_start > MAX_LIMBS - alpha) {
        return;
    }

    uint64_t weights[LIMB_Q][MAX_OUT_COLS];
    uint64_t inv[LIMB_Q];
    uint64_t out_mod[MAX_OUT_COLS], out_s[MAX_OUT_COLS], out_m[MAX_OUT_COLS];
    #pragma HLS ARRAY_PARTITION variable=weights complete dim=0
    #pragma HLS ARRAY_PARTITION variable=inv complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_mod complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_s complete dim=0
    #pragma HLS ARRAY_PARTITION variable=out_m complete dim=0

    if (fused)
        Load_HKS_Params(mem_in2, alpha, digit_start, weights, inv, out_mod, out_s, out_m);
    Load_Transform_Input(mem_in1, alpha, digit_start, fused);

    // Fused schedule per input tower: INTT then in-place SCALE, followed by
    // one BConv and all complement NTTs. SCALE and butterflies enter through
    // the same Execute_Transform call site and share its four MultMod lanes.
    const int steps = fused ? MAX_OUT_COLS + alpha : alpha;
    TRANSFORM_STEPS: for (int step = 0; step < steps; ++step) {
        #pragma HLS LOOP_TRIPCOUNT min=1 max=8 avg=5
        #pragma HLS LOOP_FLATTEN off
        const bool input_phase = fused && step < 2 * alpha;
        const bool scale_only = input_phase && (step & 1);
        const bool complement = fused && !input_phase;
        const bool is_ntt = fused ? complement : opcode == OP_NTT;
        const int q = input_phase ? (step >> 1) : step;
        int limb = digit_start + q;
        int source_limb = q;
        uint64_t scale_factor = scale_only ? inv[q] : 0;
        if (complement) {
            if (step == 2 * alpha) {
                Compute_BConv_Systolic(poly_buffer_1, weights, out_mod, out_s, out_m,
                                      MAX_OUT_COLS - alpha, alpha);
            }
            const int p = step - 2 * alpha;
            limb = p < digit_start ? p : p + alpha;
            source_limb = LIMB_Q + p;
        }
        Execute_Transform(limb, is_ntt, source_limb, scale_only, scale_factor);
    }
    Store_Transform_Output(mem_out, fused, alpha, digit_start);
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
    #pragma HLS INTERFACE m_axi port=mem_in1  offset=slave bundle=gmem0 depth=196617 max_widen_bitwidth=256
    #pragma HLS INTERFACE m_axi port=mem_in2  offset=slave bundle=gmem1 depth=196614 max_widen_bitwidth=256
    #pragma HLS INTERFACE m_axi port=mem_out  offset=slave bundle=gmem2 depth=32768 max_widen_bitwidth=256

    #pragma HLS INTERFACE s_axilite port=mem_in1          bundle=control
    #pragma HLS INTERFACE s_axilite port=mem_in2          bundle=control
    #pragma HLS INTERFACE s_axilite port=mem_out          bundle=control
    #pragma HLS INTERFACE s_axilite port=opcode           bundle=control
    #pragma HLS INTERFACE s_axilite port=num_active_limbs bundle=control
    #pragma HLS INTERFACE s_axilite port=mod_index        bundle=control
    #pragma HLS INTERFACE s_axilite port=return           bundle=control

    #pragma HLS ARRAY_PARTITION variable=poly_buffer_1 complete dim=1
    #pragma HLS ARRAY_PARTITION variable=poly_buffer_1 cyclic factor=CG_BUF_PARTITION dim=3
    #pragma HLS ARRAY_PARTITION variable=poly_buffer_2 cyclic factor=PE_PARALLEL dim=3
    #pragma HLS ARRAY_PARTITION variable=result_buffer cyclic factor=PE_PARALLEL dim=3

    #pragma HLS BIND_STORAGE variable=poly_buffer_1 type=ram_t2p impl=bram
    #pragma HLS BIND_STORAGE variable=poly_buffer_2 type=ram_2p impl=bram
    #pragma HLS BIND_STORAGE variable=result_buffer type=ram_2p impl=bram

    #pragma HLS ARRAY_PARTITION variable=NTTTwiddleFactor complete dim=2
    #pragma HLS ARRAY_PARTITION variable=NTTTwiddleFactor cyclic factor=PE_PARALLEL dim=3
    #pragma HLS ARRAY_PARTITION variable=INTTTwiddleFactor complete dim=2
    #pragma HLS ARRAY_PARTITION variable=INTTTwiddleFactor cyclic factor=PE_PARALLEL dim=3
    #pragma HLS BIND_STORAGE variable=NTTTwiddleFactor type=ram_2p impl=uram
    #pragma HLS BIND_STORAGE variable=INTTTwiddleFactor type=ram_2p impl=uram

    switch(opcode) {
#ifndef HKS_DISABLE_FUSED_OPCODE
        case OP_HKS_DIGIT:
#endif
        case OP_NTT:
        case OP_INTT:
            Execute_Transform_Operation(mem_in1, mem_in2, mem_out, opcode,
                                        num_active_limbs, mod_index);
            break;

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
            Compute_Mult(poly_buffer_1, poly_buffer_2, result_buffer, MODULUS, S, M, num_active_limbs, mod_index);
            Store(result_buffer, mem_out, num_active_limbs, mod_index);
            break;

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

            Compute_BConv_Systolic(poly_buffer_1, in_w, out_mod, out_S, out_m_barrett, sizeP, LIMB_Q);
            Store(poly_buffer_1, mem_out, sizeP, LIMB_Q);
            break;
        }

        default:
            // Reserved/unknown opcodes, including retired value 7, access no memory.
            std::cout << "[FPGA] Unknown opcode: " << unsigned(opcode) << std::endl;
            break;
    }
}
