#include "../include/top.h"
#include "../include/load.h"
#include "../include/arithmetic.h"
#include "../include/ntt_kernel.h"
#include "../include/interleave.h"
#include "../include/mod_mult_kernel.h"
#include "../include/mod_add_kernel.h"
#include "../include/mod_sub_kernel.h"
#include "../include/bconv.h"



// -------------------------
// Store the Memory--TODO：testing
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
// Store the TwiddleFactor
// ------------------------
static uint64_t NTTTwiddleFactor[MAX_LIMBS][BU_NUM][RING_DIM];
static uint64_t INTTTwiddleFactor[MAX_LIMBS][BU_NUM][RING_DIM];

void Top(
    const uint64_t *mem_in1,
    const uint64_t *mem_in2,
    uint64_t *mem_out,
    const uint8_t opcode,
    const int num_active_limbs,
    const int mod_index
){

    #pragma HLS INTERFACE m_axi port=mem_in1  offset=slave bundle=gmem0
    #pragma HLS INTERFACE m_axi port=mem_in2  offset=slave bundle=gmem1
    #pragma HLS INTERFACE m_axi port=mem_out  offset=slave bundle=gmem0

    #pragma HLS INTERFACE s_axilite port=mem_in1  bundle=control
    #pragma HLS INTERFACE s_axilite port=mem_in2  bundle=control
    #pragma HLS INTERFACE s_axilite port=mem_out  bundle=control
    #pragma HLS INTERFACE s_axilite port=opcode    bundle=control
    #pragma HLS INTERFACE s_axilite port=num_active_limbs bundle=control
    #pragma HLS INTERFACE s_axilite port=mod_index bundle=control
    #pragma HLS INTERFACE s_axilite port=return    bundle=control


    #pragma HLS ARRAY_PARTITION variable=poly_buffer_1 cyclic dim=2 factor=SQRT
    #pragma HLS ARRAY_PARTITION variable=poly_buffer_2 cyclic dim=2 factor=SQRT
    #pragma HLS ARRAY_PARTITION variable=result_buffer cyclic dim=2 factor=SQRT

    #pragma HLS BIND_STORAGE variable=poly_buffer_1 type=ram_2p impl=bram
    #pragma HLS BIND_STORAGE variable=poly_buffer_2 type=ram_2p impl=bram
    #pragma HLS BIND_STORAGE variable=result_buffer type=ram_2p impl=bram

    switch(opcode) {
        case OP_INIT: {
            std::cout << "[FPGA] Initializing Modulus Parameters..." << std::endl;
            
            // 简单布局：Q模数在索引0,1,2，P模数在索引3,4
            // 无padding，与Host端一致
            
            init_Q_Loop:
            for (int i = 0; i < LIMB_Q; i++){
                MODULUS[i] = mem_in1[i];
                K_HALF[i] = mem_in1[LIMB_Q + i];
                M[i] = mem_in1[LIMB_Q*2 + i];
                
                #ifndef __SYNTHESIS__
                std::cout << "[FPGA Init] Q[" << i << "]: MOD=" << MODULUS[i] 
                          << ", K=" << K_HALF[i] << ", M=" << M[i] << std::endl;
                #endif
            }
            init_P_Loop:
            for (int j = 0; j < LIMB_P; j++){
                // P模数从索引LIMB_Q开始，即索引3,4
                int idx = LIMB_Q + j;
                MODULUS[idx] = mem_in2[j];
                K_HALF[idx] = mem_in2[LIMB_P + j];
                M[idx] = mem_in2[LIMB_P*2 + j];
                
                #ifndef __SYNTHESIS__
                std::cout << "[FPGA Init] P[" << j << "] (idx=" << idx << "): MOD=" << MODULUS[idx] 
                          << ", K=" << K_HALF[idx] << ", M=" << M[idx] << std::endl;
                #endif
            }
            init_NTTTwiddle_Loop:
            for (int l = 0; l < LIMB_Q + LIMB_P; l++){
                int offset = LIMB_Q * 3;
                for (int b = 0; b < BU_NUM; b++){
                    for (int t = 0; t < RING_DIM; t++){
                        NTTTwiddleFactor[l][b][t] = mem_in1[offset + l*RING_DIM + t];
                    }
                }
            }
            init_INTTTwiddle_Loop:
            for (int l = 0; l < LIMB_Q + LIMB_P; l++){
                int offset = LIMB_P * 3;
                for (int b = 0; b < BU_NUM; b++){
                    for (int t = 0; t < RING_DIM; t++){
                        INTTTwiddleFactor[l][b][t] = mem_in2[offset + l*RING_DIM + t];
                    }
                }
            }
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

        case OP_MUL:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Load(mem_in2, poly_buffer_2, num_active_limbs, mod_index);
            Compute_Mult(poly_buffer_1, poly_buffer_2, result_buffer, MODULUS, K_HALF, M, num_active_limbs, mod_index);
            Store(result_buffer, mem_out, num_active_limbs, mod_index);
            break;


        case OP_NTT:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            for (int l = 0; l < MAX_LIMBS; l++){
                InterLeave(poly_buffer_1[l], true);
            }
            Compute_NTT(poly_buffer_1, NTTTwiddleFactor, INTTTwiddleFactor, MODULUS, K_HALF, M, true, num_active_limbs, mod_index);
            Store(poly_buffer_1, mem_out, num_active_limbs, mod_index);
            break;

        case OP_INTT:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            Compute_NTT(poly_buffer_1, NTTTwiddleFactor, INTTTwiddleFactor, MODULUS, K_HALF, M, false, num_active_limbs, mod_index);
            for (int l = 0; l < MAX_LIMBS; l++){
                InterLeave(poly_buffer_1[l], false);
            }
            Store(poly_buffer_1, mem_out, num_active_limbs, mod_index);
            break;

    
        case OP_BCONV:
            Load(mem_in1, poly_buffer_1, num_active_limbs, mod_index);
            static uint64_t in_w[LIMB_Q][LIMB_P];
            for (int q = 0; q < LIMB_Q; q++){
                for (int p = 0; p < LIMB_P; p++){
                    in_w[q][p] = mem_in2[q*LIMB_P + p];
                }
            }
            Compute_BConv(poly_buffer_1, in_w, MODULUS, num_active_limbs, mod_index);
            Store(poly_buffer_1, mem_out, num_active_limbs, mod_index);
            for (int l = 0; l < MAX_LIMBS; l++){
                for (int row = 0; row < SQRT; row++){
                    for (int col = 0; col < SQRT; col++){
                        std::cout << "poly_buffer_1[" << l << "][" << row << "][" << col << "] = " << poly_buffer_1[l][row][col] << std::endl;
                    }
                }
            }
            break;

        default:
            std::cout << "[FPGA] Unknown opcode: " << opcode << std::endl;
            break;
    }
}