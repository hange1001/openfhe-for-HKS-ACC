#include "../include/bconv.h"
#include <hls_stream.h>

// ------------------------------------------------
// BConv: Base Conversion using Weight-Stationary Systolic Array
// 
// Computes: out[p][row][col] = sum_{q=0}^{LIMB_Q-1} (in_x[q][row][col] * in_w[q][p]) mod MODULUS[p]
// 
// Architecture:
//   - LIMB_Q rows × LIMB_P columns of PEs
//   - X values flow horizontally (left to right)
//   - Partial sums flow vertically (top to bottom)
//   - Weights are stationary in each PE
//
// Timing:
//   - Row q input is skewed by q cycles
//   - Column p output arrives after (LIMB_Q + p) cycles latency
//   - Total cycles: LIMB_Q + RING_DIM + LIMB_P - 1
// ------------------------------------------------

// Total cycles needed: skew + data + flush
static const int TOTAL_CYCLES = LIMB_Q + RING_DIM + LIMB_P - 1;

// =================================================
// Top Level - Systolic BConv Kernel
// =================================================
void bconv_systolic(
    uint64_t in_x[LIMB_Q + LIMB_P][SQRT][SQRT], 
    const uint64_t in_w[LIMB_Q][LIMB_P], 
    const uint64_t MODULUS[LIMB_Q + LIMB_P],
    int num_active_limbs,
    int mod_idx_offset
) {
#pragma HLS INTERFACE m_axi port=in_x   bundle=gmem0
#pragma HLS INTERFACE m_axi port=in_w   bundle=gmem1
#pragma HLS INTERFACE s_axilite port=MODULUS bundle=control 
#pragma HLS INTERFACE s_axilite port=num_active_limbs bundle=control
#pragma HLS INTERFACE s_axilite port=mod_idx_offset bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    // =========================================
    // Load weights into local buffer
    // =========================================
    ap_uint<64> local_w[LIMB_Q][LIMB_P];
#pragma HLS ARRAY_PARTITION variable=local_w complete dim=0
    
    Load_W: for (int q = 0; q < LIMB_Q; ++q) {
        for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS PIPELINE II=1
            local_w[q][p] = in_w[q][p];
        }
    }

    // =========================================
    // Systolic array registers
    // x_reg[q][p]: x value between PE[q][p-1] and PE[q][p]
    // sum_reg[p][q]: partial sum between PE[q-1][p] and PE[q][p]
    // =========================================
    ap_uint<64> x_reg[LIMB_Q][LIMB_P + 1];
    ap_uint<128> sum_reg[LIMB_P][LIMB_Q + 1];
    
#pragma HLS ARRAY_PARTITION variable=x_reg complete dim=0
#pragma HLS ARRAY_PARTITION variable=sum_reg complete dim=0

    // Initialize registers to zero
    Init_X_Reg: for (int q = 0; q < LIMB_Q; ++q) {
#pragma HLS UNROLL
        for (int p = 0; p <= LIMB_P; ++p) {
#pragma HLS UNROLL
            x_reg[q][p] = 0;
        }
    }
    
    Init_Sum_Reg: for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS UNROLL
        for (int q = 0; q <= LIMB_Q; ++q) {
#pragma HLS UNROLL
            sum_reg[p][q] = 0;
        }
    }

    // Output collection counters
    int valid_count[LIMB_P];
#pragma HLS ARRAY_PARTITION variable=valid_count complete

    Init_Count: for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS UNROLL
        valid_count[p] = 0;
    }

    // =========================================
    // Main systolic loop
    // =========================================
    Systolic_Loop: for (int t = 0; t < TOTAL_CYCLES; ++t) {
#pragma HLS PIPELINE II=1
        
        // -----------------------------------------
        // Stage 1: Save current register values
        // (needed for correct systolic timing)
        // -----------------------------------------
        ap_uint<64> x_curr[LIMB_Q][LIMB_P + 1];
        ap_uint<128> sum_curr[LIMB_P][LIMB_Q + 1];
#pragma HLS ARRAY_PARTITION variable=x_curr complete dim=0
#pragma HLS ARRAY_PARTITION variable=sum_curr complete dim=0
        
        Save_X: for (int q = 0; q < LIMB_Q; ++q) {
#pragma HLS UNROLL
            for (int p = 0; p <= LIMB_P; ++p) {
#pragma HLS UNROLL
                x_curr[q][p] = x_reg[q][p];
            }
        }
        
        Save_Sum: for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS UNROLL
            for (int q = 0; q <= LIMB_Q; ++q) {
#pragma HLS UNROLL
                sum_curr[p][q] = sum_reg[p][q];
            }
        }
        
        // -----------------------------------------
        // Stage 2: Feed new inputs (with skewing)
        // Row q starts at cycle q
        // -----------------------------------------
        Feed_X: for (int q = 0; q < LIMB_Q; ++q) {
#pragma HLS UNROLL
            int data_idx = t - q;
            
            if (data_idx >= 0 && data_idx < RING_DIM) {
                // 使用row-major顺序匹配Load函数
                int row = data_idx / SQRT;
                int col = data_idx % SQRT;
                x_reg[q][0] = in_x[q][row][col];
            } else {
                x_reg[q][0] = 0;
            }
        }
        
        Init_Sum: for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS UNROLL
            sum_reg[p][0] = 0;
        }
        
        // -----------------------------------------
        // Stage 3: PE Array Computation
        // Read from curr (previous cycle values)
        // Write to reg (for next cycle and output)
        // -----------------------------------------
        PE_Row: for (int q = 0; q < LIMB_Q; ++q) {
#pragma HLS UNROLL
            PE_Col: for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS UNROLL
                // Read from previous cycle's values
                ap_uint<64> x_in = x_curr[q][p];
                ap_uint<128> sum_in = sum_curr[p][q];
                
                // MAC operation
                ap_uint<128> prod = ((ap_uint<128>)x_in * (ap_uint<128>)local_w[q][p]) % MODULUS[LIMB_Q + p];
                ap_uint<128> sum_out = (sum_in + prod) % MODULUS[LIMB_Q + p];
                
                // Write outputs (for next cycle)
                x_reg[q][p + 1] = x_in;        // Pass x to right
                sum_reg[p][q + 1] = sum_out;   // Pass sum to bottom
            }
        }
        
        // -----------------------------------------
        // Stage 4: Collect outputs (with de-skewing)
        // Column p output is ready after LIMB_Q + p cycles
        // -----------------------------------------
        Collect: for (int p = 0; p < LIMB_P; ++p) {
#pragma HLS UNROLL
            int latency = LIMB_Q + p;
            
            if (t >= latency && valid_count[p] < RING_DIM) {
                // Read result from bottom of column p
                ap_uint<128> result = sum_reg[p][LIMB_Q];
                // 使用row-major顺序匹配Store函数
                int row = valid_count[p] / SQRT;
                int col = valid_count[p] % SQRT;
                in_x[LIMB_Q + p][row][col] = (uint64_t)(result);
                
                valid_count[p]++;
            }
        }
    }
}

// =================================================
// Host Side Wrapper
// =================================================
void Compute_BConv(
    uint64_t in_x[LIMB_Q + LIMB_P][SQRT][SQRT], 
    const uint64_t in_w[LIMB_Q][LIMB_P], 
    const uint64_t MODULUS[LIMB_Q + LIMB_P],
    int num_active_limbs,
    int mod_idx_offset
) {
    bconv_systolic(
        in_x,
        in_w,
        MODULUS,
        num_active_limbs,
        mod_idx_offset
    );
}
