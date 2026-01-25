#include "../include/bconv.h"
#include <hls_stream.h>

// ------------------------------------------------
// BConv: Base Conversion using Weight-Stationary Systolic Array
// 
// Computes: out[p][ri] = sum_{q=0}^{sizeQ-1} (in_x[q][ri] * in_w[q][p]) mod out_mod[p]
// 
// Architecture:
//   - LIMB_Q rows × MAX_OUT_COLS columns of PEs (3×5)
//   - X values flow horizontally (left to right)
//   - Partial sums flow vertically (top to bottom)
//   - Weights are stationary in each PE
//   - sizeP controls how many output columns are active
//
// Timing:
//   - Row q input is skewed by q cycles
//   - Column p output arrives after (LIMB_Q + p) cycles latency
//   - Total cycles: LIMB_Q + RING_DIM + MAX_OUT_COLS - 1
// ------------------------------------------------

// Total cycles needed: skew + data + flush (use max columns)
static const int TOTAL_CYCLES = LIMB_Q + RING_DIM + MAX_OUT_COLS - 1;

// =================================================
// Top Level - Systolic BConv Kernel
// =================================================
void bconv_systolic(
    uint64_t in_x[MAX_LIMBS][SQRT][SQRT],           // Input[0..LIMB_Q-1], Output[LIMB_Q..LIMB_Q+sizeP-1]
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],      // Weights [sizeQ × sizeP]
    const uint64_t out_mod[MAX_OUT_COLS],           // Output moduli for each column
    int sizeP                                        // Number of active output columns
) {
#pragma HLS INTERFACE m_axi port=in_x   bundle=gmem0
#pragma HLS INTERFACE m_axi port=in_w   bundle=gmem1
#pragma HLS INTERFACE s_axilite port=out_mod bundle=control 
#pragma HLS INTERFACE s_axilite port=sizeP bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    // =========================================
    // Load weights into local buffer
    // =========================================
    ap_uint<64> local_w[LIMB_Q][MAX_OUT_COLS];
#pragma HLS ARRAY_PARTITION variable=local_w complete dim=0
    
    Load_W: for (int q = 0; q < LIMB_Q; ++q) {
        for (int p = 0; p < MAX_OUT_COLS; ++p) {
#pragma HLS PIPELINE II=1
            local_w[q][p] = in_w[q][p];
        }
    }

    // Load output moduli
    ap_uint<64> local_mod[MAX_OUT_COLS];
#pragma HLS ARRAY_PARTITION variable=local_mod complete
    
    Load_Mod: for (int p = 0; p < MAX_OUT_COLS; ++p) {
#pragma HLS PIPELINE II=1
        local_mod[p] = out_mod[p];
    }

    // =========================================
    // Systolic array registers
    // x_reg[q][p]: x value between PE[q][p-1] and PE[q][p]
    // sum_reg[p][q]: partial sum between PE[q-1][p] and PE[q][p]
    // =========================================
    ap_uint<64> x_reg[LIMB_Q][MAX_OUT_COLS + 1];
    ap_uint<128> sum_reg[MAX_OUT_COLS][LIMB_Q + 1];
    
#pragma HLS ARRAY_PARTITION variable=x_reg complete dim=0
#pragma HLS ARRAY_PARTITION variable=sum_reg complete dim=0

    // Initialize registers to zero
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

    // Output collection counters
    int valid_count[MAX_OUT_COLS];
#pragma HLS ARRAY_PARTITION variable=valid_count complete

    Init_Count: for (int p = 0; p < MAX_OUT_COLS; ++p) {
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
        ap_uint<64> x_curr[LIMB_Q][MAX_OUT_COLS + 1];
        ap_uint<128> sum_curr[MAX_OUT_COLS][LIMB_Q + 1];
#pragma HLS ARRAY_PARTITION variable=x_curr complete dim=0
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
        
        // -----------------------------------------
        // Stage 2: Feed new inputs (with skewing)
        // Row q starts at cycle q
        // -----------------------------------------
        Feed_X: for (int q = 0; q < LIMB_Q; ++q) {
#pragma HLS UNROLL
            int data_idx = t - q;
            
            if (data_idx >= 0 && data_idx < RING_DIM) {
                int row = data_idx / SQRT;
                int col = data_idx % SQRT;
                x_reg[q][0] = in_x[q][row][col];
            } else {
                x_reg[q][0] = 0;
            }
        }
        
        Init_Sum: for (int p = 0; p < MAX_OUT_COLS; ++p) {
#pragma HLS UNROLL
            sum_reg[p][0] = 0;
        }
        
        // -----------------------------------------
        // Stage 3: PE Array Computation (all columns)
        // Read from curr (previous cycle values)
        // Write to reg (for next cycle and output)
        // -----------------------------------------
        PE_Row: for (int q = 0; q < LIMB_Q; ++q) {
#pragma HLS UNROLL
            PE_Col: for (int p = 0; p < MAX_OUT_COLS; ++p) {
#pragma HLS UNROLL
                // Read from previous cycle's values
                ap_uint<64> x_in = x_curr[q][p];
                ap_uint<128> sum_in = sum_curr[p][q];
                
                // MAC operation using column-specific modulus
                ap_uint<64> mod_p = local_mod[p];
                ap_uint<128> prod = ((ap_uint<128>)x_in * (ap_uint<128>)local_w[q][p]) % mod_p;
                ap_uint<128> sum_out = (sum_in + prod) % mod_p;
                
                // Write outputs (for next cycle)
                x_reg[q][p + 1] = x_in;        // Pass x to right
                sum_reg[p][q + 1] = sum_out;   // Pass sum to bottom
            }
        }
        
        // -----------------------------------------
        // Stage 4: Collect outputs (with de-skewing)
        // Column p output is ready after LIMB_Q + p cycles
        // Only collect for active columns (p < sizeP)
        // -----------------------------------------
        Collect: for (int p = 0; p < MAX_OUT_COLS; ++p) {
#pragma HLS UNROLL
            int latency = LIMB_Q + p;
            
            // Only collect if this column is active
            if (p < sizeP && t >= latency && valid_count[p] < RING_DIM) {
                ap_uint<128> result = sum_reg[p][LIMB_Q];
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
    uint64_t in_x[MAX_LIMBS][SQRT][SQRT], 
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS], 
    const uint64_t out_mod[MAX_OUT_COLS],
    int sizeP
) {
    bconv_systolic(in_x, in_w, out_mod, sizeP);
}
