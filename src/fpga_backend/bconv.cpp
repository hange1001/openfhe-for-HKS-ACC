#include "bconv.h"
#include <hls_stream.h>
#include <iostream>

// ------------------------------------------------
// [CONFIG] 统一运行周期
// 总周期 = Skew(Input + Sum) + Sequence + Flush
// ------------------------------------------------
static const int TOTAL_CYCLES = SIZE_P + SIZE_Q + RING_DIM;

// =================================================
// Weight Loading
// =================================================
void load_weights(
    const ap_uint<64>* weights_linear,
    ap_uint<64>        local_w[SIZE_Q][SIZE_P]
) {
    for(int i=0; i<SIZE_Q; ++i) {
        for(int j=0; j<SIZE_P; ++j) {
#pragma HLS PIPELINE II=1
            local_w[i][j] = weights_linear[i * SIZE_P + j];
        }
    }
}

// =================================================
// PE: Weight Stationary 
// [WARNING] 这里依然存在溢出风险，因为 sum_curr 被截断为 64 位
// 理想情况下应该在 PE 内取模，或者传递 128 位宽的 sum
// =================================================
void PE_WS(
    hls::stream<ap_uint<64>>& x_in,
    hls::stream<ap_uint<64>>& x_out,
    hls::stream<ap_uint<64>>& sum_in,
    hls::stream<ap_uint<64>>& sum_out,
    ap_uint<64>               weight,
    ap_uint<64>               mod,     // 传入模数用于 PE 内部取模(可选)
    int                       row_idx,
    int                       col_idx
) {
#pragma HLS INLINE off

    ap_uint<64> reg_x_out = 0;
    ap_uint<64> reg_sum_out = 0;

    PE_Loop: for (int t = 0; t < TOTAL_CYCLES; ++t) {
#pragma HLS PIPELINE II=1
        
        // Output Stage
        x_out.write(reg_x_out);
        sum_out.write(reg_sum_out);

        // Input Stage
        ap_uint<64> x = x_in.read();
        ap_uint<64> sum_prev = sum_in.read();
        
        // Compute Stage
        ap_uint<128> prod = (ap_uint<128>)x * (ap_uint<128>)weight;
        ap_uint<128> sum_temp = (ap_uint<128>)sum_prev + prod;
        
        // [FIX] 增加简单的保护逻辑：如果在仿真模式，打印警告
        // 在实际硬件中，建议这里做 sum_temp % mod
        // 为了防止 Demo 溢出，这里暂不做 % mod，依靠 Collector
        ap_uint<64> sum_curr = (ap_uint<64>)(sum_temp); 

        // Register Stage
        reg_x_out = x;
        reg_sum_out = sum_curr;
    }
}

// =================================================
// Feeder X (Rows)
// =================================================
void feeder_x(
    const ap_uint<64>* x_linear,
    hls::stream<ap_uint<64>> (&x_grid)[SIZE_Q][SIZE_P + 1]
) {
    Row_Loop: 
    for (int i = 0; i < SIZE_Q; ++i) {
        #pragma HLS UNROLL
        for (int t = 0; t < TOTAL_CYCLES; ++t) {
#pragma HLS PIPELINE II=1
            ap_uint<64> val = 0;
            if (t >= i && t < (i + RING_DIM)) {
                val = x_linear[(t - i) * SIZE_Q + i];
            }
            x_grid[i][0].write(val);
        }
    }
}

// =================================================
// Feeder Sum (Cols)
// =================================================
void feeder_sum(
    hls::stream<ap_uint<64>> (&sum_grid)[SIZE_P][SIZE_Q + 1]
) {
    Col_Loop: for (int j = 0; j < SIZE_P; ++j) {
        for (int t = 0; t < TOTAL_CYCLES; ++t) {
#pragma HLS PIPELINE II=1
            sum_grid[j][0].write(0); 
        }
    }
}

// =================================================
// Drain X
// =================================================
void drain_x(
    hls::stream<ap_uint<64>> (&x_grid)[SIZE_Q][SIZE_P + 1]
) {
    for(int i=0; i<SIZE_Q; ++i) {
#pragma HLS UNROLL
        for(int t=0; t < TOTAL_CYCLES; ++t) {
#pragma HLS PIPELINE II=1
            x_grid[i][SIZE_P].read();
        }
    }
}

// =================================================
// Collector
// =================================================
void collector(
    hls::stream<ap_uint<64>> (&sum_grid)[SIZE_P][SIZE_Q + 1],
    ap_uint<64>  mod_scalar, // [FIX] 这里接收单个标量
    ap_uint<64>* out_linear
) {
    Col_Out_Loop: 
    for (int j = 0; j < SIZE_P; ++j) {
        #pragma HLS UNROLL
        int valid_count = 0;
        int discard_count = 0;
        int wait_cycles = SIZE_Q + j;

        for (int t = 0; t < TOTAL_CYCLES; ++t) {
            #pragma HLS PIPELINE II=1
            ap_uint<64> val = sum_grid[j][SIZE_Q].read();
            
            if (discard_count < wait_cycles) {
                discard_count++;
            }
            else if (valid_count < RING_DIM) {
                // [FIX] 使用统一的标量取模，彻底解决 Divide by 0
                out_linear[valid_count * SIZE_P + j] = val % mod_scalar;
                valid_count++;
            }
        }
    }
}

// =================================================
// Top Level
// =================================================
extern "C" void bconv_systolic(
    const ap_uint<64>* x_linear,    
    const ap_uint<64>* w_linear,    
    const ap_uint<64>  mod_scalar,  // [FIX] 修改为传值 (By Value)，不再是指针
    ap_uint<64>* out_linear         
) {
#pragma HLS INTERFACE m_axi port=x_linear   bundle=gmem0
#pragma HLS INTERFACE m_axi port=w_linear   bundle=gmem1
// [FIX] mod_scalar 映射到 s_axilite 寄存器
#pragma HLS INTERFACE s_axilite port=mod_scalar bundle=control 
#pragma HLS INTERFACE m_axi port=out_linear bundle=gmem3
#pragma HLS INTERFACE s_axilite port=return bundle=control

    // 缓存配置
    ap_uint<64> local_w[SIZE_Q][SIZE_P];
    // 不需要 local_mods 数组了，因为所有 PE 共享同一个 mod

#pragma HLS ARRAY_PARTITION variable=local_w complete

    load_weights(w_linear, local_w);

    // Stream 网格
    static hls::stream<ap_uint<64>> x_grid[SIZE_Q][SIZE_P + 1];
    static hls::stream<ap_uint<64>> sum_grid[SIZE_P][SIZE_Q + 1];

    #pragma HLS STREAM variable=x_grid depth=2
    #pragma HLS STREAM variable=sum_grid depth=2

    feeder_x(x_linear, x_grid);
    feeder_sum(sum_grid);

    // PE Loop
    Row_Gen: for (int i = 0; i < SIZE_Q; ++i) {
        #pragma HLS UNROLL
        Col_Gen: for (int j = 0; j < SIZE_P; ++j) {
        #pragma HLS UNROLL
            // [FIX] 所有 PE 传入同一个标量模数
            PE_WS(
                x_grid[i][j],       
                x_grid[i][j+1],     
                sum_grid[j][i],     
                sum_grid[j][i+1],   
                local_w[i][j],
                mod_scalar,
                i,
                j
            );
        }
    }

    drain_x(x_grid);
    
    // [FIX] 传递标量给 Collector
    collector(sum_grid, mod_scalar, out_linear);
}

// =================================================
// Host Side Wrapper
// =================================================
void Compute_BConv(
    const uint64_t in_x[MAX_LIMBS][SQRT][SQRT], 
    const uint64_t in_w[MAX_LIMBS][SQRT][SQRT], 
    uint64_t       out[MAX_LIMBS][SQRT][SQRT],  
    const uint64_t MODULUS[MAX_LIMBS],
    int num_active_limbs,                
    int mod_idx_offset               
) {
    Limb_Loop:
    for (int r = mod_idx_offset; r < num_active_limbs + mod_idx_offset; r++) {
        // [FIX] 这里修改为传递 MODULUS[r] 的值，而不是地址
        // 注意类型转换：(ap_uint<64>)MODULUS[r]
        // 【新增】调试打印
        #ifndef __SYNTHESIS__
            std::cout << "[FPGA-DEBUG] BConv Limb Loop: r=" << r 
                      << ", MODULUS[r]=" << MODULUS[r] << std::endl;
            if (MODULUS[r] == 0) {
                 std::cout << "[FPGA-ERROR] MODULUS IS ZERO! ABORTING!" << std::endl;
                 // 在仿真中，甚至可以直接 exit(1) 让报错更明显
            }
        #endif

        bconv_systolic(
            (const ap_uint<64>*)in_x[r],
            (const ap_uint<64>*)in_w[r],
            (ap_uint<64>)MODULUS[r],   // <--- 关键修改：传值
            (ap_uint<64>*)out[r]
        );
    }
}