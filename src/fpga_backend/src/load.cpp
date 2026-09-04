#include "define.h"
#include "../include/load.h"

// Fixed-length burst boundary: do not merge this loop with a runtime limb loop.
static void Load_Tower(const uint64_t* mem_in,
                       uint64_t buffer[MAX_LIMBS][SQRT][SQRT], int src, int dst) {
    #pragma HLS INLINE off
    LOAD_COEFF: for (int k = 0; k < RING_DIM; ++k) {
        #pragma HLS PIPELINE II=1
        #pragma HLS UNROLL factor=PE_PARALLEL
        buffer[dst][k / SQRT][k % SQRT] = mem_in[src * RING_DIM + k];
    }
}

static void Store_Tower(uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
                        uint64_t* mem_out, int src, int dst) {
    #pragma HLS INLINE off
    STORE_COEFF: for (int k = 0; k < RING_DIM; ++k) {
        #pragma HLS PIPELINE II=1
        #pragma HLS UNROLL factor=PE_PARALLEL
        mem_out[dst * RING_DIM + k] = buffer[src][k / SQRT][k % SQRT];
    }
}

void Load(
    const uint64_t *mem_in,
    uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
    const int num_active_limbs,
    const int mod_index // 【新增】偏移量
) {
    #pragma HLS INLINE off
    for (int l = mod_index; l < mod_index + num_active_limbs; l++) {
        Load_Tower(mem_in, buffer, l - mod_index, l);
    }
}

void Store(
    uint64_t buffer[MAX_LIMBS][SQRT][SQRT],
    uint64_t *mem_out,
    const int num_active_limbs,
    const int mod_index // 【新增】偏移量
) {
    #pragma HLS INLINE off

    for (int l = mod_index; l < mod_index + num_active_limbs; l++) {
        Store_Tower(buffer, mem_out, l, l - mod_index);
    }
}
