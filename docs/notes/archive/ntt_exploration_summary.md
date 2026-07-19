# OpenFHE FPGA Backend - NTT Implementation Exploration Report

**Date**: April 15, 2026  
**Project Root**: `/home/timhan/FHE/openfhe-for-HKS-ACC`  
**FPGA Backend**: `src/fpga_backend/`

---

## 1. PROJECT DIRECTORY STRUCTURE

```
src/fpga_backend/
├── include/              # Header files for HLS kernels
├── src/                  # HLS C++ implementation
├── testbench/            # C++ testbenches (csim)
├── Solution/             # Vitis HLS project (auto-generated)
├── Solution_ntt/         # Previous NTT experiment directory
├── Makefile              # Build automation
├── *.tcl                 # HLS synthesis/C-sim scripts
├── *.cfg                 # Device configuration files
└── HW_XCLBIN_ANALYSIS.md # Comprehensive synthesis failure analysis
```

### Key Source Files:
- **ntt_kernel.h** (154 lines) - NTT kernel function declarations
- **ntt_kernel.cpp** (424 lines) - NTT implementation with ping-pong buffering
- **ntt_kernel_tb.cpp** (708 lines) - Comprehensive HLS C-sim testbench
- **top.cpp** (259 lines) - Top-level interface, memory management, opcode dispatch
- **arithmetic.cpp** (108 lines) - Modular arithmetic (AddMod, MultMod via Barrett)
- **bconv.cpp** (146 lines) - Base conversion (systolic array)
- **interleave.cpp** - NTT data permutation/deinterleaving

---

## 2. CONSTANTS & CONFIGURATION (define.h)

```c
// Core Ring Dimensions
static const int RING_DIM = 1 << 12;        // 4096
static const int SQRT = 1 << 6;             // 64
static const int LOG_SQRT = 6;
static const int STAGE = 12;                // log2(RING_DIM)

// Processing Element Parallelism
static const int BU_NUM = 32;               // Butterfly Units (processing elements)
static const int PE_PARALLEL = 8;           // UNROLL factor / Twiddle factor copies
                                            // All cyclic partitions & UNROLL factors use this

// Multi-moduli Support
static const int LIMB_Q = 3;                // Q moduli count (indices 0, 1, 2)
static const int LIMB_P = 2;                // P moduli count (indices 3, 4)
static const int MAX_OUT_COLS = 5;          // LIMB_Q + LIMB_P
static const int MAX_LIMBS = 8;             // LIMB_Q + MAX_OUT_COLS (5+3=8)

// Opcodes
#define OP_INIT   0
#define OP_ADD    1
#define OP_SUB    2
#define OP_MULT   3
#define OP_NTT    4
#define OP_INTT   5
#define OP_BCONV  6
#define OP_AUTO   7
```

---

## 3. NTT KERNEL ARCHITECTURE

### 3.1 Main NTT Function Signature

```cpp
void NTT_Kernel(
    uint64_t in_memory[SQRT][SQRT],           // 2D input: 64×64 = 4096 elements

    const uint64_t modulus,                   // Prime modulus
    const uint64_t K_HALF,                    // Barrett parameter: ceil(log2(mod))
    const uint64_t M,                         // Barrett parameter: floor(2^(2*K_HALF) / mod)

    const uint64_t ntt_twiddle_memory[PE_PARALLEL][RING_DIM],   // NTT twiddles: 8 copies × 4096
    const uint64_t intt_twiddle_memory[PE_PARALLEL][RING_DIM],  // INTT twiddles: 8 copies × 4096

    bool is_ntt                               // true=NTT, false=INTT
);
```

### 3.2 Multi-limb Wrapper

```cpp
void Compute_NTT(
    uint64_t in_memory[MAX_LIMBS][SQRT][SQRT],
    const uint64_t ntt_twiddle_memory[MAX_LIMBS][PE_PARALLEL][RING_DIM],
    const uint64_t intt_twiddle_memory[MAX_LIMBS][PE_PARALLEL][RING_DIM],
    const uint64_t modulus[MAX_LIMBS],
    const uint64_t K_HALF[MAX_LIMBS],
    const uint64_t M[MAX_LIMBS],
    bool is_ntt,
    int num_active_limbs,
    int mod_idx_offset
);
```

---

## 4. ARCHITECTURE: PING-PONG BUFFERING

```cpp
// NTT_Kernel (lines 260-381 in ntt_kernel.cpp)

// Dual buffering to eliminate Read-After-Write (RAW) dependencies
uint64_t buf_A[SQRT][SQRT];
uint64_t buf_B[SQRT][SQRT];

#pragma HLS ARRAY_PARTITION variable=buf_A cyclic factor=PE_PARALLEL dim=2
#pragma HLS ARRAY_PARTITION variable=buf_B cyclic factor=PE_PARALLEL dim=2

// Main loop: 12 stages × 64 rows (768 iterations total)
STAGE_LOOP:
for (int j = 0; j < STAGE; j++) {
    ROW_LOOP:
    for (int k = 0; k < SQRT; k++) {
        #pragma HLS PIPELINE II=1
        #pragma HLS DEPENDENCE variable=buf_A inter false
        #pragma HLS DEPENDENCE variable=buf_B inter false
        
        // Ping-Pong Protocol:
        // Even stages (j=0,2,4,...): read buf_A, write buf_B
        // Odd stages (j=1,3,5,...): read buf_B, write buf_A
        // Different physical arrays => no RAW dependencies
        
        if ((j & 1) == 0) {
            read_data(stage_index, k, ReadData, buf_A);
            // ... processing ...
            rewrite_data(stage_index, k, RepermuteData, buf_B);
        } else {
            read_data(stage_index, k, ReadData, buf_B);
            // ... processing ...
            rewrite_data(stage_index, k, RepermuteData, buf_A);
        }
    }
}
```

**Key Insight**: STAGE=12 is even, so final stage writes to buf_A. Writeback copies from buf_A to in_memory.

---

## 5. PROCESSING ELEMENT (PE) IMPLEMENTATION

### 5.1 Configurable_PE Function

```cpp
void Configurable_PE(
    const uint64_t &input1,
    const uint64_t &input2,
    const uint64_t &twiddle_factor,
    uint64_t &res1,
    uint64_t &res2,
    const uint64_t &modulus,
    const uint64_t &K_HALF,
    const uint64_t &M,
    const bool &is_ntt
) {
    #pragma HLS INLINE
    
    if (is_ntt) {
        // NTT Butterfly (Forward)
        // res1 = (input1 + input2 * twiddle) mod p
        // res2 = (input1 - input2 * twiddle) mod p
        MultMod(input2_temp, twiddle_factor, modulus, M, K_HALF, temp);
        AddMod(input1_temp, temp, modulus, true);   // +
        res1_temp = input1_temp;
        
        input1_temp = input1;
        AddMod(input1_temp, temp, modulus, false);  // -
        res2_temp = input1_temp;
    } else {
        // INTT Butterfly (Inverse)
        // res1 = (input1 + input2) / 2
        // res2 = (input1 - input2) * twiddle / 2
        AddMod(input1_temp, input2_temp, modulus, true);
        temp1 = input1_temp;
        // Divide by 2 using bit shift + odd compensation
        res1 = (temp1 >> 1) + ((temp1 & 1) ? ((modulus + 1) >> 1) : 0);
        
        MultMod(res2_temp, twiddle_factor, modulus, M, K_HALF, temp);
        res2 = (temp >> 1) + ((temp & 1) ? ((modulus + 1) >> 1) : 0);
    }
}
```

### 5.2 PE Parallelism (compute_core)

```cpp
void compute_core(
    uint64_t PermuteData[SQRT],
    uint64_t TwiddleFactor[BU_NUM],
    uint64_t NTTData[SQRT],
    uint64_t modulus, uint64_t K_HALF, uint64_t M,
    bool is_ntt
) {
    #pragma HLS INLINE
    
    int data_pairs_total = SQRT >> 1;        // 32 pairs
    int pairs_per_pe = data_pairs_total / BU_NUM;  // 32 / 32 = 1 pair per PE
    
    for (int l = 0; l < BU_NUM; l++) {       // 32 PEs
        #pragma HLS UNROLL factor=PE_PARALLEL  // Only 8 unroll, rest serialize
        for (int m = 0; m < pairs_per_pe; m++) {
            int global_pair_index = l * pairs_per_pe + m;
            int idx1 = global_pair_index * 2;
            int idx2 = global_pair_index * 2 + 1;
            Configurable_PE(
                PermuteData[idx1],
                PermuteData[idx2],
                TwiddleFactor[l],
                NTTData[idx1],
                NTTData[idx2],
                modulus, K_HALF, M, is_ntt
            );
        }
    }
}
```

**Architecture Note**: 
- BU_NUM = 32 butterfly units in design
- PE_PARALLEL = 8 (actual unroll factor)
- Only 8 PEs unroll per iteration; remaining 24 serialize
- Each PE processes exactly 1 pair of elements per stage

---

## 6. MODULAR ARITHMETIC IMPLEMENTATIONS

### 6.1 Barrett Modular Multiplication

```cpp
// arithmetic.cpp lines 57-108

void MultMod(
    const uint64_t &a,
    const uint64_t &b,
    const uint64_t &mod,
    const uint64_t &m,        // m = floor(2^(2*k_half+2) / mod)
    const uint64_t &k_half,   // k_half = ceil(log2(mod))
    uint64_t &res_mod
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=1
    
    // Step 1: Full precision multiplication
    uint128_t res_mult = (uint128_t)a * b;
    #pragma HLS BIND_OP variable=res_mult op=mul impl=dsp latency=4
    
    // Step 2: Barrett reduction - estimate quotient q
    uint64_t res_mult_high = (uint64_t)(res_mult >> (k_half - 1));
    
    // Step 3: res_mult_high * m
    uint128_t res_mult_shift = (uint128_t)res_mult_high * m;
    #pragma HLS BIND_OP variable=res_mult_shift op=mul impl=dsp latency=4
    
    // Step 4: q = (res_mult_high * m) >> (k_half + 1)
    uint64_t q = (uint64_t)(res_mult_shift >> (k_half + 1));
    
    // Step 5: r = z - q * mod
    uint64_t q_times_mod = q * mod;
    #pragma HLS BIND_OP variable=q_times_mod op=mul impl=dsp latency=4
    uint64_t r = (uint64_t)res_mult - q_times_mod;
    
    // Step 6: Final correction (Barrett guarantees r < 3*mod)
    if (r >= mod) r -= mod;
    if (r >= mod) r -= mod;
    
    res_mod = r;
}

// Pipeline: II=1, LATENCY=4 cycles (3 DSP multiplications)
// Uses 3 × DSP48 (@ latency=4 each, scheduled in parallel)
```

### 6.2 Modular Addition

```cpp
// arithmetic.cpp lines 3-26

void AddMod(
    uint64_t &a,
    const uint64_t &b,
    const uint64_t &mod,
    const bool &is_add
) {
    #pragma HLS INLINE
    
    unsigned __int128 temp_res;
    if (is_add) {
        temp_res = (unsigned __int128)a + b;
        if (temp_res >= mod) {
            temp_res -= mod;
        }
        a = (uint64_t)temp_res;
    } else {
        // Subtraction: a - b (modulo-aware)
        if (a >= b) {
            a = a - b;
        } else {
            temp_res = (unsigned __int128)a + mod - b;
            a = (uint64_t)temp_res;
        }
    }
}
```

---

## 7. MEMORY ARCHITECTURE & ARRAY PARTITIONING

### 7.1 Top-level Static Buffers (top.cpp)

```cpp
// Lines 17-35

// Polynomial data buffers (BRAM)
static uint64_t poly_buffer_1[MAX_LIMBS][SQRT][SQRT];    // 5×64×64 × 8B = 160 KB
static uint64_t poly_buffer_2[MAX_LIMBS][SQRT][SQRT];    // Same
static uint64_t result_buffer[MAX_LIMBS][SQRT][SQRT];    // Same

#pragma HLS ARRAY_PARTITION variable=poly_buffer_1 cyclic dim=3 factor=PE_PARALLEL
#pragma HLS ARRAY_PARTITION variable=poly_buffer_2 cyclic dim=3 factor=PE_PARALLEL
#pragma HLS ARRAY_PARTITION variable=result_buffer cyclic dim=3 factor=PE_PARALLEL
#pragma HLS BIND_STORAGE variable=poly_buffer_1 type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=poly_buffer_2 type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=result_buffer type=ram_2p impl=bram

// Twiddle factors (URAM, 8 copies per limb)
static uint64_t NTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM];    // 8×8×4096 × 8B = 2.1 MB
static uint64_t INTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM];   // Same

#pragma HLS ARRAY_PARTITION variable=NTTTwiddleFactor complete dim=2    // 8 separate URAM banks
#pragma HLS ARRAY_PARTITION variable=INTTTwiddleFactor complete dim=2
#pragma HLS BIND_STORAGE variable=NTTTwiddleFactor type=ram_2p impl=uram
#pragma HLS BIND_STORAGE variable=INTTTwiddleFactor type=ram_2p impl=uram
```

### 7.2 Array Partition Strategy

| Array | Partition | Reason |
|-------|-----------|--------|
| poly_buffer_* (dim=3) | cyclic factor=8 | Matches PE_PARALLEL unroll factor |
| ntt/intt_twiddle (dim=2) | complete | 8 copies → 8 separate physical URAM banks |
| temp arrays in NTT_Kernel | complete | No bank conflicts for small arrays |
| InterLeave temp_buffer (dim=2) | cyclic factor=8 | Matches data partition |

---

## 8. HLS PRAGMAS USAGE PATTERNS

### Key Pragmas Found:

```cpp
// Dataflow/Pipelining
#pragma HLS INLINE / #pragma HLS INLINE off
#pragma HLS PIPELINE II=1 / II=2 / II=4
#pragma HLS LOOP_TRIPCOUNT min=1 max=5 avg=3

// Array Partitioning
#pragma HLS ARRAY_PARTITION variable=X cyclic factor=8 dim=2
#pragma HLS ARRAY_PARTITION variable=X complete dim=1
#pragma HLS ARRAY_PARTITION variable=X complete  // Partition all dimensions

// Loop Optimization
#pragma HLS UNROLL / #pragma HLS UNROLL factor=8
#pragma HLS DEPENDENCE variable=X inter false    // Ignore inter-iteration dependencies

// Interface/Memory
#pragma HLS INTERFACE m_axi port=mem_in1 offset=slave bundle=gmem0
#pragma HLS INTERFACE s_axilite port=opcode bundle=control
#pragma HLS BIND_STORAGE variable=X type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=X type=ram_2p impl=uram
#pragma HLS BIND_OP variable=X op=mul impl=dsp latency=4

// Resource Hints
#pragma HLS LOOP_FLATTEN off
```

---

## 9. NTT KERNEL PIPELINE STAGES (12 Stages)

Each stage processes SQRT rows (64 rows) with II=1 pipelining:

| Stage | Input | Processing | Output | Twiddle Count |
|-------|-------|------------|--------|--------------|
| 0 | SQRT elements | Bit-reversal, index generation | Permuted | BU_NUM (32) |
| 1 | SQRT elements | Butterfly 2×2 | Unpermuted | BU_NUM (32) |
| ... | ... | Cooley-Tukey decimation-in-time | ... | ... |
| 11 | SQRT elements | Final butterfly pass | Result | BU_NUM (32) |

**Total Operations Per Limb**: 12 stages × 64 rows × (32 PEs × 1 pair) = 24,576 butterfly ops

---

## 10. TWIDDLE FACTOR GENERATION & INDEXING

### 10.1 Twiddle Index Generation

```cpp
void generate_twiddle_index(int j, int k, int TwiddleIndex[BU_NUM]) {
    #pragma HLS INLINE
    
    for (int l = 0; l < BU_NUM; l++) {
        #pragma HLS UNROLL factor=PE_PARALLEL
        
        if (j < (STAGE >> 1)) {  // j < 6
            TwiddleIndex[l] = (1 << j) - 1 
                            + (k >> ((STAGE >> 1) - j)) 
                            + (l >> (STAGE - j - 1));
        } else {                 // j >= 6
            TwiddleIndex[l] = (1 << j) - 1 
                            + (k << (j - (STAGE >> 1))) 
                            + (l >> (STAGE - j - 1));
        }
    }
}

// Layout: twiddle_memory[0] = w^0 = 1 (stage 0)
//         twiddle_memory[1..2] = stage 1 roots
//         twiddle_memory[3..6] = stage 2 roots
//         ... (total 4096 roots for 12 stages)
```

### 10.2 Twiddle Retrieval with Bank Conflict Avoidance

```cpp
void permute_twiddle_factors(
    uint64_t TwiddleFactor[BU_NUM],
    const uint64_t NTTTWiddleRAM[PE_PARALLEL][RING_DIM],
    int TwiddleIndex[BU_NUM]
) {
    #pragma HLS INLINE
    
    // Critical: PE_PARALLEL (8) copies, each PE accesses its own copy
    for (int l = 0; l < BU_NUM; l++) {
        #pragma HLS UNROLL factor=PE_PARALLEL
        // Even if all indices are identical, they map to different physical RAMs
        TwiddleFactor[l] = NTTTWiddleRAM[l % PE_PARALLEL][TwiddleIndex[l]];
    }
}

// Result: Zero port conflicts despite identical indices
// (because l % 8 distributes across 8 URAM banks)
```

---

## 11. DATA PERMUTATION FUNCTIONS

### 11.1 Index Computation

```cpp
void generate_input_index(int stage, int address, int output_indices[SQRT]) {
    // Computes SQRT (64) permutation indices for input
    // Related to bit-reversal and stage-specific shuffling
    for (int i = 0; i < SQRT; i++) {
        #pragma HLS UNROLL factor=PE_PARALLEL
        int shift_val = (STAGE >> 1) - stage;    // = 6 - stage
        int mask1 = (1 << (shift_val + 1)) - 1;
        int mask2 = ~((1 << (shift_val + 1)) - 1) & ((1 << LOG_SQRT) - 1);
        
        int iwire = i;
        int temp2 = (iwire & 1) << shift_val;
        int index = ((iwire & mask2) | temp2 | ((iwire & mask1) >> 1)) + address;
        output_indices[i] = index & (SQRT - 1);
    }
}

void generate_output_index(int stage, int address, int output_indices[SQRT]) {
    // Computes inverse permutation for output
    // Guaranteed: output_indices[input_indices[i]] == i (for all i)
    // ...similar masking logic...
}
```

---

## 12. CHARACTERIZATION OF EXISTING IMPLEMENTATION

### Strengths:
1. **Ping-Pong Buffering**: Eliminates RAW dependencies, enables II=1
2. **Bank Conflict Avoidance**: 8 separate twiddle copies avoid port conflicts
3. **Complete PE Unrolling**: 8 PEs unroll per iteration (reduces stall)
4. **Barrett Modular Multiplication**: 3× DSP48 @ II=1, efficient reduction
5. **Comprehensive Testbench**: 11 unit tests validate correctness (ntt_kernel_tb.cpp)

### Critical Issues:
1. ❌ **BRAM Overflow**: Twiddle factors need 10+ MB; U55C has only 32 MB BRAM, but usage is 206% over budget
2. ❌ **Timing Violation**: Slack = −0.33 ns (requires 4.33 ns @ 250 MHz target)
3. ❌ **LUT Saturation**: BConv systolic array uses 30 hardware dividers (198K LUT)
4. ⚠️ **Serialization**: Only 8/32 PEs unroll; 24 PEs serialize per stage

---

## 13. RECOMMENDATIONS FOR CG-NTT MODULE

### Architecture Adaptations Needed:

1. **Twiddle Factor Storage**:
   - Current: All 8 copies stored in URAM (2.1 MB total)
   - **CG-NTT**: Stream twiddle from external DDR or use constant folding for smaller CG window

2. **Data Layout**:
   - Current: 2D [SQRT][SQRT] = [64][64]
   - **CG-NTT**: May benefit from different tiling (e.g., [8][512] or [16][256]) depending on CG parameters

3. **PE Configuration**:
   - Current: 32 PEs, 8 unroll
   - **CG-NTT**: May need different ratio (e.g., 16 PEs, 16 unroll) for other geometries

4. **Index Generation**:
   - Adapt `generate_input_index` / `generate_output_index` for CG-specific permutation patterns
   - May need different masking/shifting logic

5. **Twiddle Indexing**:
   - Current: Linear indexing via `generate_twiddle_index`
   - **CG-NTT**: Potential for constant-time root powers if CG uses simplified root structure

---

## 14. KEY FILES REFERENCE

| File | Lines | Purpose |
|------|-------|---------|
| `define.h` | 82 | Constants, dimensions, opcodes |
| `ntt_kernel.h` | 154 | NTT function prototypes |
| `ntt_kernel.cpp` | 424 | Core NTT implementation + ping-pong logic |
| `arithmetic.h` | 32 | Arithmetic function declarations |
| `arithmetic.cpp` | 108 | ModAdd, ModMult (Barrett), Karatsuba |
| `top.cpp` | 259 | Top-level interface, memory management |
| `ntt_kernel_tb.cpp` | 708 | 11-test validation suite |
| `bconv.cpp` | 146 | Base conversion (systolic) |
| `interleave.cpp` | ~80 | Data permutation (bit-reversal related) |
| `HW_XCLBIN_ANALYSIS.md` | 323 | Detailed synthesis failure analysis |

---

## 15. COMPILATION & TESTING

### C-Simulation (no FPGA synthesis):
```bash
cd src/fpga_backend
g++ -std=c++14 -O2 -DFPGA_STANDALONE_TEST \
    -I./include \
    testbench/ntt_kernel_tb.cpp \
    src/ntt_kernel.cpp src/arithmetic.cpp \
    -o ntt_tb && ./ntt_tb
```

### Vitis HLS C-Sim:
```bash
vitis_hls csynth.tcl -mode csim -project Solution
```

### Synthesis (current status: **FAILS** - 206% BRAM, Slack=-0.33ns):
```bash
vitis_hls csynth.tcl -mode synth -project Solution
```

---

*End of Exploration Report*
