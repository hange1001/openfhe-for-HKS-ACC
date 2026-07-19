# OpenFHE NTT Implementation - Comprehensive Exploration Report
**Date:** April 15, 2026  
**Project:** OpenFHE for HKS Accelerator (FPGA Backend + Software NTT)

---

## 1. PROJECT STRUCTURE

### Root Directory
```
/home/timhan/FHE/openfhe-for-HKS-ACC/
├── src/
│   ├── core/           # Software NTT implementations
│   ├── pke/            # Public Key Encryption schemes
│   └── fpga_backend/   # FPGA HLS kernel implementations ⭐
├── benchmark/          # Performance benchmarks
├── build/              # CMake build outputs
├── test/               # Unit tests
└── docs/               # Documentation
```

### Key Directories for NTT
- **FPGA Backend:** `src/fpga_backend/` (HLS C++, Xilinx Vitis)
- **Software NTT:** `src/core/include/math/hal/intnat/` (native integer backend)
- **Transform Headers:** `src/core/include/math/hal/*/transform*.h`

---

## 2. NTT FILE LOCATIONS & CLASSIFICATIONS

### 2.1 FPGA Backend NTT Implementation (HLS)

**Main NTT Kernel Files:**

| File | Type | Lines | Purpose |
|------|------|-------|---------|
| `src/fpga_backend/src/ntt_kernel.cpp` | Source | 424 | Core NTT/INTT butterfly operations, ping-pong buffering |
| `src/fpga_backend/include/ntt_kernel.h` | Header | 154 | NTT kernel function signatures |
| `src/fpga_backend/testbench/ntt_kernel_tb.cpp` | Testbench | 500+ | C-Sim verification of NTT functions |

**Supporting Kernel Files:**

| File | Type | Purpose |
|------|------|---------|
| `src/fpga_backend/src/arithmetic.cpp` | Source | Modular arithmetic (MultMod, AddMod, Karatsuba) |
| `src/fpga_backend/include/arithmetic.h` | Header | Arithmetic function signatures |
| `src/fpga_backend/src/bconv.cpp` | Source | Base conversion via systolic array |
| `src/fpga_backend/include/bconv.h` | Header | BConv function signatures |
| `src/fpga_backend/src/interleave.cpp` | Source | Interleaving transformation for NTT prep |
| `src/fpga_backend/include/interleave.h` | Header | Interleave function signatures |
| `src/fpga_backend/src/mod_add_kernel.cpp` | Source | Modular addition kernel |
| `src/fpga_backend/src/mod_sub_kernel.cpp` | Source | Modular subtraction kernel |
| `src/fpga_backend/src/mod_mult_kernel.cpp` | Source | Modular multiplication kernel |
| `src/fpga_backend/src/top.cpp` | Source | Top-level kernel dispatcher & memory manager |
| `src/fpga_backend/include/top.h` | Header | Top-level interfaces |
| `src/fpga_backend/src/load.cpp` | Source | DDR ↔ BRAM data transfer |
| `src/fpga_backend/src/auto.cpp` | Source | Auto-sync operations |

**Configuration & Headers:**

| File | Type | Purpose |
|------|------|---------|
| `src/fpga_backend/include/define.h` | Header | Core constants (RING_DIM, PE_PARALLEL, LIMB_Q, LIMB_P, etc.) |
| `src/fpga_backend/include/opcode.h` | Header | Operation codes (OP_NTT, OP_INTT, OP_BCONV, etc.) |
| `src/fpga_backend/include/memory.h` | Header | Memory layout definitions |

### 2.2 Software NTT Implementation (Native Backend)

| File | Type | Purpose |
|------|------|---------|
| `src/core/include/math/hal/intnat/transformnat.h` | Header | NumberTheoreticTransform declarations |
| `src/core/include/math/hal/intnat/transformnat-impl.h` | Header | NTT implementation for native backend |
| `src/core/extras/ntt1.cpp` | Example | NTT test/demo code (basic testing) |
| `src/core/extras/ntt2.cpp` | Example | Another NTT test variant |

### 2.3 Other Transform Implementations

| File | Type | Purpose |
|------|------|---------|
| `src/core/include/math/hal/bigintfxd/transformfxd-impl.h` | Header | NTT for fixed-size big integers |
| `src/core/include/math/hal/bigintdyn/transformdyn-impl.h` | Header | NTT for dynamic-size big integers |
| `src/core/include/math/hal/bigintntl/transformntl-impl.h` | Header | NTT for NTL backend |

---

## 3. KEY CONSTANTS & CONFIGURATION

### From `src/fpga_backend/include/define.h`

```cpp
// =========================================================
// Core Dimensions
// =========================================================
static const int RING_DIM = 1 << 12;           // 4096 (polynomial degree)
static const int SQRT = 1 << 6;                // 64 (spatial decomposition: 64×64 = 4096)
static const int LOG_SQRT = 6;                 // log2(64)

// =========================================================
// Processing Elements (PE) Parallelism
// =========================================================
static const int BU_NUM = 32;                  // Number of butterfly units
static const int PE_PARALLEL = 8;              // Parallelism factor (UNROLL, cyclic partition, Twiddle replicas)

// =========================================================
// Modulus Structure (CRT with Q and P moduli)
// =========================================================
static const int LIMB_Q = 3;                   // Q moduli count (indices 0, 1, 2)
static const int LIMB_P = 2;                   // P moduli count (indices 3, 4)
static const int MAX_OUT_COLS = LIMB_Q + LIMB_P;  // 5 (for base conversion)
static const int MAX_LIMBS = LIMB_Q + MAX_OUT_COLS;  // 8 (total limbs)

// =========================================================
// Transform Stages
// =========================================================
static const int STAGE = 12;                   // log2(RING_DIM) = log2(4096)

// =========================================================
// Operation Codes
// =========================================================
#define OP_INIT   0     // Initialize parameters
#define OP_ADD    1     // Modular addition
#define OP_SUB    2     // Modular subtraction
#define OP_MULT   3     // Modular multiplication
#define OP_NTT    4     // Forward NTT
#define OP_INTT   5     // Inverse NTT
#define OP_BCONV  6     // Base conversion
#define OP_AUTO   7     // Auto/reserved

// =========================================================
// Baby-Step Giant-Step (BSGS) for Automorphisms
// =========================================================
static const int K_LIST[SQRT] = {1, 2, 3, ..., 64};  // 64 baby-steps
```

### Interpretation of Constants

- **RING_DIM = 4096:** NTT operates on degree-4096 polynomials
- **SQRT = 64:** 2D spatial decomposition: 64×64 elements → efficient bank access
- **PE_PARALLEL = 8:** 8-way parallelism in:
  - Butterfly unrolling
  - Cyclic array partitioning
  - Twiddle factor memory replication
- **LIMB_Q=3, LIMB_P=2:** RNS moduli structure for FHE (3 moduli for main computation, 2 for auxiliary)
- **STAGE=12:** NTT has 12 radix-2 stages (2^12 = 4096)

---

## 4. NTT KERNEL ARCHITECTURE

### 4.1 Configurable_PE (Processing Element)

**Location:** `src/fpga_backend/src/ntt_kernel.cpp` (lines 193-243)

**Signature:**
```cpp
void Configurable_PE(
    const uint64_t &input1,          // First radix-2 element
    const uint64_t &input2,          // Second radix-2 element  
    const uint64_t &twiddle_factor,  // Twiddle factor ω
    uint64_t &res1,                  // Output element 1
    uint64_t &res2,                  // Output element 2
    const uint64_t &modulus,         // Modulus q
    const uint64_t &K_HALF,          // Barrett param (log2(modulus))
    const uint64_t &M,               // Barrett param (multiplier)
    const bool &is_ntt               // Transform direction
);
```

**Butterfly Logic (NTT Mode):**
```
if (is_ntt) {
    temp = MultMod(input2 * twiddle_factor, modulus)
    res1 = AddMod(input1 + temp, modulus)          // (x + y*ω) mod q
    res2 = SubMod(input1 - temp, modulus)          // (x - y*ω) mod q
}
```

**Butterfly Logic (INTT Mode):**
```
if (!is_ntt) {
    temp1 = AddMod(input1 + input2, modulus)
    res2 = SubMod(input1 - input2, modulus)
    res1 = RightShift(temp1, 1)                     // (x + y) / 2
    res2 = MultMod(res2 * twiddle_factor) >> 1      // (x - y) * ω / 2
}
```

**Key Features:**
- Supports both forward NTT and inverse INTT with same PE
- Uses modular arithmetic (AddMod, MultMod, SubMod)
- Twiddle factors pre-computed and fetched from dedicated memory

### 4.2 NTT_Kernel (Single Limb)

**Location:** `src/fpga_backend/src/ntt_kernel.cpp` (lines 245-381)

**Signature:**
```cpp
void NTT_Kernel(
    uint64_t in_memory[SQRT][SQRT],                 // 64×64 polynomial data
    const uint64_t modulus,                         // Prime modulus
    const uint64_t K_HALF,                          // Barrett parameter
    const uint64_t M,                               // Barrett multiplier
    const uint64_t ntt_twiddle_memory[PE_PARALLEL][RING_DIM],   // NTT twiddles (8 copies)
    const uint64_t intt_twiddle_memory[PE_PARALLEL][RING_DIM],  // INTT twiddles (8 copies)
    bool is_ntt                                     // Direction
);
```

**Algorithm Overview:**

1. **Ping-Pong Buffering:** Two 64×64 arrays (buf_A, buf_B) eliminate read-after-write (RAW) dependencies
   - Even stages: read buf_A, write buf_B
   - Odd stages: read buf_B, write buf_A
   - Final result automatically in correct buffer due to STAGE=12 (even)

2. **12 Stages × 64 Rows:**
   - Outer loop: 12 NTT stages
   - Inner loop: 64 rows (one row per pipeline cycle, II=1)
   - 32 butterfly units per row (BU_NUM=32 pairs = 64 elements)

3. **Per-Stage Operations:**
   - **compute_indices():** Calculate input/output permutation indices
   - **read_data():** Fetch data from source buffer with bit-reversal addressing
   - **permutate_data():** Reorder data using input indices
   - **generate_twiddle_index():** Calculate twiddle factor bank accesses
   - **permute_twiddle_factors():** Fetch twiddles with round-robin load distribution
   - **compute_core():** Execute 32 PEs in parallel with pipelined butterfly operations
   - **repermute_data():** Reorder outputs using output indices
   - **rewrite_data():** Write results back to destination buffer

**Pragma Directives:**
```cpp
#pragma HLS ARRAY_PARTITION variable=buf_A cyclic factor=PE_PARALLEL dim=2
#pragma HLS ARRAY_PARTITION variable=buf_B cyclic factor=PE_PARALLEL dim=2
#pragma HLS ARRAY_PARTITION variable=ntt_twiddle_memory complete dim=1    // 8 complete copies
#pragma HLS DEPENDENCE variable=buf_A inter false                           // No inter-iteration deps
#pragma HLS DEPENDENCE variable=buf_B inter false
#pragma HLS PIPELINE II=1                                                   // 1 row/cycle
```

### 4.3 Compute_NTT (Multi-Limb Wrapper)

**Location:** `src/fpga_backend/src/ntt_kernel.cpp` (lines 385-424)

**Signature:**
```cpp
void Compute_NTT(
    uint64_t in_memory[MAX_LIMBS][SQRT][SQRT],                      // 8 limbs × 64×64
    const uint64_t ntt_twiddle_memory[MAX_LIMBS][PE_PARALLEL][RING_DIM],
    const uint64_t intt_twiddle_memory[MAX_LIMBS][PE_PARALLEL][RING_DIM],
    const uint64_t modulus[MAX_LIMBS],
    const uint64_t K_HALF[MAX_LIMBS],
    const uint64_t M[MAX_LIMBS],
    bool is_ntt,
    int num_active_limbs,                                            // Which limbs to process
    int mod_idx_offset                                               // Starting limb index
);
```

**Purpose:** Loops over active limbs, calling NTT_Kernel for each RNS modulus independently.

---

## 5. MODULAR ARITHMETIC

### 5.1 Modular Multiplication (Barrett)

**Location:** `src/fpga_backend/src/arithmetic.cpp` (lines 57-108)

**Signature:**
```cpp
void MultMod(
    const uint64_t &a,              // Operand A
    const uint64_t &b,              // Operand B
    const uint64_t &mod,            // Modulus
    const uint64_t &m,              // Barrett multiplier M
    const uint64_t &k_half,         // log2(modulus)
    uint64_t &res_mod               // Result: (a*b) mod q
);
```

**Algorithm (Barrett Reduction):**
```
1. Full multiply:     z = a * b (128-bit result)
2. Extract high 64:   z_high = z >> (k_half - 1)
3. Estimate quotient: q_approx = (z_high * m) >> (k_half + 1)
4. Compute remainder: r = z - q_approx * mod
5. Correction (1-2 iterations):
   while (r >= mod) r -= mod;
```

**Pragmas:**
```cpp
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS BIND_OP variable=res_mult op=mul impl=dsp latency=4      // DSP48 for multiply
```

**Latency:** 4 cycles (matches DSP48 multiplier)

### 5.2 Modular Addition/Subtraction

**Location:** `src/fpga_backend/src/arithmetic.cpp` (lines 3-26)

**Signature:**
```cpp
void AddMod(
    uint64_t &a,
    const uint64_t &b,
    const uint64_t &mod,
    const bool &is_add              // true=add, false=subtract
);
```

**Logic:**
```
if (is_add) {
    temp = a + b
    if (temp >= mod) temp -= mod;
    a = temp;
} else {
    if (a >= b) a = a - b;
    else a = a + mod - b;
}
```

### 5.3 Karatsuba Multiplication (Utility)

**Location:** `src/fpga_backend/src/arithmetic.cpp` (lines 28-51)

**Purpose:** 64-bit × 64-bit → 128-bit multiplication (currently unused; using native `*` operator instead).

---

## 6. TWIDDLE FACTOR ARCHITECTURE

### 6.1 Twiddle Memory Organization

**In top.cpp (lines 34-35):**
```cpp
static uint64_t NTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM];   // 8×8×4096
static uint64_t INTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM];  // 8×8×4096
```

**Dimensions:**
- **MAX_LIMBS:** 8 (one set per RNS modulus)
- **PE_PARALLEL:** 8 (replicated for each PE to avoid bank conflicts)
- **RING_DIM:** 4096 (one entry per possible twiddle index)

**Pragmas (top.cpp, lines 69-72):**
```cpp
#pragma HLS ARRAY_PARTITION variable=NTTTwiddleFactor complete dim=2
#pragma HLS ARRAY_PARTITION variable=INTTTwiddleFactor complete dim=2
#pragma HLS BIND_STORAGE variable=NTTTwiddleFactor type=ram_2p impl=uram
#pragma HLS BIND_STORAGE variable=INTTTwiddleFactor type=ram_2p impl=uram
```

**Key Points:**
- 8 complete copies eliminate collision when multiple PEs access simultaneously
- Each PE gets dedicated replica → 0 bank conflicts
- Stored in URAM (U55C: 960 URAM blocks ≈ 34 MB total capacity)

### 6.2 Twiddle Factor Access Pattern

**Location:** `src/fpga_backend/src/ntt_kernel.cpp` (lines 94-118)

**generate_twiddle_index():** Computes which twiddle index to fetch based on stage & row:
```
if (stage < STAGE/2) {
    idx = (1 << stage) - 1 + (row >> (STAGE/2 - stage)) + (unit >> (STAGE - stage - 1))
} else {
    idx = (1 << stage) - 1 + (row << (stage - STAGE/2)) + (unit >> (STAGE - stage - 1))
}
```

**permute_twiddle_factors():** Maps BU_NUM butterfly units to PE_PARALLEL replicas:
```cpp
for (int l = 0; l < BU_NUM; l++) {
    #pragma HLS UNROLL factor=PE_PARALLEL
    TwiddleFactor[l] = NTTTWiddleRAM[l % PE_PARALLEL][TwiddleIndex[l]];
}
```

---

## 7. PING-PONG BUFFER (Double Buffering)

### Implementation Details

**Location:** `src/fpga_backend/src/ntt_kernel.cpp` (lines 258-361)

**Buffers:**
```cpp
uint64_t buf_A[SQRT][SQRT];     // 64×64, ping
uint64_t buf_B[SQRT][SQRT];     // 64×64, pong
```

**Partitioning:**
```cpp
#pragma HLS ARRAY_PARTITION variable=buf_A cyclic factor=PE_PARALLEL dim=2
#pragma HLS ARRAY_PARTITION variable=buf_B cyclic factor=PE_PARALLEL dim=2
```

**Protocol:**
```cpp
STAGE_LOOP: for (int j = 0; j < STAGE; j++) {
    if ((j & 1) == 0) {
        // Even stage: READ from buf_A, WRITE to buf_B
        read_data(..., buf_A);
        rewrite_data(..., buf_B);
    } else {
        // Odd stage: READ from buf_B, WRITE to buf_A
        read_data(..., buf_B);
        rewrite_data(..., buf_A);
    }
}
```

**Benefit:** Eliminates read-after-write (RAW) dependencies between iterations, enabling II=1 pipelining.

**Final Writeback:**
```cpp
if ((STAGE & 1) == 0) {
    in_memory[i][l] = buf_A[i][l];  // STAGE=12 (even), result in buf_A
} else {
    in_memory[i][l] = buf_B[i][l];
}
```

---

## 8. BASE CONVERSION (BConv)

### 8.1 Systolic Array Architecture

**Location:** `src/fpga_backend/src/bconv.cpp`

**Purpose:** Convert polynomial from Q moduli (LIMB_Q=3) to P moduli (LIMB_P=2) or arbitrary RNS conversion.

**Core Function:**
```cpp
void bconv_core(
    const uint64_t in_x[LIMB_Q][RING_DIM],          // 3×4096 (Q limbs)
    uint64_t out_x[MAX_OUT_COLS][RING_DIM],         // 5×4096 (Q+P limbs)
    const uint64_t in_w[LIMB_Q][MAX_OUT_COLS],      // 3×5 weight matrix
    const uint64_t out_mod[MAX_OUT_COLS],
    const uint64_t out_k_half[MAX_OUT_COLS],
    const uint64_t out_m_barrett[MAX_OUT_COLS],
    int sizeP                                        // Active output columns
);
```

**Algorithm:** For each coefficient n, compute all output columns in parallel:
```cpp
for (int p = 0; p < MAX_OUT_COLS; ++p) {        // 5 columns (unrolled)
    for (int q = 0; q < LIMB_Q; ++q) {          // 3 rows (unrolled)
        prod[q] = MultMod(in_x[q][n], in_w[q][p], ...);  // 3 MultMod in parallel
    }
    // Reduction tree: sum 3 products mod out_mod[p]
    s01 = (prod[0] + prod[1]) % out_mod[p];
    s = (s01 + prod[2]) % out_mod[p];
    out_x[p][n] = s;
}
```

**Parallelism:** 5 output columns × 3 input limbs × II=1 → 15 parallel MultMod operations.

---

## 9. INTERLEAVING TRANSFORMATION

### Purpose & Implementation

**Location:** `src/fpga_backend/src/interleave.cpp`

**Purpose:** Reorder polynomial coefficients before/after NTT to optimize memory access patterns.

**Function:**
```cpp
void InterLeave(
    uint64_t data[SQRT][SQRT],    // 64×64 in-place
    const bool is_right_shift     // Direction: true=right, false=left
);
```

**Algorithm:** Spatial rotation with ping-pong buffering:
1. **Stage 1:** data → temp_buffer (with circular shift per row)
2. **Stage 2:** temp_buffer → data (write-back)

**Shift Logic:**
```cpp
if (is_right_shift) {
    // For each output position k, read from source (k - row_index)
    temp_buffer[i][k] = data[i][(k - i + SQRT) & (SQRT - 1)];
} else {
    // For each output position k, read from source (k + row_index)
    temp_buffer[i][k] = data[i][(k + i) & (SQRT - 1)];
}
```

**Pragmas:**
```cpp
#pragma HLS ARRAY_PARTITION variable=temp_buffer cyclic factor=PE_PARALLEL dim=2
#pragma HLS ARRAY_PARTITION variable=data cyclic factor=PE_PARALLEL dim=2
#pragma HLS PIPELINE II=1
#pragma HLS UNROLL factor=PARALLEL_FACTOR
```

---

## 10. TOP-LEVEL KERNEL (Top.cpp)

### Main Function Signature

**Location:** `src/fpga_backend/src/top.cpp` (lines 37-44)

```cpp
void Top(
    const uint64_t *mem_in1,         // Primary input (DDR)
    const uint64_t *mem_in2,         // Secondary input (DDR)
    uint64_t *mem_out,               // Output (DDR)
    const uint8_t opcode,            // Operation code
    const int num_active_limbs,      // Number of RNS limbs to process
    const int mod_index              // Starting limb index
);
```

### Kernel Interfaces

**AXI Memory Master (DDR):**
```cpp
#pragma HLS INTERFACE m_axi port=mem_in1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=mem_in2 offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=mem_out offset=slave bundle=gmem2
```

**AXI Slave Control (Register Interface):**
```cpp
#pragma HLS INTERFACE s_axilite port=mem_in1 bundle=control
#pragma HLS INTERFACE s_axilite port=mem_in2 bundle=control
// ... rest of control signals
#pragma HLS INTERFACE s_axilite port=return bundle=control
```

### On-Chip Memory

```cpp
static uint64_t poly_buffer_1[MAX_LIMBS][SQRT][SQRT];     // 8×64×64
static uint64_t poly_buffer_2[MAX_LIMBS][SQRT][SQRT];     // 8×64×64
static uint64_t result_buffer[MAX_LIMBS][SQRT][SQRT];     // 8×64×64
static uint64_t NTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM];
static uint64_t INTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM];
```

**Storage Binding:**
```cpp
#pragma HLS BIND_STORAGE variable=poly_buffer_1 type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=poly_buffer_2 type=ram_2p impl=bram
#pragma HLS BIND_STORAGE variable=NTTTwiddleFactor type=ram_2p impl=uram
```

### Operation Dispatch

```cpp
switch(opcode) {
    case OP_INIT:   // Initialize modulus & twiddle factors
    case OP_ADD:    // Polynomial addition
    case OP_SUB:    // Polynomial subtraction
    case OP_MULT:   // Polynomial multiplication
    case OP_NTT:    // Forward NTT (with interleaving)
    case OP_INTT:   // Inverse INTT (with interleaving)
    case OP_BCONV:  // Base conversion
}
```

---

## 11. FUNCTION CALL HIERARCHY

```
Top()
├── OP_INIT
│   ├── Load Q moduli & parameters
│   ├── Load P moduli & parameters
│   ├── Initialize NTTTwiddleFactor (broadcast to 8 copies)
│   └── Initialize INTTTwiddleFactor (broadcast to 8 copies)
│
├── OP_NTT
│   ├── Load(mem_in1 → poly_buffer_1)
│   ├── InterLeave(poly_buffer_1, is_right_shift=true)
│   ├── Compute_NTT(poly_buffer_1 → result_buffer)
│   │   └── NTT_Kernel() ×num_active_limbs
│   │       ├── compute_indices()
│   │       ├── read_data()
│   │       ├── permutate_data()
│   │       ├── generate_twiddle_index()
│   │       ├── permute_twiddle_factors()
│   │       ├── compute_core()
│   │       │   └── Configurable_PE() ×BU_NUM
│   │       │       ├── MultMod()
│   │       │       ├── AddMod()
│   │       │       └── SubMod()
│   │       ├── repermute_data()
│   │       └── rewrite_data()
│   └── Store(result_buffer → mem_out)
│
├── OP_INTT
│   ├── Load(mem_in1 → poly_buffer_1)
│   ├── Compute_NTT(..., is_ntt=false)
│   │   └── NTT_Kernel(is_ntt=false)  [Similar to NTT but reverse stages & PE logic]
│   ├── InterLeave(poly_buffer_1, is_right_shift=false)
│   └── Store(poly_buffer_1 → mem_out)
│
├── OP_BCONV
│   ├── Load Q limbs (0..LIMB_Q-1)
│   ├── Compute_BConv() [Base conversion Q → P or arbitrary]
│   └── Store results
│
└── OP_MULT, OP_ADD, OP_SUB
    ├── Load both operands
    ├── Compute_Mult/Add/Sub()
    └── Store result
```

---

## 12. KEY ALGORITHM SEQUENCES

### 12.1 Forward NTT Sequence

```
OP_NTT:
  1. Load polynomial from DDR to poly_buffer_1
  2. InterLeave(right_shift=true)    → Reorder for NTT memory pattern
  3. Compute_NTT(is_ntt=true)        → Run 12 stages of radix-2 NTT
     For stage j = 0 to 11:
       For row k = 0 to 63:
         1. compute_indices(j, k, input_idx[], output_idx[])
         2. read_data from appropriate buffer (buf_A/buf_B)
         3. permutate_data(input_idx)
         4. generate_twiddle_index(j, k, tw_idx[])
         5. permute_twiddle_factors(NTT twiddles, tw_idx)
         6. compute_core() → Execute 32 PEs in parallel
         7. repermute_data(output_idx)
         8. rewrite_data to opposite buffer
  4. Store results to DDR
```

### 12.2 Inverse INTT Sequence

```
OP_INTT:
  1. Load polynomial from DDR to poly_buffer_1
  2. Compute_NTT(is_ntt=false)
     For stage j = 0 to 11:
       stage_index = 11 - j (reverse order)
       For row k = 0 to 63:
         1-7. Same as forward NTT
         6. compute_core(is_ntt=false)  → Different PE logic
  3. InterLeave(left_shift=false)     → Reorder back from NTT domain
  4. Store results to DDR
```

---

## 13. DATA STRUCTURES & MEMORY LAYOUT

### 13.1 Polynomial Storage (DDR Format)

```
Polynomial[MAX_LIMBS][SQRT][SQRT]
= [limb0, limb1, ..., limb7][row0..row63][col0..col63]

Each limb is independent RNS representation:
  poly_q(x) = sum_{i=0}^{63} sum_{j=0}^{63} P[q][i][j] * x^(i*64 + j)

Where q ∈ {0,1,2} for Q, {3,4} for P
```

### 13.2 Twiddle Factor Storage (URAM Format)

```
NTTTwiddleFactor[MAX_LIMBS][PE_PARALLEL][RING_DIM]
= [limb0..limb7][replica0..replica7][index0..index4095]

Each replica stores the same twiddle sequence:
  ω^k mod q_limb  for k = 0, 1, ..., 4095

Access pattern (round-robin):
  Butterfly unit l reads from replica[l % PE_PARALLEL][index[l]]
```

---

## 14. PERFORMANCE-CRITICAL ASPECTS

### 14.1 Pipeline II = 1

- **NTT_Kernel:** Processes 1 row per cycle (64 rows = 64 cycles per stage)
- **12 stages × 64 rows = 768 cycles** for one NTT transform

### 14.2 Bank Conflict Elimination

- **Cyclic Partitioning (dim=2):** Distributes 64 columns across 8 banks
  - Column 0,8,16,24,32,40,48,56 → Bank 0
  - Column 1,9,17,25,33,41,49,57 → Bank 1
  - ... etc (stride-8 interleaving)

- **Twiddle Replication:** 8 complete copies ensure zero collision when all PEs access simultaneously

### 14.3 Initialization Overhead

- **Twiddle loading (OP_INIT):** ~1M cycles for loading 4 sets × 8 replicas × 4096 values
- **Broadcast pattern:** Each value loaded once, replicated 8 times

---

## 15. SOFTWARE NTT (Native Backend)

### Location & Structure

**Main Implementation:**
- `src/core/include/math/hal/intnat/transformnat-impl.h` (900+ lines)

**Key Classes:**
```cpp
template <typename VecType>
class NumberTheoreticTransformNat {
public:
    static void ForwardTransformIterative(const VecType& element,
                                         const VecType& rootOfUnityTable,
                                         VecType* result);
    static void InverseTransformIterative(const VecType& element,
                                         const VecType& rootOfUnityInverseTable,
                                         VecType* result);
    // ... more methods
};

template <typename VecType>
class ChineseRemainderTransformFTTNat {
    // Cooley-Tukey FFT variant for FTT
    // Precomputed root-of-unity tables for various moduli
};

template <typename VecType>
class BluesteinFFTNat {
    // Bluestein algorithm for arbitrary-size NTT (not power-of-2)
};

template <typename VecType>
class ChineseRemainderTransformArbNat {
    // Arbitrary-size NTT using cyclotomic polynomials
};
```

### Butterfly Logic (Software)

```cpp
// From transformnat-impl.h, lines 164-186
for (logm = 1; logm <= logn; logm++) {
    for (j = 0; j < n; j += (1 << logm)) {
        for (i = 0; i < (1 << (logm-1)); i++) {
            omega = rootOfUnityTable[indexes[i]];
            indexEven = j + i;
            indexOdd = indexEven + (1 << (logm-1));
            
            oddVal = result[indexOdd];
            omegaFactor = ModMul(omega, oddVal, modulus, mu);
            evenVal = result[indexEven];
            
            result[indexEven] = (evenVal + omegaFactor) % modulus;
            result[indexOdd] = (evenVal - omegaFactor + modulus) % modulus;
        }
    }
}
```

---

## 16. TESTBENCH STRUCTURE

### 16.1 NTT Testbench (`ntt_kernel_tb.cpp`)

**Coverage:**
1. `exact_log2()` - Bit-width calculation
2. `generate_input_index/output_index()` - Index permutation logic
3. `read_data/rewrite_data()` - Memory access patterns
4. `permutate_data/repermute_data()` - Data reordering
5. `generate_twiddle_index()` - Twiddle fetch logic
6. `Configurable_PE()` - Butterfly operation (NTT + INTT)
7. `NTT_Kernel()` - Full single-limb NTT
8. `Compute_NTT()` - Multi-limb NTT with round-trip verification

**Compilation:**
```bash
g++ -std=c++14 -O2 -DFPGA_STANDALONE_TEST \
    -I../include ntt_kernel_tb.cpp \
    ../src/ntt_kernel.cpp ../src/arithmetic.cpp \
    -o ntt_tb && ./ntt_tb
```

### 16.2 BConv Testbench (`bconv_tb.cpp`)

**Tests base conversion Q → P with various weight matrices.**

### 16.3 Auto Testbench (`auto_test.cpp`)

**Tests auto-sync operations for FHE schemes.**

---

## 17. BUILD SYSTEM

### CMakeLists Configuration

**Primary Build:**
```bash
cd build
cmake ..
make
```

**HLS Synthesis:**
```bash
cd src/fpga_backend
vivado_hls -f csynth.tcl  # Synthesis script
```

**Binary Generation:**
```bash
make u55c  # Generate .xclbin for Alveo U55C FPGA
```

---

## 18. INTEGRATION WITH FPGA MANAGER

**Location:** `src/core/include/FpgaManager.h`

**Purpose:** Bridges software OpenFHE operations to hardware FPGA kernels.

**When OPENFHE_FPGA_ENABLE is defined:**
- NTT operations can offload to FPGA
- Base conversions use hardware BConv
- Data marshaling handles DDR transfers

---

## 19. SUMMARY TABLE

| Aspect | Value | Details |
|--------|-------|---------|
| **Polynomial Degree** | 4096 | Ring dimension |
| **Spatial Layout** | 64×64 | 2D decomposition |
| **PE Count** | 32 | Butterfly units per stage |
| **Parallelism Factor** | 8 | UNROLL/cyclic/replicas |
| **NTT Stages** | 12 | log2(4096) |
| **Cycles/NTT** | ~768 | 12 stages × 64 rows |
| **RNS Moduli (Q)** | 3 | LIMB_Q |
| **RNS Moduli (P)** | 2 | LIMB_P |
| **Twiddle Storage** | URAM | 8×8×4096 entries per limb |
| **Poly Storage** | BRAM | 8×64×64 entries per buffer |
| **Memory Interface** | AXI-4 | DDR via M_AXI masters |
| **Control Interface** | AXI Lite | Register-based (s_axilite) |

---

## 20. RECOMMENDED EXPLORATION NEXT STEPS

1. **Run Testbenches:**
   ```bash
   cd src/fpga_backend/testbench
   g++ -DFPGA_STANDALONE_TEST -I../include ntt_kernel_tb.cpp ../src/*.cpp -o ntt_tb
   ./ntt_tb  # Verify algorithm correctness
   ```

2. **Analyze HLS Synthesis:**
   - Review `Solution/solution1/.autopilot/db/` for synthesis logs
   - Check `vitis_hls.log` for resource utilization

3. **Study Bit Reversal Logic:**
   - `generate_input_index()` and `generate_output_index()`
   - Understand modular decomposition of permutation

4. **Explore Twiddle Generation:**
   - Look for root-of-unity computation in OP_INIT
   - Understand Barrett parameter (k_half, m) calculation

5. **Profile Memory Access:**
   - Trace `read_data()` and `rewrite_data()` patterns
   - Verify no bank conflicts with cyclic partitioning

---

**End of Report**
