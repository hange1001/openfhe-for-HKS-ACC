# NTT Implementation Exploration - Final Summary

**Status**: ✅ Complete  
**Report Generated**: April 15, 2026  
**Location**: `/home/timhan/FHE/openfhe-for-HKS-ACC/src/fpga_backend/`

---

## EXECUTIVE SUMMARY

This exploration provides a **complete architectural blueprint** of the existing NTT (Number Theoretic Transform) implementation in the OpenFHE hardware accelerator. The implementation is sophisticated with many design optimizations, but has critical resource constraints that prevent successful FPGA synthesis on U55C devices.

### Key Findings:
- ✅ **Well-architected**: Ping-pong buffering, bank conflict avoidance, modular arithmetic optimization
- ❌ **Critical issues**: 206% BRAM overflow, timing violation (Slack = −0.33 ns), LUT saturation
- 📊 **Parameters**: RING_DIM=4096, SQRT=64, STAGE=12, PE=8, BU_NUM=32
- 🔧 **Modular arithmetic**: Barrett reduction @ II=1, 3× DSP48 per multiplication

---

## ARTIFACT FILES

### Main Report
- **`ntt_exploration_summary.md`** (19 KB) - Comprehensive 15-section analysis with code examples

### Generated During Exploration
All files in `/home/timhan/FHE/openfhe-for-HKS-ACC/src/fpga_backend/`:

#### Headers (include/)
```
arithmetic.h            - ModAdd, MultMod, Karatsuba declarations
define.h               - Constants (RING_DIM=4096, PE_PARALLEL=8, etc.)
ntt_kernel.h           - NTT kernel function prototypes (154 lines)
top.h                  - Top-level interface
bconv.h                - Base conversion systolic array
interleave.h           - Data permutation functions
```

#### Implementation (src/)
```
ntt_kernel.cpp         - Core NTT + ping-pong buffering (424 lines)
  ├─ NTT_Kernel()      - Single-limb NTT with 12 stages
  ├─ Compute_NTT()     - Multi-limb wrapper
  ├─ Configurable_PE() - NTT/INTT butterfly unit
  └─ 7 helper functions (indexing, permutation, twiddle)

arithmetic.cpp         - Modular operations (108 lines)
  ├─ MultMod()         - Barrett modular multiplication
  ├─ AddMod()          - Modular addition/subtraction
  └─ Karatsuba()       - 64×64→128 multiplication

top.cpp                - Top-level interface (259 lines)
  ├─ Static buffers (poly_buffer_1, poly_buffer_2, result_buffer)
  ├─ Twiddle storage (NTTTwiddleFactor, INTTTwiddleFactor)
  └─ Opcode dispatch (OP_INIT, OP_NTT, OP_INTT, OP_BCONV, etc.)

bconv.cpp              - Base conversion (146 lines)
  ├─ bconv_core()      - Systolic computation
  └─ Compute_BConv()   - Wrapper with memory management

interleave.cpp         - NTT data permutation (~80 lines)
```

#### Testing
```
ntt_kernel_tb.cpp      - Comprehensive testbench (708 lines)
  ├─ 11 unit tests (exact_log2, index generation, PE butterfly, etc.)
  ├─ Software reference implementations
  └─ Round-trip verification (NTT→INTT)
```

---

## ARCHITECTURAL OVERVIEW

### 1. Core Dimensions
```
RING_DIM  = 4096              (2^12, degree of polynomial)
SQRT      = 64                (2^6, used for 2D buffering)
STAGE     = 12                (log2 RING_DIM, number of NTT stages)
BU_NUM    = 32                (Butterfly units, design intention)
PE_PARALLEL = 8               (Actual unroll factor in HLS)
```

### 2. Ping-Pong Buffer Architecture
```
Time: 0     1     2     3  ...  11    12 (back to buf_A)
      ───────────────────────────────────
buf_A: R     W     R     W        R
buf_B:       R     W     R     W  
                   (ping-pong alternates to eliminate RAW dependencies)
```

### 3. Data Flow Through One Stage
```
Input Data (64 elements)
    ↓
Index Generation
    ↓
Read from buffer A/B
    ↓
Input Permutation
    ↓
Generate Twiddle Indices
    ↓
Fetch Twiddle Factors (8 separate copies for zero conflicts)
    ↓
Compute Core (32 PEs, 8 unroll factor)
    ├─ 8 PE instances in parallel
    ├─ 24 PE instances serialize
    └─ Each PE: 1 butterfly = 1 MultMod + 2 AddMod
    ↓
Output Permutation
    ↓
Write to buffer B/A
    ↓
(Repeat 12 times for all stages, then writeback to in_memory)
```

### 4. Processing Element (Butterfly)
```
NTT Butterfly (is_ntt=true):
  res1 = (input1 + input2 * twiddle) mod p
  res2 = (input1 - input2 * twiddle) mod p

INTT Butterfly (is_ntt=false):
  res1 = (input1 + input2) / 2
  res2 = ((input1 - input2) * twiddle) / 2

Each butterfly requires:
  - 1 MultMod (3 DSP @ II=1, 4-cycle latency)
  - 2 AddMod (combinational logic)
```

---

## ARRAY PARTITIONING STRATEGY

| Array | Dimensions | Partition | Size | Purpose |
|-------|-----------|-----------|------|---------|
| poly_buffer_1 | [8][64][64] | cyclic dim=3 factor=8 | 160 KB | Input polynomial (BRAM) |
| poly_buffer_2 | [8][64][64] | cyclic dim=3 factor=8 | 160 KB | Input polynomial 2 (BRAM) |
| result_buffer | [8][64][64] | cyclic dim=3 factor=8 | 160 KB | Output (BRAM) |
| **NTTTwiddleFactor** | [8][8][4096] | complete dim=2 | **2.1 MB** | **8 copies for 0 conflicts** |
| **INTTTwiddleFactor** | [8][8][4096] | complete dim=2 | **2.1 MB** | **8 copies for 0 conflicts** |
| buf_A (NTT_Kernel) | [64][64] | cyclic dim=2 factor=8 | 32 KB | Ping-pong buffer A |
| buf_B (NTT_Kernel) | [64][64] | cyclic dim=2 factor=8 | 32 KB | Ping-pong buffer B |

**Critical Issue**: The 4.2 MB twiddle storage alone requires ~2560 BRAM_18K blocks. U55C has only 4032 total. With poly buffers and other structures, total usage = **8334 blocks (206% of budget)**.

---

## HLS PRAGMAS PATTERNS

### Pipelining
```cpp
#pragma HLS PIPELINE II=1                    // Target 1 iteration per cycle
#pragma HLS PIPELINE II=2 / II=4            // Relaxed II for higher latency
```

### Array Partitioning
```cpp
#pragma HLS ARRAY_PARTITION variable=X cyclic factor=8 dim=2    // Cyclic across 8 banks
#pragma HLS ARRAY_PARTITION variable=X complete dim=1           // Full unroll in dim 1
#pragma HLS ARRAY_PARTITION variable=X complete                 // Unroll all dimensions
```

### Dependency Directives
```cpp
#pragma HLS DEPENDENCE variable=buf_A inter false   // Ignore inter-iteration RAW
#pragma HLS DEPENDENCE variable=buf_B inter false   // (Allows II=1)
```

### Resource Binding
```cpp
#pragma HLS BIND_STORAGE variable=X type=ram_2p impl=bram   // Use BRAM
#pragma HLS BIND_STORAGE variable=X type=ram_2p impl=uram   // Use URAM
#pragma HLS BIND_OP variable=X op=mul impl=dsp latency=4    // Use DSP48
```

### Loop Optimization
```cpp
#pragma HLS UNROLL factor=8              // Unroll factor
#pragma HLS LOOP_TRIPCOUNT min=1 max=5   // Hint for synthesis
#pragma HLS LOOP_FLATTEN off              // Prevent flattening
#pragma HLS INLINE / #pragma HLS INLINE off
```

---

## MODULAR ARITHMETIC: BARRETT REDUCTION

### MultMod() Implementation
```cpp
// Input: a, b (< mod), twiddle factor
// Output: (a * b) mod p

// 1. Full multiplication (64×64→128, DSP48)
res_mult = a * b;                                  // DSP @ II=1, LAT=4

// 2. High bits extraction
res_mult_high = res_mult >> (k_half - 1);

// 3. Estimate quotient (another DSP @ II=1, LAT=4)
q_estimate = (res_mult_high * m) >> (k_half + 1);

// 4. Correction loop
q_times_mod = q_estimate * mod;                    // DSP @ II=1, LAT=4
r = res_mult - q_times_mod;

// 5. Final reduction
if (r >= mod) r -= mod;
if (r >= mod) r -= mod;  // At most 2 subtractions
```

**Performance**: II=1 pipeline, 4-cycle latency, 3× DSP48 per MultMod

---

## TWIDDLE FACTOR ADDRESSING

### Memory Layout
```
twiddle_memory[0]                    = w^0 = 1              (stage 0)
twiddle_memory[1..2]                 = roots for stage 1    (2 roots)
twiddle_memory[3..6]                 = roots for stage 2    (4 roots)
twiddle_memory[7..14]                = roots for stage 3    (8 roots)
...
twiddle_memory[(1<<11)-1 .. RING_DIM-1] = stage 11 roots    (2048 roots)

Total: 1 + 2 + 4 + 8 + ... + 2048 = 4096 roots
```

### Zero Bank Conflict Design
```cpp
// 32 PEs, but 8 twiddle copies (PE_PARALLEL)
for (int l = 0; l < BU_NUM; l++) {
    #pragma HLS UNROLL factor=PE_PARALLEL
    // Key: l % PE_PARALLEL distributes across 8 copies
    TwiddleFactor[l] = NTTTWiddleRAM[l % 8][TwiddleIndex[l]];
    //                                       ^^^^^^
    //                      Each PE accesses different copy
}

// Result: Even if indices are identical, no port conflicts
// (different physical URAM banks)
```

---

## CRITICAL SYNTHESIS FAILURES

### P0 Issues (Blocking)

1. **BRAM Overflow: 206%**
   - Required: 8334 BRAM_18K
   - Available: 4032 BRAM_18K (U55C)
   - Root: Twiddle factor storage (4.2 MB)

2. **Timing Violation: Slack = −0.33 ns**
   - Target: 4 ns cycle (250 MHz)
   - Actual: 4.33 ns critical path
   - Root: NTT_Kernel_Pipeline has 3 DSP multiplications in series

3. **LUT Explosion: 198K LUT for dividers alone**
   - BConv uses 30× hardware dividers (128-bit division)
   - Each divider ≈ 6607 LUT
   - 30 × 6607 = 198,210 LUT (15% of total)

### Synthesis Report Location
```
src/fpga_backend/Solution/solution1/syn/report/
├── csynth.rpt              - Main synthesis report
├── csynth_design_size.rpt  - Instruction count analysis
└── HW_XCLBIN_ANALYSIS.md   - Detailed failure breakdown
```

---

## TESTBENCH VALIDATION

### ntt_kernel_tb.cpp (708 lines)
11 comprehensive tests covering:

1. **exact_log2()** - Logarithm computation
2. **Index Generation** - Input/output permutation mutual inversion
3. **compute_indices()** - Composite function verification
4. **read_data/rewrite_data** - Ping-pong buffer stability
5. **permutate_data/repermute_data** - Data permutation correctness
6. **generate_twiddle_index()** - Index range validation
7. **Configurable_PE (NTT)** - Forward butterfly (1000× random tests)
8. **Configurable_PE (INTT)** - Inverse butterfly (1000× random tests)
9. **NTT_Kernel Round-trip** - NTT→INTT data recovery
10. **Compute_NTT Multi-limb** - Multi-limb with offset/count
11. **permute_twiddle_factors** - Twiddle access correctness

**Test Status**: All tests passing in C-simulation ✅

---

## RECOMMENDATIONS FOR CG-NTT DESIGN

### 1. Twiddle Factor Management
- **Current**: 8 complete copies in URAM (2.1 MB per limb)
- **CG-NTT**: Consider:
  - Streaming from external DDR during computation
  - Constant folding for small CG windows
  - Reduced precision intermediate twiddles
  - Compute-on-the-fly if FPGA permits

### 2. Data Layout Flexibility
- **Current**: 2D [64][64] fixed layout
- **CG-NTT**: May benefit from:
  - Different tile sizes (e.g., [8][512] for narrow PE count)
  - Rectangular buffers matching CG geometry
  - Dynamic stride calculations

### 3. Processing Element Reconfiguration
- **Current**: 32 PEs conceptual, 8 actual unroll
- **CG-NTT**: Possible adjustments:
  - Reduce PE count for smaller problems
  - Increase unroll factor for wider parallelism
  - Mixed-radix butterfly designs

### 4. Index Generation Parameterization
- **Current**: Hardcoded for SQRT=64, STAGE=12
- **CG-NTT**: Generalize:
  - Templated dimension parameters
  - Runtime dimension queries
  - CG-specific stride/offset calculations

### 5. Modular Arithmetic (Reusable)
- ✅ Barrett MultMod: Can be directly reused (II=1, flexible moduli)
- ✅ AddMod: Fully portable
- ✅ Karatsuba: Sufficient for 64-bit multiplication

---

## FILE REFERENCE QUICK LOOKUP

| Feature | File | Lines | Key Functions |
|---------|------|-------|-----------------|
| Constants | `define.h` | 1-82 | RING_DIM, SQRT, PE_PARALLEL, MAX_LIMBS |
| Arithmetic | `arithmetic.h` / `.cpp` | 32 / 108 | MultMod, AddMod, Karatsuba |
| NTT Core | `ntt_kernel.h` / `.cpp` | 154 / 424 | NTT_Kernel, Compute_NTT, Configurable_PE |
| Data Permutation | `ntt_kernel.cpp` | 20-105 | generate_input_index, compute_indices, permutate_data |
| Twiddle Access | `ntt_kernel.cpp` | 94-118 | generate_twiddle_index, permute_twiddle_factors |
| Buffering | `ntt_kernel.cpp` | 260-380 | Ping-pong logic, DEPENDENCE pragmas |
| Top Interface | `top.cpp` | 1-259 | Top(), memory init, opcode dispatch |
| Base Conv | `bconv.cpp` | 1-146 | bconv_core, Compute_BConv (systolic) |
| Data Interleave | `interleave.cpp` | ~80 | InterLeave (bit-reversal permutation) |
| Testing | `ntt_kernel_tb.cpp` | 1-708 | 11-test validation suite |

---

## NEXT STEPS FOR CG-NTT DEVELOPMENT

1. **Review Architecture**: Use Section 2-11 of `ntt_exploration_summary.md`
2. **Study Test Framework**: Run `ntt_kernel_tb.cpp` locally to understand the validation approach
3. **Clone PE Design**: `Configurable_PE()` + modular arithmetic are directly reusable
4. **Adapt Index Generation**: Modify `generate_input_index()` for CG geometry
5. **Prototype in C++**: Create CG-NTT stub, validate with testbench before synthesis
6. **Plan Memory Architecture**: Decide on twiddle storage strategy (stream vs. cached)
7. **HLS Synthesis**: Iteratively optimize pragmas for CG parameters

---

## SUMMARY STATISTICS

| Metric | Value |
|--------|-------|
| **Project LOC** | ~2,000 |
| **NTT Core LOC** | 424 |
| **Test LOC** | 708 |
| **Total BRAM (Current)** | 8334 BRAM_18K (206% budget) |
| **Total URAM (Current)** | 0 BRAM_18K (assigned but counts as BRAM for U55C) |
| **DSP (Current)** | 1010 (11% util) |
| **LUT (Current)** | 716K (54% util) |
| **FF (Current)** | 352K (13% util) |
| **NTT Stages** | 12 |
| **Butterflies per Stage** | 32 (32 PEs × 1 pair) |
| **Total Butterflies per Limb** | 24,576 (12 × 64 × 32) |
| **Twiddle Copies** | 8 (for zero conflicts) |
| **Modular Mult Latency** | 4 cycles |
| **Pipeline Target II** | 1 cycle |

---

**Report Generated**: April 15, 2026  
**Location**: `/home/timhan/FHE/openfhe-for-HKS-ACC/`  
**Files**: `ntt_exploration_summary.md`, `EXPLORATION_FINDINGS.md` (this file)

