# HKS digit fusion PoC — 2026-09-04

Scope: native Top single-digit ModUp fusion, no KeyMult, OpenFHE/XRT
integration, board run, P&R or completed RTL co-simulation.

## Conditions and reproducibility

- Vitis HLS 2023.2; xcu55c-fsvh2892-2L-e; 6 ns clock; 0.75 ns uncertainty.
- N=4096, Q=3, P=2, PE_PARALLEL=4, CG_BUF_PARTITION=8.
- Working checkout based on a313d9a0a32bf6fd43d08e1e1535029c78c4318d,
  including the user's pre-existing staged changes plus this PoC. Reports
  are NOT a claim about that clean commit.
- `make hks-csynth` and `make hks-baseline` in src/fpga_backend.
  The baseline disables only dispatch of OP_HKS_DIGIT using
  HKS_DISABLE_FUSED_OPCODE; unused helpers are eliminated.
- Older Solution/Top reports were not overwritten and are NOT the comparator.

## Functional evidence

- Existing Top before implementation: 12/12.
- Final native test suite: Top 18/18; HKS 22 valid cases, zero failures.
- Vitis C-sim: 0 errors, same suite (see vitis-csim-output.txt).
- AddressSanitizer: suite passes, no reported out-of-bounds or leaks.
- Independent scalar cyclic NTT and approximate CRT reference, plus separate
  opcode reference; 30/60-bit distinct primes, all legal digit ranges,
  random/zero/q-1 patterns, inactive metadata poison, stale scratch,
  invalid descriptors, pre-INIT rejection and output canaries.
- The reference matches the existing native Top cyclic/bit-reversed EVAL
  convention. OpenFHE negacyclic transform interoperability is NOT verified.

## Resource comparison

| Resource | Baseline | Fused | Delta |
|---|---:|---:|---:|
| BRAM_18K | 640 | 896 | +256 (+40.0%) |
| DSP | 1566 | 2088 | +522 (+33.3%) |
| FF | 69805 | 106709 | +36904 (+52.9%) |
| LUT | 181969 | 278242 | +96273 (+52.9%) |
| URAM | 96 | 96 | 0 |

Sources: baseline_top_csynth.rpt and fused_top_csynth.rpt, Utilization Estimates.
Top Memory remains 372 BRAM_18K / 96 URAM. The +256 BRAM_18K is added
instance-local scratch: BConv 128, NTT 64, INTT 64. Fused NTT/INTT instances
plus prescale account for the added 522 DSP. Existing BConv multipliers are
lifted/shared at Top in the fused design; simply summing both wrapper names
as two complete arrays would double-count these multipliers.

## Timing and model audit

- Estimated period: 5.581 ns for both designs. Budget slack = -0.331 ns
  after uncertainty. All loop constraints are NOT met. No P&R signoff.
- Top and Execute_HKS_Digit total latency: `?` (dynamic AXI/control loops).
  No total latency, throughput or measured speedup is reported.
- Fused Copy_HKS_Tower: 4098 cycles per tower; prescale/padding: 12312 cycles.
- Fused NTT/INTT: 10046 cycles for one limb. Fresh baseline: 9534. The
  wrapper's RESHAPE is 1026 vs 514 cycles; its memory interface matters.
- Both current designs use an 8499-cycle CG core. The historical Sep 3
  report used 14643. Thus the pre-synthesis prediction based on that
  historical report is invalid as a same-source performance comparator.
- Reconstructing only internal stages from the NEW reports gives
  `5*10046 + 5*4098 + 12312 + BConv(c)` = 90238 (alpha=2) / 90750 (alpha=1)
  cycles. These are derived stage sums, not measured whole-kernel latencies;
  they exclude AXI, metadata and control overhead. The initial estimates
  118364/118876 differ by more than 15% because the transform unit price was
  stale; no fitting coefficient was introduced.
- Exact useful-payload accounting remains 1056 -> 416 KiB for two digits,
  excluding INIT, metadata, dummy BO transfers and protocol overhead. That
  traffic saving does NOT establish a speedup without RTL/system measurement.

## Follow-ups, outside this PoC

1. Reduce cloned compute/local scratch and additional wrapper-copy overhead.
2. Run representative RTL co-simulation for cycle-level evidence.
3. Verify/adapt native versus OpenFHE negacyclic EVAL semantics before XRT
   integration; preserve the old dispatch path as a fallback.
4. Close timing and measure end-to-end behavior when a board is available.
