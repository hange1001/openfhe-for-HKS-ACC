# Shared NTT/INTT engine — 2026-09-04

Checkpoint before refactoring: `c073b11` (local commit; not pushed).
It includes the previous staged work and the C-model-verified OpenFHE integration.
Before committing: Top 18/18, HKS 22 cases, OpenFHE and full-library ASan passed.
One pre-existing trailing blank line in `问题回答.md` was preserved.

## Hardware change

Previously, standalone NTT/INTT and the fused HKS wrapper synthesized four
direction-specialized engines. The new Top has one call chain:

`Top → Execute_Transform_Operation → Execute_Transform → CG_Transform_Kernel`.

All three opcodes use the same `poly_buffer_2` workspace and runtime direction.
Every butterfly lane has one MultMod call, with direction-dependent operands,
addressing, twiddle selection and inverse /2 normalization. The two stage-parity
pipeline controllers share the same four MultMod instances and ping-pong memory;
they are not separate NTT/INTT engines. Compatibility template wrappers are not
reachable from Top. Neither Top ABI nor OpenFHE admission/transport changed.

The fused schedule is alpha INTTs, pre-scale/BConv, then 5-alpha NTTs. Original
EVAL towers bypass into result_buffer. Complement towers are copied into the
shared workspace before NTT and into result_buffer afterward. Compared with the
checkpoint this adds one on-chip tower copy per complement, not an AXI transfer.
There is still no inter-digit overlap or DATAFLOW.

## Same-configuration synthesis comparison

Vitis HLS 2023.2, U55C xcu55c-fsvh2892-2L-e, N=4096, PE=4, Q=3/P=2,
6 ns clock / 0.75 ns uncertainty, AXI maximum widening 512 bits.

| Resource | Checkpoint fused Top | Shared Top | Change |
|---|---:|---:|---:|
| BRAM_18K | 896 | 704 | -192 (-21.4%) |
| DSP | 2088 | 1392 | -696 (-33.3%) |
| FF | 106709 | 79656 | -27053 (-25.4%) |
| LUT | 278242 | 175140 | -103102 (-37.1%) |
| URAM | 96 | 96 | unchanged |

Three removed transform wrappers explain exactly 3×64=192 BRAM and
3×232=696 DSP saved. The shared transform including flatten/reshape has
64 BRAM / 232 DSP. Top persistent memory remains 372 BRAM / 96 URAM.
Remaining DSP: BConv 15×58=870, shared transform 4×58=232,
standalone OP_MULT 4×58=232, pre-scale 1×58=58; total 1392.
Two BConv wrapper-local scratch allocations still remain (128 BRAM each).
Forward/inverse twiddle memories are not merged by this change.

Shared core latency is 8523 cycles; including adapters, 10068 cycles per limb.
Both butterfly-loop parity bodies achieve II=1. The checkpoint forward core
was 8499 cycles; this is a small per-core scheduling cost, not an end-to-end
comparison. Additional complement copies and AXI/control must also be counted.
Top total latency is still undefined in csynth. Do not infer a speedup.

Top estimated period 5.541 ns gives budget slack 6-0.75-5.541=-0.291 ns:
**timing is not closed**. Core estimate 5.043 ns is not whole-design P&R.
All loop constraints were satisfied in this synthesis, which does not imply
that the overall clock budget is met.

## Verification

- Full native Top: 18/18; HKS 22 valid cases plus invalid inputs/canaries.
- Compatibility CG suite: 11/11, including template calls, the packed 512-bit
  wrapper across 8 limbs and a nonzero-offset subset.
- Vitis C simulation: same full regression, 0 errors.
- Real OpenFHE integration: 470 checks, 1,523,712 residues compared exactly.
  Added independent and batched transforms after INIT and every fused digit,
  including nonzero offsets and both direction transitions. Cached context is
  retained across these serialized test-only calls; no extra INIT is inserted.
- Full-library ASan with leak detection: passed the expanded integration suite.
- Generated RTL: one engine call chain, four lane MultMod instances.
- RTL smoke co-simulation: Verilog/xsim PASS. 24 transactions, 3 valid fused
  cases, 0 mismatches; includes pre-INIT/invalid-descriptor no-write checks.
  See `Top_cosim_smoke.rpt` and `vitis-cosim-smoke.txt` (C post-check PASS).
  Its aggregate min/avg/max latency mixes INIT, standalone, fused and rejected
  opcodes, so it is not a fused-digit latency or speedup measurement.
  Vitis warned that optional `zip` was unavailable; simulation and post-check
  still completed successfully with exit status 0.
- No board, XRT run, place-and-route, or measured hardware speedup.

## Reproduce

From `src/fpga_backend` in WSL/Linux (run Vitis commands sequentially):

```sh
make test-hks-digit
HKS_PROJECT_TAG=shared_v1 make hks-csim
HKS_PROJECT_TAG=shared_v1 make hks-csynth
python3 check_shared_transform.py Solution/hks_digit_csynth_shared_v1/solution1
HKS_PROJECT_TAG=shared_v1 make hks-cosim-smoke
```

The smoke test keeps real N=4096 and distinct 60-bit moduli. It compares three
digit shapes (alpha=2/start=0, alpha=1/start=2, alpha=3/start=0), standalone
transforms and fused results against the independent cyclic NTT/CRT reference.
It is not the full OpenFHE test running on RTL; OpenFHE remains C-model verified.
Tagged project directories preserve the original fused/baseline evidence.

OpenFHE and ASan build/test instructions are in the adjacent OpenFHE report and
`src/fpga_backend/HKS_DIGIT.md`. The new source is left uncommitted for review;
the original checkpoint remains intact.
