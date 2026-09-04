# OP_HKS_DIGIT prototype

当前版本已移除 FPGA AUTO：编号 7 退役，其他指令编号与 256 位接口不变。
OpenFHE CPU 自同构仍保留。删除前检查点为 `480dc91`，本轮结果以
[移除 AUTO 报告](../../docs/reports/hls/hks_no_auto_20260904/README.md) 为准；下文旧资源/时序章节保留其历史口径。

Scope: one digit of ModUp, not a complete key switch. OpenFHE integration now
exists with an opt-in HLS C-model / XRT transport. The C-model is verified;
XRT is compile-checked only. No KeyMult fusion, cross-call polynomial cache,
or throughput claim is included. Modulus/twiddle initialization is cached.

## Interface contract

- `opcode = OP_HKS_DIGIT (8)`; the Top signature is unchanged.
- The optimized Top uses three 256-bit AXI data interfaces. Device buffer base
  addresses must be 32-byte aligned; the uint64 word layout and logical lengths
  are unchanged. XRT uses device BOs, not raw CPU vector addresses. C-model
  pointers retain ordinary uint64 alignment. Narrow/unaligned metadata offsets
  remain handled by HLS's AXI adapter; no host-side ap_uint packing is required.
- `num_active_limbs = alpha` (1..LIMB_Q), `mod_index = digit_start`
  (0..LIMB_Q-alpha). Call OP_INIT first with consistent Q/P parameters.
- `mem_in1`: alpha*RING_DIM uint64 words, compact digit-local rows in the
  EVAL ordering determined by OP_INIT's twiddle tables. For OpenFHE, use the
  exact context roots and MathUtils::BuildCgTwiddle_Unified (negacyclic tables).
  No Host bit-reversal or twist is then required.
- `mem_in2`: LIMB_Q*MAX_OUT_COLS weights (row-major), followed by LIMB_Q
  QHatInv values. Only the first alpha rows and first MAX_OUT_COLS-alpha
  columns are active. Inactive entries are ignored. Weight column p denotes
  global index `p < digit_start ? p : p+alpha`; inverse row q belongs to
  global Q index digit_start+q.
- `mem_out`: MAX_OUT_COLS*RING_DIM words in Q0,Q1,Q2,P0,P1 order. Original
  digit EVAL towers are unchanged; complement towers are transformed CRT outputs.
- Residues/active constants must be canonical. Moduli must be odd, greater
  than one and below 2^62 with matching Barrett/twiddle tables from OP_INIT.
- Invalid digit ranges or uninitialized/unsupported modulus slots produce
  no output write. There is no new status port; validate descriptors on Host.

## Device sequence

1. Load into poly_buffer_2 and broadcast original EVAL into result_buffer in
   the same pass (no separate bypass copy).
2. INTT with each source tower's global modulus/twiddle index.
3. Multiply by QHatInv modulo the source qi; compact into poly_buffer_1,
   explicitly zeroing inactive BConv rows.
4. BConv into poly_buffer_1 rows LIMB_Q..LIMB_Q+complement_count-1.
5. Load each complement directly from poly_buffer_1 into the shared transform
   banks; NTT, then write directly to its global result_buffer slot.
6. Store the complete result once.

Reuses all three existing Top polynomial buffers. No additional persistent
ring-sized array. Internal BConv/NTT scratch and bank overhead still count
toward area. AXI loops remain in non-inlined helpers. This is sequential
fusion, not DATAFLOW or cross-digit output-stationary scheduling.

Since the `c073b11` checkpoint, OP_NTT, OP_INTT and OP_HKS_DIGIT all enter
`Execute_Transform_Operation`. Its single transform call uses runtime direction,
one four-lane CG engine and one set of ping-pong buffers. No flat/flat_out
intermediate arrays or Copy_HKS_Tower passes remain in Top. A template
wrapper is retained for direct-kernel tests, but is not reachable from Top.
Forward/inverse twiddle memories and the two BConv scratch allocations remain.
The standalone transform wrapper remains compatible. The optimized generated
RTL hoists four shared MultMod units to Top and also uses that pool for OP_MULT;
zero DSP in the CG submodule report therefore does NOT mean zero multipliers.
All opcodes are still sequential; PE_PARALLEL=4 and the BConv array are unchanged.

Fixed-length, non-inlined single-tower transfer helpers are intentional. Keeping
the runtime limb loop outside this boundary lets Vitis infer 256-bit bursts;
merely unrolling a variable-length flattened transfer produced 64-bit traffic.

## Board-free validation

Compare fused output against separate OP_INTT, scalar QHatInv multiplication,
OP_BCONV, OP_NTT calls AND an independent scalar cyclic-NTT/CRT reference.
Use distinct primes, all valid digit ranges, random seeds, zero/boundary data,
stale scratch, invalid descriptors and original Top regression. Keep synthesis
outputs separate from baseline. The original cyclic test alone does not establish
OpenFHE interoperability. The new integration test separately compares forward
and inverse transforms against OpenFHE, then full ModUp and encrypted rotations.

## Payload prediction (before measurement)

Let a=alpha, c=5-a and B=RING_DIM*8=32 KiB. Exclude OP_INIT, metadata,
Host memcpy, dummy inputs/BO padding, AXI protocol overhead and launch costs.
These are useful Top payload bytes, NOT measured PCIe traffic.

- Separate: (2*a + LIMB_Q + 3*c)*B (BConv pads inputs to 3 rows).
- Fused: (a+5)*B.
- alpha=2: 512 KiB -> 224 KiB; alpha=1: 544 KiB -> 192 KiB.
- Two digits: 1056 KiB -> 416 KiB, saving 640 KiB useful payload.
- One transform call per limb: 12 calls -> 2 for two digits. A limb-batched
  comparison has a different launch count (3 per digit).

Exact interface accounting is not a latency prediction. C-sim verifies
functionality only; csynth estimates are not P&R timing or board performance.

## Reproduce

From src/fpga_backend in the WSL/Linux environment:

```sh
make test-hks-digit   # g++ with the actual HLS ap_int headers, no board
make hks-csim         # Vitis HLS C simulation
make hks-csynth       # standalone experimental synthesis
make hks-baseline     # same checkout/settings, fused opcode disabled
# Optional, potentially long; not part of the completed C-sim acceptance:
make hks-cosim
# Shared-engine projects, preserving the pre-sharing experiment directories:
HKS_PROJECT_TAG=shared_v1 make hks-csim
HKS_PROJECT_TAG=shared_v1 make hks-csynth
python3 check_shared_transform.py Solution/hks_digit_csynth_shared_v1/solution1
HKS_PROJECT_TAG=shared_v1 make hks-cosim-smoke
# Optimized interface / direct-bank version (use new project tags):
HKS_PROJECT_TAG=wide256_direct_final make hks-cosim-smoke
python3 check_shared_transform.py Solution/hks_digit_cosim-smoke_wide256_direct_final/solution1 \
  --axi-width 256 --lanes 4 --total-multipliers 20
make hks-impl HKS_IMPL_PROJECT=Solution/hks_digit_cosim-smoke_wide256_direct_final
make hks-postroute HKS_IMPL_PROJECT=Solution/hks_digit_cosim-smoke_wide256_direct_final
```

Override HLS_ROOT if Vitis HLS 2023.2 is not installed under
`/tools/Xilinx/Vitis_HLS/2023.2`. HLS_PART uses the existing Makefile setting
(U55C); the experiment uses 6 ns and 0.75 ns uncertainty. HLS projects go to
`Solution/hks_digit_{csim,csynth,cosim,baseline}`, leaving `Solution/Top` untouched.
Set HKS_PROJECT_TAG to append a unique suffix to a new experiment directory.
Run these HLS targets sequentially; Vitis shares its working-directory log.

Initial verification on 2026-09-04: native g++, AddressSanitizer and Vitis
C-sim all pass. Top: 18/18; HKS: 22 valid cases with zero mismatches. Invalid
descriptors, pre-INIT calls and output guards are checked separately.

## Pre-sharing synthesis result (checkpoint c073b11, 2026-09-04)

Both fused and freshly rebuilt opcode-disabled baseline generate RTL under
Vitis 2023.2 / U55C / 6 ns / 0.75 ns uncertainty. Archived evidence and limits:
[report summary](../../docs/reports/hls/hks_digit_poc_20260904/README.md).

| Resource | Same-source baseline | Fused Top |
|---|---:|---:|
| BRAM_18K | 640 | 896 |
| DSP | 1566 | 2088 |
| FF | 69805 | 106709 |
| LUT | 181969 | 278242 |
| URAM | 96 | 96 |

Both estimated clock periods are 5.581 ns: budget slack is
6-0.75-5.581 = -0.331 ns. Loop constraints are not all satisfied. This is
**RTL generation success, not timing closure**. Total Top/digit latency is
undefined in csynth due to dynamic control/AXI loops; do not quote a total
cycle count or acceleration factor from these reports.

The extra 256 BRAM_18K are instance-local scratch (128 for BConv, 64 each for
NTT/INTT); the top-level polynomial/twiddle memory allocation is unchanged.
Reuse at the source-buffer level does not imply zero additional hardware:
HLS creates additional transform/BConv wrapper instances. Area sharing and
timing closure were follow-ups at this checkpoint. OpenFHE compatibility is tested
with its own roots and the existing Host negacyclic table generator.

## Shared-engine synthesis result (2026-09-04)

Historical result before the direct-bank / wide-interface optimization below.

Same device/clock/PE settings: Top now uses 704 BRAM_18K, 1392 DSP,
79656 FF, 175140 LUT and 96 URAM. Compared with the fused checkpoint this
removes 192 BRAM and 696 DSP: exactly three former transform instances.
The generated RTL audit confirms one CG engine with four shared MultMod lanes.
Both butterfly parity-loop bodies achieve II=1; core latency is 8523 cycles
(10068 including flatten/reshape). Whole-Top latency remains undefined in csynth.
Top estimated period is 5.541 ns, so budget slack is still -0.291 ns, not closure.
The extra on-chip complement copy must be included in any later performance model.
Expanded OpenFHE/ASan regression passes 470 checks / 1,523,712 exact residues.
Verilog/xsim smoke co-simulation also passes: 24 transactions, three valid
fused digit shapes, independent cyclic reference and no-write checks.
Full evidence and RTL co-simulation status:
[shared-engine report](../../docs/reports/hls/hks_shared_transform_20260904/README.md).

## Matched-input resource / performance experiment

`hks-digit-openfhe-test --benchmark-export fixture.txt cpu.json` measures warm
CPU ModUp (both digits, including result allocation/destruction), verifies the
C-model against CPU, and exports actual OpenFHE inputs/constants plus the CPU
oracle. This is separate from the default correctness suite.

`make hks-cosim-perf` replays exactly INIT + alpha=2/start=0 + alpha=1/start=2
against RTL. Set HKS_RTL_FIXTURE to the absolute exported file path with basename
`openfhe_rtl_fixture.txt`. Set HKS_SOURCE_ROOT to an isolated source snapshot to
compare an older revision without changing the working tree; use a unique
HKS_PROJECT_TAG. `analyze_hks_performance.py` compares the two synthesis summaries,
per-transaction RTL cycles and CPU samples. It excludes INIT from the warm
two-digit total and explicitly labels clock, PCIe and host-overhead assumptions.
Evidence and exact commands:
[matched-input performance report](../../docs/reports/hls/hks_performance_20260904/README.md).

Measured RTL-cycle result for the same OpenFHE two-digit fixture: checkpoint
239250 cycles versus shared 268170 cycles (+12.1%), both outputs bit-exact.
At the unclosed 6 ns target this is 1.4355 versus 1.6090 ms, excluding INIT.
Warm CPU median on the tested Core Ultra 5 225H/WSL2 is 0.4631 ms (OpenMP=2).
Thus area sharing is successful but this configuration does not indicate a
CPU acceleration benefit. Extra complement copies explain most of the shared
latency increase; no real PCIe/host overhead or board speedup was measured.

## OpenFHE integration (2026-09-04)

The hook is in `KeySwitchHYBRID::EvalKeySwitchPrecomputeCore`, before CPU digit
materialization. It directly packs EVAL towers, the real PartQlHatInvModq and
PartQlHatModp constants, and reconstructs five-tower DCRTPoly outputs. On failure
it never publishes partial digits. The default path is unchanged unless an
application explicitly calls `ConfigureHKSDigitBackend`.

Current admission: CKKS, DC, N=4096, **current** Q=3, P=2, alpha=1..3,
64-bit native integers, standard unsigned ApproxSwitchCRTBasis. The integration
suite specifically validates alpha=2 and the final one-limb digit. Unsupported
shape/scheme/strategy falls back to CPU before a digit launch. In FLEXIBLEAUTOEXT
at depth 1, EvalRotate uses the fused path, while EvalMult drops to Q=2 and falls
back. Do not report that multiplication as a fused execution.

Run from the repository root (no FPGA or XRT needed, Vitis headers required):

```sh
cmake -S . -B build-hks-check \
  -DOPENFHE_FPGA_ENABLE=OFF -DOPENFHE_HKS_CMODEL_TESTS=ON \
  -DGIT_SUBMOD_AUTO=OFF -DBUILD_UNITTESTS=OFF \
  -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARKS=OFF
cmake --build build-hks-check --target hks-digit-openfhe-test -j 4
ctest --test-dir build-hks-check -R hks-digit --output-on-failure
```

Set `-DHKS_HLS_ROOT=/path/to/Vitis_HLS/2023.2` if needed. This executable links
the actual HLS Top/ap_int code and actual OpenFHE libraries, not a duplicate
software implementation of the compound opcode. It is **not** XRT sw_emu or
C/RTL co-simulation. The original integration did not change HLS arithmetic;
the subsequent shared-engine refactor is independently revalidated as above.

Application API (`keyswitch/hks_digit_offload.h`, namespace lbcrypto):

```cpp
// Board-free test application: link hks_digit_cmodel and declare the Top ABI.
ConfigureHKSDigitBackend(HKSDigitBackend::CModel, Top);
// Hardware application: requires OPENFHE_FPGA_ENABLE=ON and a rebuilt opcode-8 xclbin.
ConfigureHKSDigitBackend(HKSDigitBackend::XRT);
// CPU reference; after explicit configuration the legacy per-op hooks stay disabled.
ConfigureHKSDigitBackend(HKSDigitBackend::Off);
```

These are alternative selections, not three calls to execute consecutively.
Configure only while idle. The fused context is mutex-serialized and caches the
full ordered modulus/root vectors. A context change or transport failure
invalidates the cache. Do not interleave external Top/InitModuli calls; the fused
backend owns the device context. Key generation, KeyMult, ModDown and unsupported
configurations remain CPU operations. Applications that never configure this
API retain the legacy behavior.

`GetHKSDigitStats()` distinguishes initialization, completed digit calls,
fallbacks and failures. Its byte counters are logical opcode-8 payloads only:
two digits transfer 416 KiB plus 288 metadata bytes. INIT adds 3 MiB+120 input
bytes per context initialization. The current diagnostic XRT transport additionally
sends 320 KiB of output sentinels per two digits, to detect an old/no-op bitstream;
thus its buffer-sync bytes are **not** the 416 KiB useful-payload prediction.
XRT timings go into FpgaTransferStats; C-model runtime is never labeled PCIe or
hardware timing. A future status ABI could eliminate sentinel uploads.

## Direct-bank / 256-bit interface validation

The latest optimization keeps arithmetic parallelism unchanged, removes the
17 explicit two-digit tower-copy passes and the flat adapters, and actually
generates 256-bit GMEM0/1/2 ports. Reference results, exact reproduction commands,
RTL-cycle estimates and board-free implementation limits are recorded in
[the wide-interface report](../../docs/reports/hls/hks_wide256_direct_20260904/README.md).
Software and RTL checks pass; OOC routing completes with zero routing errors.
The original 6 ns target fails (-0.779 ns WNS with 0.75 ns user setup uncertainty).
The same routed checkpoint passes a 7 ns internal-path timing scenario
(+0.221 ns WNS), but external I/O and platform clock timing are not signed off.
The repository's 6 ns target is unchanged. INIT regresses to 295063 cycles;
warm two-digit latency improves to 139734 cycles. Neither is a board measurement.

Evidence and remaining limits:
[OpenFHE integration report](../../docs/reports/hls/hks_digit_openfhe_20260904/README.md).
