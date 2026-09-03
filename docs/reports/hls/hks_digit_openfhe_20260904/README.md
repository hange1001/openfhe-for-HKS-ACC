# OpenFHE / OP_HKS_DIGIT integration — 2026-09-04

## Verified outcome

Real OpenFHE 1.4.0 libraries now call the actual HLS Top C-model from HYBRID
EvalKeySwitchPrecomputeCore. The admission-tested path is CKKS/DC, N=4096,
Q=3/P=2, digits [2,1]. The test is explicitly opt-in; normal builds do not
require HLS headers, GMP for the model, a board, or the new test executable.

The existing Host `MathUtils::BuildCgTwiddle_Unified` already generates the
OpenFHE negacyclic tables. The previous PoC supplied cyclic tables instead;
that was a test configuration, not a hardwired limitation of CG_NTT_Kernel.
No HLS arithmetic/routing or Host permutation change was needed.

## Acceptance evidence

- Integration test: **400 checks, 991232 exact residue comparisons**, PASS.
  Includes seven zero/boundary/impulse/random ModUp inputs, both digits and all
  five towers; original input immutability; repeated context reuse; A/B/A context
  changes; NTT and INTT individually compared against OpenFHE for all five moduli.
- Encrypted EvalRotate(+1/-1): output ciphertexts bit-exact against CPU, two
  fused digit calls each, plaintext max absolute error below 1e-6 (observed
  approximately 1e-13; encryption is randomized, see captured run).
- Encrypted EvalMult: ciphertext bit-exact and decrypts correctly, but **zero
  fused calls**: FLEXIBLEAUTOEXT reduces Q to 2, so this is a CPU-fallback test.
- Rejected before digit launch: Q=2, N=8192, P=1, MP/OC scheduling, coefficient
  input, unavailable XRT. Missing C-model callback throws immediately.
- A fake old bitstream that ignores opcode 8 is detected by output sentinels;
  injected failure on the second digit throws without publishing the first digit.
  Reconfiguration to the real Top recovers correctly.
- Original HLS Top regression: **18/18**, compound tests **22 valid cases, zero
  mismatches**, unchanged by this integration.
- XRT-enabled integration source passes a syntax compilation against installed
  XRT headers. This is not an XRT execution or link/runtime acceptance.
- AddressSanitizer full-library build: PASS with detect_leaks=1; all 400 checks
  and 991232 residue comparisons also pass with instrumented OpenFHE and HLS code.
- WITH_REDUCED_NOISE compilation branch: syntax-check PASS; fused offload is
  rejected there. This is not a reduced-noise numerical validation.

Captured test output: [normal run](integration-output.txt), [ASan run](asan-output.txt).

The existing `build/unittest/pke_tests --gtest_list_tests` aborts during test-data
loading with `std::invalid_argument: stoul`, from both repository and build
directories. It never runs a selected test. Full PKE suite acceptance is therefore
**not claimed**; the unrelated CSV-loading issue was not changed in this task.

## Reproduction

See [HKS_DIGIT.md](../../../../src/fpga_backend/HKS_DIGIT.md#openfhe-integration-2026-09-04)
for a clean build directory command. This run used Ubuntu 22.04/WSL, GCC 11.4,
NATIVEINT=64, MATHBACKEND=4, WITH_OPENMP=ON, WITH_NATIVEOPT=OFF, standard noise,
Vitis HLS 2023.2 ap_int headers, and OPENFHE_FPGA_ENABLE=OFF.

```sh
cmake --build build --target hks-digit-openfhe-test -j 4
ctest --test-dir build -R hks-digit -V
make -C src/fpga_backend test-hks-digit
g++ -std=c++17 -fopenmp -fsyntax-only -DOPENFHE_FPGA_ENABLE \
  -Ibuild/src/core -Isrc/core/include -Isrc/pke/include \
  -Ithird-party/cereal/include -I/opt/xilinx/xrt/include \
  src/pke/lib/keyswitch/hks_digit_offload.cpp
```

For ASan use a separate build with `-DCMAKE_CXX_FLAGS=-fsanitize=address`, the
same C-model options and disabled unrelated examples/benchmarks/unit tests;
run with `ASAN_OPTIONS=detect_leaks=1`. Never compare ASan or C-model elapsed
times to csynth cycles or board performance.

## Interface and measurement boundaries

The XRT transport uses separate input/metadata/output BO sizes, propagates
errors, and includes opcode 8 in the nine-entry FPGA timing table. CPU fallback
does not accidentally enter legacy fine-grained FPGA hooks. After Configure(),
Off explicitly selects CPU; applications that never call Configure are unchanged.

Initialization is reused only for identical ordered moduli and roots. It is
serialized with all digit launches. Only modulus/twiddle context is cached, not
ciphertexts or keys. No cross-digit KeyMult accumulation is implemented.

Logical input/output bytes: 416 KiB per two digits, plus 288 metadata bytes,
excluding INIT. **The current XRT implementation also uploads 320 KiB of output
sentinels** per two digits, so the earlier 416 KiB is not actual BO-sync traffic.
INIT input is 3 MiB+120 bytes; its dummy output adds synchronization overhead.
There is no measured bandwidth, latency improvement or speedup claim.

No physical board, RTL co-simulation, XRT sw_emu/hw_emu or P&R was run. Existing
HLS synthesis evidence still has -0.331 ns estimated budget slack; this task
does not close timing. It adds software integration, not new synthesis results.
