# HKS resource / latency recheck — 2026-09-04

Completed: checkpoint `c073b11` and the shared NTT/INTT working tree were both
re-synthesized and C/RTL co-simulated with the SAME OpenFHE-generated input,
constants, roots and CPU oracle. Both pass (40960 residues each).
The shared design saves area but takes 12.1% more RTL cycles on this workload.
Neither version indicates a CPU speedup under the tested nominal-clock model.
No FPGA board is available: RTL time below is a nominal-clock estimate,
not a measured end-to-end accelerator time.

## Results

| Resource | Checkpoint | Shared | Change |
|---|---:|---:|---:|
| BRAM_18K | 896 | 704 | -21.4% |
| DSP | 2088 | 1392 | -33.3% |
| FF | 106709 | 79656 | -25.4% |
| LUT | 278242 | 175140 | -37.1% |
| URAM | 96 | 96 | unchanged |

| Transaction / workload | Checkpoint cycles | Shared cycles | Checkpoint at 6 ns | Shared at 6 ns |
|---|---:|---:|---:|---:|
| INIT, excluded from warm total | 196759 | 196759 | 1180.554 us | 1180.554 us |
| digit alpha=2, start=0 | 121419 | 133827 | 728.514 us | 802.962 us |
| digit alpha=1, start=2 | 117831 | 134343 | 706.986 us | 806.058 us |
| Both digits, cached context | 239250 | 268170 | 1435.500 us | 1609.020 us |

Warm CPU median is 640.583 us with one OpenMP thread and 463.148 us with two.
CPU-time / shared-RTL-estimate is respectively 0.398x and 0.288x, **not >1x**.
Equivalently the shared nominal RTL estimate is about 3.47 times the two-thread
CPU median, before adding real host/device overhead. The older design is also
slower than this CPU baseline (ratio 0.323x). This does not predict every input
shape, machine or future hardware design, and is not a measured FPGA speedup.

Shared vs checkpoint latency rises by 28920 cycles / 173.520 us (+12.0878%).
The source schedule adds one full on-chip copy per complement: 3+4=7 copies,
each 4096 elements at approximately one element/cycle (about 28672 data cycles,
plus handshakes). This accounts for almost all of the increase; the CG core's
8499→8523 cycles per transform adds only a small scheduling cost. The previous
four engines were not executing this digit's transforms simultaneously.

The shared ten transforms including flatten/reshape account for about
10×10068=100680 cycles / 604.080 us. The remaining roughly 1004.940 us covers
AXI transfers, copies, pre-scale, BConv and control. This is not all arithmetic,
nor is all of it removable. Generated Top ports GMEM0/1/2 are each **64-bit**;
the configured 512-bit maximum widening did not make them 512-bit interfaces.

The next performance candidates are eliminating the shared-workspace extra
copy and increasing effective movement width/overlap, followed by measuring
transform parallelism. No such datapath changes were made in this experiment.
Area savings alone must not be presented as latency or CPU-speedup savings.

Both designs still miss the clock budget: estimated periods 5.581/5.541 ns,
budget slack -0.331/-0.291 ns with 0.75 ns uncertainty. Shared synthesis meets
loop II constraints; that does not mean whole-design timing closure.

## Workload and timing boundaries

- CKKS HYBRID/DC, N=4096, Q=3, P=2, digits alpha=[2,1], starts=[0,2].
- OpenFHE `EvalKeySwitchPrecomputeCore`: materialize both complete ModUp digits.
  This is not full key switching, rotation, EvalMult or an application.
- Deterministic input: integration test `Input(cc, 3)`; actual OpenFHE negacyclic
  roots, QHatInv and BConv weights. Export includes the CPU's expected results.
- CPU: Release/-O3, WITH_NATIVEOPT=OFF, OpenMP ON; 20 warmups + 500 samples.
  Timed scope includes result allocation/destruction; excludes context setup,
  fixture generation and validation. CPU/library caches are warm.
- Host: Intel Core Ultra 5 225H, WSL2 Ubuntu 22.04 (8 visible logical CPUs).
  This is an unpinned desktop microbenchmark; scheduling/thermal noise remains.
- RTL: same Vitis HLS 2023.2 / U55C / 6 ns / 0.75 ns uncertainty settings.
  Each run has exactly 3 transactions: INIT, digit alpha=2/start=0, digit alpha=1/start=2.
  Sum only the last two latency values for cached-context ModUp; list INIT separately.
  Includes generated kernel/AXI simulation behavior, not actual HBM/PCIe/platform.
- No XRT allocation, host packing/reconstruction, DMA, driver or PCIe latency
  is included in RTL time. Clock closure has not been achieved; 6 ns is an assumption.

## Measurement evidence

Fixture SHA-256 (the one-/two-thread exporters produced identical files):
`f8979c6688fb1ae78ab2d1f99119b32e05d16552cab26c3af33ffe2c0357e348`.

CPU observations from this run (microseconds, same complete two-digit ModUp):

| OpenMP max threads | Minimum | Median | 95th percentile |
|---|---:|---:|---:|
| 1 | 541.120 | 640.582 | 963.587 |
| 2 | 410.987 | 463.148 | 773.967 |

Raw 500-sample series are retained in `cpu-omp1.json` and `cpu-omp2.json`.
The shared re-synthesis repeats 704 BRAM / 1392 DSP / 79656 FF / 175140 LUT /
96 URAM. The shared OpenFHE RTL fixture passes all 40960 expected residues;
digit latencies are 133827 and 134343 cycles (268170 cycles total, 1609.020 us
at the assumed 6 ns clock). INIT is separate: 196759 cycles (1180.554 us).
The checkpoint independently passes the same fixture. Raw evidence:
`checkpoint/shared-top-csynth.rpt`, `checkpoint/shared-cosim.rpt`,
`checkpoint/shared-transactions.rpt`, `checkpoint/shared-vitis-output.txt` and
the machine-readable `comparison.json`. `shared-rtl-audit.json` reconfirms
one transform engine and four shared MultMod instances. Both logs include
successful C post-checks, not just the pre-simulation C-model run.
The complete default OpenFHE integration and full-library ASan regression
also passed after adding the optional benchmark mode.

## Reproduction commands

Build the existing `hks-digit-openfhe-test` target with C-model tests enabled.
From the repository root, create an ignored experiment directory and export:

```sh
mkdir -p src/fpga_backend/build/hks_perf_20260904/checkpoint
git archive c073b11 src/fpga_backend | tar -x -C src/fpga_backend/build/hks_perf_20260904/checkpoint
OMP_NUM_THREADS=2 build/bin/tests/hks-digit-openfhe-test --benchmark-export \
  src/fpga_backend/build/hks_perf_20260904/openfhe_rtl_fixture.txt \
  src/fpga_backend/build/hks_perf_20260904/cpu-omp2.json
OMP_NUM_THREADS=1 build/bin/tests/hks-digit-openfhe-test --benchmark-export \
  src/fpga_backend/build/hks_perf_20260904/omp1-fixture.txt \
  src/fpga_backend/build/hks_perf_20260904/cpu-omp1.json
cmp src/fpga_backend/build/hks_perf_20260904/omp1-fixture.txt \
    src/fpga_backend/build/hks_perf_20260904/openfhe_rtl_fixture.txt
```

From `src/fpga_backend`, run HLS sequentially (the directory log is shared):

```sh
HKS_PROJECT_TAG=shared_perf_r1 \
HKS_RTL_FIXTURE="$PWD/build/hks_perf_20260904/openfhe_rtl_fixture.txt" \
make hks-cosim-perf

HKS_PROJECT_TAG=checkpoint_perf_r1 \
HKS_SOURCE_ROOT="$PWD/build/hks_perf_20260904/checkpoint/src/fpga_backend" \
HKS_RTL_FIXTURE="$PWD/build/hks_perf_20260904/openfhe_rtl_fixture.txt" \
make hks-cosim-perf

python3 analyze_hks_performance.py \
  --checkpoint Solution/hks_digit_cosim-perf_checkpoint_perf_r1/solution1 \
  --shared Solution/hks_digit_cosim-perf_shared_perf_r1/solution1 \
  --cpu build/hks_perf_20260904/cpu-omp1.json \
  --cpu build/hks_perf_20260904/cpu-omp2.json
```

The shared working tree is not replaced or committed by this experiment.
