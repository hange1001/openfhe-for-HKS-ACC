#!/bin/bash
# 编译 test-fpga-modules host 应用
# 宏由 build/src/core/config_core.h 提供 (WITH_BE4/MATHBACKEND=4/NATIVEINT=64/OPENMP)
set -euo pipefail
ROOT=/home/timhan/FHE/openfhe-for-HKS-ACC
cd "$ROOT"

OUT=build/test-fpga-modules

g++ -std=c++17 -O2 -fopenmp \
  -DOPENFHE_FPGA_ENABLE \
  -I /home/timhan/FHE/openfhe-for-HKS-ACC/build/src/core \
  -I "$ROOT/src/core/include" \
  -I "$ROOT/src/pke/include" \
  -I "$ROOT/src/binfhe/include" \
  -I "$ROOT/third-party/cereal/include" \
  -I /opt/xilinx/xrt/include \
  "$ROOT/src/pke/examples/test-fpga-modules.cpp" \
  -L"$ROOT/build/lib" \
  -lOPENFHEpke -lOPENFHEcore -lOPENFHEbinfhe \
  -L/opt/xilinx/xrt/lib \
  -lxrt_coreutil -lxrt_swemu \
  -lm -ldl -lpthread \
  -o "$OUT"

echo "=== BUILD OK ==="
ls -la "$OUT" 2>/dev/null || true
