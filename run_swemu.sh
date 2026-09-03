#!/bin/bash
# sw_emu 仿真运行脚本：生成 emconfig -> source Vitis -> 运行 host 应用
set -euo pipefail
source /tools/Xilinx/Vitis/2023.2/settings64.sh >/dev/null 2>&1

ROOT=/home/timhan/FHE/openfhe-for-HKS-ACC
PLATFORM=xilinx_u55c_gen3x16_xdma_3_202210_1
FB=$ROOT/src/fpga_backend

# 1) 确保 emconfig.json 在 xclbin 同目录
cd "$FB/build/sw_emu"
if [ ! -f emconfig.json ]; then
  echo ">>> Generating emconfig.json (platform=$PLATFORM) ..."
  emconfigutil --platform "$PLATFORM" --nd 1 2>&1 | tail -5
  ls -la emconfig.json
else
  echo ">>> emconfig.json already present"
fi

# 2) 设置 sw_emu 仿真环境
export XCL_EMULATION_MODE=sw_emu
export LD_LIBRARY_PATH=$ROOT/build/lib:/opt/xilinx/xrt/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PATH=/opt/xilinx/xrt/bin:$PATH

# 3) 运行 host 应用（xclbin 路径由 GetXclbinPath() 自动定位到 build/sw_emu）
echo ">>> Running test-fpga-modules in sw_emu mode"
echo "    XCL_EMULATION_MODE=$XCL_EMULATION_MODE"
cd "$ROOT"
RESULT_DIR="$FB/build/sw_emu"
mkdir -p "$RESULT_DIR"
RESULT_LOG="$RESULT_DIR/simulation_result.log"
echo ">>> Simulation result will be saved to: $RESULT_LOG"
timeout 600 "$ROOT/build/test-fpga-modules" 2>&1 | tee "$RESULT_LOG"
echo ">>> EXIT_CODE=${PIPESTATUS[0]}"
