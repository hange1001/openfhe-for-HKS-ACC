# OpenFHE HKS Accelerator

基于 [OpenFHE](https://github.com/openfheorg/openfhe-development) 的研究型 fork，面向
RNS-CKKS Hybrid Key Switching（HKS）的软硬件协同加速。项目包含 Xilinx Alveo
U55C 上的 Vitis HLS kernel、OpenFHE Host/XRT 集成，以及用于比较 Digit-Centric
（DC）、Max-Parallel（MP）和 Output-Centric（OC）数据流的软件仿真器。

> 本仓库不是 OpenFHE 官方发行版。当前实现是研究原型，接口参数和 FPGA 构建流程
> 仍针对特定实验配置。

## 项目内容

- **HKS FPGA kernel**：实现 NTT/INTT、模运算、BConv 和融合的
  `OP_HKS_DIGIT` 数据路径。
- **Host 集成**：通过 XRT 管理设备、buffer、kernel 调用和 CPU fallback。
- **共享 ABI**：Host 与 HLS 使用同一份协议常量，启动时检查 ABI 版本和硬件能力。
- **三种 HKS 策略**：DC、MP、OC 软件路径及操作数、内存峰值和耗时统计。
- **事件驱动仿真器**：在相同资源约束下比较数据流、output tiling、容量和调度策略。
- **功能验证**：包含 HLS testbench、OpenFHE 集成测试和真 OC reference test。

当前 FPGA 融合边界是单个 digit 的 ModUp：

```text
OpenFHE EvalRotate / KeySwitch
  -> HKS precompute
  -> TryHKSDigitOffload
  -> XRT
  -> Top(OP_HKS_DIGIT)
  -> INTT -> pre-scale -> BConv -> NTT
```

Automorphism 仍在 CPU 路径执行；仓库目前没有宣称完整 HKS 全链路都已卸载到 FPGA。

## 代码结构

| 路径 | 作用 |
|---|---|
| [`src/fpga_backend/`](src/fpga_backend/) | HLS kernel、testbench、Tcl 和 Vitis/Make 构建入口 |
| [`src/fpga_backend/src/top.cpp`](src/fpga_backend/src/top.cpp) | FPGA 顶层函数 `Top` 和 opcode 分发 |
| [`src/core/include/FpgaManager.h`](src/core/include/FpgaManager.h) | OpenFHE Host/XRT 桥接层 |
| [`src/core/include/fpga_hks_protocol.h`](src/core/include/fpga_hks_protocol.h) | Host/HLS 共享 ABI |
| [`src/pke/include/keyswitch/hks_strategy.h`](src/pke/include/keyswitch/hks_strategy.h) | DC/MP/OC 选择、统计和内存跟踪 |
| [`src/pke/lib/keyswitch/hks_digit_offload.cpp`](src/pke/lib/keyswitch/hks_digit_offload.cpp) | 融合 digit offload 入口 |
| [`src/pke/examples/hks-benchmark.cpp`](src/pke/examples/hks-benchmark.cpp) | HKS 策略 benchmark |
| [`src/pke/unittest/UnitTestHKSTrueOC.cpp`](src/pke/unittest/UnitTestHKSTrueOC.cpp) | 真 OC 功能 reference 和一致性测试 |
| [`tools/hks_flow_sim/`](tools/hks_flow_sim/) | DC/MP/OC 资源约束仿真器 |

## 环境要求

CPU 构建：

- Linux
- CMake 3.16.3 或更高版本
- GCC 9+ 或 Clang 10+
- 支持 C++17 的工具链

FPGA 流程额外需要：

- Xilinx Vitis HLS 2023.2
- XRT 2.16
- Alveo U55C 平台文件（默认
  `xilinx_u55c_gen3x16_xdma_3_202210_1`）

软件仿真器使用 Python 3；读取仓库中的 YAML 配置需要
[`PyYAML`](https://pyyaml.org/)。

## 快速开始

### 1. CPU-only 构建

```bash
cmake -S . -B build \
  -DOPENFHE_FPGA_ENABLE=OFF \
  -DBUILD_EXAMPLES=ON \
  -DBUILD_UNITTESTS=ON

cmake --build build \
  --target hks-benchmark hks_true_oc_tests \
  -j"$(nproc)"
```

运行三种策略中的任意一种：

```bash
./build/bin/examples/pke/hks-benchmark --strategy DC --iters 10
./build/bin/examples/pke/hks-benchmark --strategy MP --iters 10
./build/bin/examples/pke/hks-benchmark --strategy OC --iters 10
```

在 CPU-only 构建中，benchmark 会自动使用软件路径。

### 2. 功能测试

```bash
./build/unittest/hks_true_oc_tests
```

该 target 验证 DC/MP/OC 结果一致性、单目标塔 BConv 与全宽 BConv 的一致性，以及
端到端 CKKS 解密正确性。

### 3. 数据流仿真器

```bash
python3 -m pip install pyyaml

cd tools/hks_flow_sim
python3 hks_sim.py --config configs/p4_reproduction.yaml
python3 -m unittest discover -s tests -t .
```

组合基线和 KeyMult 配置进行扫描：

```bash
python3 hks_sim.py \
  --config configs/p4_reproduction.yaml \
  --config configs/keymult_nominal.yaml \
  --boundary full_hks \
  --invocation per_hks \
  --oc-tile-sweep
```

仿真器报告的是资源约束下的周期模型。默认结果不包含 PCIe、驱动、buffer 分配和
Host 调度开销，因此不能直接当作板卡端到端实测时间。

## FPGA 构建

先完成一次 OpenFHE Host 构建，以生成 `build/src/core/config_core.h`，然后加载
Vitis HLS 环境：

```bash
source /tools/Xilinx/Vitis_HLS/2023.2/settings64.sh

make -C src/fpga_backend csim \
  MODULE=Top \
  OPENFHE_ROOT="$PWD/"

make -C src/fpga_backend csynth \
  MODULE=Top \
  OPENFHE_ROOT="$PWD/"
```

可用的 HLS 顶层模块包括：

- `Top`
- `Compute_CG_NTT`
- `NTT_Kernel`
- `Compute_NTT`
- `Compute_BConv`
- `Compute_BConv_Naive`
- `Compute_BConv_Systolic`

生成 Vitis 构建产物：

```bash
make -C src/fpga_backend sw_emu OPENFHE_ROOT="$PWD/"
make -C src/fpga_backend hw_emu OPENFHE_ROOT="$PWD/"
make -C src/fpga_backend hw     OPENFHE_ROOT="$PWD/"
```

`OPENFHE_FPGA_ENABLE` 默认开启；如果 CMake 在
`/opt/xilinx/xrt/{include,lib}` 找不到 XRT，会给出警告并自动关闭 FPGA 支持。

## 当前约束

- 共享 ABI 当前固定为 `N=4096`、`Q=3`、`P=2`，修改参数时必须同步验证
  Host、HLS kernel 和 xclbin。
- 默认器件是 `xcu55c-fsvh2892-2L-e`，目标时钟周期为 6 ns。
- FPGA 运行需要与源码和 ABI 匹配的本地 xclbin；仓库不包含预编译 bitstream。
- `OP_HKS_DIGIT` 聚焦融合 ModUp，尚未覆盖完整 HKS 和 Automorphism。
- 软件仿真用于架构探索，不替代真实板卡上的 wall-clock 测量。

## 上游项目与许可证

本项目派生自 OpenFHE。OpenFHE 的安装、API 和算法文档请参考：

- [OpenFHE 官方文档](https://openfhe-development.readthedocs.io/)
- [OpenFHE GitHub 仓库](https://github.com/openfheorg/openfhe-development)
- [OpenFHE 设计论文](https://eprint.iacr.org/2022/915)

代码沿用仓库中的 [BSD 2-Clause License](LICENSE)。使用或引用 OpenFHE 时，请同时
遵循上游项目的许可证与引用要求。
