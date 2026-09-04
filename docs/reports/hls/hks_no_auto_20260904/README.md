# HKS 加速器移除 FPGA AUTO — 2026-09-04

## 范围与检查点

按用户要求，先保存删除前版本为 Git 提交 `480dc91`，再移除 FPGA AUTO。
该提交保留单套 NTT/INTT、256 位接口、去整塔拷贝、中文报告及原始验证证据，没有推送。
旧版功能验证通过，但 6 ns 物理时序未闭合；其 7 ns 结果不能直接套用到新版。

本轮不增加并行度，不改变 N=4096、Q=3、P=2 和 PE=4；保留 INIT、模加/减/乘、
NTT/INTT、BConv 与 OP_HKS_DIGIT。只融合 HKS 的单 digit ModUp，不声称已覆盖完整 HKS。

## 删除与兼容边界

- 移除 Top 的 AUTO 分支、Load_Auto_Meta、硬件实现 `auto.cpp`、头文件 `auto.h`、
  专用 HLS 测试 `auto_test.cpp`、FPGA 自同构诊断示例 `test-ckks-auto.cpp`。
- Makefile/CMake 不再编译 AUTO；主机 AutoOffload API 与 OpenFHE 中对应卸载分支一并移除。
- 保留 OpenFHE CPU 自同构和旋转功能；现有系数域算法未重写。
- 指令 7 退役且不复用，其余编号不变。Top 对旧指令不读写外部内存、不修改上下文；
  沿用无状态返回的 ABI。主机通用 Execute 对 7 在访问设备或缓冲区前抛出异常。
- 旧 AUTO 功能测试改为退役指令保护测试，随后继续执行模运算和 HKS，以检查状态保持。
  保留所有原 HKS 正确性覆盖；不是简单删除失败用例。
- 结构审计新增 `--no-auto`，要求可达 RTL 层级中无 AUTO/其元数据加载模块。

删除文件和旧接口均可从 `480dc91` 恢复；历史证据没有覆盖。新版改动暂未提交。

## 当前验证结果

- 原生 Top 18/18 通过（其中一项变为旧指令保护），HKS 22 个合法用例与边界检查通过。
- OpenFHE Release：472 项检查、1,523,712 个模数元素精确一致，含正反向 EvalRotate
  解密、EvalMult 回退，以及新增的主机旧指令拒绝检查。
- 全库 ASan（含泄漏检查）：同样 472 项检查、1,523,712 个模数元素精确一致。
- XRT 启用分支编译检查通过（退出码 0、无诊断输出），未运行 XRT 或连接板卡。
- 真实 OpenFHE fixture 的 RTL 对拍：3 次调用（INIT + 两个 digit），40,960 个元素精确一致。
- 扩展 RTL smoke：35 次调用、4 种 digit，通过退役指令、非零 tower 偏移、模运算、
  变换、HKS 与非法参数/输出哨兵检查。
- perf 和 smoke 的整个生成 Verilog 目录逐文件一致；结构审计通过，旧含 AUTO RTL
  在启用 `--no-auto` 后被正确拒绝，在未启用时仍通过原共享/端口检查。
- 新物理实现及独立 post-route 检查已完成；全部网络布通。默认 6 ns 时序通过，
  但加入 0.75 ns setup uncertainty 后未闭合；7 ns 情景通过。

同条件 HLS 对比（Vitis 2023.2、U55C、6 ns、0.75 ns 不确定度）：

| 指标 | 检查点（含 AUTO） | 移除 AUTO 后 |
|---|---:|---:|
| BRAM_18K | 708 | 688 |
| DSP | 1160 | 1160 |
| FF | 78528 | 78323 |
| LUT | 169991 | 147680 |
| URAM | 96 | 96 |
| AXI GMEM0/1/2 位宽 | 256/256/256 | 256/256/256 |
| 变换延迟（含加载/写回），周期 | 8526 | 8526 |

结构审计通过：无 AUTO，一套双向变换、四路共享模乘、整机 20 个 MultMod，
两个蝶形循环 II 均为 1。HLS 估计周期仍为 5.250 ns，不代表物理时序闭合。

### 资源到底省在哪里

总 LUT 减少 22,311（约 13.1%），并非仅删掉 `Compute_Auto` 这一个模块的 2,884 LUT：

- 顶层选择逻辑：20,045→6,837 LUT，减少 13,208。
- 全部子模块：149,236→139,998 LUT，减少 9,238；其中原 AUTO 与其元数据加载
  分别占 2,884 和 72 LUT，其余变化来自移除该访存使用者后的模块/存储接口重新综合。
- 顶层存储逻辑 LUT 增加 145，表达式 LUT 减少 10；合计正好减少 22,311。
- BRAM 减少的 20 个均来自顶层共享存储映射（372→352），不是删掉一个 20 BRAM 的
  AUTO 私有缓存。HLS 原 AUTO 模块本身 BRAM 为 0。全部子模块 BRAM 仍为 336。
- DSP 不变：AUTO 不是模乘阵列，移除它没有减少 HKS 的算术并行度。

这是 HLS 估计，不与 Vivado 的 CLB LUT 或 BRAM Tile 数混为一个口径。

布局布线后的物理资源对比如下；百分比均相对于删除前物理结果：

| 物理资源 | 检查点（含 AUTO） | 移除 AUTO 后 | 变化 |
|---|---:|---:|---:|
| CLB LUT | 107471 | 91520 | -15951（-14.84%） |
| CLB Register | 64062 | 61076 | -2986（-4.66%） |
| BRAM Tile | 359.5 | 348 | -11.5（-3.20%） |
| DSP | 1160 | 1160 | 0 |
| URAM | 96 | 96 | 0 |

可见删除 AUTO 后不只是 HLS 估计下降，物理 LUT、寄存器和 BRAM 也下降；DSP 不变，
因为 HKS 的模乘阵列没有减少。物理报告中 RAMB36E2=340、RAMB18E2=16。

## 周期与性能边界

沿用删除前同一份 OpenFHE 输入，RTL 结果如下：

| 调用 | 含 AUTO 检查点，周期 | 移除 AUTO 后，周期 |
|---|---:|---:|
| INIT（冷启动） | 295063 | 295063 |
| digit：alpha=2，start=0 | 70131 | 70131 |
| digit：alpha=1，start=2 | 69603 | 69603 |
| 两个 digit 暖态合计 | 139734 | 139734 |

当前 ModUp 从未调用 AUTO，因此移除后周期不变，INIT 退化也没有被本轮修复。
139,734 周期在 6 ns 时钟下为 0.838404 ms，但 6 ns 加 0.75 ns setup uncertainty 未闭合；
通过严格检查的 7 ns 情景为 0.978138 ms。两者都只是 RTL 周期换算，不含 PCIe、驱动和
主机打包，不能当作上板墙钟时间。本轮没有在 CAD 高负载下重新测 CPU。

### 运行时间、加速比与 slack 汇总（本轮补齐）

统一工作负载为两个 ModUp digit（alpha=2/1），不是完整 HKS，也不是 EvalRotate。
下表时间均由 RTL 周期换算，不含 PCIe/驱动；严格裕量下应看 7 ns 一列。

| 指标 | 删除前（480dc91） | 删除 AUTO 后 |
|---|---:|---:|
| 暖态两个 digit，RTL 周期 | 139734 | 139734 |
| 暖态，假定 6 ns | 0.838404 ms | 0.838404 ms |
| INIT，假定 6 ns | 1.770378 ms | 1.770378 ms |
| 冷启动 INIT + 两个 digit，假定 6 ns | 2.608782 ms | 2.608782 ms |
| 相对删除前同频加速比 | 1.000x（基准） | 1.000x |
| CPU/RTL 暖态时间比，假定 6 ns | 0.5524x | 0.5524x |
| HLS 预算余量：6 - 0.75 - 5.25 | 0.000 ns | 0.000 ns |
| 布线后默认 6 ns WNS | -0.029 ns | +0.078 ns |
| 布线后 WNS，6 ns + setup 0.75 ns | -0.779 ns | -0.672 ns |
| 布线后 WNS，7 ns + setup 0.75 ns | +0.221 ns | +0.328 ns |

CPU/RTL = 0.4631475/0.838404 = 0.5524，小于 1；即该名义 FPGA 时间约为 CPU 的
1.81 倍，不是“加速了 0.55 倍”。CPU 使用先前空闲环境的两线程基线，不是本轮重测。
按严格通过的 7 ns 换算，FPGA 为 CPU 两线程中位数的 2.11 倍。HLS 预算余量不是布线后
WNS；本轮真实布局布线证明默认 6 ns 仅勉强通过，加入 0.75 ns 工程余量即失败。

### CPU 与 FPGA 运行数据

CPU 数据来自相同 OpenFHE 参数和两个 digit 的暖态 `Precompute`，20 次预热、500 次采样；
FPGA 是真实 OpenFHE fixture 的 RTL co-sim 周期，不含板级传输。两者尚不是严格的端到端
公平加速比，只能用于定位计算瓶颈。

| 路径 | min | median | p95 | max | 备注 |
|---|---:|---:|---:|---:|---|
| CPU，OMP=1 | 0.541120 ms | 0.640583 ms | 0.963587 ms | 2.312079 ms | 实测墙钟 |
| CPU，OMP=2 | 0.410987 ms | 0.463148 ms | 0.773967 ms | 1.699541 ms | 实测墙钟 |
| FPGA RTL，6 ns | — | 0.838404 ms | — | — | 名义周期换算；严格裕量失败 |
| FPGA RTL，7 ns | — | 0.978138 ms | — | — | 严格裕量通过；仍不含主机/PCIe |

FPGA 冷启动另有 INIT 295063 周期；连同两个 digit 共 434797 周期，即 6 ns 下
2.608782 ms、7 ns 下 3.043579 ms。一次两-digit 调用的逻辑外部数据量约 416.28 KiB
（输入 96 KiB、metadata 288 B、输出 320 KiB），co-sim 没有模拟实际 PCIe 带宽与固定开销。

### FPGA 慢在哪里

下表按 HLS 子模块延迟和 RTL 总周期做归属，和为 139734 周期；其中 1353 周期作为
顶层控制/调用余量。它是定位模型，不是波形逐拍 profiler。

| 阶段 | 周期 | 占比 | 6 ns | 7 ns |
|---|---:|---:|---:|---:|
| 10 次变换，共 3 INTT + 7 NTT | 85260 | 61.02% | 0.511560 ms | 0.596820 ms |
| ├─ 蝶形计算 | 64690 | 46.30% | 0.388140 ms | 0.452830 ms |
| └─ 每次变换的本地 bank 加载/写回 | 20520 | 14.69% | 0.123120 ms | 0.143640 ms |
| QHatInv 预缩放（两个 digit） | 24624 | 17.62% | 0.147744 ms | 0.172368 ms |
| 两次 BConv | 14924 | 10.68% | 0.089544 ms | 0.104468 ms |
| 顶层输入加载 | 3117 | 2.23% | 0.018702 ms | 0.021819 ms |
| 顶层结果写回 | 10370 | 7.42% | 0.062220 ms | 0.072590 ms |
| 参数、控制与调用余量 | 1439 | 1.03% | 0.008634 ms | 0.010073 ms |
| **总计** | **139734** | **100%** | **0.838404 ms** | **0.978138 ms** |

这里的“变换”不是 AUTO。alpha=2 digit 做 2 次 INTT + 3 次 NTT，alpha=1 digit 做
1 次 INTT + 4 次 NTT，所以总计 3 INTT + 7 NTT。AUTO 已删除，而且原先也不在
`OP_HKS_DIGIT` 调用链中。

最直接的结论是：仅 10 次变换在 6 ns 下就需 0.511560 ms，已经超过 CPU 两线程整段
中位数 0.463148 ms；因此当前设计即使把 BConv 完全变成零开销，也不能超过该 CPU 基线。
第一瓶颈是变换次数与其局部 bank 搬运，第二是 QHatInv 预缩放，第三才是 BConv。
仅删掉变换内部 20520 周期的复制，模型仍约 0.715284 ms；再把预缩放四路化并融合，
约 0.604692 ms，仍慢于 CPU 两线程。这说明要得到真正加速，需进一步减少变换/访存边界、
让数据跨后续 KeyMult/ModDown 驻留，而不是只优化 BConv。

## 物理验证状态

2026-09-04 已对本轮通过对拍的 perf 工程完成 `hks-impl` 和 `hks-postroute`，目标器件
U55C，239992 条网络全部布通，无未布线网络。默认 6 ns：WNS=+0.078 ns、TNS=0、
WHS=+0.008 ns；显式加入 0.75 ns setup uncertainty 后 WNS=-0.672 ns、
TNS=-332.331 ns、1638 个失败端点；同一已布线设计改按 7 ns 检查时
WNS=+0.328 ns、TNS=0。因此当前可保守采用的仿真换算是 7 ns，不把 6 ns 宣称为有裕量闭合。

最差 setup 路径从顶层 `poly_buffer_1` 的 RAMB36E2 到 BConv 的 `local_in_x` RAMB36E2，
数据路径 5.618 ns，其中逻辑 0.907 ns（16.14%）、布线 4.711 ns（83.86%），逻辑级数为0。
这说明当前时序主要受跨 BRAM 物理距离/连线影响，而不是 `/64` 或 `%64` 的算术。
拥塞报告没有 level>5 的有效拥塞窗口、没有跨 SLR 网络；因此 256 位 AXI 接口已经可布通，
但内部共享 buffer 到 BConv 本地 buffer 的宽搬运仍是时序风险点。

原始实现、资源、默认/严格时序、关键路径、拥塞与布通报告均以 `physical-*` 名称归档在本目录。

## 预乘模乘器复用：改动前预测，尚未实施

用户本轮要求估算与数据说明，没有要求开始实现。以下为基于现有代码和报告的条件预测，
不能计入验证结果；保持四路变换硬件与 BConv 3×5 阵列不变。

### 当前可归属成本

`Prepare_HKS_BConv_Input` 每 digit 固定遍历 3×4096=12288 个位置，外加 24 周期流水/
控制开销，共 12312 周期。有效行乘 QHatInv，无效行写零，因此 alpha=1/2 时间相同。
两个 digit 共 24624 周期，占总周期 17.622%。模块 HLS 资源为 DSP58、LUT3948、FF4049、
BRAM/URAM0（使用共享缓存）。其中单个 MultMod 为 DSP58、LUT3283、FF1757；不能假定
共享后整个预乘模块的 LUT/FF 都能无代价删除。

### 两种“复用”必须分清

- **只复用一路模乘器、维持原遍历调度**：目标 DSP1160→1102（-5%），但预乘仍约
  12312 周期/digit，整核周期大致不变，且可能增加少量选择/调度开销。
- **复用现有四路模乘器，并合并 INTT 写回与预乘/compact 写入**：不新增物理算术通道，
  但让预乘利用变换阶段结束后空闲的四路资源。目标总周期约 116000～123000，
  比当前减少约 12%～17%，同频约 1.14～1.20x。这里预乘有效吞吐由一路变四路，
  不等于保持预乘吞吐不变；它保持的是硬件计算通道数量不增加。

模型依据：

```text
当前总周期                         C0 = 139734
当前独立预乘成本                   P0 = 2 × 12312 = 24624
四路但仍保留独立三行遍历（模型）   P4 ≈ 2 × (3 × 4096 / 4 + 24) = 6192
对应总周期                             ≈ C0 - P0 + P4 = 121302
预乘融合到已有四字/拍 INTT 写回：      ≈ C0 - P0 + 少量流水/调度 = 115110 + overhead
```

24 周期仅沿用当前预乘的填充量估算候选调度，并非新 RTL 测量；融合方案必须保证共享
端口和模乘均能维持 II=1。已有 INTT 写回每 tower 约 N/4 次迭代，不再重复计算整遍
写回的节省。无效行须在 BConv 装载时显式补零/屏蔽，不能直接删清零造成旧数据污染。
若必须保留独立清零或额外 bank 搬运，则更接近上面的 121k 周期情景。

仅消除当前独立预乘、其他步骤完全不变时，模型上限为 139734/115110=1.214x。
因此本项单独不会带来 2x；超过该范围的收益必须另有明确来源，不能算到“共享”名下。
同样按 6 ns 假设，工程预测暖态约 0.696～0.738 ms；对旧 CPU 基线仍只有约
0.63～0.67x 的 CPU/RTL 时间比，依然不能宣称快于 CPU。

### 资源与时序的预测边界

- DSP 的明确目标是去掉一个物理 MultMod：1160→1102。必须用生成 RTL 审计确认
  整机 MultMod 从20变19；若 HLS 复制了四路共享池，说明实现没有达到目标。
- LUT 只作低置信度预算：删除乘法器本体约3283 LUT，预留约0.5～2k LUT 给局部选择/
  调度后，预期净减少约1～3k（整核约0.7%～2%）。这不是综合值；共享接口扇出或存储
  重新映射可能吃掉收益，甚至导致 LUT 增长。FF 不给未经综合的硬指标。
- BRAM688/URAM96 暂按不变预算，不把尚未设计的缓存合并收益提前计入。
- 不预测 WNS 必然改善：删除乘法器省面积，新操作数 MUX 又会增加路径/扇出。
  性能预测假设可维持同频；例如周期减少15%而频率下降15%，时间收益恰好抵消。
- 当前 INTT 是每级蝶形做模二分，没有一遍独立的尾部 N逆元模乘可以直接“免费合并”。
  需要新增共享执行阶段/控制，不得通过简单再写四个 MultMod 调用就声称已复用。

预测依据原始报告已归档：`prepare-bconv-csynth.rpt`、`transform-csynth.rpt`；
本轮没有修改 HLS/Host 算术实现，没有启动复用版本综合，也没有重新提交 Git。

## 地址除余核查：本例没有除法器

用户指出 `Execute_Transform` 中 `k / SQRT`、`k % SQRT` 可能产生除法器。
已核对 `define.h:25–27`：SQRT 是编译期常量64，LOG_SQRT=6，k由非负循环生成且小于4096。
因此本例 `/64` 等价于 `>>6`（逻辑上的 k[11:6]），`%64` 等价于 `&63`（k[5:0]）。
这里是固定接线/位选择，不是一般的动态除法和模约减。STAGE=12也使最终A/B选择固定为A。

进一步核对本轮生成的 `Top_Execute_Transform_Pipeline_TRANSFORM_LOAD.v` 和
`Top_Execute_Transform_Pipeline_TRANSFORM_STORE.v`：二维索引被展平并映射到bank后，
地址直接使用 i[9:1]、i[9:0]、i[0] 和 limb 字段拼接。例如 STORE：

```verilog
assign lshr_ln_fu_356_p4 = {{ap_sig_allocacmp_i_4[9:1]}};
assign tmp_s_fu_399_p3 = {{limb}, {trunc_ln170_fu_396_p1}};
assign trunc_ln170_fu_396_p1 = i_4_reg_464[9:0];
```

这两个地址生成模块未见除法/取余算子或其子模块；LOAD/STORE 的 HLS 统计分别为
0DSP/628LUT/24FF 和 0DSP/327LUT/25FF。不能仅凭0DSP断言没有除法器，本结论还核查了
实际RTL的地址驱动与子模块。真正保留的代价包括complement来源MUX、bank选择与写使能，
以及每次加载1026周期、写回1026周期；它们是片上访存，不是PCIe传输。

源码随后已显式改为 `row = k >> LOG_SQRT`、`col = k & (SQRT - 1)`，并增加
`static_assert(SQRT == (1 << LOG_SQRT))`；原生 Top/HKS 回归仍为18/18、HKS 22 cases通过。
该改写帮助阅读，但不能因此宣称新省了除法器或加速。
变量除数、非二次幂常数或可能为负的有符号索引需要另外核验，不能无条件套用上述改写。
上述 P&R 使用的是改写前源码；由于原生成 RTL 已证明常量除余被实现为相同位选择，
物理数据可用于判断当前架构，但没有把文本改写冒充为重新综合后的证据。

## 复现

从 `src/fpga_backend` 运行，测试数据沿用检查点使用的同一份 OpenFHE 导出结果：

```sh
make test-hks-digit
HKS_PROJECT_TAG=no_auto_r1 \
HKS_RTL_FIXTURE="$PWD/build/hks_perf_20260904/openfhe_rtl_fixture.txt" \
make hks-cosim-perf
HKS_PROJECT_TAG=no_auto_r1 make hks-cosim-smoke
python3 check_shared_transform.py \
  Solution/hks_digit_cosim-perf_no_auto_r1/solution1 \
  --axi-width 256 --lanes 4 --total-multipliers 20 --no-auto
diff -qr Solution/hks_digit_cosim-perf_no_auto_r1/solution1/syn/verilog \
         Solution/hks_digit_cosim-smoke_no_auto_r1/solution1/syn/verilog
python3 analyze_hks_performance.py \
  --checkpoint Solution/hks_digit_cosim-perf_wide256_direct_final/solution1 \
  --shared Solution/hks_digit_cosim-perf_no_auto_r1/solution1 \
  --cpu build/hks_perf_20260904/cpu-omp2.json \
  --output ../../docs/reports/hls/hks_no_auto_20260904/comparison.json
make hks-impl HKS_IMPL_PROJECT=Solution/hks_digit_cosim-perf_no_auto_r1
make hks-postroute HKS_IMPL_PROJECT=Solution/hks_digit_cosim-perf_no_auto_r1
```

OpenFHE / ASan 从仓库根目录运行（沿用已配置的 Release 和 ASan 构建目录）：

```sh
cmake --build build --target hks-digit-openfhe-test test-fpga-modules -j2
OMP_NUM_THREADS=2 build/bin/tests/hks-digit-openfhe-test
cmake --build build-hks-asan --target hks-digit-openfhe-test -j2
OMP_NUM_THREADS=2 ASAN_OPTIONS=detect_leaks=1 \
  build-hks-asan/bin/tests/hks-digit-openfhe-test
g++ -std=c++17 -fsyntax-only -fopenmp -DBUILTIN_INFO_AVAILABLE -DPARALLEL \
  -DMATHBACKEND=4 -DOPENFHE_FPGA_ENABLE -I/opt/xilinx/xrt/include \
  -Isrc/core/include -Isrc/binfhe/include -Isrc/pke/include \
  -Ibuild/src/core -Ithird-party/cereal/include src/pke/examples/test-fpga-modules.cpp
```

## 原始证据与版本指纹

- `native-test.txt`：原生 HKS/Top；`openfhe-test.txt` / `asan-test.txt`：集成精确对拍。
- `openfhe-build.txt` / `asan-build.txt`：重编译记录；`xrt-compile-output.txt`：启用 XRT 的
  语法编译输出（成功时为空，完整命令见上文；不代表硬件运行）。
- `top-csynth.rpt` / `rtl-audit.json`：新资源与可达层级；`audit-old-version.json` /
  `audit-old-rejected.txt`：旧版通过原审计、被无 AUTO 审计拒绝的反例。
- `perf-cosim.rpt` / `perf-transactions.rpt` / `perf-vitis-output.txt`：真实输入 RTL；
  `smoke-*`：扩展 smoke 原始证据；`comparison.json`：同配置周期/资源对比。
- 两个新工程 `Top.v` SHA-256 均为
  `1ad4015311c2a0c31d454198e6e818f25e3d7631717d651d8c79b7b219f3a13b`。
- 沿用的 OpenFHE fixture SHA-256：
  `f8979c6688fb1ae78ab2d1f99119b32e05d16552cab26c3af33ffe2c0357e348`。
- 删除前检查点 `480dc91`；历史原始报告保留在 `../hks_wide256_direct_20260904/`，未覆盖。

没有板卡，RTL 周期、HLS 估计与 OOC 时序均不能替代包含 HBM/PCIe/驱动的端到端实测。
本目录已归档完成的验证证据；物理实现仍在运行，其结果尚未归档。
