# HKS 直接访问片上存储与 256 位接口优化报告 — 2026-09-04

验证负载：OpenFHE CKKS/DC 的 ModUp，N=4096、Q=3、P=2，两个 digit 的大小为 [2,1]。

实现与无板卡验证已完成：软件和 RTL 功能检查通过，布局布线成功，但原定的 6 ns
时序目标未通过。同一份布线结果在 7 ns 时钟周期、0.75 ns 用户建立时间不确定度下，
内部路径时序通过。该结果不等于整个平台时序签核，也不是端到端 FPGA 实测性能。

本文正文使用中文；代码标识、复现命令、哈希及工具生成的原始报告保留原文，便于核查。

## 一、改动与保持不变的条件

- PE_PARALLEL=4、四路蝶形和 3×5 BConv 阵列保持不变。
- 独立 NTT/INTT 指令与融合指令继续共用一套运行时可切换方向的变换引擎。
- 初次 AXI 加载时，将原始 EVAL 数据同时写入旁路结果位置，避免后续再单独拷贝。
- 变换加载直接从输入缓存或 BConv 输出读入 bank A；写回直接到系数结果或最终 EVAL
  结果位置。仅在局部选择两个固定的数据源/目标，没有新增通用交叉开关。
- 移除两个 digit 中全部 17 次显式 Copy_HKS_Tower，不只是最初共享工作区新增的 7 次。
- 移除 flat/flat_out 中间数组。NTT A/B 双缓冲和 BConv 行/列临时缓存仍保留，
  用于满足既有访存带宽和数据顺序要求；并不是所有数据搬运都消失了。
- 生成 RTL 中的 GMEM0/1/2 均为 256 位。uint64 ABI 和字节排列不变；设备端缓冲区
  基地址须按 32 字节对齐。
- 当前 Vitis 版本需要明确的固定长度单塔传输函数，才能推导出自动突发拓宽。
  两个仅调整循环的中间版本未达到目标位宽，不计为最终验收通过的接口实现。

## 二、HLS 资源对比

对比条件一致：Vitis HLS 2023.2、U55C、目标周期 6 ns、不确定度 0.75 ns。
基线是上一版共享 NTT/INTT 引擎的方案，不是更早的多套变换引擎提交。

| 指标 | 上一版共享引擎 | 本轮优化后 |
|---|---:|---:|
| BRAM_18K | 704 | 708 |
| DSP | 1392 | 1160 |
| FF | 79656 | 78528 |
| LUT | 175140 | 169991 |
| URAM | 96 | 96 |
| AXI GMEM0/1/2 | 64/64/64 | 256/256/256 |
| 变换延迟（含加载/写回适配），周期 | 10068 | 8526 |
| HLS 估计时钟周期，ns | 5.541 | 5.250 |
| HLS 预算时序裕量，ns | -0.291 | 0.000 |

资源变化来源：变换缓冲从 64 减到 32 BRAM，三个 AXI 适配器合计从 12 增到 48 BRAM，
所以净变化为 -32+36=+4 BRAM。Top 的常驻存储仍为 372 BRAM 和 96 URAM，两份 BConv
临时缓存也仍保留。此段 BRAM 均沿用上表的 BRAM_18K 口径。

Vitis 还把四个变换模乘单元提升到 Top，与独立 OP_MULT 共享，减少 4 个 MultMod，
对应减少 232 DSP；不是提高了算术并行度。整个 Top 有 20 个 MultMod：BConv 15 个、
NTT/INTT/MULT 共用 4 个、预缩放 1 个。CG 子模块报告显示 DSP=0，是因为这些物理
乘法器被统计在上层，不能理解为变换不需要乘法器。

## 三、软件与 RTL 结构验证

- 原生测试：Top 18/18；HKS 22 个合法用例，以及非法描述符与边界哨兵检查均通过。
- 新增 start=1、连续两个塔的 ADD/SUB/MULT/AUTO 回归：使用各通道不同、接近模数的
  数据，并检查输出保护区，覆盖三个拓宽后的 AXI 接口。
- OpenFHE Release 和 ASan（含 detect_leaks=1）各通过 470 项检查，
  1,523,712 个模数元素逐项精确一致。
- 最终 Vitis 全量 C-sim：0 errors，覆盖 Top/HKS/宽端口回归。
- 旧独立 CG 包装器：11/11 通过。
- RTL 结构：一套变换引擎、四路接口对应四个物理共享模乘单元；端口均为 256 位，
  不再存在 Copy_HKS_Tower 模块。
- 原生程序和 ASan 的运行时间只用于正确性验证，不能作为 FPGA 时延。
- 扩展 RTL smoke：PASS，共 35 次调用、4 个融合 digit 用例，包括 alpha=2/start=1、
  全部合法 alpha 值、方向切换、60 位素数、非法调用、边界哨兵以及上述四种端口回归。
  其中独立标量参考使用循环变换；真实 OpenFHE 的负循环测试数据另行对拍，二者不混用。
- 两个奇偶阶段蝶形循环的 II 均为 1，每次调用延迟 537 周期。

## 四、同输入 RTL 性能对比（功能对拍通过）

最终性能协同仿真使用真实 OpenFHE CPU 结果作参考，逐项核对 40,960 个模数元素。
重放的仍是此前导出的同一份测试数据，其 SHA256 为：
f8979c6688fb1ae78ab2d1f99119b32e05d16552cab26c3af33ffe2c0357e348。

| 调用 | 上一版周期 | 优化后周期 | 优化后按 6 ns 换算 |
|---|---:|---:|---:|
| INIT（不计入暖态总计） | 196759 | 295063 | 1770.378 us |
| alpha=2, start=0 | 133827 | 70131 | 420.786 us |
| alpha=1, start=2 | 134343 | 69603 | 417.618 us |
| 两个暖态 digit 合计 | 268170 | 139734 | 838.404 us |

暖态内核周期减少 47.8935%，相对上一版 FPGA 设计的同频周期比为 1.919 倍。
此前 CPU 双线程基线的中位时间为 463.1475 us；CPU 时间除以上表名义 RTL 时间为
0.5524 倍，仍未达到相对 CPU 的加速。以上不含 PCIe、XRT、真实 HBM 和主机打包开销。
由于 6 ns 目标未通过布局布线时序，上表仅用于同频换算，不能当作已验证的工作点。

INIT 退化了 98,304 周期，按 6 ns 为 589.824 us。该冷启动成本单列，没有藏入或忽略在
暖态指标中；本轮尚未优化初始化循环的突发传输适配。若只初始化一次、只计算一组两个
digit，则名义总时间为 2.608782 ms，原来为 2.789574 ms，仅减少约 6.48%。因此是否复用
已初始化的上下文，会明显影响收益。

10 次变换连同直接加载/写回共需 85,260 周期；其余 54,474 周期包含预缩放、BConv、
外部 AXI 与控制开销，不能全当作冗余拷贝。计算并行度保持不变。

## 五、最终独立模块物理验证（OOC）

Vivado 2023.2 已在 U55C 上完成同一份生成 RTL 的布局布线。239,992 条可布线网络
全部布通，路由错误为 0。注意：即便时序不通过，实现工具进程仍可能返回 0，
不能仅凭进程退出码判断时序验收是否通过。

布线后资源为：CLB LUT 107,471（8.24%）、FF 64,062（2.46%）、DSP 1,160（12.85%）、
BRAM 359.5 个存储块（17.83%；358 个 RAMB36 加 3 个 RAMB18，等效 719 个 BRAM_18K），
以及 URAM 96（10%）。这些是 Vivado 物理实现统计，不是前文的 HLS 估算。
尚未对上一版进行同条件 P&R，因此不声称得到旧版与新版的物理资源差值。

| 已布线设计的时序情景 | 建立 WNS，ns | 建立 TNS，ns | 建立时间失败端点数 | 保持 WHS，ns |
|---|---:|---:|---:|---:|
| 6 ns，导出时的默认约束 | -0.029 | -0.029 | 1 | +0.019 |
| 6 ns，用户建立时间不确定度 0.75 ns | -0.779 | -1574.646 | 6083 | +0.019 |
| 7 ns，用户建立时间不确定度 0.75 ns | +0.221 | 0.000 | 0 | +0.019 |

默认导出约束包含工具的抖动估计，但缺少显式用户不确定度。补上 0.75 ns 后，最差
建立路径的总不确定度为 0.785 ns。三种情景均没有保持时间违例。7 ns 情景仅在读取
同一个布线检查点后调整时钟约束，没有修改 RTL、重新综合、重新布线，也没有修改
仓库配置中的 6 ns 目标。

在已检查的 7 ns 内部时序情景（142.857 MHz）下，两个暖态 digit 按周期换算为
0.978138 ms，INIT 为 2.065441 ms。相对此前 0.4631475 ms 的 CPU 基线，CPU/RTL
暖态时间比为 0.4735 倍：尚未计入传输开销就已经慢于 CPU。不能把旧版 6 ns 与新版
7 ns 的时间比叫作同频加速比；同频周期降幅仍是 47.8935%。

最差建立路径是 poly_buffer_1 BRAM，经 Compute_Auto 的存储体选择和条件模取负，
写入 result_buffer BRAM，对应 `src/fpga_backend/src/auto.cpp` 第 80–82 行。
数据路径延迟为 5.681 ns，其中逻辑 1.544 ns、布线 4.137 ns（72.821%）。
它不是 AXI 端口自身的关键路径，但这并不能证明接口拓宽对布局没有间接影响，
后者需要旧版的同条件物理对照实验。后续若保持并行度不变，可检查该路径的流水
寄存器边界以及共享缓存选择/控制信号的扇出。

拥塞报告保留了布局阶段 level 5、路由初始阶段 level 5/6 的热点，部分局部存储资源
窗口占用达到 100%。这些是早期阶段的局部指标，不是全芯片利用率，也不是最终路由
失败数。高扇出信号包括复位（2830 个负载）和共享缓存地址/控制（约 1158 个负载）。
该 OOC 设计没有跨 SLR 网络；当前拓宽版能够布通，但接入平台后拥塞情况仍可能改变。

物理验证的限制如下：

- 内部未约束端点为 0，但 580 个输入端口、526 个输出端口缺少外部延迟约束；
  ap_clk 未通过 HD.CLK_SRC 绑定平台时钟缓冲器。因此时钟源、偏斜和接口时序尚未
  完成整个平台级别的签核。
- DRC 保留 200 条 DSP MREG 流水寄存器告警及 1 条“无可布线负载”告警；没有 DRC
  错误，但没有隐藏或豁免这些告警。
- 未链接 U55C 平台 shell/HBM/PCIe，未上板，也未进行布线后门级仿真。
  功能证据是 RTL 协同仿真，物理证据是独立静态时序分析，不能替代端到端实测。

原始证据保存在本目录：`implementation.log`、`implementation-utilization*.rpt`、
`implementation-drc.rpt`、`implementation-methodology.rpt`、`exported-Top.xdc`、
`postroute.log` 及 `postroute-*.rpt`。首次诊断因 Vivado 不支持 Tcl 的 redirect 命令
而失败，失败日志 `postroute-first-attempt.log` 已保留。修正脚本并去掉已弃用的
`-nets` 后，完整重跑成功；没有改动数据通路。

已布线检查点的 SHA256 为：
4a50d902cf776fecd6ad3ef48e4b4e2fe3c205eb57025f034fedfefbf9dbd552。

## 六、AUTO 与当前 HKS 复合算子的边界

**AUTO 是自同构（Automorphism），用于旋转相关的多项式置换；它不属于当前
OP_HKS_DIGIT 的 ModUp 数据通路。** 同一个通用 Top 同时保留了这两个独立指令。

1. **算法职责不同。** 自同构将 `a(X)` 变成 `a(X^k) mod (X^N+1)`，合适的 k 对应
   CKKS 槽位旋转。它也改变密文所关联的秘密密钥，因此完整旋转要配合密钥切换。
   HKS 指混合密钥切换，不是自同构本身，也并非每种 HKS 使用场景都需要 AUTO。
2. **本仓库的旋转顺序。** `base-leveledshe.cpp` 的 EvalAutomorphism 先执行
   KeySwitchInPlace，再对两个密文分量执行 AutomorphismTransform；密钥生成时使用
   逆自同构后的秘密密钥，所以不能把这里描述为固定的“先 AUTO 再 HKS”。
   EvalFastRotation 则先消费已预计算的 digits 做快速密钥切换，再执行自同构。
3. **本轮只融合 HKS 的预处理。** `keyswitch-hybrid.cpp` 在
   EvalKeySwitchPrecomputeCore 中调用 TryHKSDigitOffload。OP_HKS_DIGIT 内部依次执行
   INTT、QHatInv 预缩放、BConv、NTT，并保留原 EVAL digit 塔；不调用 Compute_Auto，
   也没有包含后续评估密钥乘加和 ModDown。因此这里的周期不是完整旋转或完整 HKS。
4. **已有 AUTO 电路的具体做法。** `auto.cpp` 的 Compute_Auto 在系数域按输出地址
   反查输入位置，使用 `kinv = k^(-1) mod 2N`；绕过 N 的项需要模取负，因为
   `X^N = -1`。这不是把 EVAL 数组直接循环移位。仓库现有带 precomp 参数的
   AutomorphismTransform 包装路径采用 INTT、系数域 AUTO、NTT；同文件不带 precomp
   的重载另有 EVAL 域直接置换实现，不能据此认为系数域转换是自同构数学上的必需步骤。
5. **当前卸载与物理边界。** ConfigureHKSDigitBackend 会设置 HksDigitOnly，旧单算子
   卸载入口随之禁用，AUTO 留在 CPU 路径。尽管本次融合调用不执行 OP_AUTO，通用 Top
   仍包含它的硬件，且同用 ap_clk；整核静态时序检查仍要覆盖这条合法指令路径。
   所以 AUTO 能限制 Top 的可用频率，却不会在两次 OP_HKS_DIGIT 的周期统计中凭空
   增加 AUTO 执行周期。不能仅因这次没有调用，就将其设置为 false path。

核对位置：`src/fpga_backend/src/top.cpp` 第 266、370、429 行；
`src/pke/lib/keyswitch/keyswitch-hybrid.cpp` 第 332 行；
`src/pke/lib/schemebase/base-leveledshe.cpp` 第 348、389、440 行；
`src/core/include/lattice/hal/default/poly-impl.h` 第 314、370 行；
`src/pke/lib/keyswitch/hks_digit_offload.cpp` 第 60 行；
`src/core/include/FpgaManager.h` 第 324 行。
本次仅补充说明，没有删除 AUTO、修改算术逻辑或增加时序例外。

## 七、复现方法与证据一致性

在 `src/fpga_backend` 目录下，使用此前导出的真实 OpenFHE 测试数据运行：

```sh
make test-hks-digit
HKS_PROJECT_TAG=wide256_direct_r3 make hks-csynth
HKS_PROJECT_TAG=wide256_direct_final \
HKS_RTL_FIXTURE="$PWD/build/hks_perf_20260904/openfhe_rtl_fixture.txt" \
make hks-cosim-perf
HKS_PROJECT_TAG=wide256_direct_final make hks-cosim-smoke
python3 check_shared_transform.py \
  Solution/hks_digit_cosim-perf_wide256_direct_final/solution1 \
  --axi-width 256 --lanes 4 --total-multipliers 20
make hks-impl HKS_IMPL_PROJECT=Solution/hks_digit_csynth_wide256_direct_r3
make hks-postroute HKS_IMPL_PROJECT=Solution/hks_digit_csynth_wide256_direct_r3
```

物理实现使用同一最终数据通路的 r3 综合快照；性能和 smoke 分别使用独立工程。
三个完整的 `syn/verilog` 目录已逐文件、逐字节比较，无差异。Top.v 的 SHA256 为：
21249a3901fe1b8299c384e0ddd272d9a02b9c6543f7bb0c6feadde9e3708345。

实现命令在独立日志目录运行，执行的是 OOC IP 实现，不是平台 shell/HBM/PCIe
链接或上板测试。导出的 XDC 只有 6 ns 时钟，没有显式 0.75 ns 用户不确定度和外部
I/O 延迟；postroute 命令在已布线检查点上补充检查 0.75 ns 建立时间不确定度。
OOC 时钟源/偏斜和未约束外部接口仍构成整个平台时序签核的限制。

对比已归档的共享版基线时，使用 `analyze_hks_performance.py` 并显式传入
`--allow-width-change`，因为位宽正是本轮有意改变的变量。保持器件、时钟和负载
一致，保留原始调用周期与综合报告。CPU 基线采用此前机器空闲时的暖态微基准，
不要在综合、实现同时争用 CPU 时重测并混入对比。
