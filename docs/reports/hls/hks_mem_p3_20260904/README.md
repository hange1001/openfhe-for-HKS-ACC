# P3：共享工作区与原位预乘

状态：C-model、综合结构、RTL smoke/perf、随机 AXI stall 与 OOC 布线验收均通过。
源码检查点为 `53cf757`；布线报告在检查点生成后补入本目录。

## 改动范围

基于 P2 `32c04b5`，保持 PE=4、Q3/P2、N4096、三个 AXI256、指令编号不变。
`poly_buffer_1` 就是 A_work；不是再新建一套 A_work。输入物理槽为 `q`，对应全局模数
`digit_start+q`；BConv 补基物理槽为 `3+p`，对应模数 `p<digit_start?p:p+alpha`。

CG 引擎直接将这个工作区与一条 scratch 交替读写，12 级后结果回到工作区。
删除原 `TRANSFORM_LOAD`、`TRANSFORM_STORE`、独立 bank_a/bank_b；预乘读取并写回
同一个工作系数。最后在 AXI 写回边界选择旁路 EVAL 或工作区结果，没有另做整塔复制。
保留 `poly_buffer_2` 作为通用算子的第二操作数，保留 `result_buffer` 作为通用结果/EVAL旁路。

## 综合结果（不是上板测量）

Vitis HLS 2023.2，U55C，6 ns + 0.75 ns uncertainty。

| 指标 | P2 | P3 r3 |
|---|---:|---:|
| BRAM_18K | 440 | 424 |
| DSP | 1,160 | 1,160 |
| FF | 76,770 | 81,549 |
| LUT | 154,852 | 180,505 |
| URAM | 96 | 96 |
| 单变换核周期 | 6,469 | 6,481 |
| 单变换含外围周期 | 8,526 | 6,483 |
| 两个 butterfly parity 循环 II | 1 / 1 | 1 / 1 |

BRAM 净减16：原两条局部变换塔32 blocks变成一条scratch16 blocks，原有工作区复用。
LUT 增25,653（约16.57%）：不再用局部A bank隔离多塔工作区，增加运行时塔选择、读MUX及写使能。
这是真实代价，不宣称所有资源都下降。核心多12周期来自每级增加1拍的调度开销；
含外围净省2,043周期/变换，去掉拷贝仍有明显收益。

## 真实 OpenFHE RTL 周期

| 调用 | P2 | P3 | 变化 |
|---|---:|---:|---:|
| INIT（冷启动） | 295,063 | 295,063 | 0 |
| alpha=2, start=0 | 62,979 | 52,779 | -10,200 |
| alpha=1, start=2 | 57,843 | 47,643 | -10,200 |
| 暖态两个digit | 120,822 | 100,422 | -20,400（-16.88%） |

同频核内加速约1.203×。每digit的5次变换净省10,215周期，其他边界控制净多15周期，
合计净省10,200。实测100,422相对事前95,000预测偏差+5.71%，位于±15%区间。
6ns名义换算0.602532ms；7ns假设换算0.702954ms，时钟可采用性以本版布线后报告为准。

## 结构与正确性依据

- 可达层次只有一个 `Execute_Transform_Operation → Execute_Transform → CG_Transform_Work`。
- 4 路变换模乘、整机20路Barrett模乘；没有将 II 下降后工具串行化的版本误当作资源优化。
- 64个工作区T2P bank、8个scratch T2P bank；检查实际两个可写端口及128条工作区写使能。
- 端口预算检查器穷举8槽位、双方向、12级，所有系数每级访问一次、每bank每拍至多2访问。
  它证明算法端口预算，不能单独替代RTL调度和仿真。
- r1不解除依赖时II=5；r2仅解除经证明的跨迭代WAW后II=3；r3保留方向分支内的静态bank地址，恢复II=1。
  没有沿用不加论证的整套 `DEPENDENCE false`；RAW和intra依赖及stage排空仍保留。
- native新增完全重叠和两种部分重叠的输入输出测试；HLS多AXI主机cosim不模拟地址别名，二者口径分开。

## 验收状态

- Top18/18、HKS22、独立BConv9/9通过（包含无效行毒化）。
- OpenFHE Release / ASan各472项检查、1,523,712个余数精确一致；ASan检测泄漏开启。
- 最终完整C-sim通过；RTL smoke 35次调用/4种digit通过，perf 40,960余数精确一致。
  smoke与perf的生成Verilog目录逐文件一致。
- 使用 `config_cosim -random_stall` 重放 perf fixture 通过，覆盖 AXI 随机背压下的
  数据完整性；它验证协议/流水线停顿，不改变无 stall 的周期基线。
- 无板卡，任何周期×时钟换算仅是核内时间，不含PCIe/驱动/主机准备；INIT单列。

## 布线后结果

Vivado 2023.2，U55C `xcu55c-fsvh2892-2L-e`。全部 235,225 条可路由 net 均完成
路由，routing error 为 0。

| 时序情景 | WNS | TNS | setup 失败端点 | WHS | THS | hold 失败端点 | 判定 |
|---|---:|---:|---:|---:|---:|---:|---|
| 默认 6 ns | +0.029 ns | 0 | 0 | +0.007 ns | 0 | 0 | 通过 |
| 6 ns + 0.75 ns uncertainty | -0.721 ns | -595.871 ns | 3,140 | +0.007 ns | 0 | 0 | 不通过 |
| 7 ns + 0.75 ns uncertainty | +0.279 ns | 0 | 0 | +0.007 ns | 0 | 0 | 通过 |

因此 P3 可以按保守的 7 ns 工程情景签核，不能宣称带 0.75 ns 裕量的 6 ns 已通过。
最差 setup 路径位于 BConv `Load_W` 内部：data path 5.951 ns，其中 routing
5.436 ns（91.346%），logic 0.515 ns（4 levels）；旧的整塔工作区到局部 bank
路径已消失。该结果说明当前临界项主要是布局/走线，而不是新增深逻辑链。

| 布线后资源 | 用量 | 器件占比 |
|---|---:|---:|
| CLB LUT | 108,392 | 8.31% |
| CLB Register | 64,058 | 2.46% |
| Block RAM Tile | 272 | 13.49% |
| DSP | 1,160 | 12.85% |
| URAM | 96 | 10.00% |

HLS 资源估计与 Vivado 布线后利用率属于不同阶段、不同统计口径，不能直接相减；
两组数字均保留，前者用于版本间快速审计，后者用于器件容量判断。原始证据见
`postroute_margin/`，PE 扫参与布线分析见 `doc/PE扫参与P3布线结果_2026-09-05.md`。

完整试验过程见 `doc/P3_P4_Shoup_实施日志.md`。失败项目r1/r2保留在本地Solution下，不纳入通过证据。
