# HKS 片上存储 P4：预乘复用变换模乘 lane

## 结论

P4 将 `ApproxSwitchCRTBasis` 的 `QHatInv` 预乘从独立的
`Prepare_HKS_BConv_Input` 硬件移入唯一的 `CG_Transform_Work`。每个输入塔按
`INTT → SCALE` 顺序执行，随后进行一次 BConv 和补集塔 NTT。三类操作通过同一个
`Execute_Transform` 调用点进入同一物理变换引擎，保持 `PE_PARALLEL=4`。

生成 RTL 的结构审计确认：

- NTT、INTT、SCALE 共用一组 4 路 `MultMod`，不是三组各 4 路；
- Top 一共 19 个 `MultMod`：BConv 15 路，变换引擎 4 路；P3 为 20 路；
- 两个蝶形奇偶客户端和一个 SCALE 客户端均请求 4 lane，三者互斥调度；
- SCALE、两个蝶形循环均为 `II=1`；独立预乘模块已从可达 RTL 中消失；
- 直接工作区、单 scratch、8 个 T2P bank、256-bit 三组 AXI 接口均保持不变。

P4 通过 OpenFHE、原生、ASan、HLS C-sim 和两轮 OpenFHE fixture 的 RTL 精确比对。
板卡不可用，本文只报告 RTL 周期下界和 OOC 物理实现，不称为上板端到端时间。

## 性能与资源

同一 OpenFHE fixture（`INIT + alpha=2 + alpha=1`）的无随机 stall RTL co-sim：

| 指标 | P3 | P4 | 变化 |
|---|---:|---:|---:|
| INIT | 295,063 cycles | 295,063 cycles | 0 |
| alpha=2 digit | 52,779 | 46,671 | -6,108 |
| alpha=1 digit | 47,643 | 44,571 | -3,072 |
| 暖态两 digit | 100,422 | 91,242 | -9,180（-9.141%） |
| 同频架构加速 | 1.000× | 1.1006× | +10.06% |

测前预测为减少约 9,216 周期；实测少 36 周期，偏差 0.39%。P4 暖态时间换算为
0.547452 ms（名义 6 ns）或 0.638694 ms（保守 7 ns 场景），均不含 PCIe、驱动、
buffer 分配和 Host 调度，只是 kernel RTL 周期下界。

| HLS 口径 | P3 | P4 | 变化 |
|---|---:|---:|---:|
| BRAM_18K | 424 | 424 | 0 |
| DSP | 1,160 | 1,102 | -58（-5.00%） |
| FF | 81,549 | 80,186 | -1,363 |
| LUT | 180,505 | 178,865 | -1,640 |
| URAM | 96 | 96 | 0 |

DSP 恰好减少一套 Barrett `MultMod`（58 DSP），与“删掉独立预乘、复用现有四路”的
结构解释一致。BRAM/URAM 不变，因为 P4 没有再删除整塔存储。

## 调度探索及拒绝项

P4 不是仅凭 C++ 函数名判断复用，先后保留了以下实验：

1. r1 把 SCALE 与蝶形访问揉进同一循环，HLS RAM 端口调度退化为 `II=2`，拒绝。
2. r2 分开 SCALE/蝶形循环但放在同一父引擎，三类客户端共用四路模乘，全部 `II=1`；
   RTL 40,960 个余数精确通过。
3. r3 删除多余 WAW 提示，仅保留 SCALE 的窄 `inter RAW=false`，结果与周期不变；
   再次通过 40,960 个余数精确比对，作为最终证据工程。
4. r4 完全删除 SCALE RAW 提示，Vitis 将原位读写误判为距离 1 相关，`II=22`，拒绝。
5. r5 拆成低/高固定 bank 两段，Vitis 仍误判距离 8 相关，只能 `II=3`，拒绝。

最终源码使用 r3 调度。独立 bank 检查器穷举 4,096 个系数，证明每个流水迭代地址
不重叠、每个活动 bank 每拍只有一次读和一次写；两轮 RTL 精确比对进一步验证依赖提示
没有改变语义。xsim 运行时仍会对该明确覆盖的提示给出依赖监视告警，不能把告警数量当成
算术失败；最终余数和哨兵才是功能判据。

## 验证结果

- OpenFHE Release：`hks-digit-openfhe` 通过；
- OpenFHE ASan：通过，零内存错误；
- 原生 Top：18/18；HKS：22 个有效形状，非法 descriptor 与 canary 通过；
- HLS C-sim：0 errors；
- RTL OpenFHE fixture：两次均精确通过 40,960 个余数；
- 结构审计：一套变换、4 路共享模乘、Top 共 19 路、AXI256、无 AUTO、无独立 SCALE；
- OOC 布线：226,368/226,368 条可路由 net 全部布通，routing error 为 0；
- 默认 6 ns：WNS/TNS = +0.122/0 ns，WHS/THS = +0.019/0 ns；
- 6 ns + 0.75 ns setup uncertainty：WNS/TNS = -0.628/-474.418 ns，不通过；
- 7 ns + 0.75 ns setup uncertainty：WNS/TNS = +0.372/0 ns，通过。

布线后资源如下。HLS 与 Vivado 的 RAM/LUT 统计口径不同，不能跨行直接相减。

| Vivado post-route | P3 | P4 | 变化 |
|---|---:|---:|---:|
| CLB LUT | 108,392 | 106,098 | -2,294 |
| CLB register | 64,058 | 61,700 | -2,358 |
| DSP | 1,160 | 1,102 | -58 |
| BRAM tile | 272 | 272 | 0 |
| URAM | 96 | 96 | 0 |

默认 6 ns 最差 setup 路径从 `MultMod` 的流水寄存器到 `poly_buffer_1` 的 BRAM
写口，data path 为 5.622 ns，其中逻辑 1.649 ns、走线 3.973 ns（70.667%）。
P3 最差路径位于 BConv `Load_W`；P4 后关键路径已经迁移，但仍以走线为主。

## 证据文件

- `rtl-audit.json`：生成 RTL 层次、实例数、端口、II 与 HLS 资源；
- `top-csynth.rpt`：Top 综合报告；
- `transform-work-csynth.rpt`：共享变换父引擎报告；
- `scale-loop-csynth.rpt`：SCALE 循环 `II=1` 报告；
- `perf-cosim.rpt`：RTL co-sim 总结；
- `perf-transactions.rpt`：INIT 和两个 digit 的逐事务周期。
- `postroute/route_status.rpt`：最终可布通性；
- `postroute/utilization.rpt`：Vivado 布线后资源；
- `postroute/timing_default.rpt`、`timing_setup075.rpt`、
  `timing_period7_setup075.rpt`：三种时序口径；
- `postroute/paths_setup075.rpt`、`high_fanout.rpt`、`congestion.rpt`：
  关键路径、扇出和拥塞证据。

## 当前边界

- 固定参数仍为 N=4096、Q=3、P=2，本节点没有扩展工业参数形状；
- 本节点不做 digit 并行，也不改变 PE；
- BConv 15 路模乘仍是 Barrett，Shoup 接入属于下一个独立提交；
- OOC 的 AXI 外部端口没有板级 input/output delay，但内部无未约束 endpoint；
- 没有板卡、平台 shell 链接或真实 PCIe 测量，不能由本报告声称系统加速比。
