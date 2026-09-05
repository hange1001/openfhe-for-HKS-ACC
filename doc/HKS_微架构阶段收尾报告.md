# HKS 加速器微架构阶段收尾报告

日期：2026-09-05

稳定节点：`e948a69`（P4）

范围：CKKS/DC、N=4096、Q=3、P=2、单 digit ModUp 融合；PE_PARALLEL=4。

## 1. 收尾结论

本轮微架构优化在 P4 收口，不继续接入 Shoup，也不实施 P5 剩余的输出旁路改写。
P0～P4 已完成数据生命周期、片上存储、变换复用和预乘算术复用，并通过 OpenFHE、
native、ASan、HLS、RTL 精确比对及 OOC 布线验收。

这不是“已经没有优化空间”，而是剩余方向具有不同前提：

- `OP_INIT` 优化主要改善冷启动/上下文频繁切换，不改善同上下文暖态；
- BConv 两个副本是控制/寄存器包装副本，不是两套 15 路算术阵列；
- Shoup 主要释放 DSP，当前参数下几乎不缩短 II=1 的 BConv 主循环；
- 工业参数下 NTT 容量、BConv 分块和参数化会改变当前最优结构，应在工业参数确定后重新设计。

因此，P4 是一个功能、资源、周期和物理时序都闭环的稳定微架构基线。后续工作转为
“工业参数参数化”和“板级系统测量”，不在本节点继续堆叠未经验证的局部改动。

## 2. 最终数据流与存储所有权

```text
OP_INIT（上下文变化时）
  ├─ gmem0 → Q 参数 + NTT twiddle  ┐
  └─ gmem1 → P 参数 + INTT twiddle ┴─> 片上上下文缓存

每个 HKS digit（同上下文重复）
  AXI load → A_work/Q slots
      ├─ INTT：A_work ↔ 单塔 scratch
      ├─ SCALE：原位读/写 A_work，复用同一组 4 路 MultMod
      └─ BConv：直接读 Q slots，直接写 P slots
                    └─ NTT：P slot ↔ 同一单塔 scratch
                              └─ AXI store（EVAL 塔旁路 + 补基塔）
```

最终结构要点：

- NTT/INTT/SCALE 只有一个运行时调用入口和一组 4 路 Barrett `MultMod`；
- 变换只有一个共享 scratch tower，不再在函数边界复制完整塔；
- BConv 不再有 `local_in_x/local_out_x`，直接消费/产生共享工作区；
- 三个 AXI master 均为 256 bit；独立固定长度单塔 helper 保持 burst；
- AUTO 已删除；保留独立 INIT、ADD/SUB/MULT、NTT/INTT、BCONV 和融合 HKS_DIGIT。

## 3. P0 到 P4 的结果

### 3.1 周期

| 同一 OpenFHE fixture | P0 | P4 | 变化 |
|---|---:|---:|---:|
| alpha=2 digit | 70,131 | 46,671 | -23,460 |
| alpha=1 digit | 69,603 | 44,571 | -25,032 |
| 暖态两 digit | 139,734 | 91,242 | -48,492（-34.70%） |
| 同频架构加速 | 1.000× | 1.531× | +53.1% |
| OP_INIT | 295,063 | 295,063 | 不变 |

P4 暖态时间下界为 0.547452 ms（名义 6 ns）或 0.638694 ms（保守 7 ns）。
这些数字不含 PCIe、驱动、buffer 分配和 Host 调度，不能称为上板端到端加速比。

### 3.2 资源与物理结果

| HLS 口径 | P0 | P4 | 变化 |
|---|---:|---:|---:|
| BRAM_18K | 688 | 424 | -264（-38.37%） |
| DSP | 1,160 | 1,102 | -58（-5.00%） |
| FF | 78,323 | 80,186 | +1,863 |
| LUT | 147,680 | 178,865 | +31,185 |
| URAM | 96 | 96 | 0 |

| Vivado post-route | P0 | P4 | 变化 |
|---|---:|---:|---:|
| CLB LUT | 91,520 | 106,098 | +14,578（+15.93%） |
| register | 61,076 | 61,700 | +624 |
| DSP | 1,160 | 1,102 | -58 |
| BRAM tile | 348 | 272 | -76（-21.84%） |
| URAM | 96 | 96 | 0 |

P4 的 226,368 条可路由 net 全部布通。默认 6 ns WNS/TNS 为 +0.122/0 ns；
6 ns + 0.75 ns setup uncertainty 为 -0.628/-474.418 ns；保守 7 ns + 0.75 ns
为 +0.372/0 ns，WHS/THS 为 +0.019/0 ns。最终关键路径从 `MultMod` 流水寄存器
到 `poly_buffer_1` BRAM 写口，5.622 ns 中走线占 70.667%。

这组结果说明存储和周期下降不是免费的：统一工作区的动态寻址增加 LUT，但 post-route
仍然布通并提高了默认/保守时序余量。该取舍在当前固定参数下可接受。

## 4. Amdahl 视角：OP_INIT 是什么性质

把一次“两个 digit 的暖态 ModUp”记为一个重复工作单元，则当前 RTL 周期模型为：

`T(K) = 295063 + K × 91242`

| 同一上下文内重复次数 K | 总周期 | INIT 占比 |
|---:|---:|---:|
| 1 | 386,305 | 76.38% |
| 10 | 1,207,483 | 24.44% |
| 100 | 9,419,263 | 3.13% |

因此 `OP_INIT` 是应用层的固定串行段：digit 并行不能覆盖它，而且冷启动时占比极高；
但它不是硬件内部的一条完全串行长链。两张 twiddle 表分别从 gmem0/gmem1 读取，
各有 196,608 项、II=1，HLS 父模块估算为 196,674 周期，说明两条独立路径发生重叠。
完整 RTL 事务为 295,063 周期，额外差值包含 Top/AXI 事务调度，不能只用子模块估算替代。

当前 twiddle 第一维按 `MAX_LIMBS=8` 装载，但数学上只有 Q+P=5 个模数有效。每个方向
有 3×12×2048=73,728 个无效 64-bit 项，两方向合计 1,179,648 byte，占表项 37.5%。
未来可把 twiddle limb 数从工作区跨度中解耦，只加载/存储 5 个有效模数，并评估 256-bit
四项并行装载。若粗略按 5/8 线性缩放，冷启动一轮的理想预测约为 1.40×；这是未综合、
未 RTL 测量的上界预测，不能写成已实现收益。暖态上下文复用时该优化收益为零。

## 5. 两个 BConv 副本的准确含义

P4 可达 RTL 中确实存在两个包装模块：

| 包装器 | 用途 | latency | LUT | FF | DSP |
|---|---|---:|---:|---:|---:|
| `Compute_BConv_Systolic` | 独立 OP_BCONV | 4,145 | 5,260 | 3,949 | 0 |
| `Compute_BConv_Systolic_24_25` | 融合 HKS_DIGIT | 4,145 | 5,520 | 3,942 | 0 |

这两个模块各有 Load_W/Load_Mod、阵列控制和流水寄存器，但 15 路 `MultMod` 被 HLS
提升到 Top 后共享。整核只有 19 路 Barrett：BConv 15 路 + 变换 4 路，即 1,102 DSP；
如果 BConv 算术也复制一份，DSP 不可能仍为 19×58。

统一包装器的资源收益上限约为删除其中一份，即约 5.3k LUT / 3.9k FF，而不是再省
870 DSP。保留两个包装器的原因是独立 OP_BCONV 和嵌套融合路径处于不同调用层次。
未来有两种方向：HKS-only kernel 删除独立入口；或重构为 Top 的单一运行时 BConv
调用点并保留兼容。两者都需要重新检查控制 MUX、高扇出、II 和 post-route，不作为本轮收尾改动。

## 6. 为什么 P5 不再实施

原计划 P5 包含“减少 result 中转”和“最终物理收尾”。后半部分已经由 P4 完成：
AXI256、burst、非法调用零写入、全部网络布通、7 ns + 0.75 ns 正 slack、暖态相对
P0 降低 34.70%、BRAM 降低且 DSP 不增加，均已满足。

剩余旁路改写没有清晰的资源回收目标：

- 原 EVAL 的保存已融合在第一次 AXI load 中，不存在额外的整塔读取 pass；
- 补基塔已经从 P_work 直接写回 AXI；
- `result_buffer` 仍被 ADD/SUB/MULT 使用，仅删除 HKS 对它的访问不会删除这组物理 RAM；
- 改为最后重新读取 `mem_in1` 会增加 DDR 访问，提前直写 `mem_out` 又会收紧输入输出别名约定。

因此 P5 剩余改动的确定收益接近零，风险是 ABI/别名和走线回归。本轮以“审计后不实施”
关闭 P5，而不是为了完成编号制造一次无收益改动。

## 7. Shoup 延后方向

独立原语已经验证，但没有接入 Top/OpenFHE ABI：

| 单路模乘 | DSP | FF | LUT | latency | II |
|---|---:|---:|---:|---:|---:|
| Barrett `MultMod` | 58 | 1,692 | 3,269 | 19 | 1 |
| Shoup（独立实测） | 36 | 898 | 948 | 12 | 1 |

若未来只替换 BConv 的 15 路，整核 DSP 预测为 `4×58 + 15×36 = 772`，相对 P4
减少 330（-29.95%），BConv 模乘池减少 37.93%。但是主循环已有 II=1、tripcount=4103，
单路 latency 缩短 7 拍预计只让两个 digit 合计减少约十几拍（低于 0.02%），实际值必须
由整机 RTL 验证。Shoup 更适合工业参数下 DSP 成为扩展上限时再接入，而不是当前的速度收尾项。

仓库中的 `shoup.h/.cpp` 是未接入 KERNEL_SRC 的独立 WIP，不影响 `e948a69` 的可达 RTL。
未来接入必须同时完成 ABI 版本保护、OpenFHE `PrepModMulConst` 元数据、独立/融合两个入口、
旧 bitstream fail-closed、RTL 实例审计和重新布线，不能只替换一行乘法调用。

## 8. 微架构阶段交付边界

已闭环：

- P0～P4 逐节点源码与 Git 检查点；
- OpenFHE 真实参数/旋转解密、native/ASan、HLS C-sim；
- 两轮真实 fixture RTL 精确比对，每轮 40,960 个余数；
- RTL 周期、HLS 资源、可达实例、AXI 位宽和 OOC post-route 证据；
- 中文实施日志与每阶段报告。

未闭环且不冒充完成：

- 无 U55C 板卡，未测 PCIe/HBM/XRT 墙钟和系统加速比；
- OOC 外部 AXI 端口没有平台 shell 的 input/output delay；
- 参数仍是 N=4096/Q3/P2，不代表工业参数资源与瓶颈；
- 未做完整 HKS 的 KeyMult/ModDown、跨 digit 流水或多 digit 并行；
- 6 ns + 0.75 ns 进取目标未闭合，保守 7 ns 场景闭合。

下一阶段建议先冻结工业参数与板级测量口径，再重新做容量/带宽/算力模型。否则继续在
当前小参数上优化，容易把局部收益误当成工业参数结论。

## 9. 证据入口

- P4 中文报告：`docs/reports/hls/hks_mem_p4_20260905/README.md`
- P4 实施日志：`doc/P3_P4_Shoup_实施日志.md`
- 存储计划：`doc/HKS_片上存储与数据搬运优化实施计划.md`
- 融合接口说明：`src/fpga_backend/HKS_DIGIT.md`
- P4 Git：`e948a69`
