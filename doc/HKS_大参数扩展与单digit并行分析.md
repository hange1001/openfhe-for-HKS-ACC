# HKS 大参数扩展与单 digit 并行分析

日期：2026-09-04。代码基线：32c04b5（P2），P1 为 2181989。
补充澄清：本文讨论的是 digit 内部并行；用户实际询问不同 digit 之间的并行，见同目录 `HKS_digit间并行与BConv资源拆解.md`，其中也记录了 P2 BRAM 归因勘误。
本次为代码与已有证据复核、架构预算；未修改硬件源码或重跑综合。
范围为当前 OP_HKS_DIGIT 的 ModUp，不包含完整 HKS 的 KeyMult/ModDown。

## 1. 当前状态

P1 已跳过无效塔的预乘和装载，Feed_X 对无效行注入 0；P2 已删除 BConv 的 local_in_x/local_out_x 和整块 Load_X/Store_X，直接读写 poly_buffer_1。
P3 尚未实现：Execute_Transform 仍分配 bank_a/bank_b 并执行 TRANSFORM_LOAD/STORE；预乘仍从 poly_buffer_2 写入 poly_buffer_1。

| 项目 | P2 结果 | 口径 |
|---|---:|---|
| 两 digit 暖态 | 62979+57843=120822 周期 | RTL co-sim，比 P0 下降 13.53% |
| 10 次 NTT/INTT，含适配搬运 | 85260 周期，约 70.57% | 每次 8526，综合层级周期归因 |
| 其中 CG 核心 | 64690 周期，约 53.54% | 每次 6469 |
| 两次 BConv | 约 8290 周期，约 6.86% | 2×wrapper 综合延迟 4145，非独立波形计时 |
| BConv 模乘 | 15×58=870 DSP | 占整核 1160 DSP 的 75% |
| 整核资源 | 440 BRAM_18K / 1160 DSP / 96 URAM / 154852 LUT | HLS 估计 |

原始 perf 输出确认 40960 个 residue 与 CPU oracle 精确一致；结构审计确认一套双向变换、4 路蝶形、整机 20 个 MultMod、三个 256-bit AXI。
P2 的 5.250 ns 是 HLS 估算，不能作为布线后频率；历史无 AUTO 版本的 7 ns 签核不能直接沿用。
P2 功能与 HLS/RTL 阶段证据齐全，物理验收仍待后续完成。
证据目录：`docs/reports/hls/hks_mem_p2_20260904/`，重点为 top-csynth.rpt、rtl-audit.json、perf-vitis-output.txt。

## 2. 算法规模与物理规模

定义 L 为当前层 Q 塔数、K 为特殊基 P 塔数、alpha 为当前 digit 有效输入塔数。
补基输出塔数 beta=L+K-alpha，包含其他 Q 塔和特殊 P 塔，不只是特殊基 P。
预乘后：`y[p,k] = sum(q=0..alpha-1, x_scaled[q,k]*W[q,p]) mod modulus[p]`。
单 digit BConv 有 N×alpha×beta 个乘积，变换有 alpha 次 INTT 和 beta 次 NTT。

当前 bconv_systolic.cpp 的 PE_Row/PE_Col 完全展开，物理规模却为 LIMB_Q×MAX_OUT_COLS=L×(L+K)。
active_q 只控制有效数据，不会在运行时消除物理 PE。
当前 alpha=2/beta=3 和 alpha=1/beta=4 的有效乘积槽位占比分别为 6/15、4/15，不表示其余单元已时钟门控。
扩展时应独立定义 BCONV_ROWS/BCONV_COLS/COEFF_LANES，和逻辑塔数分开。

下表大规模形状仅用于预算，未生成生产参数或通过安全检查。实际应从目标 OpenFHE context 导出 N、各层 L、K、各 digit alpha 和模数位宽。

| L/K/alpha/beta | 原代码全展开 L×(L+K) | 仅展开实际 digit alpha×beta | 仅 BConv DSP 估计，前者/后者 |
|---|---:|---:|---:|
| 3/2/2/3 | 15 PE | 6 PE | 870 / 348 |
| 12/4/4/12 | 192 PE | 48 PE | 11136 / 2784 |
| 24/8/8/24 | 768 PE | 192 PE | 44544 / 11136 |

预算沿用当前相同 64-bit Barrett 实现的 58 DSP/PE。xcu55c 整器件有 9024 DSP；平台 shell 占用、布线和频率需另验。
固定 3×5 为 870 DSP；固定 4×8 约 1856 DSP。4×8 若保留当前 4 路变换和 1 路预乘，模乘相关 DSP 约 (32+4+1)×58=2146。
这说明固定阵列存在资源空间，不构成整核可实现承诺。

固定 R×C、每周期一个系数位置时：`T_BConv ≈ N×ceil(alpha/R)×ceil(beta/C)`。
还需加上权重切换、填充排空、部分和访问、bank 冲突和外存等待。
例如 alpha=8/beta=24，4×8 的主体约 6N，3×5 约 15N。跨输入行块必须保留部分和，不能重新清零覆盖。

只增大 N 时 BConv 约按 N 增长，NTT 约按 N log2(N) 增长；BConv 的相对压力不一定上升。
当 L、alpha、beta 增大而物理阵列固定时，BConv 才可能成为主要瓶颈。
同频、忽略尾块与带宽时，设 H 为总蝶形吞吐、M 为 BConv 模乘吞吐：
`T_BConv/T_transform ≈ 2H alpha beta / (M(alpha+beta)log2(N))`。

## 3. 存储边界比单纯 DSP 预算更紧迫

当前正反向 twiddle 各为 `[MAX_LIMBS][STAGE][N/2]`，MAX_LIMBS=2L+K。
双方向展开表原始数据量为 `8×MAX_LIMBS×N×log2(N)` 字节。
以 N=65536、L=24、K=8 为预算例：当前布局需 448 MiB；只保留真实 32 模数仍需 256 MiB。
三份 `[56][N]` 多项式数组还需 84 MiB。
整器件 2016 BRAM36 Tile+960 URAM 合计仅约 42.61 MiB 原始位容量（含校验位），可用数据容量更低，端口也不能任意互换。
容量来源：`docs/reports/hls/hks_no_auto_20260904/physical-utilization.rpt` 的 Available 列。

因此要同时改变 twiddle 驻留和工作区生命周期：比较按塔/stage 缓存、压缩重复因子索引、片上生成等方案，并计入带宽或生成器资源。
只把表搬到 HBM 不会消除供数需求，4 个蝶形每周期仍需 4 个 twiddle。
N=65536 时每塔为 512 KiB；可以评估 alpha 塔 Q_work 常驻、补基输出按 C 塔分批。
为 BConv/NTT 重叠准备两个输出 tile 时，需 2C 个塔的容量，scratch 和原 EVAL 旁路另算。
只有核对输入/输出别名契约后，才可提前回写 EVAL 并回收存储。

## 4. 单 digit 内并行

### 独立模数塔的变换并行

同一 digit 的 alpha 个 INTT 独立，补基后的 beta 个 NTT 也独立。
可增加第二套双向引擎，两套都能执行 NTT/INTT，保留各自内部方向复用。
第二套 4 蝶形引擎的模乘约增加 232 DSP，还需要第二份 scratch、独立塔端口、twiddle 供数和调度。
应比较 1 套×8 蝶形与 2 套×4 蝶形：总算力相同，前者塔尾批浪费少，后者存储带宽可以局部化。
当前两 digit 的变换批数由 5+5 变为 3+3。假设适配搬运也能完全并发且每批仍为 8526 拍，乐观模型为 120822-4×8526=86718 拍，约 1.39×。
这个数字未计端口/AXI 竞争、额外控制或频率变化，不是实测；P3 消除拷贝后必须重新比较。

### BConv 与 NTT 按输出塔批次重叠

固定阵列处理大 beta 时，可先保留一套双向变换引擎：

1. 完成输入塔 INTT 和预乘，得到 Q_work。
2. BConv 完成第 0 批输出塔的全部 N 个系数和全部 alpha 项累加，发布 ready。
3. 变换引擎对第 0 批执行 NTT，同时 BConv 在另一存储区域生成第 1 批。
4. NTT 消费完成发布 free，后续批次复用；通过背压处理两阶段速率不匹配。

需要输出 tile 所有权和互不冲突的物理 bank；只添加 DATAFLOW 或共用数组参数不能保证并发。
当前全列阵列输出塔只相差少量拍到齐，直接重叠尾部的收益很小。分批输出后才有长重叠区间。
即使隐藏当前全部约 8290 个 BConv 周期，整核理想收益也仅约 1.074×；此方案主要服务大参数扩展。

### 系数维并行的边界

BConv 的不同系数 k 独立，可增加系数 lane，但需更多 PE 或分摊已有 PE，并增加带宽与部分和存储。
NTT 已有 4 个蝶形 lane，每拍需 8 次系数读、8 次系数写、4 次 twiddle 读。两套 4 lane 或一套 8 lane 的总需求都会翻倍。
y[p,k] 依赖全部 alpha 个输入项；当前 CG 又有整塔跨 stage 访问，不能将小系数块分别做小 NTT 后直接拼接。
更细粒度 INTT→BConv→NTT 流式融合，需要另行证明次序、依赖、缓冲和速率匹配。

## 5. 建议顺序

先完成 P3 的工作区直接变换，并导出真实目标 context 的容量/带宽预算。
随后让 BConv 的逻辑尺寸与固定 R/C 分离，先验证 3×5 分块和非整倍数尾块，再评估 4×8。
降低当前单 digit 延迟时比较 1×8 与 2×4 变换；面向大参数时评估 BConv/NTT 输出 tile 重叠。
这些是建议，不更改已确认的 P3/P4 范围。P4 复用变换 lane 做预乘期间不能再把同一 lane 计作并发变换算力。

长期可研究 BConv 宽乘积累加后约减，CPU 的 HAVE_INT128 分支已有类似形式。
硬件须证明至少 b_x+b_w+ceil(log2(alpha)) 的累加位宽和 Barrett 输入范围，不能直接删除逐项约减。

## 6. 外部参考

- [OpenFHE 参数选择](https://openfhe.discourse.group/t/help-with-parameter-selection/1833/2)：模数链、HYBRID 扩展模数和安全环维度之间的关系。
- [OpenFHE NumLargeDigits](https://openfhe.discourse.group/t/what-does-setnumlargedigits-do/1890)：分解、扩展模数及参数约束，不能只为阵列尺寸任意调整。
- [AMD U55C 用户指南](https://docs.amd.com/r/en-US/ug1469-alveo-u55c/Card-Features)：XCU55C 与 HBM 平台特性。
