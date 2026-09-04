# HKS digit 间并行与 BConv 资源拆解

日期：2026-09-04。用户澄清：本次讨论不同 digit 同时执行，不是同一个 digit 内的塔并行。
资源证据基线是 P2 提交 32c04b5、mem_p2_r1。复核期间曾观察到 PE_PARALLEL 从 4 改为 32；用户说明正在扫参，取值随后继续变化。本次不修改扫参配置，所有归档数字均按四路版本注明。
本次只做源码、生成 RTL、原始综合报告分析及日志记录，没有修改硬件、重新综合或提交 Git。

## 1. Digit 间并行是否可行、是否必要

算法上可行。不同 digit 的输入 Q 子集、预乘参数和 BConv 权重独立，ModUp 输出可以分别生成，没有 digit j 必须先消费 digit j-1 结果的依赖。
完整 HKS 中，每个 digit 与自己的评估密钥相乘，再把各 digit 对两路输出的贡献做模加归并；不能让多个任务无仲裁写同一累加器。
当前 OP_HKS_DIGIT 仅实现 ModUp，KeyMult/ModDown 还在核外，不能把 ModUp 并行收益称作完整 HKS 加速比。

源码依据：

- `src/pke/lib/keyswitch/hks_digit_offload.cpp`：TryHKSDigitOffload 有状态锁，按 part 循环调用 Transfer。
- `src/core/include/FpgaManager.h`：对应传输路径 run.wait()，完成当前调用再返回。
- `src/pke/lib/keyswitch/keyswitch-hybrid.cpp`：各 partsCt 独立构造；EvalFastKeySwitchCoreExt 按 digit 做 cj×aj、cj×bj 并累加。
- FPGA Top 使用全局可变 poly_buffer/metadata 和单套变换资源；当前接口并非可重入多上下文调度器。

只在 Host 多开线程不会自动生成并行硬件；需要多 CU，或核内多个 digit 上下文和资源调度。

| 方案 | 作用与限制 |
|---|---|
| 完整复制两套 digit 执行器 | 两 digit 的计算可并发，最直观；会重复 BConv/控制和在途工作存储，提高 twiddle 与 AXI 带宽要求 |
| 两套变换引擎服务多个 digit，共用 BConv | 将新增算力放到变换上；digit 在独立上下文中排队使用 BConv，需要存储所有权和公平/依赖调度 |
| 一套变换引擎，跨 digit 交错 BConv/预乘/I/O | 可隐藏其他阶段，但所有 NTT/INTT 仍占用同一个引擎，不会使变换算力翻倍 |

P2 两 digit 是 62979、57843 拍，合计 120822。假设完整复制、同频、带宽无竞争，理想并发时间为 max=62979 拍，约 1.918×；仅暖态 ModUp，不含 INIT、Host/PCIe 或完整 HKS 的其他阶段。
完整复制两套当前模乘资源约 2320 DSP。若两套四路变换共用 15 路 BConv 和一路预乘，则模乘资源约 (8+15+1)×58=1392 DSP；共享排队使其时延不能直接套用完整复制的 1.918×。
可共享只读 twiddle/参数内容，但要解决并发端口，不等于缓存完全不用增加。

必要性取决于目标：若要求单次 HKS 的低延迟且有多个 digit、存储带宽充足，digit 并行有价值；若面积/存储优先，则不是必选。
实际 digit 数会受层级和分解参数影响，可能有尾 digit 或只有一个 digit。固定多执行器会有负载不均或闲置；不得只为填满硬件任意增加 NumLargeDigits。
当前建议先完成共享存储及新的变换并行度评估，再决定是否加多 digit 上下文。当前 BConv 已占 P2 的 870/1160 DSP，而执行时间约占 6.9%，优先整套复制并不经济。
PE_PARALLEL=32 后该占比会变，必须重做整核测量；当前 BConv 行列常量仍为 3×5。

## 2. 870 DSP 的逐项来源

每个 PE 调用一次 MultMod，但这个函数并非一个普通乘法器。
令 z=x×w，m 为 Barrett 预计算因子，qhat 为近似商：

| 操作 | 原始综合映射 | 单 PE DSP | 15 PE DSP |
|---|---|---:|---:|
| z=x×w | 64×64→128 | 16 | 240 |
| high64(z)×m | 64×64→128 | 16 | 240 |
| low64(z)×m | 64×64→128 | 16 | 240 |
| qhat×mod | 有效低 64-bit 乘法，64×64→64 | 10 | 150 |
| 剩余减法、比较、移位、校正、阵列模加 | fabric | 0 | 0 |
| 合计 | 3 个全乘积 + 1 个截断乘积 | 58 | 870 |

其中真实 x×w 为 240 DSP（27.59%），Barrett 约减中的乘法为 630 DSP（72.41%）。
qhat×mod 的源变量虽然声明 128 位，后续只消费低位，报告/RTL 已裁剪为 64 位有效结果，不能按源类型估算成完整宽乘法。
MultMod 延迟 19 拍、II=1，表示可每拍接收新输入；不能再把 58 DSP 乘以 19。
权重和模数在一次调用内稳定，但由运行时加载，所以并非综合时常量，不会自动变成常数乘法电路。

证据：`src/fpga_backend/Solution/hks_digit_csynth_mem_p2_r1/solution1/syn/report/MultMod_csynth.rpt`，及归档 `docs/reports/hls/hks_mem_p2_20260904/top-csynth.rpt` 的 Bind Op Report。

## 3. LUT、FF 和共享资源的归属

单个 MultMod：3269 LUT、1692 FF、58 DSP。
LUT 分类为 Expression 2280、子实例 797、寄存器实现辅助 192；FF 为子实例 408、Register 1284。
其中三个运行时桶形移位各 423 LUT，合计 1269 LUT/PE；15 PE 合计 19035 LUT，只是模乘内部的一部分。
因此 DSP=0 的移位/校正仍可能占较多 LUT；流水线状态/参数对齐消耗 FF。

| 当前可达硬件 | DSP | LUT | FF | BRAM_18K |
|---|---:|---:|---:|---:|
| 15 个共享 BConv MultMod | 870 | 49035 | 25380 | 0 |
| HKS BConv 包装层及阵列控制 | 0 | 5517 | 3942 | 0 |
| 独立 OP_BCONV 包装层及阵列控制 | 0 | 5260 | 3949 | 0 |
| 上述三项小计 | 870 | 59812 | 33271 | 0 |
| poly_buffer_1 工作存储 | 0 | 选择/连接逻辑位于上层，不单独归因 | 同左 | 128 |

两份包装层真实存在：Top.v 实例化 Top_Compute_BConv_Systolic，Top_Execute_Transform_Operation.v 实例化 Top_Compute_BConv_Systolic_32_33。
Top.v 只有 15 个 Top_MultMod，包装层通过外提端口共享该池；因此不能将报告两个调用树里的 15 路各算一次。
包装层的零 DSP 也不表示 BConv 没有模乘器。
小计不含共享工作区、AXI、外部参数加载和顶层连接逻辑，不能视为可独立裁剪的完整 BConv 面积；也不是布局布线后资源。
HKS 包装层内部：阵列循环 4846 LUT/1991 FF，权重加载 213 LUT/971 FF，模数参数加载 125 LUT/965 FF，余下约 333 LUT/15 FF 为上层控制。
阵列循环包括模加/条件减、地址、选通、输入与部分和流水状态；x_curr/sum_curr 的源码临时数组不能机械地当成另一份物理寄存器。
poly_buffer_1 同时服务其他算子，BRAM 不能全部记为 BConv 专属开销。

## 4. P1→P2 BRAM 归因勘误

原 P2 README 将净减 248 解释成 local 128 + 工作区重映射节省 120；该解释与原始存储报告、生成 RTL 不符。

| 项目 | P1 | P2 | 变化，BRAM_18K |
|---|---:|---:|---:|
| HKS BConv local_in_x/local_out_x | 128 | 0 | -128 |
| 独立 OP_BCONV local_in_x/local_out_x | 128 | 0 | -128 |
| poly_buffer_1 | 8×15=120 | 64×2=128 | +8 |
| 其余资源 | 312 | 312 | 0 |
| 整核 BRAM_18K | 688 | 440 | -248 |

正确核算：688 - 2×128 + 8 = 440。
每份旧局部缓冲为 local_in_x 48 + local_out_x 80；存在两份包装层缓存，共享的是算术池。
P2 成功删除了两份局部缓存；但“省去工作区 120 BRAM”及“只有一份 BConv wrapper”的旧文字应以本次原始证据复核为准。
原始证据：P1 top-csynth.rpt Storage Report 的两条 wrapper（1319/1401 行），P2 的 poly_buffer_1 64 个独立 RAM 条目，以及两版本对应的生成 RTL。
本次记录勘误，未覆盖历史报告原文。

## 5. 哪些优化值得先研究

优先核查重复的约减成本：研究每个输出列对宽乘积求和，再在列末做约减，或对固定长度输入行块累加后约减。
只作算术结构比较，如果 15 个原始乘法各 16 DSP、5 个列末约减仍能各用 42 DSP，则为 15×16+5×42=450 DSP，比 870 少 420。
这个 450 不是实现预测：新的累加位宽、近似商误差、校正次数、归约吞吐和时序都需要重新证明/综合，现有 MultMod 的两次校正保证不能直接沿用。
无符号累加至少需 b_x+b_w+ceil(log2(alpha)) 位。还应核查输入 x 属于源模数、输出属于目标模数的范围关系。

另一方向是缩小物理 PE 数后分块复用，以更多周期换 DSP；决定多少资源分给变换，要依据新的 PE32 综合/仿真结果。
LUT 优化可研究收紧合法模数/shift 范围或专用约减结构；运行时变量写成若干 switch 分支不保证移位器消失。
统一两份包装层可减少部分重复控制/寄存器，但不会再次省下 15×58 DSP，因为这部分已经共享。

## 6. 用户后续重点：替换 BConv 约减方式

用户说明并行度正在扫参，本轮暂停该方向的设计讨论，专门评估模约减。建议第一候选为固定乘数预计算的 Shoup 模乘；第二候选为宽累加后分组/列末约减。

### 6.1 为什么 Shoup 匹配当前数据流

每个 BConv 权重 w=W[q,p] 在一个 digit 内对 N 个系数重复使用。令 B=2^64，目标模数为 p；Host 对每个权重预计算 w_pre=floor(w×B/p)，在 FPGA 上重复使用。
这是一项与系数无关的预计算，除法留在 Host，不生成 FPGA 除法器；权重可以运行时切换，不必是综合常数。

FPGA 计算：

1. qhat=high64(x×w_pre)。
2. r=low64(x×w)-low64(qhat×p)，按无符号 64-bit 环绕相减。
3. 若 r>=p，则 r-=p，得到标准余数。

适用范围是 0<=x<B、0<=w<p、p<B/2；当前 HKS 要求 p<2^62，满足该条件。w=0 也成立。
并不要求 x<p：BConv 输入来自源模数 q，可能大于目标 p，这里仍正确。
因为 0<=wB/p-w_pre<1 且 0<=x*w_pre/B-qhat<1，可得 0<=x*w-qhat*p<xp/B+p<2p<B。
因此低位乘积与环绕减法能还原该范围内的真实差值，最后一次条件减足够。
无效输出列的 p=0 必须绕过预计算并维持现有零注入语义，不得对零模数计算 w_pre。

该输入范围由 Harvey 论文 Algorithm 2 及第 3 节明确支持：
[Faster arithmetic for number-theoretic transforms](https://arxiv.org/pdf/1205.2926)。
本项目 `src/core/include/math/hal/intnat/ubintnat.h` 的 PrepModMulConst/ModMulFastConst 已使用同类结构；其实现采用 qhat+1 和符号判断补加 p 的等价校正形式。
不能将它概括成“CPU BConv 内积已经逐项使用 Shoup”：当前 HAVE_INT128 的内积路径采用宽乘积累加后约减，预乘使用 ModMulFastConst。

### 6.2 面积预算与时延边界

| 项目 | 当前 Barrett | Shoup 候选 |
|---|---|---|
| 乘法结构 | 3 个完整 64×64 乘积 + 1 个低半乘积 | 1 个高半乘积 + 2 个低半乘积 |
| 单 PE DSP | 已测 58 | 按现有映射约 16+10+10=36 |
| 15 PE DSP | 已测 870 | 纸面约 540，减少 330（37.93%） |
| 通用动态 shift | 当前三个桶形移位 | 仅固定 high/low 位选取，无该组动态 shift |
| 末端校正 | 两次条件减 | 一次条件减 |
| 额外预计算 | 当前按模数的 m/S | 每个权重一个 64-bit w_pre，15 个为 120 字节 |

高半乘积预算保守按完整乘法 16 DSP，具体 HLS 可能采用不同实现；低半乘积也需重新综合确认，不能把 540 写成已验证结果。
额外 120 字节是预计算表原始大小；若删除 BConv 专用 m/S 寄存器，净寄存器增长可能更少，不能混入仍供 NTT 等使用的全局参数。
Host/FPGA 的权重元数据协议需一致调整并考虑版本兼容，不能直接改变既有缓冲区长度后沿用旧 Host。

速度应区分三个指标：单模乘 latency、循环 II、布局布线后的时钟。
当前 MultMod=19 拍、II=1；BConv 主循环同样 II=1、主体约 N 拍。Shoup 缩短乘法依赖链并去掉动态 shift，有机会降低单次 latency、排空开销和关键路径。
但相同时钟、相同阵列和 II=1 下，约 4096 拍的主体并不会因 DSP 少 38% 而少 38%；频率提升必须用布局布线证实，不能根据资源直接换算加速比。

### 6.3 其他路线与优先级

Montgomery 可以作为对照。对固定权重，可预先将 w 转为 wR mod p，再让普通域 x 与该权重做 Montgomery 乘法，得到普通域 xw mod p；不必为了这一次 BConv 把整条数据流转换到 Montgomery 域。
其实际 DSP/latency/II 取决于字宽和实现，尚无本项目测量，不能直接断言比 Shoup 更省。
特殊形状素数约减需要模数具有特定结构，不能假设 OpenFHE 当前给出的所有 RNS 素数都适用。
宽累加后约减能减少约减器数量，但须重做累加范围和约减精度设计，不能直接把几个 Shoup 残值相加后只减一次 p。

推荐先做独立 Shoup 算术模块对比，保持扫参工作不受影响，验证范围与现有 BConv 一致；随后再决定接入。这是实施建议，本次没有实现或启动综合。
验收应包含：x>=目标p、0/最大权重、接近2^62的模数、跨模数输入、无效列；模乘/阵列 II 不退化；综合资源与 latency；完整 BConv/HKS/OpenFHE residue 回归；最终物理时序。

### 6.4 本次数学检查

使用大整数精确运算，对上述 low/high-word Shoup 公式进行 100288 组边界与确定性随机抽查，其中 87485 组 x>=目标p。
包含 x=0、p-1、p、p+1、2^64-1，w=0/1/p-1，以及小模数至接近 2^62 的奇数模数；全部满足余数等于精确 xw mod p，且校正前 r<2p。
这只是数学模型检查，不是新 HLS C-sim、RTL 仿真或硬件性能测量。

## 7. OpenFHE 接入与 P3/P4/Shoup 实施顺序

日期：2026-09-04。用户询问接入方式与顺序，本次仅核对接口并给出建议，未开始实现。扫参配置不作修改。

### 7.1 OpenFHE 数学与数据格式保持一致

Shoup 仍计算标准域 x×w mod p，返回 [0,p) 的精确 residue；不改变 RNS 素数、权重数学定义、NTT/EVAL 格式或 digit 排列。
因此主要改动集中在现有 FPGA 适配层、参数协议和 BConv 计算单元，CPU 的原始数学路径保留作独立参考及回退。
当前入口为 KeySwitchHYBRID::EvalKeySwitchPrecomputeCore → TryHKSDigitOffload → Transfer → Top。

Host 已从 GetPartQlHatModp(level,part) 取得权重，并从 GetParamsComplPartQ(level,part) 取得目标模数。
每个有效权重生成 `weights[i][j].PrepModMulConst(target_modulus[j])`，即 floor(w×2^64/p)，转换为 64-bit 元数据。
必须使用这一输出列的目标模数，不能误用输入 q_i；每个 digit 的补基列还包含其他 Q 塔。
本地接口依据：`src/core/include/math/hal/intnat/ubintnat.h` 的 PrepModMulConst。当前 offload 已限制 NATIVEINT=64。
这种计算只依赖参数，不依赖密文系数，可按准确的有序源/目标模数、权重及 digit/level 身份缓存；context/参数变化与后端切换时需失效，不能只以 N 或目标模数做缓存键。
统计中应将首次预计算与暖态缓存命中分开，仍计入真实 Host 端开销。

### 7.2 参数协议与修改点

当前 HKS payload 是 18 个 uint64_t：15 个 weights + 3 个 QHatInv。
建议保留原 payload 前缀，再追加 15 个 weights_precon，得到 33 个有效数据字（144→264 字节，未计新增协议头）。
QHatInv 的每系数预乘仍由 FPGA 完成；追加的参数预计算不等于把这个预乘移回 CPU。

| 位置 | 所需修改 |
|---|---|
| src/pke/lib/keyswitch/hks_digit_offload.cpp | 生成/缓存每权重预计算值，更新 kMeta、打包偏移和字节统计 |
| src/core/include/FpgaManager.h | HksDigitTransfer 的固定 in2.size()==18 检查及协议匹配；BO 分配目前按 vector.size()，会随新长度调整 |
| src/fpga_backend/include/hks_digit.h 与共享 ABI 定义 | 统一 payload 长度、weights/precon/inv 偏移及版本，避免 Host/HLS 各自硬编码 |
| src/fpga_backend/src/top.cpp | Load_HKS_Params 读取 precon，按 alpha/补基列有效性保护；向 BConv 传递对应权重、precon、目标模数 |
| src/fpga_backend/src/bconv_systolic.cpp 及接口 | PE 改调独立的 Shoup 模乘；其他变换/预乘继续使用原有 Barrett 模乘，不全局替换 MultMod |
| OpenFHE fixture 导出、RTL fixture 解析和集成测试 | 更新固定18字长度、fixture版本和字节统计，重新导出真实输入及CPU oracle |

版本兼容必须明确：33 是 payload 数据字数，不是已定最终协议总长；版本/长度信息及 Host 与 kernel 能力匹配需一起设计。
不能让新 kernel 按更长布局读取旧 Host 的 18 字缓冲；也不能把旧 kernel 忽略新增字段后的正常结果误认为新硬件已启用。
保留 Top 的指针传输形式，不要求改变密文输入与结果数组的格式。零模数无效列不做除法，保持无效行/列保护语义。

独立 OP_BCONV 同样调用现有 systolic 核心，必须同步迁移：FpgaManager::BConvOffload 当前打包 15 weights + 5 mod + 5 S + 5 m 共 30 字，新布局需携带 weights_precon；Load_BConv_Params、测试台同步更新。
否则只修改 HKS 入口而保留可达的旧 BConv 算术，可能使 Shoup 和 Barrett 两套 BConv 同时存在，无法兑现资源节省。
同时迁移两入口后仍需 RTL 审计，确认共享同一组 BConv Shoup 单元；两份包装控制是否进一步统一另行记录。
NTT/SCALE 所需的全局 Barrett m/S 仍保留，不能随 BConv 的参数迁移一起删除。

### 7.3 推荐的整核推进顺序

**P3 → P4 → BConv Shoup 接入**；Shoup 独立算术模块的验证可以提前开展，但不要求改动正在扫参的 Top。

| 顺序 | 核心问题 | 选择依据 |
|---|---|---|
| P3 | 数据存在哪里、怎么直接访问 | 工作区直接作变换 A bank、共享 scratch、预乘原位；先固定真实存储与端口契约，消除主要复制开销 |
| P4 | 预乘使用哪一组物理计算单元 | 在 P3 的原位存储上复用变换模乘器，排空/valid/写地址对齐；先完成原计划的存储与预乘性能收益 |
| BConv Shoup | BConv 的单次模乘如何实现 | 独立改变 BConv 算术和 Host 参数协议，方便定位正确性问题、单独归因 DSP/LUT/latency 变化 |

严格依赖是 P3 的存储契约先于既定 P4 原位多lane调度；BConv Shoup 与 P4 没有数学上的前置依赖。
选择排在 P4 后是工程优先级：避免在一次改动里同时重写存储、资源共享和参数协议，也便于给每阶段建立可复现差分。
本轮 Shoup 仅限 BConv；若以后要让 NTT/预乘也换 Shoup，才需要重新审查 P4 共享数据通路和额外 twiddle 预计算存储。
原计划 P4 的四路及 20→19/1160→1102 是 PE=4 基线数字。扫参后使用哪个版本，应以该版冻结基线为准，不在本轮替用户决定并行度。

### 7.4 接入验收

先独立比较 Shoup 与精确宽乘积取模，覆盖 x>=p、最大权重、近2^62模数和无效列；随后检查模乘综合 DSP/LUT、latency 与 II。
再验证独立 BConv、整 digit、OpenFHE Release/ASan、真实 fixture RTL co-sim，逐 residue 精确一致，不能只看 CKKS 解密误差容限。
回归至少保留既有 472 checks/1523712 residues 与 40960 residue RTL fixture 覆盖，并更新协议错配/长度/参数缓存失效用例。
资源验收按物理可达实例去重计数，确认两入口共用 Shoup、旧 BConv Barrett 池不再保留，变换/预乘共享关系不退化。
周期、INIT/预计算冷启动、Host打包开销分别报告；最终频率和性能用新布局布线结果，不能套用纸面 DSP 估计。
