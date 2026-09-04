# HKS 片上存储与数据搬运优化实施计划

日期：2026-09-04。状态：P0/P1/P2 已完成并通过验收（暖态 139734→120822，报告见 docs/reports/hls/hks_mem_pX_20260904/）；P3/P4 留到下一轮。

## 0. 已确认决策（2026-09-04 与用户澄清）

1. 对比锚点：以历史基线 139734 周期为准；P0 重跑当前源码所得数字仅用于解释与历史基线的差异，差异未解释前不开始优化。各阶段实际节省按同源码前后差值（P0 重跑→P1→P2）计量。
2. 期限：约两天，本轮只推进 P0/P1/P2；P3（统一 A_work 共享存储/bank 拓扑改造）、P4（预乘复用 4 路模乘器）留到下一轮。
3. 通用 OP_BCONV 同样迁入统一 A_work（AXI→A_work→systolic→A_work→AXI），其 local 数组与 HKS 路径一并删除，共用同一运行时执行入口与 active_q 接口。
4. A_work 按 8 塔地址跨度预留（塔0..6 + scratch 1 塔 = 256KiB），实际同时占用 ≤6 塔；P_work 物理塔号固定从 LIMB_Q 起，不做地址压缩。
5. 基线检查点：2026-09-04 经用户授权，将验证过的无 AUTO 源码（含地址位选择改写）及其证据报告提交 Git，作为防改错检查点（提交 4b7eabf）；后续阶段确需新检查点仍先向用户申请授权。

## 1. 先明确两项设计含义

### 1.1 预乘采用原位计算，但不是所有算子都原位

预乘的目标是对有效输入塔执行：

`Q_work[q][k] ← Q_work[q][k] × QHatInv[q] mod MODULUS[digit_start + q]`。

同一个地址先读、经模乘流水线、再写回；不再为 compact/预乘结果额外分配完整数组。
这只消除中间副本，不消除必要的读、模乘、写，也不保证与 INTT 零开销重叠。
首版保留“INTT 完成 → 原位预乘 → BConv”顺序，降低数学与调度同时变化的风险。

- NTT/INTT 仍需 A/B 两份工作区进行级间 ping-pong，不改成危险的单数组覆盖。
- BConv 读取 Q 区、写入不重叠的补基区；不覆盖尚未消费的输入。
- 原 EVAL 输出旁路必须保留；若 INTT 覆盖了输入，就必须另存或安全地提前写到最终输出。

### 1.2 Q_work 同时是共享工作存储和 NTT/INTT 的 A bank

“共享”描述多个阶段分时访问同一存储；“NTT buffer”描述该存储在变换阶段的角色。
两者并不矛盾。目标不是 `Q_work → 独立NTT A buffer`，而是直接把某个 Q_work 塔作为 A。
P_work 同理：BConv 写完之后，NTT 直接在该塔与共享 scratch 之间迭代。

Q_work/P_work 是拟议名称，当前源码尚未实现。分时共享也不意味着多个阶段可无仲裁地并发访问。

## 2. 范围、目标与非目标

### 固定范围

- CKKS/DC，N=4096、当前 Q=3、P=2，单 digit ModUp；性能比较固定两个 digit：alpha=2/start=0 与 alpha=1/start=2。
- 保持一套双向 transform、4 个蝶形 lane、BConv 3×5 阵列、三个 256-bit AXI 接口。
- 不恢复 AUTO；指令7继续退役。默认保留 INIT、ADD/SUB/MULT、NTT/INTT、BCONV、HKS_DIGIT 的现有接口和语义。
- 不修改 OpenFHE 数学语义、单位根、EVAL 位序、模数范围与非法描述符无写入约定。

### 必达目标

1. 取消无效输入行的整塔清零与装载；在消费者入口显式注入0。
2. 消除 HKS BConv 的 `Load_X/Store_X` 全量中转及其重复 local 数组。
3. 让工作区直接作为 transform A bank，消除 HKS 的独立 `TRANSFORM_LOAD/STORE` 中转。
4. 预乘仅处理有效塔并原位执行；进一步分时复用现有4路模乘器，不新增物理算术 lane。
5. 功能、接口、RTL结构、周期、资源及物理时序均有可追溯证据。

### 非目标与边界

- 本轮不要求完整 HKS、KeyMult/ModDown 融合、Q=2 降层硬件支持或跨调用驻留。
- 不将“减少复制”写成“零访存”，不以浅 FIFO 取代必须保存的完整塔。
- 不默认删除独立算子指令。若兼容路径使资源无法回收，单列 HKS-only kernel 方案供用户选择，不擅自改 ABI。
- 无板卡：不能承诺或验收 PCIe/HBM/XRT 端到端加速比；保留 INIT 成本单列。

## 3. 基线与证据

以 `docs/reports/hls/hks_no_auto_20260904/` 为历史基线证据目录。
`480dc91` 是删除 AUTO 前检查点，不能冒充当前无 AUTO 源码检查点。
当前工作区还包含未提交的 AUTO 删除和地址位操作改写；旧 P&R 对应位操作改写前源码。
实施前必须冻结当前快照并重新建立同源码基线，不能混用源码、RTL和物理结果。
（已确认：对比锚点仍为历史 139734；当前源码重跑数字仅用于解释差异，差异未解释前不开始优化。冻结快照按用户要求改为 Git 提交检查点，另行记录提交哈希。）

| 项目 | 历史基线 |
|---|---:|
| 两个 digit RTL 周期 | 139734（70131 + 69603） |
| INIT RTL 周期 | 295063 |
| 10次变换外层加载/写回 | 20520周期 |
| 两次BConv加载/写回 | 6664周期 |
| 两次预乘阶段，含无效行写0 | 24624周期 |
| HLS BRAM_18K / DSP / FF / LUT / URAM | 688 / 1160 / 78323 / 147680 / 96 |
| 物理 BRAM Tile / DSP / Register / CLB LUT / URAM | 348 / 1160 / 61076 / 91520 / 96 |
| 默认6ns WNS | +0.078ns |
| 6ns + 0.75ns setup uncertainty WNS | -0.672ns |
| 同布线7ns + 0.75ns setup uncertainty WNS | +0.328ns |

周期分项来自 HLS 层级归属，不是逐拍波形 profiler；最终总周期以同 fixture 的 RTL co-sim 为准。
默认6ns通过不等于带工程裕量通过；物理结果是 OOC 内部时序，不是平台 shell/板级签核。

## 4. 先制定计算访问契约

| 使用者 | 稳态读需求 | 稳态写需求 | 布局/生命周期要求 |
|---|---|---|---|
| NTT/INTT，4蝶形lane | 8系数 + 4 twiddle/拍 | 8系数/拍 | 源、目的区分离；当前8-bank T2P；12级反复使用一塔 |
| 预乘，首版 | 1有效系数/拍 | 同地址1系数/拍 | 使用源qi与digit-local逆元，写地址必须随流水线延迟 |
| 预乘，共享4路版 | 4有效系数/拍 | 同地址4系数/拍 | 与蝶形互斥执行，共用同一组4个模乘器 |
| BConv | 每q一项，最多3项/拍 | 每有效p一项，HKS为2～4项/拍 | Q/P地址不重叠；输入、输出按systolic时序错位 |
| AXI输入/输出 | 每方向理想4系数/拍 | 每方向理想4系数/拍 | 256-bit burst；真实吞吐含握手与背压 |

NTT读 `[i,i+2048]`、写 `[2i,2i+1]`；INTT反向。每块A/B按系数低3位划为8bank，
NTT读侧/INTT写侧可能同bank两访问，因此不能把 T2P 随意替换为单读单写 RAM。
INTT末级每拍输出8系数，而共享4路预乘每拍只能处理4系数，且蝶形阶段模乘器正在工作；
不允许把“末级直连预乘”宣传成无额外周期。

## 5. 目标存储组织与所有权

### 首选组织：统一 A_work 地址空间 + 一塔 scratch

为减少动态切换 Q/P 两组物理数组带来的接口特化，优先实验：

- `A_work[MAX_LIMBS][RING_DIM]`：系数维 cyclic=8、T2P；Q_work/P_work是其中不重叠的逻辑区。
- HKS输入 Q_work 使用物理塔0..alpha-1；源模数/twiddle索引仍为 digit_start+q。
- HKS补基 P_work 使用物理塔 LIMB_Q+p；目标模数/twiddle索引为 `p<digit_start ? p : p+alpha`。
- `scratch[RING_DIM]`：cyclic=8、T2P，只服务当前单塔transform，所有塔分时复用。
- 初期保留输出旁路存储及兼容算子必需的第二操作数存储，逐阶段做生命周期分析后回收。

Transform控制中分开保存“物理塔号”和“全局模数号”，禁止把compact后的塔号当作modulus索引。
偶数stage读A写scratch、奇数stage反向；当前STAGE=12，最终结果回到A。
增加编译期偶数级约束，或明确实现奇数级结果所有权切换，不用隐含假设掩盖未来参数变化。

| 数据 | 产生阶段 | 最后消费者 | 可覆盖时点 |
|---|---|---|---|
| 原digit EVAL | AXI输入 | 最终输出旁路 | 已安全写出且不存在输入别名危险后 |
| Q_work | 输入、INTT、预乘 | BConv | BConv全部输入消费完毕后 |
| P_work | BConv | NTT、最终输出 | 对应塔写出完成后 |
| scratch | 当前变换stage | 下一个stage | 相关流水线排空后，可供下一塔使用 |
| weights/inv/modulus | INIT或digit metadata | 当前digit所有计算 | digit结束或显式重新初始化后 |

INTT逐塔完成，但BConv同时消费多个Q流；BConv同时生成多个P流，但NTT逐塔处理。
因此Q/P塔的保存是实际需求，而复制到另一套同容量数组不是必需条件。

如果统一 A_work 的宽MUX/深BRAM导致时序不佳，再比较独立Q/P局部bank方案，不能只以源码更简洁选型。
纯HKS数据容量可以讨论Q3塔+补基最多4塔+scratch1塔=256KiB；这不是保持全部通用指令时的整核BRAM预测。
（已确认：按此 8 塔地址跨度预留，实际同时占用 ≤6 塔；P_work 物理塔号固定从 LIMB_Q 起，不做压缩。）
独立OP_BCONV最多5个输出塔，兼容模式不得误裁为4个。

## 6. 分阶段执行

### P0：冻结基线与访问契约

**目标**：建立可复现、同源码的前后比较。

**具体操作**：

1. 检查工作区差异，记录HEAD、源码diff/哈希、fixture哈希、工具版本、器件、时钟与测试参数。
2. 创建独立实验标签，如 `mem_p0_r1`（Vitis HKS_PROJECT_TAG，非 Git 标签）；不得覆盖 `no_auto_r1` 和旧报告。
3. 重新运行当前源码native/C-sim、结构审计及真实OpenFHE RTL fixture，确认位操作改写后的周期。
4. 列出每块数组、可达实例、bank数、RAM类型、生产者/消费者；核实BConv wrapper BRAM（已核实：单实例 Compute_BConv_Systolic=128 BRAM_18K，local_in_x 48+local_out_x 80，OP_BCONV与HKS两个调用点共用；不是“两套各128”）。
5. 需要Git检查点时仅提交明确相关文件并先确认授权；2026-09-04 已按用户要求完成基线提交（见第 0 节决策 5）。

**验收**：测试全部通过；源码/RTL/报告映射唯一；当前与历史139734周期若有差异必须先解释，未解释前不开始优化。

### P1：消除无效塔清零与装载

**目标**：不增加并行度，消除alpha以外的无用整塔处理。

**实现方法与操作**：

1. 为BConv内部接口增加 `active_q`；HKS传alpha，独立OP_BCONV传LIMB_Q，更新全部调用和测试。
2. `Prepare_HKS_BConv_Input` 只循环q<alpha，暂保留现有compact布局。
3. 暂时保留wrapper时，Load_X也仅装载active_q行。
4. 在Feed_X中先判断q<active_q及地址有效，再读取RAM；无效行直接注入0，不能靠“旧值乘零”替代有效性约定。

**验收**：alpha交替、毒化旧缓存、非法参数测试通过；RTL不再整塔清零/装载无效行；II不退化；DSP不增加。
独立预测：两digit主体少约12288次写零+1536次装载迭代≈13824周期（约9.9%）。
实际节省须显著为正；与模型偏差超过10%必须逐项解释，不按预测填报告。

### P2：BConv直接读写工作区，取消局部中转

**目标**：去掉 `poly_buffer_1→local_in_x→local_out_x→poly_buffer_1` 的两次整块复制。

**实现方法与操作**：

1. 把systolic核心接口改为工作区+输入/输出塔映射；直接从Q区Feed_X，直接向补基区Collect。
2. 删除大数组local_in_x/local_out_x及Load_X/Store_X；保留小型weights/modulus寄存器与阵列流水寄存器。
3. 通用BCONV与HKS使用同一个运行时执行入口，避免调用点特化重新生成两套缓存/计算控制。
   （已确认：通用 OP_BCONV 同样迁入 A_work 直接读写，其 local 数组与 HKS 路径一并删除。）
4. 编写bank访问检查器：遍历输入/输出索引、stage/方向，检查地址范围、每bank逐周期R/W数量。
5. 逻辑循环中读地址为t-q、写地址为t-(3+p)，但HLS流水线延迟会改变实际同拍相位；
   必须用调度报告/RTL实际RAM使能复核，不能只靠这组下标的模8错位证明无冲突。
6. 若推不出稳定II=1，先修正bank/相位或加局部寄存器切分，不恢复全量复制来掩盖问题。

**验收**：systolic II=1；输入输出精确对拍；两类调用不再有完整local_in/out存储；
模乘阵列不复制，DSP≤基线；HLS BRAM严格下降，BConv总周期下降。
旧wrapper复制归属6664周期是P0基线的独立上界项，与P1已经删掉的无效装载不得重复相加。
P&R检查旧BRAM→local_in_x关键路径确已消失；新路径仍需单独验收。

### P3：工作区直接承担NTT A bank，预乘原位化

**目标**：取消transform外围整塔LOAD/STORE，并取消预乘compact副本。

**实现方法与操作**：

1. 将相关数据区统一为A_work逻辑视图；AXI输入直接进入compact Q区，模数索引保持独立。
2. 修改 `Execute_Transform` / `CG_Transform_Banks`：接受工作塔号与一份scratch，不再定义独立bank_a。
3. 让每级直接在A_work当前塔与scratch之间ping-pong，移除HKS可达的TRANSFORM_LOAD/STORE循环。
4. INTT后，预乘先按现有单lane实现对Q_work原位读乘写；仅访问q<alpha。
5. 写地址、模数号、valid与计算延迟对齐；保留真实RAW依赖，不能沿用旧代码的DEPENDENCE false来掩盖新冲突。
6. BConv直接消费原位预乘后的Q_work；补基NTT直接处理P_work。
7. 梳理ADD/SUB/MULT和独立NTT/INTT所需存储，复用现有空间；不得新增A_work后把三套旧poly/result数组全留着。

**验收**：单套transform、一个共享scratch；正反向蝶形II均为1；HKS路径无外层整塔复制；
所有合法alpha/start、独立算子偏移与原EVAL旁路正确；BRAM和总周期较P2下降。
P0归属20520周期为模型参照，不包含仍必须进行的AXI读写或原位预乘。
若动态塔MUX导致时序恶化，比较更局部的bank组织，不能以“少了数组”直接判定成功。

### P4：预乘复用现有4路模乘器

**目标**：原位预乘达到4有效系数/拍，并消除独立预乘模乘器。

**实现方法与操作**：

1. 为共享计算引擎增加互斥的BUTTERFLY/SCALE阶段，单一物理4路乘法资源池。
2. SCALE阶段每lane处理一个系数，输出写回原地址；不得简单新增四处MultMod调用就称为共享。
3. 各阶段切换前排空乘法流水线；对齐4路地址/valid、检查尾部和模数切换。
4. 保持INTT→SCALE→BConv串行；把预乘移到INTT前或融合末级作为后续可选实验，不混入首版。
5. 若以后尝试前移，先独立验证每塔标量乘与INTT交换关系，并保持原EVAL旁路不被缩放。

**验收**：SCALE II=1、4系数/拍；RTL追踪同一4个模乘实例在两个阶段被使用；
整核MultMod预期20→19、DSP预期1160→1102，映射差异必须解释；不得以吞吐提高换来新增4路硬件。
有效预乘主体从3N≈12288次串行迭代降到3N/4=3072拍，另计填充和切换；不宣称零周期。

### P5：最终输出/旁路融合与物理收尾

**目标**：减少result中转，保留burst，并完成可采用时钟下的整体验收。

**实现方法与操作**：

1. 先审计Host/C-model/API的mem_in1、mem_in2、mem_out别名约定。
   提前写出可能覆盖尚未读完的输入；在证明不重叠前保留原旁路buffer，不能隐式收紧ABI。
2. 安全情况下，输入加载时把原EVAL写入最终全局塔位置；补基NTT完成后直接从P_work做4-wide AXI drain。
3. 各输出塔恰好写一次，非法调用仍零写入；检查32-byte设备对齐与元数据非对齐偏移。
4. 保留固定长度单塔burst helper等已验证机制，检查是否退化为64-bit AXI。
5. 若残留长距离RAM连接，加入有证据的寄存器/FIFO切分并核对II；不能靠普通临时变量保证物理流水寄存器存在。
6. 运行OOC P&R和严格post-route检查，整理最终源代码、RTL、报告、时钟、周期与资源清单。

**验收**：三个AXI实际256bit、burst无退化；功能/结构全通过；所有网络布通；
7ns+0.75ns setup uncertainty下WNS≥0、TNS=0、WHS≥0、THS=0，且时钟/约束覆盖完整。
6ns+0.75ns闭合为进取目标；未达到时必须明确标为未达到，不能用默认6ns结果代替。
最终同fixture暖态周期至少比P0减少10%，物理BRAM低于基线、DSP不增加；
同一通过时序检查的时钟下耗时必须下降。未达目标则报告缺口并保留上一通过阶段，不标“完成优化”。

## 7. 验证矩阵与操作入口

### 每阶段必要验证

| 层级 | 检查内容 | 通过标准 |
|---|---|---|
| 独立数学参考 | NTT/INTT、预乘、CRT、全digit | 所有模数元素精确一致，不只检验往返 |
| native/C-sim | 原Top18/18、HKS22及新增用例 | 原覆盖不得减少；随机、零、q-1、distinct/大素数全部通过 |
| 状态/保护 | alpha1/2/3与所有start、交替调用、毒化缓存、重复INIT | 零旧数据污染、canary不变、非法描述符零写入、退役opcode7不变 |
| OpenFHE集成 | Release和ASan、真实twiddle/fixture、旋转解密与降层回退 | 不低于现有472 checks/1523712 residues；ASan含泄漏零错误 |
| RTL co-sim | 真实fixture40960 residues、扩展smoke35调用/4种digit及新增用例 | 零差异；记录INIT和每digit周期；不能只以C-sim代替 |
| 结构/调度 | bank数/类型、RAM R/W使能、实例数、II、AXI位宽 | 实现与本阶段契约一致，无隐藏复制/算术特化 |
| 物理实现 | 实际资源、布通、严格setup/hold、关键路径 | 按P5物理标准；OOC与板级边界清楚 |

新增必测：原位写地址与乘法延迟对齐、输出旁路精确不变、独立BCONV sizeP=1..5、
Q/P物理塔号与全局modulus号错配反例、stage切换排空、AXI随机背压与别名处理。
复用的scalar参考不得调用被测的bank寻址helper，否则可能同错同过。

### 已有命令（在WSL的 src/fpga_backend 下串行执行）

以下以P1为例；后续各阶段换唯一标签，不覆盖既有项目。

```sh
make test-hks-digit
HKS_PROJECT_TAG=mem_p1_r1 make hks-csim
HKS_PROJECT_TAG=mem_p1_r1 make hks-csynth
HKS_PROJECT_TAG=mem_p1_r1 make hks-cosim-smoke
HKS_PROJECT_TAG=mem_p1_r1 \
HKS_RTL_FIXTURE="$PWD/build/hks_perf_20260904/openfhe_rtl_fixture.txt" \
make hks-cosim-perf
python3 check_shared_transform.py \
  Solution/hks_digit_cosim-perf_mem_p1_r1/solution1 \
  --axi-width 256 --lanes 4 --total-multipliers 20 --no-auto
make hks-impl HKS_IMPL_PROJECT=Solution/hks_digit_cosim-perf_mem_p1_r1
make hks-postroute HKS_IMPL_PROJECT=Solution/hks_digit_cosim-perf_mem_p1_r1
```

P4成功共享后审计数量改为19，并新增同实例时分使用检查；现有脚本可能依赖旧层级名，
应同步更新结构断言，不能删断言“让测试通过”。bank逐拍检查器是待新增项目，当前没有现成命令。
OpenFHE/ASan/XRT编译沿用既有报告中的构建入口，执行前核对当前CMake选项与测试二进制路径。

## 8. 改动文件与交付物

| 文件/目录 | 计划操作 |
|---|---|
| src/fpga_backend/src/top.cpp | 生命周期调度、物理/逻辑塔映射、原位预乘、输入输出交接 |
| src/fpga_backend/src/cg_ntt.cpp、include/cg_ntt.h | 共享工作bank接口、scratch及SCALE阶段 |
| src/fpga_backend/src/bconv_systolic.cpp、include/bconv_systolic.h | active_q、直接Feed/Collect、删除整块local副本 |
| src/fpga_backend/src/load.cpp及相关头文件 | 必要时适配工作区/旁路，保持burst |
| src/fpga_backend/testbench/ | 全部兼容回归、bank/lifetime/原位新增测试 |
| src/fpga_backend/check_shared_transform.py | 对新层级做严格可达实例/位宽/共享审计，保留旧设计反例 |
| src/fpga_backend/HKS_DIGIT.md | 仅在实现验证后更新真实内存和设备执行约定 |
| docs/reports/hls/hks_mem_pX_日期/ | 每阶段中文报告、原始证据、源代码/fixture哈希与预测偏差 |
| AI_Cowork/log.md | 记录每阶段方法、试错、验证、未通过项；不代写个人学习笔记 |

当前用户指定计划放在 `doc/`（单数），因此新建该目录；不移动原有 `docs/` 报告。

## 9. 执行节奏与停止条件

- P0→P1→P2→P3→P4→P5，阶段独立验证，失败时停在本阶段修正，禁止把多处未经验证改动一起综合。
- 先跑native/小规模访问检查，再HLS，再RTL；P2去长连线、P3改bank拓扑、P5最终版本应做P&R。
- 已确认本轮约两天：仅推进P0/P1/P2，P3/P4留到下一轮；不能以期限替代验收。
- 周期减少但频率下降、存储减少但LUT/MUX大增、指令兼容性丢失，都不视为无条件成功。
- 需要删除兼容指令、改变别名规则、增加物理lane或扩大到完整HKS时，先请求用户确认。
- 交付中文总结必须区分：完成/未完成、实测/模型、HLS/物理资源、RTL时间/上板墙钟；禁止合并收益重复记账。

最终判定：同一份工作数据被不同阶段接续使用，而不是每进入一个C++函数就重新复制一份。
