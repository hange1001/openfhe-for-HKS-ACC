# AI 协作日志

> 记录每次 AI Agent 会话的执行过程、结果和注意事项。
> 每次会话追加一条记录，包含日期、任务、AI 执行摘要和验证结果。

---

## 会话记录模板

```markdown
### [日期] 任务简述

**Agent**: Claude Code (Opus 4.x)

**任务**: 简要描述

**执行步骤**:
1. 步骤 1 - 结果
2. 步骤 2 - 结果

**修改文件**:
- `path/to/file` - 修改说明

**验证结果**:
- [ ] 编译通过
- [ ] 测试通过
- [ ] 代码审查通过

**注意事项 / 踩坑记录**:
- 发现的问题和解决方案
```

---

## 会话记录

<!-- 新会话记录追加在下方 -->

### [2026-09-04] 编写共享工作存储与原位预乘实施计划

**Agent**: Codex。按用户要求在新建的 `doc/` 目录保存中文实施计划，未改动硬件实现。
澄清预乘是同地址读乘写而非零访存；Q_work既是共享工作区，也直接承担NTT/INTT的A bank。
计划包含P0基线、P1无效行、P2直接BConv、P3工作bank/原位化、P4模乘复用、P5输出和物理验收，
明确接口兼容、别名风险、实际流水线bank冲突、周期收益不能重复相加和无板卡边界。
交付：`doc/HKS_片上存储与数据搬运优化实施计划.md`；本次未运行综合、未提交Git。

### [2026-09-04] 从计算单元访问契约反推 HKS 存储与搬运

**Agent**: Codex。用户指出不能继续按软件函数边界从Top复制到子函数，要求解释现有搬运、
计算单元需求及可优化空间。本轮只做代码/报告审计和架构推导，不修改硬件源码。

**当前真实链路**:
- AXI输入同时写 `poly_buffer_2` 与 `result_buffer`；每次变换显式执行
  `poly_buffer→bank_a→A/B十二级ping-pong→poly/result_buffer`。
- INTT 后 `Prepare_HKS_BConv_Input` 将全局limb压紧、乘QHatInv并写 `poly_buffer_1`；
  每digit固定扫3塔，不足alpha的行也整塔写0。
- BConv又执行 `poly_buffer_1→local_in_x→systolic→local_out_x→poly_buffer_1`，
  随后每个补基塔再复制进transform bank做NTT。
- 两digit可明确归属的纯片上LOAD/STORE为transform 20520cycles、BConv 6664cycles，
  合计27184cycles（19.45%，6ns为163.104us）；这是优化候选，不等于可直接全额扣除。

**计算访问契约**:
- NTT/INTT：4蝶形/拍，需8个系数读、8个系数写和4个twiddle读；12级之间必须保留
  A/B ping-pong或等价双工作区。可删的是外层整塔装入/倒出，不是级间工作集。
- prescale：当前1系数/拍、1个模乘；两digit固定写6塔但有效仅3塔，另12288次为写0。
- BConv core：每q行每拍1个输入，共3个；每有效p列每拍1个输出，HKS为3/4个。
  `LOAD_PAR=8`仅服务额外整塔复制，不是core吞吐需求。
- AXI：256bit=4系数/拍，与transform内部512bit访问不同；需要bank解耦，但无需每个函数
  都再拥有一份完整数组。

**资源与物理证据**:
- 通用OP_BCONV与HKS调用点生成两个BConv wrapper；各自 `local_in_x/local_out_x`
  占128 BRAM_18K，总计两套256 BRAM_18K。算术被上提共享，不能从wrapper DSP=0误判无计算。
- 当前最差路径就是 `poly_buffer_1 BRAM→HKS BConv local_in_x BRAM`：5.618ns，
  route占83.86%。说明为复制而造的512bit片内连接同时损害BRAM与时序。

**建议的数据所有权方向**:
- 由HKS工作区而非C++函数拥有数据：`Q_work[3][N]` 保存INTT/BConv输入，
  `P_work[5][N]` 保存BConv/NTT输出，另保留一塔共享transform scratch。
- transform直接把Q/P work当A bank，scratch当B bank；12级为偶数时结果回到A，
  从而取消transform外层LOAD/STORE。
- BConv直接读Q_work、写P_work；先用地址skew证明bank冲突，再由HLS II报告和P&R验收。
  不能只依靠C++形参或pragma猜端口。
- 原EVAL旁路是输出语义，不能删值；可在AXI输入时直接写最终输出区域。外部ModUp输出仍
  必须传输，只有继续融合KeyMult/ModDown并片上驻留才能消掉。
- 不用浅FIFO替代全部塔：BConv需多Q塔同时存在，单套NTT又依次消费多个P塔；这些生命周期
  决定至少要保留一份完整Q/P工作集。

**实施顺序建议**: 先去无效行清零/装载并加alpha有效门控；再去BConv local副本；
再让work bank直接参与NTT；最后处理AXI旁路/更大范围算子融合。每步都检查正确性、
MultMod实例数、BRAM、周期、II与post-route WNS。

### [2026-09-04] 去 AUTO 物理实现完成，并定位 CPU/FPGA 性能瓶颈

**Agent**: Codex。用户要求列出 CPU/FPGA 运行数据并分析慢点，同时追问“变换”是否仍是 AUTO。
继续按 mentor 区分实测、RTL 换算和预测；只更新工程记录，不改个人学习笔记。

**执行与结果**:
- 完成 `hks-impl` 及独立 `hks-postroute`。239992 条网络全部布通，未布线网络为0；
  物理资源为 CLB LUT 91520、Register 61076、BRAM Tile 348、DSP 1160、URAM 96。
- 相对含 AUTO 检查点，物理 LUT -15951（-14.84%）、Register -2986（-4.66%）、
  BRAM Tile -11.5（-3.20%），DSP/URAM 不变。
- 默认6ns：WNS=+0.078ns、TNS=0；6ns加0.75ns setup uncertainty：
  WNS=-0.672ns、TNS=-332.331ns、1638个失败端点；7ns加同样裕量：WNS=+0.328ns。
- 最差 setup 路径为顶层 poly_buffer_1 RAMB36E2 到 BConv local_in_x RAMB36E2，
  数据路径5.618ns，其中route 4.711ns（83.86%）、logic 0.907ns（16.14%）、0逻辑级。
  无level>5有效拥塞窗口、无跨SLR网络；256位端口可以布通，但内部BRAM宽搬运仍有走线风险。
- CPU同工作负载暖态500次：OMP1 median=0.640583ms，OMP2 median=0.463148ms。
  FPGA两个digit为139734cycles：6ns名义0.838404ms，严格裕量通过的7ns为0.978138ms；
  均不含主机/PCIe，因此不是已上板端到端加速比。
- FPGA周期归属：3次INTT+7次NTT共85260（61.02%），QHatInv预缩放24624（17.62%），
  两次BConv 14924（10.68%），顶层结果写回10370（7.42%），其余输入/参数/控制约3.26%。
  这里的“变换”不是AUTO；AUTO已删除，原本也不在OP_HKS_DIGIT调用链中。
- 将 `top.cpp` 的 `/SQRT`、`%SQRT` 显式改为 `>>LOG_SQRT`、`&(SQRT-1)`，增加幂次
  `static_assert`；原生Top/HKS回归18/18、HKS 22 cases通过。P&R使用改写前但综合语义等价的
  源码，没有把该文本改写冒充为重新物理实现。

**结论**:
- 第一瓶颈是10次NTT/INTT及其局部bank加载/写回，第二是预缩放，第三才是BConv。
  仅变换阶段在6ns下已需0.511560ms，超过CPU两线程整段0.463148ms；只优化BConv不足以反超。
- 原始物理证据以 `physical-*` 归档到 `docs/reports/hls/hks_no_auto_20260904/`；中文报告已回填。
- 删除AUTO后的功能改动与本次位寻址可读性修改仍未提交或推送。

### [2026-09-04] 补齐去 AUTO 性能口径，并记录模乘器复用的改动前预测

**Agent**: Codex。用户要求预估复用收益，并追问去 AUTO 的运行时间、加速比和 slack。
按 mentor 区分测量与预测；按 mindmap-learning 只更新工程报告/日志，不代写学习笔记。
本轮不实现复用，不改变正在物理实现的源码，不提交或推送。

**核验结果**:
- 去 AUTO 的真实 OpenFHE RTL 周期仍为139734；6ns条件换算暖态0.838404ms，
  INIT1.770378ms，INIT+两digit2.608782ms；相对删除前同频1.000x。
- 旧CPU暖态基线0.4631475ms，CPU/名义RTL时间比0.5524x（FPGA约慢1.81倍）；
  不是完整HKS/旋转性能，未含PCIe/驱动，不把旧CPU基线称为本轮重测。
- HLS预算余量为0.000ns；新P&R仍存活并进入Detail Placement，routed DCP尚未生成。
  新布线WNS/TNS缺项是尚未产出，不能用旧6ns+0.75ns的-0.779ns、旧7ns的+0.221ns替代。
- 当前预乘模块：12312cycles/digit、DSP58/LUT3948/FF4049；MultMod本体DSP58/LUT3283/FF1757。
  两digit的预乘占17.622%，每次遍历固定3N含无效行清零。

**改动前预测（不是实现结果）**:
- 仅共享一路、维持原调度：DSP目标1160→1102；整核周期大致不变，可能略增控制开销。
- 复用现有四路并合并INTT写回：独立四路遍历模型约121302cycles；融合写回理想模型
  115110+overhead。工程预测116000～123000cycles，即减少12%～17%、同频1.14～1.20x。
- 明确DSP目标-58（-5%）；LUT净减1～3k为低置信度预算，可能被新增MUX抵消；
  BRAM/URAM先按不变，FF与slack不报未经验证的硬指标。模乘器物理数应20→19。
- 模型要求共享后II=1、四字/拍访存、BConv无效行显式补零；不假设新增算术并行度。
  当前INTT逐级模二分，不存在可直接免费合并的尾部N逆元模乘。
- 6ns条件下预测0.696～0.738ms，对旧CPU仍只有0.63～0.67x；不声称单项优化后超过CPU。

**记录位置**: `docs/reports/hls/hks_no_auto_20260904/README.md` 新增完整运行时间/比值/slack表
及复用模型，归档预乘与变换原始综合报告；仅工程文档更新。

**同轮追加：回应地址除余的硬件代价质疑**:
- SQRT为编译期64，k范围0～4095；逻辑上`/64`是高位、`%64`是低6位。
- 实际生成LOAD/STORE RTL经展平、bank映射后使用i[9:1]/i[9:0]/i[0]和limb拼接；
  核对地址赋值与模块层级，没有该处的除法/取余器。不是仅凭0DSP判断。
- LOAD=1026cycles/628LUT/24FF，STORE=1026cycles/327LUT/25FF，均0DSP；
  这里真正的代价是片上搬运和选择逻辑。没有以性能优化名义做等价语法改写。
- P&R后续已进入Initial Net Routing，尚无最终routed DCP；中间未布通网络数不代表最终失败。

### [2026-09-04] 先保存检查点，再移除 FPGA AUTO（功能验证完成，物理验证运行中）

**Agent**: Codex。按用户明确要求执行；工程报告使用中文，保留个人学习笔记与学习路线，
不改 OB/inbox、不推送。沿用 ai-cowork 工程记录，遵守 mindmap-learning 的字段边界。

**执行顺序与范围**:
1. 删除前重新运行 native Top/HKS，通过后提交 `480dc91`：
   `feat(fpga): checkpoint shared 256-bit HKS before AUTO removal`。
   含宽接口/共享引擎/去拷贝、中文报告及原始证据；80 个文件，提交后工作区干净。
2. 删除 `src/fpga_backend/src/auto.cpp`、`include/auto.h`、`testbench/auto_test.cpp`
   和 `src/pke/examples/test-ckks-auto.cpp`；均可从检查点恢复。
3. 移除 Top AUTO/元数据分支、构建源项、AutoOffload 和 OpenFHE 卸载分支。
   保留 CPU 自同构、旋转、HKS 算子；不改 HKS 数据通路、PE=4 和 256 位接口。
4. 指令 7 退役，其余编号不变。主机在访问设备前拒绝 7；硬件不访问外部内存，
   保持上下文；旧 AUTO 测试替换为退役保护，然后继续合法算术/HKS 调用。
5. 新增 `--no-auto` 可达 RTL 审计；同时用旧版作反例，确保不是恒通过检查。

**验证结果（本轮新运行）**:
- [x] native Top 18/18、HKS 22 合法用例和非法参数/边界检查。
- [x] OpenFHE Release 与全库 ASan（含泄漏检查）各 472 checks、1,523,712 个元素精确一致；
  正反向 EvalRotate 解密正确，EvalMult 降层 CPU 回退仍通过。
- [x] XRT 启用主机分支语法编译通过，未运行 XRT 或连接板卡。
- [x] 真实 OpenFHE fixture RTL 40,960 个元素精确一致；扩展 smoke 35 次调用、4 种 digit 通过。
- [x] perf/smoke 的整个生成 Verilog 目录逐文件一致；无 AUTO、一套双向 CG、四路共享模乘、
  整机 20 个 MultMod、三个 256 位端口、两个蝶形 II=1；旧 RTL 被 `--no-auto` 正确拒绝。
- [x] HLS BRAM_18K 708→688，DSP 1160 不变，FF 78528→78323，LUT 169991→147680，URAM 96 不变。
- [x] 同 fixture 周期不变：INIT 295063；digit[2,0] 70131，digit[1,2] 69603，合计 139734。
- [ ] 新 OOC P&R 已启动：`Solution/hks_digit_cosim-perf_no_auto_r1`，
  日志 `src/fpga_backend/build/hks_digit/no-auto-implementation.log`。物理结果未计入通过。

**试错与限制**:
- 首次提交的 whitespace 检查遇原始工具报告行尾空格；仅对源码/手写文档检查，保留原始证据原文。
- 审计反例最初的 shell 匹配串写成 `Retired AUTO`，实际诊断为 `Retired automorphism hardware`；
  修正匹配后反例通过，审计自身正确拒绝旧设计。比较脚本最初输出路径多了一层 `..`，随后修正。
- 删除 AUTO 不会减少 ModUp 的运算周期；资源下降不能直接称为加速。新布局布线完成前，
  不套用旧 7ns 时序，不声明新版 6ns 闭合。当前融合范围仍只是单 digit ModUp，不是完整 HKS。
- 没有板卡，未验证 PCIe/HBM/驱动端到端性能。删除后改动尚未再次提交。
- 最终 `git diff --check` 通过，暂存区为空，`docs/notes` 无本轮改动；已归档软件/RTL
  完整证据，核实测试程序实际位于 `build/bin/tests/`，不是 examples 目录。
- 中文报告与原始证据：`docs/reports/hls/hks_no_auto_20260904/README.md`。

### [2026-09-04] 报告中文化，澄清 AUTO 与 HKS 的边界

**Agent**: Codex。用户要求报告使用中文，并询问 AUTO 在 HKS 中的作用。

**执行结果**:
- 将宽接口优化报告正文完整改为中文，保留原始数值、命令、哈希及验收限制；
  工具生成的 `.rpt` / 日志保留原文，不覆盖或翻译原始证据。
- 报告新增 AUTO 边界说明：它是通用 Top 的独立自同构指令，不是 OP_HKS_DIGIT
  的 ModUp 步骤。本仓库 EvalAutomorphism 先 KeySwitchInPlace 再对两个分量自同构；
  当前融合卸载只接管 EvalKeySwitchPrecomputeCore 中的 digit 预处理。
- 核对 HksDigitOnly：融合模式禁用旧单算子卸载，AUTO 留在 CPU；但通用 Top 中
  仍有 AUTO 电路，因此它参与整核静态时序检查，不参与本次两个融合 digit 的执行周期。
- 记录用户偏好：本项目报告正文使用中文，代码标识和工具原始证据保留原文。

**验证与范围**: 核对 Top 指令分发、OpenFHE 旋转/密钥生成、融合入口及 FpgaManager
的实际代码；仅修改报告和工程记录，没有改动 RTL/算术逻辑、删除 AUTO 或增加 false path，
未重跑硬件测试，前轮测试结果仍对应原设计。按 mentor 分清算法和实现层次；
按 mindmap-learning 只维护工程记录，不代写个人学习笔记，不同步 OB，不提交 Git。

### [2026-09-04] 保持并行度，拓宽 AXI 与优化整塔拷贝（完成；6ns时序未闭合）

**Agent**: Codex。用户明确要求「并行度先不提高，拓宽接口，优化拷贝塔，并软硬件验证」，
并补充「记得写日志」。沿用 ai-cowork 工程记录；按 mindmap-learning 边界保留学习路线、
个人笔记和收件箱，不同步 OB，不新增 Git 提交或推送。

**实现与不变量**:
- PE_PARALLEL=4、NTT 四路蝶形和 BConv 3×5 阵列不变，仍只有一个运行时双向 CG 引擎。
- 加载时广播原 EVAL 到旁路结果；BConv 输出直接进入 CG A/B bank，变换后直接写目标。
  去掉两个 digit 的全部 17 次独立 Copy_HKS_Tower，以及 flat/flat_out 中间数组。
  保留 NTT 双缓冲和 BConv 多行/多列缓存；未增加算术并行度、未加入跨 digit DATAFLOW。
- 三个 AXI 数据端口实际从 64 位变成 256 位。Top uint64 ABI 和数据排列不变；
  设备地址要求 32 字节对齐，元数据偏移仍由 AXI adapter 处理。

**试验与弯路**:
1. r1：四字并行循环加 32 字节对齐，生成端口仅 64/128/64；OpenFHE RTL 40960 residues
   通过，但不接受为拓宽完成。r2 改连续循环、禁止展平，仍为 64/128/64。
2. r3：将运行时 limb 循环与固定 4096 元素传输分离到非内联单塔 helper。
   HLS 才推导出 256/256/256；不能仅凭 max_widen_bitwidth=256 宣称拓宽。
3. 移走内部适配数组后，Vitis 把四个 CG MultMod 提升到 Top，并与 OP_MULT 共享。
   原审计只数 CG 子模块，误报 0；已增加四路外部接口与顶层物理池双重检查。
   整机实际 20 个 MultMod，而非算术消失或并行度增加。
4. OOC实现正常退出不代表时序通过：日志明确6ns失败。首次布线后诊断脚本误用
   Vivado不支持的redirect，已保留失败日志，改为直接输出诊断并完整重跑成功；
   同时去掉已弃用的report_timing -nets，未修改数据通路。

**已完成验证**:
- 最终 native：Top 18/18、HKS 22 valid cases、非法参数与 canary；新增非零塔偏移
  ADD/SUB/MULT/AUTO 宽端口回归；旧 CG 包装器 11/11。
- 最终 OpenFHE Release / ASan（补跑detect_leaks=1）：各470checks、1523712residues一致。
- Vitis全量C-sim 0 errors；扩展RTL smoke 35次调用、4个融合case、0 failures/PASS，
  覆盖alpha=1/2/3、start=1、方向切换、60bit模数、非法调用及ADD/SUB/MULT/AUTO端口。
- 最终性能 RTL：真实 OpenFHE 同一 fixture，INIT+两个 digit，40960 residues 精确一致，PASS。
- 结构：一个 CG 引擎、4 个共享模乘、总计20个 MultMod，三端口256位，无 Copy_HKS_Tower。
- 性能RTL、smoke与P&R输入的三个syn/verilog目录逐文件完全一致；两个蝶形循环均II=1。
- 同配置 HLS：BRAM 704→708、DSP 1392→1160、FF 79656→78528、LUT 175140→169991，
  URAM96不变。NTT适配/缓冲 -32 BRAM、三 AXI adapter 合计 +36 BRAM，净 +4。
- 两 digit 周期 268170→139734（-47.8935%），按6ns为1.609020→0.838404ms。
  CPU双线程旧空闲环境中位0.4631475ms，CPU/RTL估计比0.5524，仍不宣称CPU加速。
- 冷启动 INIT 196759→295063周期（1.180554→1.770378ms），退化必须单列；暖态排除INIT。

**最终物理验证与限制**:
- Vivado OOC P&R已完成：239992/239992条可布线网络全部布通，路由错误0。
  物理资源LUT107471、FF64062、DSP1160、BRAM359.5个36K等效块（719个18K）、URAM96。
  此处是Vivado口径，不能与旧HLS资源直接相减；未做旧版对照P&R。
- 默认6ns：WNS -0.029ns、TNS -0.029ns、1个失败端点；补0.75ns用户setup裕量后：
  WNS -0.779ns、TNS -1574.646ns、6083个失败端点，明确不通过。
- 同一布线检查点只改时钟约束为7ns并保留0.75ns裕量：WNS +0.221ns、TNS 0、
  失败端点0；三种情景hold最差均+0.019ns。源码/HLS的6ns目标未修改。
- 最差路径：poly_buffer_1 BRAM→AUTO选bank/模取负→result_buffer BRAM；
  数据路径5.681ns，其中布线4.137ns（72.821%），不是AXI端口本身的关键路径。
  初始路由有level5/6局部拥塞，最后布通；无旧版物理对照，不能排除拓宽对布局的间接影响。
- 7ns内部时序情景下暖态两digit估计0.978138ms，INIT 2.065441ms；仍慢于CPU0.4631475ms。
  原0.838404ms仅是未达成6ns下的同频换算，不能当实际工作频率下的性能。
- 580个输入、526个输出缺外部延迟，ap_clk缺平台HD.CLK_SRC；内部未约束端点0。
  DRC保留200条DSP流水寄存器告警和1条无可布线负载告警；无DRC错误。
  没有板卡/XRT运行/真实HBM/PCIe/平台shell链接，也未做布线后门级仿真。
- HLS估计5.250ns不等于P&R达频；RTL功能通过+独立静态时序检查不等于整个平台签核。
  所有时间估计不含主机打包、驱动、真实DMA，仅两digit ModUp。
- 历史原始报告未覆盖；源码与本轮报告仍未提交。

**主要文件与证据**: `src/fpga_backend/src/{top,load,cg_ntt}.cpp`、对应头文件、
`testbench/hks_digit_tb.cpp`、`hks_digit.tcl`、`hks_impl.tcl`、结构/性能审计脚本；
复现与报告：`docs/reports/hls/hks_wide256_direct_20260904/README.md`。

### [2026-09-04] 先保存检查点，再统一 NTT/INTT 引擎（完成）

**Agent**: Codex。按用户明确请求实施；ai-cowork / mindmap-learning 仅用于维护工程记录，
保留学习路线和个人笔记。未推送、未重置已有改动。

**检查点**: `c073b11`。重跑原 Top/HKS、OpenFHE、全库 ASan 通过后提交，包含原暂存优化。
WSL 未配置身份，使用已存在的 Windows Git 身份作本次提交，不修改全局配置。
原 `问题回答.md` 尾部空行原样保留；正确性范围是功能模型，不是硬件验收。

**实现**: 独立 OP_NTT / OP_INTT 与 OP_HKS_DIGIT 统一到同一调度循环和工作缓冲。
运行时方向选择地址、twiddle、蝶形前后处理；每路仅一个 MultMod 调用。
保留原接口与测试模板包装；三份冗余变换实例移除，新增每补集塔一次片上拷贝。

**验收**:
- 原生 Top 18/18、HKS 22 cases、CG 11/11；Vitis C-sim 0 errors。
- OpenFHE/全库 ASan：470 checks、1523712 residues 精确一致；包括多 limb、非零偏移、
  独立 NTT/INTT 在每个融合 digit 后交错调用与 context cache 保持。
- 生成 RTL 的结构审计确认一个 CG 引擎、4 个共享模乘；旧多引擎报告被反例检查拒绝。
- Verilog/xsim smoke：24 transactions、3 valid fused cases、0 failures，C/RTL PASS。
- 同配置综合：BRAM 896→704、DSP 2088→1392、FF 106709→79656、
  LUT 278242→175140、URAM 96→96；蝶形循环 II=1。

**限制**: 估计 slack -0.291ns，时序未闭合；RTL smoke 不等于完整 OpenFHE RTL 仿真。
无 XRT 运行/P&R/上板/加速比。BConv 870 DSP、两份局部缓冲和双方向 twiddle 仍保留。
新共享版暂未提交，供用户审阅；证据与复现见
`docs/reports/hls/hks_shared_transform_20260904/README.md`。

### [2026-09-04] 复合算子接入 OpenFHE 并验证（完成）

**Agent**: Codex。按用户「接入openfhe并验证」实施；保留原暂存区与学习路线，不写个人学习笔记。

**实现**:
- 在 HYBRID EvalKeySwitchPrecomputeCore 的 CPU digit materialization 前加入显式 opt-in 接口。
  直接使用 OpenFHE 的 EVAL、单位根、PartQlHatInvModq 和 PartQlHatModp；返回完整 QlP digit。
- 新增 HLS C-model / XRT / CPU 选择、按完整 modulus/root 顺序的初始化缓存、上下文互斥、
  参数/元数据检查及失败时不发布半成品。显式配置后旧逐算子 FPGA hooks 保持关闭。
- XRT 入口按三个不同 BO 大小打包，opcode 8 纳入统计。默认未配置的应用保持原行为。
- 新增可选 CMake/CTest 目标，把真实 OpenFHE 库与真实 Top/ap_int 代码链接起来；无板卡。

**验证**:
- 集成测试 400 checks、991232 residues 精确对比通过。五个模数独立 NTT/INTT 对照；
  7 类 ModUp 输入；两个 digit 全部五塔；缓存复用和 A/B/A 参数切换；输入不变。
- EvalRotate(+1/-1) 均执行两个融合 digit；与 CPU 密文 bit-exact，解密误差远低于 1e-6。
- EvalMult 因 FLEXIBLEAUTOEXT 降到 Q=2，验证的是 CPU 回退正确性，不是融合乘法。
- Q=2/N=8192/P=1/MP/OC/COEFFICIENT/无XRT回退；缺少C-model、旧bitstream无opcode8、
  第二digit传输失败均有检查；失败不交付部分结果，重新配置后正常恢复。
- 全库 AddressSanitizer（detect_leaks=1）通过同一套集成测试。
- 原 Top 18/18、HKS 22 cases 回归通过；XRT-enabled 及 WITH_REDUCED_NOISE 分支语法编译通过。
- 旧 build/unittest/pke_tests 在列测试阶段已因 stoul 异常退出，根目录/构建目录均如此；
  不宣称全量 PKE suite 通过，不在本任务修改无关 CSV 用例加载逻辑。

**被验证修正的假设**: 上轮 cyclic oracle 不代表硬件只能 cyclic。现有 Host
BuildCgTwiddle_Unified 已支持 OpenFHE negacyclic 表，本轮用真实参数对拍确认无须改 HLS
或增加 Host bit-reversal/twist。只缓存初始化参数，没有实现跨调用密文/密钥驻留。

**测量边界**: useful payload 416 KiB + 288 metadata bytes / 两digit；当前诊断 XRT
实现另上传 320 KiB 输出 sentinel，INIT 输入另有 3 MiB+120 bytes，故不等于真实 PCIe 流量。
XRT 仅语法编译，未运行 sw_emu/hw_emu/RTL co-sim/上板/P&R。旧综合时序缺口未解决，
不声称加速比。复现与原始输出在 `docs/reports/hls/hks_digit_openfhe_20260904/`。

### [2026-09-04] OP_HKS_DIGIT 无板卡 PoC（完成）

**Agent**: Codex。用户明确要求开始实现复合算子；按本轮授权实施，不修改学习笔记。

**范围**: 单 digit 原始 EVAL bypass + INTT + QHatInv 预缩放 + BConv + complement NTT；
不接入 XRT/OpenFHE 主链，不含 KeyMult。接口约定见 `src/fpga_backend/HKS_DIGIT.md`。

**综合前预测**（2026-09-04，非上板墙钟）:
- useful payload 口径：两 digit 1056 KiB → 416 KiB；排除 INIT、metadata、dummy BO、协议开销。
- 参考本地 `Solution/Top/solution1/syn/report/` 2026-09-03 报告：同为 Vitis 2023.2、U55C、
  6 ns/uncertainty 0.75 ns；NTT/INTT 包装每 limb 15678 cycles；BConv c=3/4 分别
  由区间推算 7206/7718 cycles（每增加输出列约 512 cycles）。
- 新算子纯片内阶段预测：5*15678 + 5*4096(copy) + 3*4096(scale/pad) + BConv，
  alpha=2: 118364 cycles ±15%；alpha=1: 118876 cycles ±15%。假设 copy/scale II=1、
  transform 单价不因新调用点改变；排除 Load/Store AXI 等待、metadata、控制握手，
  不与 Top 墙钟直接比较。综合若显示不同 II/函数克隆，应分项解释而非调系数。
- 基线资源（历史同约束参考，非本轮重跑）：BRAM_18K=640、DSP=1334、FF=61659、
  LUT=186193、URAM=96。没有新增源代码层面的整环持久数组，不保证 HLS 不克隆子模块。

**已验证**: 修改前 Top 12/12；修改后本机 C++ 与 Vitis C-sim 均为 Top 18/18，
其中 HKS 22 个有效 case 与独立 scalar cyclic-NTT/CRT、分步 opcode 路径逐系数一致；
覆盖 distinct 30/60-bit primes、全部合法 digit 区间、零/边界输入、无效描述符与输出 canary。

**限制**: 现有 Top native EVAL/cyclic NTT 的等价性不等于 OpenFHE negacyclic NTT 接口验收。
本轮保留此前全部 staged changes，新增改动保持 unstaged。

**最终验证/对账**:
- native、AddressSanitizer、Vitis C-sim 通过；新增一键 test-hks-digit / hks-csim /
  hks-csynth / hks-baseline（及尚未运行的 hks-cosim）入口。
- fused 与本轮新跑 baseline 均生成 RTL；比较器只关闭 opcode 8，不回退用户已有优化。
- BRAM_18K 640→896、DSP 1566→2088、FF 69805→106709、LUT 181969→278242、URAM 96→96。
  顶层 Memory 分配不变，新增 BRAM 来自函数克隆/局部 scratch，不是新增顶层数据数组。
- 时钟估计均 5.581ns；6ns - 0.75ns uncertainty - 5.581ns = -0.331ns。
  有 II 警告，不能声明时序收敛。总 Top/HKS latency 是 `?`，不能报总周期或加速比。
- 预测偏差 >15% 已定位：旧报告 CG core 为 14643 cycles，当前源码两种设计均为 8499；
  新调用点另有 reshape 1026 vs 原调用点 514 cycles。当前片内阶段加总为
  90238/90750 cycles（非完整 kernel 延迟）；保留原预测以便追溯，不调系数拟合。
- 归档在 `docs/reports/hls/hks_digit_poc_20260904/`；接口与复现说明见 HKS_DIGIT.md。
  task.yaml YAML 校验、git diff --check 通过。未修改/提交/清理已有暂存区。

### [2026-08-26] D3 出题；查出 `0.32/0.89/2.1` 是无出处数字

**Agent**: Claude Code (Opus 5)　**方式**: 用户读完 D3 三份材料，AI 出题（角色翻转：你画我改）

**任务**: 用户报「D3 已看完，做相关任务」。按 MAP §2 第 3 项，本轮 AI 的职责是出题 + 批改，不产出图。

**过程**:
1. 读 README / MAP / task.yaml / 对表 / 推导v1 / 符号表，确认 D3 的验收题与交付物
2. 收集给定值时逐条核出处 → 撞到下面这个发现
3. 按 step 1.1 `exercise_order` + step 1.4 `warmup_exercise` 的规格出题（给定/映射/假设/过关线四件套）

**关键发现——`B0/B1/B4 = 0.32/0.89/2.1 ops/byte` 是拍脑袋常数**:
- 它标注的出处「推导v1 §四」是**失效引用**：§4 的标题是「口径对账（F2）」，全文 grep 这三个数计数为 **0**
- 全仓库四处引用（task.yaml step 1.3 / 符号表 §9.1 / 对表 §1.4·§3.3 / 本 log 08-11）**互相指认**，
  无一处含推导。git 溯源：2026-08-04 首次写下（676dd65），此后只被复制，从未被推过
- 附带：**`B2` / `B3` 两个边界在全仓库从未被定义过**，只有 B0/B1/B4 以 AI 值的形式露过面。
  「B0–B4 五点上图」这个交付要求，其实有两个点连名字都没有
- 同类：`open_questions` q1 的机器平衡点 **0.22 / 25.6** 也只有结论没有推导，且未声明分母是哪条带宽
- **影响**：D3 验收题 1「把 B0/B1/B4 放上图」不能是查表，必须重推。
  这与 step 1.3 `output_requirement`（必须填符号式）本来就是同一件事 → 已写成练习的第 0/1 步

**另一处口径隐患（画 HBM 屋顶前必须先解决）**:
- `u55c.cfg` 把三个 AXI 端口绑到 `DDR[0..2]`，但 **U55C 是 HBM-only 卡、无 DDR**
- `connectivity.ini` 绑的是 `HBM[0..2]`，但 `Makefile:72` 把 `--config connectivity.ini` **注释掉了**
- → 两份配置互相矛盾，且哪一份都没生效。「HBM 屋顶」当前没有可信的高度值

**修改文件**:
- `AI_Cowork/task.yaml` - step 1.3 新增 `exercise_roofline`（四步题，56 行）与
  `givens_for_roofline`（给定值清单，每条带出处，38 行）；
  `changes` 里那条失效引用改成核查标注；`global_validation.traceability` 存量问题 **+2 条**
- `docs/notes/L0_符号表.md` - §9.1 表中三个数标 ⚠️【出处失效，待重推】，并加一段核查说明与禁用范围
- `AI_Cowork/MAP.md` - §1 改写为 08-26 现状；§2 第 2 项打勾、第 3 项展开成四步
- `问题回答.md` - **新建**。学习收件箱（project-learning Phase 1），
  `sync-to-ob.sh` 的 SYNC_FILES 早已列它但文件一直不存在，同步必然漏一半

**验证结果**:
- [x] task.yaml YAML 解析通过：7 顶层键 / 4 phase / 17 step / 8 open_question / 无 tab
- [x] `0.32|0.89|2.1` 在 推导v1 中的出现次数 = 0（结论可复现）

**注意事项 / 踩坑记录**:
- **失效引用比错误数字更难发现**：错数字会被下游算出矛盾，失效引用只会被复制。
  这三个数活了 22 天、扩散到四个文件，靠的就是「每一处都注明了出处」这个外观
- **规矩**：`doc_ref_format` 只规定了引用的**写法**，没规定引用要被**验证**。
  建议加一条——新增「别名 §章节」引用时，当场 grep 一次目标节确认数字在场
- `AI_Cowork/scripts/sync-to-ob.sh` 的 PROJECT_DIR（`/f/project/...`）与 OB_VAULT（`/f/Obsidian/...`）
  是旧机器路径，与 doc_refs.external 里的 `C:/Users/20521/...` 是同一类漂移，用前须改

### [2026-08-14] 步 1.4 warmup_exercise 闭卷结果：A 3/5，B 2/8

**Agent**: Claude Code (Opus 5)　**方式**: 用户闭卷推导，AI 批改

**结果**: 相比 08-13 quiz 的 0.5/5 是质变——**上次是没动手，这次是动手但不复核**。
Amdahl 一般式 `1/(1-f+f/n)`、上界 `1/(1-f)` 完全写对；A4（逆用：端到端快 2 倍需计算快几倍）
方法完美；B3 写出了 max() 形式。**概念层面已通，卡在两个习惯。**

**逐题**:
- A1 ✅ 2748 / 65-14-21 全对
- A2 🟡 一般式正确但未代入（应得 1/0.35 = 2.86x）
- A3 ❌ 178.5、1141.5 都对，但 2748/1141.5 算成 1.32（真值 2.41）
- A4 ✅ 方法完美，1374-963=411 → 1785/411=4.34（写成 4.43，数字颠倒）
- A5 🟡 抓到「一边省计算一边加控制」这个关键，但未核 225 us 净减，未答出「8% 是诊断结果」
- B0 ❌ 答「传输时间为 0」，方向不对（**本组最重要的一题**）
- B1 ✅ 182.8　B2/B3/B4/B5/B6 🟡（结构对，全部停在公式未出数）　B7 ❌ 答「控制」，应为传输

**两个系统性问题**:
1. **算完不复核**：A3 的正确答案就在他自己 A2 写的公式里（1/(1-0.65+0.065)=2.41），
   手握验算工具没用；B 全程停在「加速比就是前后总时间之比」等于没写预测
   → 与 `predict_before_measure` 要求的「逼出一个能被打脸的数字」相违
2. **跨层数字直接相除**：B0 得 0。2.55x 是 csynth 周期比、123.0us 是上板墙钟，
   直接相除需两个条件（频率不变；csynth→上板修正系数 1.35 不变）。
   **这正是本项目历次翻车的同一根源**（表 5-4 两列不同源 / 8.229ms 不完整加总 / 三个 AI 跨三条边界）。
   讽刺的是用户前一天刚抓到 AI 出题时犯的这个错——**「认得出」尚未到「用得上」**

**AI 的账（第二次规格缺陷）**:
- 练习 A 用「计算/传输/控制」（Precompute 全区间 2748us），
  练习 B 用「H2D/kernel/D2H」（单次 call 182.8us），**两套不同分解，出题时没说明**。
  用户把 A 的「控制」项搬进 B，导致 B3/B4/B7 连锁答歪（B 中 T_fixed 藏在 40.6 的 H2D 里，不是独立项）
- 前一天同样因规格不清返工两次（「闭卷推导什么」「B2 优化多少」）
  → **教训：出题必须写全「给定 / 映射 / 假设 / 过关线」四件套**，step 1.1 的 exercise_order 是正确范式

**下一步**: MAP §2 第 1 项完成。进 D3（Roofline + TPU + CS149 Locality），角色翻转：用户画图，AI 批改。

### [2026-08-13·续3] 建 MAP.md 导航页；task.yaml 定位改为「数据库」

**Agent**: Claude Code (Fable 5)

**起因**: 用户反馈「之前做的都不清楚了，一头雾水；写在哪里了？应该搞个 map；
task.yaml 顺序不线性，可读性差」。诊断：6 份笔记 + 722 行 yaml，没有一个「你在这里」。

**修改文件**:
- `AI_Cowork/MAP.md` - **新建**，一页六段：§1 你在这里 / §2 线性执行队列（10 项）/
  §3 已完成时间线（9 行，产出可点开）/ §4 核心结论速查 5 条 / §5 文件地图（想知道 X → 打开 Y）/
  §6 task.yaml 的正确用法。维护规则：每次会话末尾只更新 §1/§2
- `AI_Cowork/task.yaml` - 头部注释改为「数据库不是阅读材料，执行顺序看 MAP」；
  **两处漏改的 status 修正**：task.status `not_started→in_progress`、
  phase1.status `pending→in_progress`（1.1 早已 done，状态一直没跟上——
  这也是「一头雾水」的一个来源：文件自己都不知道自己做到哪了）
- `AI_Cowork/README.md` - 顶部加 MAP 指针，定位收窄为「AI 冷启动背景页」；
  顺手修陈旧值 `PE_PARALLEL = 8 → 4`（对齐 define.h:33，标注 q7 的 pragma 硬编码问题）

**职责划分（此后生效）**:
  MAP.md = 人看的导航（在哪/下一步）；task.yaml = 数据库（搜 id 查询）；
  README = AI 冷启动背景；对表 = 冻结参考；log = 会话流水

**验证结果**:
- [x] task.yaml 727 行、无 tab、7 顶层键、35 id，结构核对通过
- [x] MAP §2 的执行序与 step 1.3 order_note 一致（1.3 → 1.5/1.6 → 1.2 → 1.4）

**踩坑记录**:
- **状态漂移无声发生**：task.status 停在 not_started 六周，phase1 停在 pending——
  每次只改 step 级 status，上层状态没人管。规矩：改 step status 时顺手核对上层
- 「一头雾水」不是记性问题，是**导航缺失**：产物都在（§3 能列出 9 行），
  但没有任何一页把「时间线 + 位置 + 下一步」放在一起。MAP 的维护成本
  （每会话更新两段）远低于每次重新解释的成本

**Agent**: Claude Code (Opus 5)

**任务**: 用户提出「task.yaml 要写好引用文档位置」。扫描全文引用状况后建表 + 立规 + 清理。

**扫描结果**（69 处引用）: md 文件 25 / 章节号 17 / 源码 file:line 14 / ADR 13

**关键发现——三份被反复引用的文档【不在仓库里】**:
- `HKS-ACC-Portfolio.md` **不存在**。它是 L0-L5 金字塔、R1 相图、R2 批深度的定义源，
  被 task.yaml `context` 与 推导v1 §5 依赖。疑在 OB vault，
  但 `scripts/sync-to-ob.sh` 的 `SYNC_FILES` 里**没有它**（只同步两个学习笔记）
- **HERA / CiFlow 的 PDF 不在仓库**。正文引用过 HERA §3.1/§3.4、CiFlow §IV.D/表 II 共 13+ 次，
  全部无法被下一个人（或下一次 AI 会话）核对
- ✅ 已核实存在：`AI_Cowork/decisions.md` 含 ADR-001~010；`docs/reports/bconv/` 含 10 份报告

**其余三类格式问题**（各已修）:
- 裸章节号：「与 §4 已查实…」——哪个文件的 §4？（实为 推导v1 §4）
- 别名不在任何表里：「笔记 §4.5」「gap 分析 §4」——靠上下文猜
- 正文写全路径：`docs/notes/L1_优化术语对表.md §3.1`——路径变了要改 N 处

**修改文件**:
- `AI_Cowork/task.yaml`
  - 新增顶层 **`doc_refs`**：`in_repo` 12 项（别名→路径）+ `external` 3 项（带 ⚠️ 与 todo）
  - 新增 `global_validation.doc_ref_format` 规则：一律「别名 §章节」，
    并规定**引用 external 文档前必须确认本地能打开——拿不到的文献不得作为推导依据**
    （traceability 的延伸：出处不可核对 = 拍脑袋常数）
  - `reference_docs` 改为指向 doc_refs，不再重复列举
  - 清理 15 处引用为别名格式
- `AI_Cowork/log.md` - 本条

**验证结果**:
- [x] 722 行、无 tab、7 个顶层键（新增 doc_refs）、35 个 id 顺序正确
- [x] 脚本扫描残留违规 = 0（跳过 doc_refs 表与规则条文自身）
- [x] **external 三项当场落实**（用户提供）：Portfolio 在
      `C:/Users/20521/Desktop/career-plan/`；HERA/CiFlow 在 Zotero。
      三条路径已写死并逐一验证存在
- [x] **doc_refs 全部 15 条路径脚本验证通过**（12 in_repo + 3 external，全 ✅）

**落实 external 时的额外发现**:
- **决定不入库**：Portfolio 是个人职业材料；论文 PDF 有 IEEE/ACM 版权，不应提交进 git。
  写死绝对路径 + 记录 Zotero key 失效时的重取方法（右键 → Show File）
- **Portfolio 与本仓库无同步机制**，`sync-to-ob.sh` 的 SYNC_FILES 里没有它 →
  两边可能各自漂移，引用其结论（L0-L5 金字塔 / R1 / R2）前须先确认未漂移。已写进 doc_refs 备注
- **顺带核查 D7 书单在 Zotero 的可得性**：BTS/ARK/SHARP ✅ 已有；
  **F1(MICRO'21)、CraterLake(ISCA'22) ❌ 缺失**——而 F1 是 D7 的 ★，需先获取。
  另发现 **Theodosian (2025)「memory-hierarchy-centric FHE acceleration」**，
  主题正是「系统框图 + memory hierarchy 找瓶颈」，比 F1 新四年，是否补进 D7 待用户定

**踩坑记录**:
- **「引用格式不统一」只是表症，真问题是引用指向了拿不到的东西**。
  扫描前以为是排版问题，扫完发现 13+ 处带章节号的文献引用无法核对——
  **这是 traceability 规则的一个盲区**：规则只管「数字有没有出处」，
  没管「出处本身能不能打开」。已把后者写进 doc_ref_format
- 同一类错误第四次出现：**组织形式决定内容**——
  reference_docs 是个平铺列表，没有「在不在仓库」这一维，
  于是仓库内外的文档并排列着，缺失了三份也没人发现

### [2026-08-13·续] §4 排空 + 对表冻结：待办归一到 task.yaml

**Agent**: Claude Code (Fable 5)

**任务**: 用户问「task.yaml 和对表到底按哪个走」→ 答：主线只有 task.yaml，对表是 step 1.3 的脚手架；
随后把对表 §4 剩余 4 项排空并冻结对表。

**修改文件**:
- `AI_Cowork/task.yaml` - 1.3 `changes` 增分层 Roofline 交付物；1.6 `validation` 增 N 维 tiling 问题
  （并连上猜想：合法切法是否恰是 HERA §3.1 的 NTT 分解——若是则 1.6 两问同一答案）；
  3.3 新增 `risk_load_imbalance`（前瞻记账：顺序执行时为 0，叠加异步双缓冲时兑现）
- `docs/notes/L0_符号表.md` - 新增 §9.1 AI 三条口径规定；**含一个新决策：分子统一 MAC 计 1**
  （理由：与脉动阵列 PE 计数直接对应；文献 ops 口径 ×2 换算）
- `docs/notes/L1_优化术语对表.md` - 头部加冻结声明；§4 标注已排空

**规矩（此后生效）**: 待办只住在 task.yaml；对表只在查术语（§1）、查口径（§3）、读材料（§5）时打开。
材料与 step 的对应：D3→1.3 roofline / D4-D5→1.3 符号表+1.6 / D7→1.3 框图 / D6→phase2/4（可推迟）。
**读材料的顺序就是做 step 1.3 的干活顺序，两者不是两条线。**

**验证结果**:
- [x] task.yaml 结构核对（缩进级，本机仍无 pyyaml）
- [x] 当日闭卷自测 5 题：**0.5 / 5**（Q1 半对，Q2–Q5 未答出）
  - 判定：并发度=1 的事实在（Q1），但 Amdahl 心算、超乘性机制、
    predict_before_measure 的反事实推理、同构错误扫查 全部没拿到
  - 两个真信号：Q1 反问「并发 HKS 还是算子？」恰是 q2 早已写下的 R2 区分；
    Q3 反驳「微架构优化可能减少读写量」对一般情况成立（B1 即例证）——
    直觉在，但没连上已有结论
  - **诊断**：本周模式漂移——step 1.1 是「用户推导、AI 检查」（成功，见 08-04/08-07 记录），
    本周变成「AI 推导、用户旁观」。文件越厚，手越空。
  - **处方**：① 明晨 30 分钟闭卷重推超乘性四行表 + Amdahl 三项记账，对 1.4 model_form 对账；
    ② 自 D3 起角色翻转：用户读论文、画 roofline、答验收题，AI 只批改不供成品

### [2026-08-13] D1–D2 收尾：查证 offload 阻塞性 + 材料清单补漏

**Agent**: Claude Code (Opus 5)

**任务**: 用户读完 CS149 第 1–2 讲，收尾 D1–D2 的三个验收问题；并核查材料清单为何比聊天里给的少

**过程**:
1. 用户提出「接下来应该分析每个微架构的 throughput 优化」→ 用 Little's Law 指出**并发度=1 时
   吞吐 = 1/延迟，两者是同一个量，「优化 throughput」当前不成立**
2. 用户指出「不是按笔记 §5 推进吗」——**确实跑偏了**，自作主张排了新的三步，
   其中第三步（II blocker / 关键路径）是 D6 内容，跳过了 D3/D4/D5。已归位
3. **查证 `FpgaManager` 的 offload 阻塞性**（D1–D2 验收问题 1 的后半段）
4. 用户核对材料清单，发现笔记 §5 比聊天里少——查实**只少了 CS149 的两讲**，其余 11 项都在

**修改文件**:
- `docs/notes/L1_优化术语对表.md` - §1.3 表更新 / §3.1 大幅增补源码结论与超乘性计算 /
  §5 补 CS149 两讲并加修订说明 / §4 待回写表增 4 行 / §6 速查改写
- `AI_Cowork/task.yaml` - 新增 `global_validation.latency_numbers_are_lower_bounds`；
  新增 `open_questions.q8`（已答）；新增 **步 3.4「Buffer 复用」**；
  步 3.0 加 `candidate_set`（3 项→4 项）；步 1.4 `model_form` 增补重叠版/超乘性/Amdahl 记账；
  步 1.2 新增**实验 D**（含先写后测的预测值）

**关键产出**:
- **结论 1：offload 完全阻塞。** 三条路径同一模式 `bo.sync ×2 → run.wait() → bo.sync`，
  无 async、无多 run 在飞、无双缓冲 → **并发度 = 1**。加号模型成立，
  但成立的原因是**实现如此，不是物理限制**（ADR-008 的推论）
- **结论 2（更重要）：182.8 μs 不是墙钟，是【下界】。**
  计时窗口开在 `Execute:415`，未计入 `xrt::bo ×3 分配`（:402-407）、`bo.write() ×2`、
  `bo_out.read()`、`bo ×3 析构` → `8.229 ms` 是**不完整的记账加总**
  → 已升级为 `global_validation` 规则：引用这些数字必须带「下界」二字
- **结论 3：新候选项「Buffer 复用」（步 3.4）。** `xrt::bo` 每次调用内部构造，
  45 次调用 = 135 次分配 + 135 次释放。**纯 Host 一个文件，是所有候选项里最便宜的**
- **结论 4：微架构优化与传输重叠是【超乘性】的。**
  1.69× / 1.49× 单独做，一起做 **3.06× > 乘积 2.52×**（`+` 变 `max()` 的后果）
  → 步 3.0 排序时**这两项必须捆绑评估**。且做完后瓶颈换人：PCIe 接管
- **Amdahl 收口**：计算清零上界 2.85×；修 OC 1.09×；OC+驻留 1.41×。
  **8% 本身就是诊断结果**——说明修 OC 没往主要矛盾上使劲

**踩坑记录**:
- **AI 判断方向搞反过一次**：对表 §3.1 初稿写「若有重叠则 8.229ms 是**上界**」，
  实际是**无重叠但漏计时 → 下界**。**教训：相反方向的偏差来源必须分开查**
  （重叠 → 高估；漏计时 → 低估），只想到一个就下结论，方向会反。已写进笔记留档
- **AI 排材料时把 CS149 从三讲砍成两讲**，砍掉的正是
  《Work Distribution and Scheduling》与《Locality, Communication and Contention》——
  **恰好覆盖对表圈出的三个真空白中的两个**。
  根因：按「天」组织日程，那两讲不属于 D1–D2 的主题又没有自己的一天，排期时被挤掉。
  **这与表 5-4「按列记账让两条不同源路径冒充同源」是同一类错误——组织形式反过来吃掉内容**
- **AI 一度脱离笔记自排计划**，被用户拉回。**笔记 §5 的顺序是有依据的**：
  D3(Roofline) 排在 D6(HLS) 前，因为**分层 roofline 决定 D6 该分析哪个模块**

**验证结果**:
- [x] 三条 offload 路径全部核对（`Execute:386` / `ModOpBatchOffload:566` / `BConvOffload:614`）
- [x] task.yaml 结构核对：640 行、无 tab、6 个顶层键、34 个 id 顺序正确（q8 已归位至 q7 后）
- [x] D1–D2 三个验收问题全部答完
- [ ] 实验 D 待上板（预测已按 predict_before_measure 规则写入 step 1.2）
- [ ] 下一步：D3（Roofline + TPU + CS149 Locality 讲）

### [2026-08-11] Step A：导师优化清单 ↔ 项目结论对表 + 材料学习计划

**Agent**: Claude Code (Opus 5)　**方式**: 先核查证据锚点，再贴标签

**任务**: 导师给出算子优化通用检查清单（并行性 / 数据依赖 / 延迟vs吞吐 / 计算vs带宽 /
tiling / reduction·同步·负载均衡 / Amdahl·Roofline，然后「从系统框图加 memory access 找瓶颈」）。
判断该怎么开展学习。

**过程**:
1. 读 task.yaml 全文 + `L1L2_推导v1.md` 全文 + `L0_符号表.md` + L0 闭卷推导章节结构，逐条比对
2. **核查结论：这不是一门要从头学的课，是一张检查清单，8 条里 5 条已答**——
   且 `L1L2_推导v1.md` §2.2 已经明确点名 Eyeriss stationary 分类，§2.1 已经算过 BConv 的
   AI 并做了 I/O-bound 交叉验证。比开会前的估计做得更多
3. 逐条标注 ✅有实有名 / 🟡有实无名 / ❌空白，每条挂证据锚点（文件+章节）
4. 圈出 3 个真空白，各写出「要答的问题」与「影响哪个 step」
5. 给 7 天材料计划，每份材料绑定一个项目内验收问题

**修改文件**:
- `docs/notes/L1_优化术语对表.md` - **新建**。主表 18 行（7 组）+ 三个空白展开 + 回写清单 + 材料计划
- `AI_Cowork/task.yaml` - step 1.3 加 `prereq_done_2026_08_11`；reference_docs 加对表

**关键产出**:
- **统计：✅×9 / 🟡×6 / ❌×3**。真工作量约 1.5 天，不是一学期
- **空白 ①「延迟 vs 吞吐」是最重要的**：L2 模型 `T = 计算+传输+控制` 用加号，
  等于假设三项不重叠——**该假设从未被声明过**。更糟的是 8.229ms = 1829+5537+863
  是**分项加总**，与 §4 已查实的 CPU 路径「总延迟是记账加总不是独立计时」是**同一个方法论问题
  在 FPGA 路径上重演**，但上次只点名了 CPU 路径
  → 下一步：查 `FpgaManager` 的 offload 是否阻塞。若非阻塞，现有全部延迟数字口径存疑
- **空白 ③ 发现 AI 口径混用**：`BConv 0.23 MAC/Byte`（on-chip AXI 边界）与
  `B0/B1/B4 = 0.32/0.89/2.1 ops/byte`（PCIe 边界）**不是同一条屋顶**。
  §2.1「BConv 卡 I/O」与 §四「AI=0.32」两句话都对，但笔记里读起来像同一个瓶颈
  → 必须画**分层 roofline**（PCIe / HBM / on-chip 三条），并在符号表 §9 立三条口径规定
  （分子 op 计法 / 分母跨哪条边界 / 是否含 evk）
- **空白 ② 负载均衡先厘清了一个陷阱**：`α_d + |C_d| = ℓ+K = 5` 对所有 d 恒等，
  所以**总量是均衡的**；不均衡只在**分算子**时出现（INTT 2:1，BConv/NTT 3:4）。
  → 于是它在 B0 边界下不构成问题，B1（每 digit 一次）下才构成 → 是步 3.3 的风险项
- **导师第 5 条对本项目的特殊含义**：瓶颈随边界而变。
  HERA 假设密文已在 HBM，图上没有 PCIe 屋顶 → 它看见的瓶颈是 NTT。
  把 PCIe 屋顶画出来才是本项目的 Fig. 1。**画图前先声明 scope，否则会得到一张
  和 HERA 一样的图，然后得出一个对本系统不成立的结论**

**踩坑记录**:
- **给已有结论贴标签这件事本身有价值**：「归约只沿 i」外部读者听不懂，
  「BConv 的 reduction 维是输入 limb，j/n 两维 embarrassingly parallel」一句就过。
  推导深度够但词汇不通用 = 外部读者一条也看不见
- **口径问题会跨路径复发**：CPU 路径的「分项加总冒充墙钟」在 2026-08-04 已查实并写进笔记，
  但 FPGA 路径的同一问题过了一周才被发现——**查出一个口径问题后应立即扫一遍所有同构位置**

**验证结果**:
- [x] 对表每一行都挂了可跳转的证据锚点（文件 + 章节号），全部核对过
- [x] task.yaml 结构核对（新键与同级键均在 col 8、无 tab、顶层键完整）。
      **注：本机无 pyyaml，未做真正的解析校验**——下次有环境时补一次 `yaml.safe_load`
- [ ] 材料计划（7 天）待执行，验收问题见对表 §5

### [2026-08-07] phase1 step 1.1 收口：题 2 + 题 3 + 符号表

**Agent**: Claude Code (Opus 5)　**方式**: 苏格拉底式（含两次纠偏）

**任务**: 完成 step 1.1 的题 2（ModUp 目标基）与题 3（ModDown），审稿用户 08-06 写的闭卷推导并修订

**过程**:
1. **审稿用户的 `L0_HKS_闭卷推导.md`**，查出一处真错误：§2.1 的 |v| 界用未约减版本，与 §2.3 的 α/2 差 2^60，中间缺了「先 mod a_i」那一步
2. **题 2**：改用「从 P4 倒推」的路径（原路径太绕），用户四步全对。结论：目标基由逐 tower 点积逼出，`|C_j| = l+K-alpha_j`
3. 用户中途纠正 AI 一处错误：**evk 存的是全 Q·P 基，不是 Q_l·P**（用户自己 §3.2 就写对了）
4. **题 3**：用户连问三次同一困惑（「取模后不就没了」「MD 里没有除法」「噪声呢」），根因是**一直在残数层面找答案**。用数字例子（59/7 mod 105）打通
5. **发现 AI 自己的符号撞车**：`s` 同时表示密钥和 `c̃ mod P`；`m` 同时表示"想要的结果"和"输入基素数个数"——是用户困惑的实质来源之一
6. 按用户批准的大纲改写 §4，新建符号表与最小可验证例子附录

**修改文件**:
- `docs/notes/L0_符号表.md` - **新建**，全项目统一符号口径（含 FPGA/建模侧，供 step 1.3/1.4 直接引用）
- `docs/notes/L0_HKS_闭卷推导.md` - §2 修订（约减那一步）、§3.1 更正（d_j ≠ c1）、§4 全面改写（4.0 动机 / 4.1 倒推四步 / 4.2 rho 的真实作用 / 4.4 噪声账 / 4.5 方法论）、新增附录 A 最小可验证例子
- `docs/notes/L0_HKS算法手推.md` - 踩坑记录从 4 条增至 6 条（新增"残数层 vs 整数层"、"符号撞车"）
- `AI_Cowork/task.yaml` - step 1.1 → **done**，derived_results 按题 1/2/3 分类重写

**关键产出（后续 step 直接消费）**:
- `|C_j| = l + K - alpha_j`（3/4）——只扩到 P 则仅 K=2，**多做 75% 是为了和 evk 对上 tower**
- **ModDown 算子账**：每 poly = K INTT-limb + 1 BConv + l NTT-limb = 6，两 poly 共 **12 次**
  → 直接命中 step 1.2「30 次未归属 kernel 调用」中的 12 次
- **噪声唯一源头是 evk 的 epsilon；BConv 不产生噪声，只放大 (1+alpha) 倍**（epsilon=0 时 BConv 误差贡献为 0）
  → alpha 因子的准确出处，闭合 `K = alpha (+1)`
- `rho` 让除法合法（不是消噪声）；`×P^{-1} mod q_i` 是否等价于除法只取决于被乘数是否 P 的整数倍

**踩坑记录**:
- **用户**：闭卷稿 §2 漏了「先 mod a_i」→ 误差界差 2^60；§3.1 把 digit 承载的整数写成 c1（应为 d_j = c1 mod Q_j）
- **AI**：符号撞车两处；题 2 推完后没看 validation 就继续加问题，约 80% 是重复劳动
- **方法论（已写入笔记）**：写代码在残数层，想问题在整数层；跨子问题推导前先建符号表；
  **验收标准写在任务文件里就是为了知道什么时候停**

### [2026-08-04] phase1 step 1.1 题 1：闭卷手推 BConv 与 CRT

**Agent**: Claude Code (Opus 5)　**方式**: 苏格拉底式，AI 出题 + 检查，用户推导

**任务**: task.yaml phase1/step1.1 的题 1——从 CRT 推出 BConv 及其近似误差界

**过程**:
1. 用户闭卷推导，AI 只给提示不给公式；中途 AI 过度苏格拉底化（把该讲的原理当成该问的问题），用户反馈后改为直接讲原理、保留推导给用户
2. 用户独立推出：X < αQ、e ≤ α−1、完整 BConv 公式、依赖表第 2/3 行
3. 用户独立发现关键观察：**误差界与 c_i 是不是逆元无关**——正确性与误差界是两条独立约束
4. 闭卷自测 11 题：6 题全对；4 题概念对但表达不精确；1 题（Q11）答反因果

**修改文件**:
- `docs/notes/L0_HKS算法手推.md` - 新建，11 节；标 ⭐（自己推出）/ ○（讲解后理解）便于复习时定位
- `AI_Cowork/task.yaml` - step 1.1 加 deliverable_path / derived_results，题 1 标完成
- 本条日志

**关键产出（后续 step 直接消费）**:
- **依赖表**：`X[i][n]` 不依赖 `j`／`W[i][j]` 不依赖 `n`／`μ_j` 只依赖 `j`
- 由依赖表推出：BConv 是 GEMM；`W` 复用 4096 次 vs `X` 复用 5 次 → **weight-stationary 是算出来的不是选出来的**；阵列 `α×K` = `LIMB_Q×MAX_OUT_COLS`
- **single-tower BConv 省输出不省输入**（因 `X` 不依赖 `j`）→ 支撑 step 1.3 真 OC 的输入重取结论
- **`K = α (+1)` 从观测升级为推导**（evk 噪声 + BConv 误差两条同量级约束）→ 支撑 step 1.6 参数外推
- 恒等式 `模数 × 补模数 = 总模数` 在单素数层与 digit 层各用一次：前者保证约减合法，后者歼灭 BConv 误差

**踩坑记录（已写入笔记 §0）**:
- 混淆 `a` 与求和式 `X`，误推 `e=0`——`mod` 展开式 `A = B − e·Q` 里 `e` 属于 `B` 不属于 `A`
- 误以为 CRT 重建是残数直接相加（等价于把基向量取成全 1）
- 算出 `e = −4/15` 非整数却未察觉——类型自检应当即时做
- 因果反向：一度认为"误差界与 P 有关"，实为"误差界决定 P 的下界"；且**两天前自己推出的最佳观察被自己遗忘**——这正是写笔记的理由

**方法论教训**:
- 苏格拉底式只在"对方有零件但没组装"时有效；缺基础概念时（CRT 重建）应直接讲，否则只是制造挫败
- 每个量先问"按定义它该是什么类型"，类型不符立刻回头，不要等推到最后

### [2026-08-04] twiddle 架构核查 + 目标转向全栈能力 + task.yaml 改组

**Agent**: Claude Code (Opus 5)

**任务**: 读 CiFlow 与 HERA 原文；确认 twiddle 是否压缩存储；按「全栈能力训练」（而非发论文）重排路线并写入 task.yaml

**执行步骤**:
1. 读 CiFlow(ISPASS'24) 与 HERA(FPGA'26) 全文（IEEE/ACM PDF 有 owner password，用 pdftotext 提取）
2. 核查 OpenFHE 的 sizeP 计算（rns-cryptoparameters.cpp:129-136）与 ComplPartQ 构造（:240-272）
3. 追查 PACK_RATIO / PACKED_TW_SIZE / local_twiddle，对照 csynth 报告的 URAM 用量
4. 追查 cg_ntt_reorder 与 NttForwardOffload 的位序处理
5. task.yaml 改组为四阶段路线，新增 step 1.5/1.6 与 open_questions q6/q7

**修改文件**:
- `AI_Cowork/task.yaml` - 大改：单一 OC 任务 → 四阶段全栈路线（13 step、7 open_questions），原 OC 步骤保留在 phase3
- `AI_Cowork/log.md` - 本条

**验证结果**:
- [x] task.yaml YAML 语法校验通过（413 行）
- [x] twiddle 问题已定性（见下）
- [ ] 位序问题待跑 test-cg-ntt 确认（列为 step 1.5）

**注意事项 / 踩坑记录**:
- **PACK_RATIO 不是压缩**：= 512/64 = 8 是 AXI 打包比（ADR-004），字节数不变。Portfolio 的「压缩 twiddle 320KB」口径错，实际是 192KB 且不是压缩
- **发现两套并存的 twiddle 架构**：Compute_CG_NTT 每 limb 流式（192KB）vs top.cpp 全 limb 常驻（3.0MB）。**csynth 基准 15,701cyc 来自前者，上板走后者——微基准与部署核不是同一设计**（前者 42% 周期花在 twiddle 加载，后者为零）
- **ADR-003 在工业参数下失效**：CG-NTT 存 STAGE×(N/2) 个 twiddle，标准 NTT 只需 N/2 → 冗余 log N 倍。N=2^16 时全驻留要 67MB（U55C 154%）。这正是 HERA §3.1 先做 NTT 分解的原因
- **cg_ntt_reorder 是孤儿函数**：只被 testbench 调用；ADR-005 描述的 host 侧 bit-reverse 在 FpgaManager.h:469 里也不存在。位序大概率被 BuildCGTwiddle 吸收，但需验证
- **PE_PARALLEL 四处不一致**：define.h=4 / cg_ntt.cpp:191 硬编码 8 / :233 注释 16 / csynth 报告对应 8。ADR-002「统一引用常量」已破坏，**docs/reports/ 全部数字属旧配置口径**
- **B1 复合 opcode 比原估计容易得多**：poly_buffer[8][SQRT][SQRT]=256KB 已够（B1 只需 160KB）；且全片上时 coefficient-wise 访问无转置代价（第一维就是 limb）。降级为「串三个已有子函数」
- **CiFlow §IV.D 原文**：算子总数与 dataflow 无关 → OC 的收益纯粹是 off-chip 流量/AI，不是计算量。gap 分析 §4 的「等价计算量」表概念上不成立
- **HERA 假设密文已在 HBM**（host 只做 "limb transfer from HBM"）→ host↔加速器这层边界仍无人研究

### [2026-08-04] 执行 L1/L2 推导 v1 + F2 口径对账

**Agent**: Claude Code (Opus 5)

**任务**: 按 2026-08-03 建立的指南实际做推导（上次只建了框架未做推导），并结合 Portfolio 的 L0–L5 金字塔补全 L1/L2 空白

**执行步骤**:
1. 读 AI_Cowork/ 全部 + PROJECT_STATUS + OC gap 分析 + 实验规划补 4-2/4-3/表 5-4 + `keyswitch-hybrid.cpp:331-590` + `define.h`
2. L1：把 BConv 还原为 GEMM (K=α_d, M=|C_d|, batch=N)，从代码数出 DC/MP/伪OC/真OC/真OC+驻留 五个点的调用数与字节数
3. 交叉验证：算出的 offload 调用数 15/15/14 与表 5-4 的 45/45/44 差恒为 30 → 映射关系坐实
4. L2：三模型代入实测系数，预测 vs 折算实测误差 <2%（但两边共享 182.8μs 常数，属自洽非独立验证）
5. F2 口径对账：确认 8.0–8.5ms 与 4.47–5.99ms 是 FPGA/CPU 两条路径，且前者 = n_call × 单次均值

**修改文件**:
- `docs/notes/L1L2_推导v1.md` - 新建，推导产物
- `docs/notes/L1L2_算子划分与性能模型指南.md` - 顶部加推导产物指针 + 三处错误修正说明
- `AI_Cowork/README.md` - 进度表 L1/L2 行改为已完成

**验证结果**:
- [x] 算子账与代码对拍（OC 的 ntt_limb=7 = 3 (p==0 Q补基) + 4 (2p×2part)，与表 5-4 一致）
- [x] 调用数模型三策略全部对上（差恒为 30）
- [ ] 独立的 Precompute 区间上板计时（F2 必做，未做）
- [ ] bytes 扫描实验标定真 BW / T_fixed（未做，两点拟合自由度为零）

**注意事项 / 踩坑记录**:
- **表 5-4 内部不自洽**：45 × 182.8μs = 8.23ms > 报的 4.47ms 总延迟；根因是延迟列是 CPU 路径、kernel 列是 FPGA 路径，两列不同源
- **答辩材料的 "OC 8.04ms 最优" = 44 × 182.7μs**，是算出来的控制模型输出而非实测，领先 2.3% 落在噪声内
- **峰值缓冲 64/128/32KB 是插桩虚构**：`MemoryTracker` 用 `sizeP`(=2) 当单位，实际 BConv 输出基是补基 |C_d|=3或4；伪OC 的 `peak_p_towers=1` 是硬编码常量而代码持有完整 `fullCompl`。伪OC 相对 DC 的峰值优势为 0 → **ADR-005/R2 的定量入口当前不成立**
- **gap 分析 §4 的真 OC BConv 次数 9 算错**，正确是 7（Section2 被数了两遍）
- **gap 分析 Gap#5 "跨 digit 累加 partial sum" 字面描述有误**：不同 digit 乘不同 evk，累加的必须是 `Σ_d ĉ_d[p]·evk_d[p]`，即需要 ModUp⊗KeyMult 融合
- **最重要的结论**：ADR-008（无状态 opcode RPC）导致每次调用重取输入，output-stationary 的复用收益被 PCIe 吃掉 → 真 OC 修好也慢（17 调用 > DC 15）。CP-005 单独做天花板 8%，加 "X 驻留复合 opcode" 是 29%

### [2026-08-03] 建立 L1 算子划分 + L2 三性能模型建模指南

**Agent**: Qoder

**任务**: 结合姜伟雄《AI Core 整理》方法论（三模型：计算/存储/控制；金字塔原理）与 HKS 项目代码库，新建算子划分推导与性能建模指南

**执行步骤**:
1. 阅读 AI_Cowork/ 全部文件、HKS-ACC-Portfolio.md、姜伟雄 PDF（8 页幻灯片，提炼三模型与金字塔原理 Pro）
2. 核对三策略代码锚点（keyswitch-hybrid.cpp DC:427-458 / MP:405-465 / OC:501-587）与表 5-4 实测数据
3. 新建指南：Part A（L1 循环嵌套形式化 + 设计空间三自由度）、Part B（L2 三模型总公式与单价表/带宽模型/控制模型）、Part C（验收清单）、Part D（与 F1-F4/R2 任务映射）、Part E（反模式）
4. 更新 README.md 进度表 + 本条日志

**修改文件**:
- `docs/notes/L1L2_算子划分与性能模型指南.md` - 新建，168 行
- `AI_Cowork/README.md` - 进度表新增 L1/L2 建模行

**验证结果**:
- [x] 指南内所有实测数字均来自实验规划.md 表 5-4（无拍脑袋常数）
- [x] 三策略行为描述已与 keyswitch-hybrid.cpp 当前代码对拍
- [ ] 指南中的推导本身待开展（本会话只建指南，未做推导）

**注意事项 / 踩坑记录**:
- define.h 当前 PE_PARALLEL=4，与 ADR-002/README 记录的 8 不一致——参数漂移无人察觉的实例，建模时需先对账（已在指南 B1 备注）
- 指南遵循"模型先于修复"原则：修 OC（CP-005）前应先完成模型 v0


### [2026-06-29] AI_Cowork 初始化

**Agent**: Claude Code (Opus 4.7)

**任务**: 在项目根目录建立 AI_Cowork/ 协作基础设施

**执行步骤**:
1. 检查 AI_Cowork/ 不存在 → 可创建
2. mkdir AI_Cowork/scripts/
3. 从 ai-cowork skill templates/ 复制 5 个模板，PROJECT_NAME 替换为 openfhe-for-HKS-ACC
4. README.md 用项目当前进度填充，避免空 placeholder

**修改文件**:
- `AI_Cowork/README.md` - 项目入口（链 docs/notes/PROJECT_STATUS.md）
- `AI_Cowork/task.yaml` - 通用空模板
- `AI_Cowork/decisions.md` - 空 ADR 模板（待逐步积累）
- `AI_Cowork/log.md` - 本会话首条记录
- `AI_Cowork/scripts/sync-to-ob.sh` - 待用户填 OB_VAULT 路径

**验证结果**:
- [x] 文件结构创建完成
- [ ] sync-to-ob.sh 待用户填实际路径

**注意事项 / 踩坑记录**:
- README.md 没用空 placeholder，直接复用 PROJECT_STATUS.md 的项目背景，让 AI 接手时无需再读 git 历史
- decisions.md 保持空模板（按 skill 文档意图，让用户主动积累，不预先填充）

### [2026-09-04] 片上存储优化 P0：基线冻结与同源码重跑

**Agent**: Qoder。用户先要求阅读 `doc/HKS_片上存储与数据搬运优化实施计划.md` 并提问，
经 4 项决策确认（锚点历史 139734、期限约两天只做 P0-P2、通用 BCONV 也迁 A_work、
8 塔跨度预留）后启动 P0。

**执行步骤**:
1. 取证发现工作树被 `git stash -u`（stash@{0}）打包重置到含 AUTO 的 480dc91，
   计划文档与 no_auto 报告目录一并被 stash，工作树与文档描述不符 → 询问用户后 pop 恢复
2. 从 stash 第三提交恢复 `docs/reports/hls/hks_no_auto_20260904/`（29 个证据文件）
   与原始计划文档；重建的文档与原版逐行比对，仅含已确认的决策性修改
3. 经用户授权提交基线检查点 `4b7eabf`（无 AUTO 源码 + 位选择改写 + 证据报告）
4. 以 `mem_p0_r1` 标签重跑：native 18/18+HKS22、csim、csynth、结构审计、smoke、
   perf co-sim（真实 fixture）、OpenFHE Release/ASan/XRT 语法检查

**修改文件**:
- `doc/HKS_片上存储与数据搬运优化实施计划.md` - 写入 4 项决策、基线提交哈希、BConv 单实例修正
- `docs/reports/hls/hks_mem_p0_20260904/` - 新增 P0 报告与全部证据
- Git 提交 `4b7eabf` - 基线检查点（20 M + 4 D + 28 A 报告文件）

**验证结果**:
- [x] 全部重跑结果与历史基线逐项一致：资源 688/1160/78323/147680/96、估算 5.250ns、
  结构审计 JSON 逐字段一致、OpenFHE 472/1523712、ASan 零泄漏
- [x] perf 周期零差异：INIT 295063、digit 70131/69603、暖态 139734、总 434797
- [x] 修正计划文档事实错误：BConv 为单实例 wrapper 128 BRAM_18K（local_in_x 48 +
  local_out_x 80），OP_BCONV 与 HKS 两个调用点共用，不是“两套各 128”

**注意事项 / 踩坑记录**:
- `git stash -u` 会把未跟踪的计划/证据文件一起打包移出工作树，恢复时若同名文件已存在
  会报 “already exists” 且 stash 保留；应从 stash 第三提交（stash@{0}^3）单独核对/恢复
- 计划文档此前由 Write 工具重建为 CRLF，与原版 LF 不一致；内容经 tr -d '\\r' 归一后
  diff 确认只差 6 处决策性修改
- 预乘 13824 周期预测与历史 README 数据自洽：12288 写零 + 1536 装载（LOAD_PAR=8）
- P1 阶段应避免重蹈：先跑 native，再 csim/csynth，再 cosim；BConv active_q 接口
  需同时更新 OP_BCONV 调用点与全部测试

### [2026-09-04] 片上存储优化 P1：消除无效塔清零与装载

**Agent**: Qoder。按计划 P1 实施 `active_q`：BConv 接口新增有效行数，
HKS 传 alpha、独立 OP_BCONV 传 LIMB_Q；Prepare 只循环 q<alpha，Load_X 只装载
有效行，Feed_X 对无效行显式注入 0（不读未装载的旧缓存）。

**执行步骤**:
1. 修改 bconv_systolic.h/cpp：签名加 active_q，Feed_X 注入 0，Load_X 按 active_q 装载
2. 修改 top.cpp：Prepare 循环界 q<alpha（去写零分支），两个调用点分别传 alpha/LIMB_Q
3. 测试台新增 TC-7/8 毒化用例：无效行填垃圾，验证不读旧值；golden 模型支持 active_q
4. mem_p1_r1 全链重跑：native、csim、csynth、结构审计、smoke、perf、OpenFHE Release/ASan

**修改文件**:
- `src/fpga_backend/include/bconv_systolic.h` / `src/bconv_systolic.cpp` - active_q 接口与注入 0
- `src/fpga_backend/src/top.cpp` - Prepare 只处理有效行，调用点传 alpha/LIMB_Q
- `src/fpga_backend/testbench/bconv_systolic_tb.cpp` - active_q golden 与毒化用例
- `docs/reports/hls/hks_mem_p1_20260904/` - P1 报告与证据

**验证结果**:
- [x] RTL perf：暖态 139734→125910，节省 13824 周期（9.89%），与模型偏差 0%；INIT 不变
- [x] 资源 BRAM 688/DSP 1160/URAM 96 不变，FF -1339、LUT +236，Fmax 190.48MHz 不变
- [x] native 18/18+HKS22、BConv 9/9（含毒化）、smoke 35/35、perf 40960 residues、
  OpenFHE 472/1523712、ASan 零泄漏、结构审计 PASS（20 MultMod、II=1）

**注意事项 / 踩坑记录**:
- Prepare 的 HLS 报告延迟仍是 12312：这是 tripcount max=3 的保守估计，实际按 alpha
  缩短；不能拿报告延迟当实测周期
- 后台任务共用一个 terminal_id 时，第二个任务可能不执行；长链任务用重定向日志 +
  ps 确认进程真实存在，不要只看命令回显
- 毒化用例若把无效行毒化值与有效值相同，会测不出"读了旧值"；毒化必须改变数值
