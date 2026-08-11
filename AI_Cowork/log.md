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
