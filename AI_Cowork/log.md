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
