# openfhe-for-HKS-ACC 学习笔记

> 项目：基于 OpenFHE 的 CKKS Hybrid Key-Switching 软硬协同加速器（FPGA backend + Host 调度）。
> 本笔记是**理解驱动**的草稿，按知识点组织（不是按时间、不是复制代码注释）。
> 施工现场笔记——不做 `[[]]` 链接，链接留到搬进 Obsidian 时做。
>
> 配套文档：项目进度看 [docs/notes/PROJECT_STATUS.md](docs/notes/PROJECT_STATUS.md)；
> 架构决策看 [AI_Cowork/decisions.md](AI_Cowork/decisions.md)；OC 策略细节看 [docs/notes/OC_strategy_gap_analysis.md](docs/notes/OC_strategy_gap_analysis.md)。

---

## 目录

1. [CKKS 与 RNS：为什么密文是一张 N×L 的矩阵](#1)
2. [HKS 混合密钥切换：ModUp / ModDown 与数位分解](#2)
3. [NTT → CG-NTT：把"变几何"改成"恒定几何"](#3)
4. [BConv 脉动阵列 + Barrett 模约减](#4)
5. [HLS 设计陷阱合集（本项目踩过的所有坑）](#5)
6. [CPU-FPGA 异构与通信瓶颈](#6)
7. [数据流调度：DC / MP / OC 三策略](#7)

---

<a id="1"></a>
## 1. CKKS 与 RNS：为什么密文是一张 N×L 的矩阵

> **一句话理解**：CKKS 把"实数向量"塞进一个大整数多项式；RNS 再把"大整数"拆成一排 64 位小整数，于是一条密文多项式变成 **N 行（多项式系数）× L 列（RNS limb）的矩阵**，每一列可以独立并行算——这就是硬件能加速的根。

### 1.1 CKKS 在算什么

- CKKS 是**近似**同态加密：明文是实数/复数向量 `m`，编码成多项式环 `R_Q = Z_Q[X]/(X^N + 1)` 上的元素。
- `N` 是 2 的幂（本项目 `N = RING_DIM = 4096`，见 [define.h:25](src/fpga_backend/include/define.h)）。真实 FHE 里 `N` 常到 2^16 ~ 2^17，本项目取小参数便于综合。
- 密文是一对多项式 `(c0, c1)`。乘法、旋转会让密文"脱离"原密钥，必须用 **key-switching** 拉回来。

### 1.2 RNS：大整数拆成小整数

- 模数 `Q` 有几百到几千 bit，CPU/FPGA 的原生字长只有 64 bit，硬算不了。
- RNS（残数数系统）用 CRT 把 `Q = q_0 · q_1 · … · q_{L-1}` 拆成一堆互质的小模数 `q_i`（每个 36~64 bit）。
- 一个大整数 `x` 变成一组余数 `(x mod q_0, x mod q_1, …)`。**加减乘都能逐 limb 独立做**，limb 之间零依赖。
- 本项目：`LIMB_Q = 3`（主模数），`LIMB_P = 2`（辅助模数），`MAX_LIMBS = 8`，见 [define.h:37-49](src/fpga_backend/include/define.h)。

### 1.3 心智模型：N×L 矩阵

把密文的一个多项式想象成矩阵：
```
        limb0  limb1  limb2  ...  (L 列 = RNS moduli)
系数0 [  ·      ·      ·   ]
系数1 [  ·      ·      ·   ]
 ...
系数N-1[ ·      ·      ·   ]
 (N 行 = 多项式系数)
```
- **行方向（N）**：NTT/INTT 在这个方向做变换（系数域 ↔ 求值域）。
- **列方向（L）**：BConv 在这个方向做基转换（一组模数 → 另一组模数），每列独立并行。
- HKS 全程就是这张矩阵在**变形状**（列数在 ModUp 里从 α 扩到 α+K，在 ModDown 里缩回来）。

### 1.Y 易错点

- **`X^N + 1` 的 "+1" 不能丢** —— 这是负循环卷积（negacyclic），不是普通循环卷积。它决定了 NTT 要用 **2N 次单位根 ψ** 而非 N 次单位根 ω（多一层 ψ^i 的前后加权）。忘了这点会导致 NTT 结果整体错位。
- **RNS 不是"精度损失"** —— CRT 是精确可逆的（只要 `x < Q`）。RNS 的小模数和 CKKS 的"近似"是两回事：近似来自 CKKS 的编码/rescale，不来自 RNS。
- **limb 数会变** —— 密文的 L 不是常数。每次乘法后 rescale 会掉一个 limb（模数链变短）。代码里到处是 `sizeQl`（当前 Q limb 数）而不是固定值，就是这个原因。

### 1.Z 自测

1. 为什么 FHE 要用 RNS，而不是直接用大整数库（如 GMP）算？

<details><summary><b>答案</b></summary>
因为硬件（CPU/FPGA/ASIC）的原生算术单元是 64 位的。RNS 把几百上千 bit 的大整数模运算拆成一排 64 位小模数上的独立运算，每个都能用硬件原生乘法器/DSP 完成，且 limb 间无依赖 → 天然并行。大整数库要做多字长进位链，串行且慢，硬件不友好。
</details>

2. 本项目 `N=4096`、`LIMB_Q=3`，那么一个 Q 侧密文多项式在内存里占多少字节？（一个系数一个 limb 用 uint64）

<details><summary><b>答案</b></summary>
N × LIMB_Q × 8 = 4096 × 3 × 8 = 98,304 字节 ≈ 96 KB。单个 limb（一列）= 4096 × 8 = 32 KB —— 这个 32 KB 在后面 FPGA 传输/峰值缓冲的讨论里反复出现，记住它。
</details>

3. CKKS 的"近似"来自哪里？RNS 分解会引入近似吗？

<details><summary><b>答案</b></summary>
近似来自 CKKS 本身：编码时把实数量化进多项式、rescale 时丢弃低位。RNS 分解（CRT）是**精确可逆**的，只要被表示的整数 < Q 就无损。把两者混为一谈是常见误解。
</details>

---

<a id="2"></a>
## 2. HKS 混合密钥切换：ModUp / ModDown 与数位分解

> **一句话理解**：key-switching 是"给密文换锁"，占了 FHE 70% 的时间；HKS 把它拆成 **ModUp（把模数从 Q 撑大到 PQ）→ 乘 evk → ModDown（缩回 Q）**，核心算子链是 **INTT → BConv → NTT**，反复出现。

### 2.1 为什么需要 key-switching

- 密文乘法/旋转后，能解密它的密钥从 `s` 变成了 `s'`（比如 `s^2` 或旋转过的 `s`）。
- 要继续算，必须把密文"重加密"回原密钥 `s` 能解的形式——这就是 key-switching，靠乘一个特殊的 **evaluation key (evk)**。
- evk 很大（论文里 100MB~400MB），且每次 HKS 只用一次（无复用）。这是 CiFlow 论文关注"片外带宽"的原因。

### 2.2 HKS 的两个阶段

HKS = **ModUp**（P1-P5）+ **ModDown**（P1-P4）。参考 CiFlow (Neda et al., ISPASS 2024) §III：

| 阶段 | 步骤 | 干什么 |
|------|------|--------|
| ModUp | P1 INTT | 求值域 → 系数域（每 tower 一次，O(N log N)） |
| | P2 BConv | 基转换：把 digit 的模数从 α 撑到 α+K 列 |
| | P3 NTT | 系数域 → 求值域（撑大后的每列） |
| | P4 Apply Key | 逐点乘 evk（在大模数 PQ 下，噪声增长小） |
| | P5 Reduce | dnum 个 digit 的结果求和成一个矩阵 |
| ModDown | P1-P4 | 反过来：INTT → BConv → NTT → 缩回模数 Q |

**关键观察**：INTT → BConv → NTT 这条链在 ModUp 和 ModDown 里都出现，是整个 HKS 的计算热点。所以硬件加速 = 加速这条链的三个算子。

### 2.3 数位分解（digit decomposition）

- 把 `L+1` 个输入 tower 分成 `dnum` 个 **digit**，每个 digit 约 `alpha = ⌈(L+1)/dnum⌉` 个 tower。
- 见 [keyswitch-hybrid.cpp](src/pke/lib/keyswitch/keyswitch-hybrid.cpp) 里的 `numPartQl`（= dnum）、`alpha`。
- **权衡**：dnum 越大，每个 digit 越小 → evk 越小、噪声越低，但 BConv/NTT 调用次数越多。这是 HKS 相比朴素 key-switching 的核心可调旋钮。

### 2.Y 易错点

- **INTT 在前，NTT 在后，不是反过来** —— ModUp P1 先做 INTT（回系数域）才能做 BConv（基转换必须在系数域，因为它是 CRT 重建，逐系数独立），BConv 完再 NTT 回求值域。搞反了整个数据流就错了。
- **BConv ≠ NTT** —— 两个都在"变换"，但方向正交：NTT 沿**行（N 系数）**变换，BConv 沿**列（L limb）**变换。BConv 本质是矩阵-向量乘（见 §4），跟 FFT 一点关系没有。
- **digit 数 dnum ≠ limb 数 L** —— dnum 是"把 L 个 limb 分成几组"，是分组数。代码里 `numPartQl` 是 dnum，`sizeQl` 是 L，`alpha` 是每组大小，三者别混。
- **最后一个 digit 可能不满 alpha 个 tower** —— 当 `L+1` 不能被 dnum 整除时，最后一个 digit 的 `sizePartQl < alpha`。写 OC 的 bypass 逻辑时必须用 `partsCt[d].GetNumOfElements()` 而不是硬编码 alpha（这是 OC 修复的一个已知坑，见 task.yaml open_questions）。

### 2.Z 自测

1. HKS 里 INTT / BConv / NTT 三个算子的调用顺序是什么？为什么 BConv 必须夹在 INTT 和 NTT 中间？

<details><summary><b>答案</b></summary>
顺序：INTT → BConv → NTT。因为 BConv（基转换）是逐系数的 CRT 重建，必须在**系数域**做；而密文平时存在**求值域**（方便逐点乘法）。所以要先 INTT 把数据拉回系数域，做完 BConv 再 NTT 送回求值域。
</details>

2. dnum（数位数）增大，对 evk 大小、噪声、计算量分别有什么影响？

<details><summary><b>答案</b></summary>
dnum ↑ → 每个 digit 更小 → evk 更小、噪声增长更低（好）；但 BConv/NTT 要按 digit 逐个做，调用次数 ↑ → 计算量和数据搬运 ↑（坏）。这是 HKS 用"计算量换噪声/密钥大小"的核心权衡旋钮。
</details>

3. 为什么说 key-switching 是 FHE 的性能瓶颈？（给出量级）

<details><summary><b>答案</b></summary>
它占 HE 整体执行时间约 70%（CiFlow / ResNet-20 私有推理数据）。原因：每次密文乘法和旋转后都要调用它，bootstrapping 里更是海量调用；单次 HKS 涉及数百次 NTT、数百 MB 输入/输出、近 500MB evk、最多 1GB 中间数据，计算与访存双重密集。
</details>

---

<a id="3"></a>
## 3. NTT → CG-NTT：把"变几何"改成"恒定几何"

> **一句话理解**：标准 NTT 每一层蝶形的访存跨度（stride）都在变，硬件要一堆 MUX 去选地址；CG-NTT（恒定几何）靠"每层固定读 `[i, i+N/2]`、写 `[2i, 2i+1]` + 完美洗牌"把地址生成变成硬连线，代价是输出顺序被打乱成 bit-reversed。

### 3.1 NTT 是什么

- NTT = 有限域上的 FFT。把多项式乘法从 O(N²) 降到 O(N log N)。
- 核心是**蝶形运算（butterfly）**：`(a, b) → (a + ω·b, a − ω·b) mod q`，`ω` 是旋转因子（twiddle factor）。
- N=4096 → `STAGE = log2(4096) = 12` 层，每层 N/2 = 2048 个蝶形。见 [define.h:44](src/fpga_backend/include/define.h)。

### 3.2 变几何 vs 恒定几何

| | 标准变几何 NTT | CG-NTT（本项目主力） |
|---|---|---|
| 每层访存跨度 | 随 stage 指数变化（1,2,4,…,2048） | **恒定**：永远读 `[i, i+N/2]` |
| 写回 | 原位 | **完美洗牌**：写 `[2i, 2i+1]` |
| 地址生成 | 运行时算，一堆 MUX | 编译期定死，硬连线 |
| 硬件代价 | MUX 爆炸、布线拥挤 | MUX 极小 |
| 输出顺序 | 标准顺序 | **bit-reversed**（要额外重排） |

- CG-NTT 单 limb 延迟 15,701 cycles，比变几何的 112,671 快 **7.17×**，且 DSP -75% / FF -89% / LUT -57%（见 [docs/papers/实验规划.md](docs/papers/实验规划.md) 表 §1b）。
- 代码：[cg_ntt.h:44](src/fpga_backend/include/cg_ntt.h) 的 `CG_NTT_Kernel` 已改成 `template<bool IS_NTT>`，编译期消除 `is_ntt` 分支，让蝶形 II 从 4 路降到 2 路（这是 ntt_path.md 建议的模板特化优化，已落地）。

### 3.3 完美洗牌与 bit-reversed 输出

- 每层写 `[2i, 2i+1]` 等价于一次 perfect shuffle（完美洗牌）置换。
- 12 层洗牌累积下来，输出数据处于 **bit-reversed 顺序**（下标二进制位翻转）。
- OpenFHE 期望标准顺序，所以本项目在 **Host 侧**做补丁：`NttForwardOffload` 出口做 bit-reverse 还原，`NttInverseOffload` 入口做 bit-reverse 预打乱（见 [docs/notes/CG-NTT-Migration-Report.md](docs/notes/CG-NTT-Migration-Report.md) §5、AI_Cowork/decisions.md ADR-005）。
- 为什么放 Host 不放硬件：bit-reversal 是自反操作（做两次=恒等），Host 拦截只需 ~10μs，远小于 PCIe 延迟，且不用重综合 HLS。

### 3.Y 易错点

- **ψ（2N 次根）vs ω（N 次根）** —— negacyclic NTT（`X^N+1`）要用 2N 次本原根 ψ，蝶形前后还要乘 ψ^i / ψ^{-i} 加权。直接套普通 FFT 的 ω 会错。
- **CG-NTT 的输出不能直接用** —— 它是 bit-reversed 的。单元测试 `top_tb.cpp` 曾经"通过"是因为测试平台偷偷调了 `cg_ntt_reorder()`，但真实部署路径没调 → CKKS 端到端全错。**"测试通过 ≠ 部署正确"**，这是本项目一个真实血泪教训（见 Migration-Report §5.2）。
- **NTT 沿行、BConv 沿列** —— 再强调一次，别把两个变换方向搞混。
- **`CG_PE_NUM` 注释写 8，实际是 4** —— [cg_ntt.h:31](src/fpga_backend/include/cg_ntt.h) 写 `// 8`，但 `PE_PARALLEL` 在 [define.h:33](src/fpga_backend/include/define.h) 是 **4**。归档的旧 NTT 报告也写 8。**以 define.h 为准**，注释和旧文档是过期的。

### 3.Z 自测

1. CG-NTT 相比标准变几何 NTT，最本质的硬件优势是什么？代价是什么？

<details><summary><b>答案</b></summary>
优势：每层访存跨度恒定（永远读 [i, i+N/2]、写 [2i, 2i+1]），地址生成从"运行时算 + 一堆 MUX"变成"编译期定死 + 硬连线"，消除 MUX 爆炸和布线拥挤 → Fmax 更高、资源更省。代价：STAGE 次完美洗牌后输出是 bit-reversed 顺序，需要额外的重排步骤还原。
</details>

2. 为什么 bit-reverse 重排放在 Host 侧而不是 FPGA 内核里做？

<details><summary><b>答案</b></summary>
①bit-reversal 是自反操作，Host 端一次 memcpy 重排即可，成本 ~10μs 远小于 PCIe 传输延迟（>40μs），可忽略；②放硬件要重综合 HLS（小时级），影响其它优化迭代；③保留了"未来在 STORE_OUT 阶段零延迟内嵌 reorder"作为 L3 优化选项。是典型的"软件补丁换开发速度"决策。
</details>

3. `cg_ntt.h` 里 `CG_NTT_Kernel` 为什么要写成 `template<bool IS_NTT>` 而不是传一个 `bool is_ntt` 参数？

<details><summary><b>答案</b></summary>
模板参数是**编译期常量**，HLS 综合时会为 true/false 各生成一份独立硬件，`if(IS_NTT)` 分支被直接消除。而运行时 `bool` 参数会让蝶形循环里保留 4 路 MUX（NTT 和 INTT 两套读写模式并存），II 被拖到 4。模板特化后退化成纯 2 路 ping-pong，II 降到 2。
</details>

### 3.末 和当前项目的关系

| 子任务 | 用到的概念 |
|--------|-----------|
| `cg_ntt.cpp` 蝶形核 | 恒定几何、perfect shuffle、ping-pong buffer |
| `FpgaManager::NttForwardOffload` | bit-reverse 输出补丁 |
| L3 优化（待做） | BUTTERFLY_LOOP II=1、URAM 预加载 twiddle |

---

<a id="4"></a>
## 4. BConv 脉动阵列 + Barrett 模约减

> **一句话理解**：BConv 本质是"输出每一列 = 输入所有列的加权求和"，一个矩阵-向量乘；用 **PE 脉动阵列**并行做，用 **Barrett 约减**把取模的除法换成乘加移位，避免 HLS 综合出巨大的硬件除法器。

### 4.1 BConv 在算什么

- 基转换：把多项式从一组 RNS 模数 `{q_i}` 转到另一组 `{p_j}`。
- 数学：`out[j][n] = (Σ_i (x[i][n] · QHatInvModq[i] mod q_i) · QHatModp[i][j]) mod p_j`。
- 看 [dcrtpoly-impl.h:1233-1254](src/core/include/lattice/hal/default/dcrtpoly-impl.h) 的三层循环：外层 `ri`（系数）、中层 `i`（输入 limb）、内层 `j`（输出 limb）。**内层 j 之间完全独立**——这是能并行成脉动阵列的根，也是 OC 单-tower 优化可行的根。

### 4.2 脉动阵列（systolic array）

- 本项目 BConv 阵列维度 = `LIMB_Q 行 × MAX_OUT_COLS 列 = 3×5`，见 [define.h:40-42](src/fpga_backend/include/define.h)。
- **权值驻留（weight-stationary）**：BConv 权值常数（`QHatModp` 等）预加载进 PE 保持不变，输入系数流过阵列。
- 实测：Systolic 8,232 cycles vs Naive 53,613 cycles，**6.5× 加速**，代价是 DSP 从 48 涨到 870（15 个 MultMod PE，每个 58 DSP）。

### 4.3 Barrett 模约减

- 取模 `a mod p` 如果直接写 `%`，HLS 会综合出**硬件除法器**（128-bit `urem`，每个 ~6607 LUT + 8651 FF）。
- Barrett：预计算 `m ≈ 2^S / p`，把除法换成"乘 m + 移位 + 最多两次条件减"。全流水线 II=1。
- 见 [FpgaManager.h:659-666](src/core/include/FpgaManager.h)：`S = pbits + 62`，`m = 2^S / p`，Host 端算好传给 FPGA。

### 4.Y 易错点

- **`%` 在 HLS 里是灾难** —— 本项目最初 `bconv.cpp` 用 128-bit `%` 综合出 **30 个除法器**，吃掉 ~198K LUT（全片 15%），是时序违例主因之一。这是"软件思维写硬件"的经典反例：软件里 `%` 一个指令，硬件里是一整块除法器。
- **Barrett 的移位量 S 要够大** —— 早期版本用截断式双步乘法导致小模数溢出（见 git log `fix(bconv): 用全精度 Barrett 替换截断式双步乘法`）。全精度 `S = pbits + 62` 才安全。
- **权值驻留 ≠ 权值不用传** —— weight-stationary 是说权值在 PE 里驻留不动，但初始化时还是要通过 AXI 从 Host 加载一次。
- **脉动阵列的延迟对齐** —— `MULTMOD_LAT=4` 硬编码了 MultMod 的流水线深度做移位补偿。如果 HLS 改变了 MultMod 实际延迟，对齐会破坏。所以 `arithmetic.cpp` 里给 MultMod 加了 `LATENCY min=4 max=4` 锁死（见 decisions.md 相关，综合后修改方案 FIX-C）。

### 4.Z 自测

1. 为什么 BConv 能做成脉动阵列并行，而不是必须串行？（从数据依赖角度）

<details><summary><b>答案</b></summary>
因为输出的每一列 j 只依赖：所有输入列的 `x[i]·QHatInvModq[i]`（各 j 共享）、第 j 列的权值 `QHatModp[i][j]`、第 j 列的输出模数 `p_j`。**不同 j 之间零依赖**，所以 5 个输出列可以同时算。dcrtpoly-impl.h 内层 `for j` 循环就是这个并行维度。
</details>

2. 把 `(x*w) % p` 直接写进 HLS 会发生什么？正确做法是什么？

<details><summary><b>答案</b></summary>
HLS 会为每个 `%` 综合一个 128 位硬件除法器（urem），每个约 6607 LUT + 8651 FF，极其昂贵且时序路径很长。本项目 BConv 曾因此产生 30 个除法器吃掉 198K LUT。正确做法：用 Barrett 约减，预计算 `m ≈ 2^S/p`，把除法换成乘 + 移位 + 条件减，全流水 II=1。
</details>

3. 【闭卷】画出 BConv 3×5 脉动阵列一个 PE 在一拍里做的事（输入、权值、累加）。

<details><summary><b>答案</b></summary>
一个 PE 在一拍里：拿输入系数 `x[i][n]`（从上游流入）、驻留的权值 `QHatModp[i][j]`，做一次 Barrett 模乘 `x·w mod p_j`，累加到本列的局部部分和 `sum[j]`，把 `x` 往下游 PE 传。3 行 PE 纵向累加完成一列输出，5 列横向并行 → 每拍推进一个系数 n，II=1。
</details>

---

<a id="5"></a>
## 5. HLS 设计陷阱合集（本项目踩过的所有坑）

> **一句话理解**：HLS 让你用 C++ 写硬件，但"能编译"和"能综合出好硬件"是两回事——本项目集齐了 BRAM 超标、MUX 爆炸、bank 冲突、II 违例、除法器爆炸五大典型坑，每个都有明确根因和修法。

### 5.1 五大坑速查表

| 坑 | 现象 | 根因 | 修法 |
|----|------|------|------|
| BRAM 超标 206% | Place 直接失败 | twiddle `[MAX_LIMBS][BU_NUM][RING_DIM]`，BU_NUM=32 维**全冗余**（32 份一样的数据） | 消除 BU 维度 / CG-NTT 重构 |
| 时序违例 -0.33ns | Slack 负 | `switch(opcode)` 各 case 裸写 AXI for 循环 → 顶层 AXI master 巨型 MUX | 提取子函数 + `INLINE off` |
| Bank 冲突 | II=3 而非 1 | 8 个 PE 同拍读同一个 twiddle bank（双端口最多 2 读） | 建 PE_PARALLEL 个物理副本，各读各的 |
| 除法器爆炸 | LUT 54% | 128-bit `%` → 30 个 urem | Barrett 约减 |
| MUX 爆炸 | LUT 爆、综合卡死 | `ARRAY_PARTITION complete` 把 64 列全展开 → 64:1 MUX | 改 `cyclic factor=PE_PARALLEL` |

### 5.2 为什么这些坑反复出现

三个 HLS 反直觉的地方：
1. **数组维度可能是冗余的** —— `[8][32][4096]` 看着合理，但如果 32 这一维每份数据都一样（广播），就是 8× 浪费 BRAM。软件里数组不占"物理端口"，硬件里占。
2. **顶层函数越"大"越糟** —— 软件里一个大 `switch` 很正常；HLS 里顶层同时调度多个 case 的 AXI 访问会生成跨 case 巨型 MUX，扇出/布线爆炸。解法是把每个 case 的逻辑塞进 `INLINE off` 子函数，顶层退化成纯 dispatcher（decisions.md ADR-007）。
3. **`complete` 分区不是免费的** —— `ARRAY_PARTITION complete` 把数组彻底打散成寄存器/独立 BRAM，维度大了就是 N:1 MUX。要按并行度用 `cyclic factor=PE_PARALLEL` 而不是无脑 complete。

### 5.3 一个关键教训：CSynth 的 -0.33ns 未必是真违例

- HLS 的 C-Synthesis 只是**高层估计器**，面对十几个子模块共享 BRAM/AXI 时，它用最悲观的公式把所有 MUX 延迟相加，容易在顶层报假违例。
- Vivado 后端的 Implementation（P&R）有物理级 retiming / logic replication，-0.33ns 级的纯控制/路由延迟常被吸收。
- **所以别死磕 CSynth，该跑 Implementation 就跑**（见 archive/fix.md 结尾）。但这不是让你忽略所有违例——算法路径上的真违例（如 MultMod 关键路径）还是要修。

### 5.Y 易错点

- **"能 sw_emu 跑通"≠"能 hw 综合"** —— `OP_MUL` 拼错成未定义宏，sw_emu 容忍了，hw 综合直接编译报错。宽松环境会掩盖问题。
- **`#pragma HLS LATENCY min=2` 是下限不是上限** —— archive/fix.md 记录过一次把它理解反的乌龙。min 是"至少几拍"，压不了关键路径。
- **`cyclic factor` 要和 `UNROLL factor` 对齐** —— 分区 8 但展开 4，或反过来，会产生 crossbar 重排开销。本项目统一用 `PE_PARALLEL` 一个常量控制所有地方（decisions.md ADR-002），就是为了避免这种不一致。
- **`addr++` 优于 `base + i*stride + j`** —— AXI 地址如果用复杂乘加公式，HLS 要在一拍里实例化 DSP 乘法器 + 级联加法器 + 64bit 地址加法，路径极长。改成自增计数器 `addr++` 直接退化成 +1 累加寄存器（archive/fix.md 的"终极解决方案"）。

### 5.Z 自测

1. twiddle factor 数组 `[MAX_LIMBS][BU_NUM][RING_DIM]` 导致 BRAM 206% 超标，根因是什么？两种修法各是什么？

<details><summary><b>答案</b></summary>
根因：BU_NUM=32 这一维在 OP_INIT 里是**广播写入**——32 份数据完全相同，纯冗余，白占 32× BRAM。修法一（治标）：消除 BU 维度降到 `[MAX_LIMBS][RING_DIM]`，但会引发 bank 冲突，需再建 PE_PARALLEL 个副本放 URAM。修法二（治本）：直接换 CG-NTT，旋转因子布局变成 `[MAX_LIMBS][STAGE][CG_HALF_N]`，从根上不需要那么多片上存储。
</details>

2. 顶层 `Top` 报 Slack=-0.33ns，为什么"把 case 里的 for 循环提取成 `INLINE off` 子函数"能解决？

<details><summary><b>答案</b></summary>
因为原来各 case 直接裸写 AXI 读写，HLS 要在顶层为共享的 AXI master（gmem）生成跨所有 case 的巨型地址 MUX，扇出和布线延迟爆炸，压到顶层关键路径上。提取成子函数 + INLINE off 后，AXI 握手逻辑被封装进各子模块内部，顶层只剩 ap_start/ap_done 的纯函数调用分发，MUX 消失，路由延迟骤降。
</details>

3. 为什么本项目坚持所有 `ARRAY_PARTITION cyclic factor` 和 `UNROLL factor` 都引用同一个 `PE_PARALLEL` 常量？

<details><summary><b>答案</b></summary>
①一致性：分区因子和展开因子不匹配会在函数边界产生 crossbar 重排开销；统一后全链路 partition 一致，无 crossbar。②可维护性：调整并行度只改 define.h 一处，不用在 10+ 个文件里逐个改。这是 decisions.md ADR-002 的核心。
</details>

---

<a id="6"></a>
## 6. CPU-FPGA 异构与通信瓶颈

> **一句话理解**：粗粒度 FHE 算子（NTT/BConv）计算量远大于 CPU 调度开销，所以用 **CPU 逐个发 opcode（RPC 模型）** 就够了；但真正的瓶颈变成了 **PCIe 小粒度传输**——每次只发 1 个 limb，DMA 固定开销摊不开。

### 6.1 任务划分：谁在 CPU、谁在 FPGA

- **FPGA**：计算密集且规则并行的算子——NTT/INTT、BConv、逐元素模乘。
- **CPU（Host）**：控制密集/数据依赖复杂的——密文管理、数位分解、rescale、密文加法、调度。
- 桥梁：`FpgaManager` 单例 + opcode。8 个 opcode 见 [define.h:57-64](src/fpga_backend/include/define.h)。

### 6.2 opcode-RPC 调度模型

- 判据：**开销比 = CPU 调度开销 / 单次计算周期**。
  - CPU 发一次 op 约 2000~3000 cycles 纯开销（AXI-lite 写寄存器 + 中断 + 上下文切换）。
  - NTT ~300K cycles、BConv ~500K cycles → 开销比 < 1%，CPU 调度无感。
  - 单次模乘 ~5 cycles → 开销比 60000%，必须合并/FSM。
- 结论：本项目粗粒度算子用 opcode-RPC 就是最优，不需要 FPGA 内部自主 FSM（decisions.md ADR-008，详见 [docs/papers/fsm_vs_cpu_scheduling.md](docs/papers/fsm_vs_cpu_scheduling.md)）。

### 6.3 真正的瓶颈：PCIe 小粒度传输

- 上板实测：H2D 均值 40.6μs/call，单 limb 32KB → 折算带宽仅 **787 MB/s**，是 PCIe Gen3 x16 理论值（~12GB/s）的 6.5%。
- 根因**不是** PCIe 带宽不足，而是 **DMA 固定建立开销（~20-30μs/call）在 32KB 小传输下占了主导**。
- 实测加速比（INTT 3.08× / BConv 2.88× / NTT 4.83×，答辩 PPT Slide 16）远低于 csynth 的 7.17× / 6.5×，就是因为 Host 每次只 emit 1 个 limb，把 csynth 报的 5-limb 批量延迟拆成了 5 次小调用，固定开销 ×5。
- 优化方向：**Batch HKS** —— 多个 HKS 共享一次 kernel launch / DMA / buffer sync，提高单次 offload 的计算密度（答辩 PPT Slide 17）。

### 6.Y 易错点

- **"传输慢"≠"带宽不够"** —— 这是本项目排查上板性能时最关键的认知。带宽（GB/s）和固定延迟（μs/call）是两个独立指标。小数据量下是固定延迟主导，加带宽没用，要减少调用次数（批处理）。
- **csynth 延迟 ≠ 上板延迟** —— csynth 只反映 kernel 理想执行延迟，不含 runtime 调度、buffer 同步、DMA 建立。上板一定更慢，差值就是通信/调度开销。
- **串行三段式掩盖了带宽真值** —— 当前 Load/Compute/Store 串行无 overlap，"PCIe 带宽"其实是 (Transfer+Kernel+Transfer)/Transfer，被 kernel 时间稀释。不做 double-buffering 就没法单独证伪带宽。
- **mem_in1 和 mem_out 别共享 gmem bank** —— 共享同一 HBM bank 会读写争用带宽减半。要用 `.cfg` 的 `sp=` 指令分配到不同 bank。

### 6.Z 自测

1. 上板 H2D 折算带宽只有 787 MB/s，远低于 PCIe 理论 12GB/s。这说明什么？该怎么优化？

<details><summary><b>答案</b></summary>
说明瓶颈不是带宽而是 **DMA 固定建立开销**（~20-30μs/call）。单 limb 只有 32KB，固定开销占了绝大部分时间，折算带宽自然低。优化不是"换更快的 PCIe"，而是**增大单次传输粒度 / 减少调用次数**——把多个 limb 或多个 HKS 攒成一次传输（Batch HKS），让固定开销摊薄。
</details>

2. 什么样的算子适合 CPU 逐个发 opcode 调度，什么样的必须放进 FPGA 内部 FSM？判据是什么？

<details><summary><b>答案</b></summary>
判据是开销比 = CPU 调度开销(~2-3K cycles) / 单次计算周期。开销比 <1%（如 NTT 300K、BConv 500K cycles）→ CPU 调度无感，用 opcode-RPC 最灵活。开销比 >10%（如单次模乘 5 cycles）→ 调度开销压倒计算，必须合并成批量 op 或在 FPGA 内部用 FSM 循环，把循环边界从"CPU-FPGA 通信边界"拉到"FPGA 内部"。
</details>

3. 为什么上板实测的 NTT 加速比（4.83×）比 csynth 报的（7.17×）低？

<details><summary><b>答案</b></summary>
csynth 的 7.17× 是 kernel 内部理想执行延迟的对比，不含通信。上板的 4.83× 包含了 runtime 调度、buffer 同步、DMA 建立开销，而且 Host 每次只发 1 个 limb（csynth 报的是 5-limb 批量），固定开销被乘以 5，稀释了 kernel 的加速。两者口径不同，差值就是系统级开销。
</details>

---

<a id="7"></a>
## 7. 数据流调度：DC / MP / OC 三策略

> **一句话理解**：同一套 FPGA 算子，Host 端**发射顺序**不同就得到三种调度——MP（全宽并行，最快但最吃缓冲）、DC（逐 digit 深度优先，中庸）、OC（逐输出 tower，缓冲最小）。三者是"存储 ↔ 延迟 ↔ 带宽"的权衡。

### 7.1 三策略速览（源自 CiFlow, ISPASS 2024）

| 策略 | 调度逻辑 | 峰值缓冲 | 特点 |
|------|---------|---------|------|
| **MP** (Max-Parallel) | 广度优先：所有 digit 同时推进（INTT全部 → BConv全部 → NTT全部） | 最大（128KB） | 并行度高，但中间数据全部同时存在 |
| **DC** (Digit-Centric) | 深度优先：一个 digit 走完 INTT→BConv→NTT 再下一个 | 中（64KB） | OpenFHE 原生模式，实现最简单 |
| **OC** (Output-Centric) | 按输出 tower：一次只算一个 output tower | 最小（32KB） | 中间数据最少，CiFlow 主推 |

- 代码：三分支在 [keyswitch-hybrid.cpp:408-587](src/pke/lib/keyswitch/keyswitch-hybrid.cpp)，MP phase1 在 :408，DC/MP BConv 在 :427，OC 在 :501-587。
- 切换靠 `SetHKSStrategy()`，走同一 FPGA bitstream（不重综合）——这正是 opcode-RPC 模型的好处（decisions.md ADR-001）。

### 7.2 实测结果与反常现象

实测（N=4096, sizeQl=3, sizeP=2, 2 digits）：

| | DC | MP | OC |
|---|---|---|---|
| Precompute 延迟 | **4.47ms** | 5.17ms | **5.99ms** |
| 峰值缓冲 | 64KB | 128KB | **32KB** |

- **反常**：OC 峰值缓冲最小（对），但**延迟最高**（5.99ms > DC 4.47ms），与 CiFlow 论文"OC 延迟最优"矛盾。

### 7.3 OC 反常的根因（本项目核心待修问题）

- 当前 OC 实现是"**DC + sizeP 倍冗余 BConv**"：[keyswitch-hybrid.cpp:534-545](src/pke/lib/keyswitch/keyswitch-hybrid.cpp) 调 `ApproxSwitchCRTBasis()` 算**全 sizeP** 个 tower，然后 `GetElementAtIndex(complIdx)` 只挑 1 个用，其余丢弃 → 计算量白白乘了 sizeP 倍。
- CiFlow 真正的 OC 有 5 个机制，当前只做对了 1 个（INTT 输出复用）。详见 [docs/notes/OC_strategy_gap_analysis.md](docs/notes/OC_strategy_gap_analysis.md)。
- 修复顺序（decisions.md ADR-010）：Section1 bypass（纯 Host）→ single-tower BConv（改 OpenFHE）→ on-chip 累加。

### 7.Y 易错点

- **单 kernel 下 DC ≈ MP** —— 因为只有一个 Top kernel 实例，所有操作串行执行，MP 的"并行"发挥不出来。MP 的优势要**多 kernel 实例**才体现。别以为 MP 天生更快。
- **OC "输出粒度更细"只是表象** —— CiFlow OC 的关键不只是"一次算一个 tower"，还有 **Section1 bypass**（mod Q 输出时对应 digit 的 tower 在原始密文里已经存在，直接复用，零计算）。只理解成"粒度细"会做出当前这个"更慢的 OC"。
- **峰值缓冲小 ≠ 延迟低** —— OC 用调用次数换缓冲，两个指标可以背离。CiFlow 在大 SRAM + 向量处理器上 OC 双优，但本项目小参数 + PCIe-bound 场景下，缓冲优势保留、延迟优势要靠修复才能兑现。
- **论文/PPT 论述可能领先于实测** —— 答辩 PPT Slide 19 写"OC 延迟略优"是**预期**不是当前实测。修复前这是材料里的隐性 bug，要注意别在答辩里说错。

### 7.Z 自测

1. 为什么当前实现里 OC 的延迟反而比 DC 高？（一句话根因）

<details><summary><b>答案</b></summary>
因为当前 OC 每次算一个 output tower 时，仍调用完整的 `ApproxSwitchCRTBasis()` 算出全 sizeP 个 tower 然后只留 1 个丢弃其余，等价于把 BConv 计算量乘了 sizeP 倍。它是"DC + sizeP 倍冗余 BConv"，不是真 OC。
</details>

2. 单 kernel 架构下 DC 和 MP 实测延迟几乎一样，为什么？MP 的优势什么时候才体现？

<details><summary><b>答案</b></summary>
因为只有一个 Top kernel 实例，所有操作（不管哪个策略）都串行执行，MP 的"广度优先并行"在单 kernel 上退化成和 DC 一样的串行序列。MP 的优势要**实例化多个 kernel**（比如 2 个 NTT kernel）才能体现——那时 MP 阶段内的多次调用可并行分发到不同 kernel。
</details>

3. CiFlow OC 的 "Section1 bypass" 是什么？为什么它能省计算？

<details><summary><b>答案</b></summary>
在 ModUp 计算 mod Q 的输出 tower 时，第 `p/alpha` 个 digit 对应的那个 tower **在原始密文的求值域系数里已经存在**（它没被 BConv 扩展过）。所以直接 `partsCtExt[bypass_digit][p] = 原始 EVAL 系数` 即可，这个 digit 完全不需要过 BConv 和 NTT，省下 1/dnum 的计算量。当前实现漏了这个机制，对所有 digit 一视同仁地过完整 BConv。
</details>

### 7.末 和当前项目的关系

| 子任务 | 用到的概念 |
|--------|-----------|
| `hks_strategy.h` 枚举 + `SetHKSStrategy` | 三策略切换、同 bitstream |
| `MemoryTracker` | 峰值缓冲 CSV 追踪、DC/MP/OC 对比 |
| OC 修复（待做） | Section1 bypass、single-tower BConv、on-chip 累加 |

---

## 附：知识点依赖关系（速查）

```
CKKS/RNS (§1) ──┬──> HKS (§2) ──┬──> NTT/CG-NTT (§3) ──> HLS 陷阱 (§5)
                │               └──> BConv/Barrett (§4) ──┘
                │
                └──> 数据流调度 DC/MP/OC (§7)
                          ↑
        CPU-FPGA 异构/通信 (§6) ──┘
```

- 想理解**硬件为什么这么设计** → §1（RNS 并行性）+ §3/§4（算子）+ §5（HLS 约束）
- 想理解**系统为什么没达到理论加速** → §6（PCIe 小粒度）+ §7（OC 调度反常）
- 想**改代码** → 先读 [AI_Cowork/task.yaml](AI_Cowork/task.yaml)（OC 修复任务）+ [AI_Cowork/decisions.md](AI_Cowork/decisions.md)（历史决策）
