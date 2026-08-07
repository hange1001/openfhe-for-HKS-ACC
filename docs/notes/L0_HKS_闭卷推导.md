# L0 闭卷推导：HKS（Hybrid Key Switching）

> step 1.1 交付物。方法：先闭卷从 CRT 第一性原理推出 ModUp / ModDown，再与代码锚点对账。
> 日期：2026-08-06。项目参数口径：l = sizeQl = 3，α = 2，K = sizeP = 2，dnum = 2，N = 4096。
> 数据形状链图系自画，不抄 CiFlow Fig.1。

符号约定：
- Q_l = ∏_{i∈active} q_i，active tower 数 l；P = ∏_{j=1..K} p_j；所有模数两两互素。
- 每个 q_i 约 2^60 bit；「tower」= 一个 Z_{q_i}[X]/(X^N+1) 上的多项式，N 个系数。
- 系数均取中心剩余 (−q/2, q/2]，记 |·| 为中心剩余绝对值。

---

## 1. 问题设定：key switching 要算什么

CKKS 密文 c = (c0, c1)，c0 + c1·s ≡ m + e (mod Q_l)。
换密钥 s → s′，需要构造 c′ 使 c′0 + c′1·s′ ≡ c0 + c1·s (mod Q_l)，噪声增幅可控。

思路：把 c1 拆成小份（digit），每份乘一份 evk。evk_j 加密的是 P·s′
（P 是辅助模数积，P 加权嵌在 evk 里而不是乘在 digit 上——这是 OpenFHE/eprint 2021/204
附录变体；Han-Ki 原文是带权分解，两者等价，见指南 §Step 1 修正记录）。于是乘出来的是
P·c1·s′，最后 ModDown 除以 P 还原。

digit 分解（tower 不相交切分）：active tower 集合 {0..l−1} 按 α 个一组切成 dnum 组：

    Q_j = ∏_{i ∈ digit j} q_i ,  c_{1,j} = c1 在 digit j 的 tower 上的投影

恒等式（整系数多项式层面）：

    c1 ≡ Σ_j CRT_j(c_{1,j}) · (Q_l / Q_j) · (Q_l / Q_j)^{-1}  (mod Q_l)      …… 形式上

但 HKS 不显式重构大整数，而是把「重构 + 换基」合并为 ModUp。噪声预算：digit 越小
（dnum 越大、α 越小）噪声越好，代价是 BConv 次数和 evk 份数增加——这是 L1 层的设计空间，
L0 只负责把每一步的代数恒等式和误差界钉死。

---

## 2. 从 CRT 推出 BConv（题 1 核心）

### 2.1 精确 CRT 重构

设基 A = {a_1,…,a_m} 两两互素，Â = ∏ a_i，Â_i = Â/a_i。
x 在基 A 上的中心剩余为 x_i ∈ (−a_i/2, a_i/2]。CRT 重构：

    v = Σ_i x_i · Â_i · (Â_i^{-1} mod a_i)  ≡ x (mod Â)

即 v = x + u·Â，u ∈ Z。u 的界由求和界给出：

    |v| ≤ Σ_i (a_i/2) · Â_i · |Â_i^{-1} mod a_i| < (1/2)·Â·Σ_i a_i

### 2.2 近似基转换（BConv）

要把 x 表示到新基 B = {b_1,…,b_k}（与 A 互素）。精确做法：先算大整数 v，再取 mod b_j。
BConv 的近似是：把 2.1 的求和式直接「延迟取模」——中间不做大整数归约，只在最后对每个
b_j 取模：

    BConv(x)_j = [ Σ_i x_i · (Â_i^{-1} mod a_i) · (Â_i mod b_j) ] mod b_j

**BConv 求和式**（这就是要徒手写出的式子；每个系数位置独立成立，以下省略系数下标）：

    ŷ_j = [ Σ_{i=1..m} x_i · QHatInv_i · QHatModp_{i,j} ] mod b_j
    其中 QHatInv_i = (Â_i^{-1 mod a_i)，QHatModp_{i,j} = Â_i mod b_j

### 2.3 「近似」误差界

因为 Σ 内部对 Â 取模的值与真重构 v 同余（mod Â），BConv 输出与精确值之差是 Â 的倍数：

    ŷ_j ≡ x + u·Â (mod b_j)，  |u| ≤ (1/2)·Σ_i a_i ≤ m/2 · max a_i

对 ModUp 的应用（A = Q_j，m = α 个 tower，q_i ≈ 2^60）：

    **每系数误差 = u·Q_j，|u| ≤ α/2 ≈ O(α)**
    即「BConv 溢出误差 ~ α·Q_j」——task.yaml step 1.1 提示里 a·Q_j 的 a 就是 α（digit 宽度）。

误差的量级由【输入基的积】决定，与输出基无关：这是后面「P ≥ max digit」的种子。

### 2.4 求和式里哪些量跨 output tower 共享、哪些按列独立（validation 1）

把 BConv 看成 GEMM：M = k（输出 tower），K' = m（输入 tower），batch = N（系数）：

    ŷ[n][j] = Σ_i x[n][i] · QHatInv_i · QHatModp[i][j] (mod p_j)

- **跨 output tower 共享**（只算一次，j 循环外）：
  x_i = c[i]·QHatInv_i mod q_i —— 每个输入 tower 每系数一个标量乘法；
  以及 Barrett 约约前的 128-bit 部分和结构。
- **按列独立**（每个输出 tower 一套）：
  QHatModp[i][j]（权重）、p_j（输出模数）、Barrett mu_j = ⌊2^128/p_j⌋。

→ 计算复杂度账：共享量把每系数成本从 O(m·k) 次乘降到 m 次共享乘 + m·k 次乘累加；
这正是 FPGA 脉动阵列 3×5 的 K 维复用结构，也是 OC 策略「算全列再丢弃」浪费的定量来源。

**对账**：dcrtpoly-impl.h:1233-1254（CPU 主路径）——外层 ri（系数）、中层 i 算
xQHatInvModqi（共享量）、内层 j 做 Mul128 累加（列独立权重），收尾用 BarrettUint128ModUint64
按列的 mu[j] 归约。三层循环次序与上面完全一致。

---

## 3. ModUp：digit j 的目标基为什么是 (Q_l \ Q_j) ∪ P（validation 2）

### 3.1 目标基的推导

P4（apply key）的运算是：

    c̃_{0,1}[t] = Σ_j ĉ_j[t] · evk_j[t]  (逐 tower t ∈ Q_l ∪ P 基上的点积)

两个事实决定目标基：

1. **evk 存在 Q_l ∪ P 全基上**（keygen 时按全 Q 加密 P·s′）。tower-wise 点积要求
   ĉ_j 在【与 evk 相同的全基】上有值 → digit j 必须被扩展到 (Q_l \ Q_j) ∪ P，
   加上自己原有的 Q_j，拼起来才覆盖全部 l+K 个 tower。
   **只扩到 P 是不够的**：Q_l 里 Q_j 之外的 tower（其余 digit 的地盘）上 ĉ_j 没有值，
   P4 的点积在这些 tower 上缺项；而这些位置的值数学上不是 0——c1 作为整数在这些模数下
   的剩余恰由 BConv 给出（带 u·Q_j 误差）。
2. **为什么用近似值而不是精确 CRT lift**：精确 lift 需要对 Q_l 全基做带误差修正的精确
   基转换（SwitchCRTBasis 多一步 αQModp 修正项），而近似误差 u·Q_j 在 ModDown 除以 P 后
   变成 u·Q_j/P ≤ u（当 P ≥ Q_j），是可接受的小噪声。用 α·K 个标量换掉一次精确 lift，
   这就是 HKS 相对 GHS 的核心交换。

**形状**：digit j 的输出基 C_j = (Q_l \ Q_j) ∪ P，|C_j| = l − a_j + K
（a_j = digit j 实际 tower 数；最后一个 digit 残缺时 a_j < α，见 q5）。
注意输出基的 tower 顺序按 ComplPartQ 的构造序（先其余 Q tower 后 P tower），
不是 QlP 序——拼成 QlP 序是 assemble 步的事。

**对账**：rns-cryptoparameters.cpp:240-272——sizeComplPartQj = (l+1) − sizePartQj + sizeP；
前 (l+1)−sizePartQj 个模数取「currDigit >= j 时自增跳过自己」的其余 digit 的 tower，
后 sizeP 个取 moduliP。与推导的 C_j 逐元素一致。

### 3.2 ModUp 各步数据形状（自画形状链，P1-P5 + 拼装）

记 N×t = t 个 tower、每 tower N 个系数。当前参数 l=3, α=2, K=2, dnum=2, a_0=2, a_1=1。

```
输入 c1（EVAL）:                     2 × (N×l) 中取 c1        →  N×3
│
│ P1a digit 分解（tower 切片，零填充，无加权）
│     keyswitch-hybrid.cpp:353-380
▼
dnum × (N×α) = 2 × (N×2)             （最后 digit 残缺 → 2×(N×2) 与 1×(N×1)，共 3 tower）
│
│ P1 INTT（EVAL → COEFF，BConv 只在系数表示下是逐系数运算）
▼
2 × (N×α) COEFF
│
│ P2 BConv：基 Q_j → C_j = (Q_l\Q_j) ∪ P，每系数加 u·Q_j 误差
│     ApproxSwitchCRTBasis；输出基大小 l+K−a_j
▼
dnum × (N×(l+K−α)) = 2 × (N×(5−a_j)) → 当前参数：2 × (N×3) + 1 × (N×4)
│
│ P3 NTT（COEFF → EVAL）
▼
同上形状，EVAL
│
│ P4(组装) assemble：Q_j 段放回原值，其余位置填 BConv 输出
│     keyswitch-hybrid.cpp:471-499 / OC 分支 :507-526
▼
dnum × (N×(l+K)) = 2 × (N×5)          ← 每个 digit 都在全基 Q_l∪P 上
│
│ P4(apply key) c̃_{0,1} = Σ_j ĉ_j ⊙ evk_j，逐 tower MAC
│     EvalFastKeySwitchCoreExt:634-657（P 段 evk 索引偏移 sizeQ，evk 存全 Q 基）
▼
2 × (N×(l+K)) = 2 × (N×5)             （c̃0、c̃1 各一份）
│
│ P5 ModDown：÷P，基 Q_l∪P → Q_l，详见 §4
▼
2 × (N×l) = 2 × (N×3)                 ← 输出密文，形状回到输入
```

总形状链（task.yaml 指定格式）：

    (N×l) → dnum×(N×α) → dnum×(N×(l+K−α)) → dnum×(N×(l+K)) → 2×(N×(l+K)) → 2×(N×l)

（中间 2×(N×(l+K)) 是 P4 输出；task.yaml 写的链把 digit 拼装与 P4 合并记了一步。）

**每个算子的 NTT 域要求**（决定 INTT/NTT 调用数，L1 算子账的根）：
- BConv：输入输出都在 COEFF（求和式是系数逐位的）；
- P4 点积：EVAL（点积即 NTT 域乘法）；
- 因此 ModUp 每 digit 必有 INTT(α) + NTT(l+K−α) 一组往返；
  其中 Q_j 段的 EVAL 版本在 assemble 时可直接复用 c1 原值（DC 分支 :473-478 还是
  多做了一次 NTT——这是 gap 分析里可抠的点）。

### 3.3 参数代入核对

l=3, α=2 → dnum=2（digit0: tower{0,1}，digit1: tower{2}，a_1=1）；
|C_0| = 3−2+2 = 3，|C_1| = 3−1+2 = 4。BConv 共 dnum = 2 次，输出 tower 数 3+4 = 7，
与 L1L2_推导v1.md 算子账（ntt_limb=7 中的 BConv 输出部分）及表 5-4 一致。

---

## 4. ModDown：推导 + 为什么不需要 digit 分解（题 3）

### 4.1 推导

输入 c̃ 在基 Q_l ∪ P 上。要算 round(c̃/P) mod Q_l。
把 c̃ 看成整系数多项式（其代表元在 Q_l·P 意义下）：

    c̃ = c̃_Q-part + c̃_P-part 的 CRT 组合

步骤（对应 ApproxModDown，dcrtpoly-impl.h:1315-1354）：

1. 取 P 侧 K 个 tower，INTT 回 COEFF；
2. BConv：P 基 → Q_l 基，得到 ŝ ≡ s + w·P (mod Q_l)，s 是 c̃_P 部分的 P-CRT 重构，
   **误差 |w| ≤ K/2**（§2.3 公式，输入基 = P，m = K）；
3. NTT 回 EVAL，与 Q 侧逐 tower 相减并乘 P^{-1} mod q_i：

       out_i = ( c̃ mod q_i − ŝ mod q_i ) · P^{-1} (mod q_i)

4. 代数验证：由 CRT，c̃ mod q_i − s mod q_i ≡ c̃ − s ≡ 0 (mod P)（因为 s 就是 c̃
   在 P 基上的重构，c̃ ≡ s mod P），所以差可被 P 整除，除以 P 后得：

       out ≡ c̃/P + w (mod Q_l)

   ModDown 自身引入的噪声 = w，**每系数 |w| ≤ K/2 = O(K)**——与 Q 无关。
5. round 的半单位舍入 ‖c̃‖/(2P) 项在实现里被吸进近似（「Approx」的另一半含义），
   同样是小常数级。

### 4.2 为什么 ModDown 不需要 digit 分解（validation 3）

三层理由：

1. **基的大小方向相反**。BConv 的误差是输入基的积的倍数。ModDown 的输入基是 P
   （K 个素数，K 很小），输出基是 Q_l（大）。从【小基】转到【大基】，误差 = K·P，
   除以 P 后只剩 O(K)。digit 分解的作用是把【大】输入基切小以压误差——这里输入基
   已经天然很小，没有可切的东西。
2. **P 本身就是一个不可再分（也不需要再分）的辅助基**。P 的 K 个素数就是「digit」
   的全部；再分只会让误差项从 K·P 变成 Σ(子基积)，而任何划分的子基积之和 ≥ K·P^{...}
   不会更优，反而增加 BConv 次数。
3. **误差预算已经够**。ModDown 噪声 O(K)（每系数 ~2·2^0）远小于 ModUp 带进来的
   α·Q_j/P ≈ O(α) 噪声和密文本底噪声，不是瓶颈，不值得为它加复杂度。

**对账**：ApproxModDown 里对 partP 只调一次 ApproxSwitchCRTBasis（P→Q 方向，
PHatInvModp / PHatModq 预计算，rns-cryptoparameters.cpp:203-216），全程无
numPartQ 概念。

### 4.3 全链路噪声预算（把题 1 的误差界用掉）

P4 输出（整系数视角）：

    c̃ = Σ_j ĉ_j · evk_j = c1·P·s′ + Σ_j u_j·Q_j·evk_j + (evk 本底噪声项)

ModDown 后：

    out ≡ c1·s′ + Σ_j u_j·Q_j·evk_j / P + O(K) (mod Q_l)

**关键噪声项 ‖Σ_j u_j·Q_j·evk_j/P‖ ≈ β · ‖c1‖ · Q_j / P**（β = digit 范数界）。
要让这一项有界（≤ β·‖c1‖），必须：

    **P ≥ Q_j = max digit**（对所有 j；最后一个残缺 digit 更小，故取 max）

此时 ModUp 近似误差被 P 完全吃掉，每系数残差 O(α)。若 P < Q_j，噪声按
2^{60(α−K)} 爆炸——正确性直接崩。这就是「P ≥ max digit」与「BConv 是近似」之间的
定量关系：**BConv 溢出误差 ~ α·Q_j 是因，ModDown ÷P 是消因的手段，P ≥ Q_j 是消干净的条件**。

**对账**：rns-cryptoparameters.cpp:129-136——maxBits = 单个 q 的最大位宽（注：实现取的是
素数位宽，故 P ≈ 2^{K·auxBits} ≥ 2^{α·qbits} 的满足由 dnum 划分保证，即每个 digit 的
总位宽 ≤ K·auxBits），sizeP = ⌈maxBits/auxBits⌉。K 由 α 和位宽推出，**不是自由参数**。
（严格说代码里 maxBits 是单素数位宽，digit 总位宽 α·qbits 的约束落在
numPartQ/α 的选取上：选 α 使得 α·qbits ≤ K·auxBits 成立——两条约束联立才闭合
P ≥ max digit。此处与代码完全一致，推导无矛盾。）

---

## 5. 三条 validation 汇总（自检清单）

1. **BConv 求和式**（§2.2）：ŷ_j = Σ_i x_i·QHatInv_i·QHatModp[i][j] mod p_j。
   共享量：x_i = c[i]·QHatInv_i（j 循环外算一次）；列独立：QHatModp[i][j]、p_j、
   Barrett mu_j。对账 dcrtpoly-impl.h:1233-1254 ✓
2. **目标基**（§3.1）：(Q_l \ Q_j) ∪ P，因为 P4 逐 tower 点积要求 ĉ_j 与 evk 同在
   Q_l∪P 全基上有值；只扩 P 则 Q 侧其余 tower 缺项。近似代替精确 lift 的代价
   u·Q_j 由 ModDown 消化。对账 rns-cryptoparameters.cpp:240-272 ✓
3. **P ≥ max digit**（§4.3）：BConv 溢出误差 ~ α·Q_j，ModDown 除以 P 后残差
   α·Q_j/P，P ≥ Q_j 时 ≤ O(α)；K = sizeP 由 ceil(maxBits/auxBits) 定死，非自由参数。
   对账 rns-cryptoparameters.cpp:129-136 ✓

---

## 6. 推导中浮出的新问题（供 task.yaml 收录）

- **OQ-1**：DC 分支 assemble 时对 Q_j 段多做一次 NTT（:473-478），而 c1 进函数时
  本来就是 EVAL——这次往返是否可省（CiFlow 的 INTT 复用 gap #1 的另一半）？
- **OQ-2**：§4.3 指出「P ≥ max digit」实际由两条约束联立（sizeP 公式 + α 划分），
  L1L2_推导v1.md 的参数账可以补一行闭合验证：α·qbits ≤ K·auxBits 是否在当前参数成立
  （2×60 ≤ 2×60，恰好压线——压线意味着 u 的最坏界没有余量，值得记一笔）。
