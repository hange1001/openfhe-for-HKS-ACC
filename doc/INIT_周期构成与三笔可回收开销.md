# INIT 阶段 295,063 周期的构成与三笔可回收开销

日期：2026-09-05。基线：P3（两 digit 暖态 100,422 周期）。
方法：只读 HLS 报告与源码，未做任何改动，未跑综合。
证据：`docs/reports/hls/hks_mem_p3_20260904/{top-csynth.rpt, perf-transactions.rpt}`、
`src/fpga_backend/src/top.cpp`、`src/fpga_backend/include/define.h`。

---

## 0. 为什么现在看 INIT

P1/P2/P3 三轮优化把两个 digit 从 139,734 压到 **100,422** 周期（−28.1%），
而 **INIT 始终是 295,063 周期，一拍未动**。

| | 周期 | 相对两 digit |
|---|---:|---:|
| INIT（冷启动） | 295,063 | **2.94×** |
| 两个 digit 暖态（P3） | 100,422 | 1.00× |

**初始化的代价已经接近有效工作的三倍**，而且随着 digit 越优化，这个比例只会更难看。
本文把这 295,063 拆开，指出三笔可回收开销，并说明其中两笔与 PE 并行度上限
（ADR-013 的 URAM 墙）是同一个根因。

---

## 1. 现状拆解

`top-csynth.rpt` 中：

```
grp_Load_Init_Params_fu_1805 | Load_Init_Params | 196674 | 1.180 ms
```

`Load_Init_Params` 占 **196,674 周期**，是 INIT 的 67%。它的主体是两个 twiddle 加载循环
（`top.cpp` `init_NTTTwiddle_Loop` / `init_INTTTwiddle_Loop`）：

```c
for (int l = 0; l < MAX_LIMBS; l++)        // 8
    for (int s = 0; s < STAGE; s++)        // 12
        for (int t = 0; t < CG_HALF_N; t++) {   // 2048
            #pragma HLS PIPELINE II=1
            NTTTwiddleFactor[l][s][t] = mem_in1[...];
        }
```

每个循环 `8 × 12 × 2048 = 196,608` 次迭代，II=1。
NTT 表走 `gmem0`、INTT 表走 `gmem1`，属不同 AXI bundle，HLS 可并行；
`Load_Init_Params` 报 196,674 ≈ 单个循环的量，说明两者**基本完全重叠**。

INIT 事务总计 295,063，减去 196,674 尚有 **98,389 周期未归属**
（≈ 196,608/2，与「两循环约 50% 重叠」一致，但未逐拍证实）。
按项目 traceability 规矩，此项标记为**未归属**，不写入任何加速比推导。

---

## 2. 三笔可回收开销

### ① 加载 8 个 limb，只用得到 5 个

```c
#define MAX_LIMBS (LIMB_Q + MAX_OUT_COLS)      // 3 + 5 = 8
static uint64_t NTTTwiddleFactor[MAX_LIMBS][STAGE][CG_HALF_N];
```

`MAX_LIMBS = 8` 是为 `poly_buffer` 定的——BConv 需要同时持有 3 个输入塔与 5 个输出塔。
但 **twiddle 是按模数索引的**，而模数只初始化了 5 个：

```c
for (int i = 0; i < LIMB_Q; i++)  MODULUS[i]           = mem_in1[i];      // 0,1,2
for (int j = 0; j < LIMB_P; j++)  MODULUS[LIMB_Q + j]  = mem_in2[j];      // 3,4
```

变换侧的 limb 下标也只到 4：`limb = digit_start + step`（0..2）与
`limb = p < digit_start ? p : p + alpha`（0..4）。

**结论：limb 5/6/7 的 twiddle 被从 PCIe 搬进片上、占着 URAM，然后永不被读——37.5% 是纯垃圾。**

> ⚠️ 落地前需复核：确认所有调用路径（含独立 OP_NTT/OP_INTT 与未来复合算子）
> 的 twiddle 下标都不超过 `MAX_OUT_COLS - 1 = 4`。本文只核了融合 digit 路径。

### ② 内层未展开，256-bit AXI 只用了 1/4

内层循环只有 `PIPELINE II=1`，**没有 `UNROLL`**，因此每拍搬 1 个 64-bit 字；
而 AXI 已经是 256-bit（ADR-011），每拍可搬 4 个。

写侧 `NTTTwiddleFactor` 是 `cyclic factor=PE_PARALLEL dim=3`，**4 个 bank 已经存在**，
不需要额外划分。对照同文件的 `Load_Transform_Tower` 已写了 `UNROLL factor=PE_PARALLEL`，
这两个循环是**漏了**，不是有意为之。

与 P1 修掉的 `Prepare_HKS_BConv_Input` 属同一类：pragma 遗漏。

### ③ `log N = 12×` 的 twiddle 冗余

CG-NTT 存 `STAGE × (N/2)` 个旋转因子，标准 NTT 只需 `N/2`——冗余 `log N = 12` 倍。
这是 ADR-003 拿存储换掉变 stride MUX 的**代价**，属设计决策而非缺陷；
去掉它需要重新引入按 stage 的寻址逻辑。

同一冗余已在 `open_questions q6` 中被指出会让 `N=2^16` 需要 67 MB 片上（154%，不可行）。

---

## 3. 收益量化

twiddle 数据量：`8 limb × 12 stage × 2048 × 8 B × 2 表 = 3.0 MB`。

| 措施 | 改动量 | INIT 周期 | URAM | PCIe 流量 |
|---|---|---:|---:|---:|
| 现状 | — | **295,063** | 96 (10%) | 3.0 MB |
| ① limb 8→5 | 循环界与数组维度改 `MAX_OUT_COLS` | ~184,000 | **60** | 1.9 MB |
| ② 内层 UNROLL 4 | **两行 pragma** | ~74,000 | 96 | 3.0 MB |
| **①+②** | 以上两处 | **~46,000（−84%）** | **60** | 1.9 MB |
| +③ 去冗余 | 重做 twiddle 寻址 | ~4,000 | **~8** | 160 KB |

> 表中 ①②的周期为按迭代数线性外推，**未经综合验证**；③为容量推算。
> 落地时须按 `predict_before_measure` 规矩先写预测再实测对账。

---

## 4. 与 PE 并行度上限是同一个根因

ADR-013 由扫参实测得出 `URAM(PE) = 24 × PE`，PE=32 时 80%，PE=64 时 160% 出界。
**这 96 块 URAM 全部来自这两张 twiddle 表**（`3.0 MB / 288 Kbit ≈ 85 块`，
加划分开销得 96，与实测一致）。

| | URAM @PE=4 | `URAM(PE)` | PE=64 需要 |
|---|---:|---|---:|
| 现状 | 96 | 24×PE | 1,536（**160% ✗**） |
| ① 之后 | 60 | 15×PE | 960（100%，勉强） |
| ①+③ 之后 | ~8 | ~2×PE | 128（**13% ✓**） |

**去掉 twiddle 冗余，URAM 墙从 PE=32 一路退到 PE=64 以外。**

而 ADR-013 已经算出：只提 PE 的周期上限是 2.46×（PE=32），
按 CPU 双线程 463 μs 对比**最好只能打平**。若 URAM 墙退开，PE=64 可达 2.75×——
仍然不够，但它把「并行度」这条路的天花板抬高了，值得与 P4 捆绑评估
（ADR-013 已记：三项收益不是相加，是分母一段段挖掉）。

---

## 5. 一个更大的结论

同一个 `log N` 倍 twiddle 冗余，在两个维度上都是硬约束：

- **小参数（N=4096）**：吃掉 URAM，把蝶形并行度锁在 PE=32
- **大参数（N=2^16）**：需 67 MB 片上，直接不可行（q6）

**ADR-003 用 `log N` 倍存储换掉变 stride MUX，这笔交易的代价可以在并行度和参数规模
两个维度上分别量化。** 这是一个可写进论文的 claim——它不是「我做了个 CG-NTT 快 7 倍」，
而是「CG-NTT 的收益有明确的适用边界，边界在哪我能算出来」。

对照 HERA (FPGA'26)：它把「NTT 分解」放在所有优化最前面，正是为了消掉这个冗余。

---

## 6. 建议优先级

| | 改动量 | 收益 | 风险 |
|---|---|---|---|
| **② UNROLL** | 两行 pragma | INIT −75% | 低。与 P1 同类，bank 已就位 |
| **① limb 8→5** | 循环界 + 数组维度 | INIT 再 −37.5%，URAM −36 | 中。**须先复核所有 twiddle 下标不超过 4** |
| ③ 去冗余 | 重做 twiddle 寻址 | URAM −90%，解锁 PE=64，同时救 N=2^16 | 高。触及 ADR-003 的核心权衡 |

②可立即做。①在复核下标后做。③应作为独立课题，与 q6 和 phase4 的 NTT 分解一并考虑。

---

## 7. 未决与边界

- 98,389 周期未归属，需逐拍波形或 HLS 调度报告确认两循环的实际重叠度。
- 本文所有周期收益均为**外推**，无一经过综合验证。
- INIT 是否值得优化取决于使用场景：一次初始化 + 多次 HKS 调用时它被摊薄；
  但当前答辩材料按「冷启动 INIT + 两个 digit」口径报数，此时它占 74.6%。
- 未评估 ① 对 `poly_buffer` 与 `MODULUS` 等其他 `MAX_LIMBS` 使用者的影响；
  这些数组确实需要 8，**不能一并改**。
