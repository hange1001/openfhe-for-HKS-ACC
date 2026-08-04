# OC 策略 vs CiFlow OC：差距分析与修复路径

> 参考文献：Neda et al., *CiFlow: Dataflow Analysis and Optimization of Key Switching for Homomorphic Encryption*, ISPASS 2024.
> 分析时间：2026-06-29

## 1. 一句话结论

当前 [keyswitch-hybrid.cpp:501-587](../../src/pke/lib/keyswitch/keyswitch-hybrid.cpp) 的 OC 分支**不是 CiFlow OC**，是"DC + sizeP 倍冗余 BConv"。这就是实测 OC 5.99ms > DC 4.47ms 的根因，详见 [实验规划.md §补5-1 表5-4](../papers/实验规划.md)。

## 2. CiFlow OC 的五个机制 vs 当前实现

| # | CiFlow OC（论文 §IV.C） | 当前实现 | 状态 |
|---|------------------------|----------|------|
| 1 | INTT 输出 on-chip 复用（所有 output tower 共享一份 INTT） | `keyswitch-hybrid.cpp:506-526` 预先 INTT 所有 digit | ✅ **唯一贴合的一项** |
| 2 | BConv 真正只算 1 个 output tower（粒度切到 tower） | `:534-545` 调用 `ApproxSwitchCRTBasis()` 算全 sizeP 然后丢弃 | ❌ **核心 gap** |
| 3 | **Section1 bypass**：算 mod Q output 时，第 `p/alpha` 个 digit **不过 BConv 也不过 NTT**，直接复用原始 EVAL 系数 | 所有 `part` 都过完整 BConv | ❌ 丢失 ~1/dnum 计算量 |
| 4 | Section1 / Section2 分两段（mod Q output 用 dnum-1 digit，mod P output 用全 dnum digit） | 单层 `for p in sizeP` 不区分 | ❌ |
| 5 | ModUp P5 partial sum **on-chip 累加**，只回写最终 tower | 直接写到 `partsCtExt[part][extPIdx]`，reduce 留给后面的 `EvalFastKeySwitchCoreExt` | ⚠️ 次要 |

## 3. 数据流对比图

```
CiFlow OC（论文 Fig.1 红色 tower 轨迹）：
  INTT 所有 digit                                          ← 已实现 ✓
  for p in [0, sizeQl + sizeP):
      if p < sizeQl:    # Section1 (mod Q output)
          bypass_digit = p / alpha
          for d in digits except bypass_digit:
              BConv_singleTower(d, target=p) → partial    ← 未实现：仍 BConv 全 sizeP
              partial → 累加器                            ← 未实现：直接写 partsCtExt
          partsCtExt[bypass_digit][p] = partsCt[bypass_digit][p_local]
                                                          ← 未实现：所有 digit 一视同仁
      else:             # Section2 (mod P output)
          for d in all_digits:
              BConv_singleTower(d, target=p) → partial
              partial → 累加器
      flush 累加器 → output[p]

当前"OC"：
  INTT 所有 digit                                          ✓
  for p in sizeP:                                          ← 只跑 sizeP 不跑 sizeQl+sizeP
      for d in all_digits:
          fullCompl = BConv(d, allTowers)                  ← 算 sizeP 倍
          pick one tower from fullCompl
          NTT(one tower) → 写 partsCtExt[d][sizeQl+p]
```

## 4. 性能算术

**当前参数**：N=4096, sizeQl=3, sizeP=2, alpha=2, numPartQl=2

| 项 | CiFlow OC | 当前"OC" | 当前 DC |
|----|-----------|----------|---------|
| BConv 次数 | `(sizeQl + sizeP) × (numPartQl - 1) + sizeP × numPartQl` ≈ **9 次 single-tower** | `sizeP × numPartQl = 4` 次 **full-sizeP** ≈ **20 次 single-tower 等价计算** | `numPartQl = 2` 次 full-sizeP ≈ **10 次 single-tower 等价计算** |
| 等价计算量 | 9 | 20 | 10 |
| 实测延迟 | — | 5.99 ms | 4.47 ms |

**这就是 OC 比 DC 还慢的算术解释**：BConv 计算量直接乘了 sizeP（=2 here，5 in BTS 参数）倍。

## 5. 修复可行性（A 部分的研究结论）

### 5.1 CPU 路径

读了 [dcrtpoly-impl.h:1150-1281](../../src/core/include/lattice/hal/default/dcrtpoly-impl.h)，`ApproxSwitchCRTBasis` 的快速路径（`HAVE_INT128 && NATIVEINT==64`）核心循环：

```cpp
#pragma omp parallel for firstprivate(sum)
for (uint32_t ri = 0; ri < ringDim; ++ri) {
    std::fill(sum.begin(), sum.end(), 0);
    for (uint32_t i = 0; i < sizeQ; ++i) {
        const auto xQHatInvModqi = m_vectors[i][ri].ModMulFastConst(...);
        for (uint32_t j = 0; j < sizeP; ++j) {          // ← 内层 j 完全独立
            sum[j] += Mul128(xQHatInvModqi, QHatModp[i][j].ConvertToInt<uint64_t>());
        }
    }
    for (uint32_t j = 0; j < sizeP; ++j) {
        ans.m_vectors[j][ri] = BarrettUint128ModUint64(sum[j], pj, modpBarrettMu[j]);
    }
}
```

**关键观察**：output tower `j` 的计算只依赖 `xQHatInvModqi`（所有 j 共享）、`QHatModp[i][j]` 第 j 列、`modpBarrettMu[j]`、`paramsP[j]`。**j 之间零依赖**。

**single-tower 重载的最小骨架**（不写实现，只论证可行性）：

```cpp
// 提议新增到 dcrtpoly.h / dcrtpoly-impl.h
PolyType ApproxSwitchCRTBasisSingleTower(
    const std::shared_ptr<Params>& paramsQ,
    const std::shared_ptr<NativeInteger>& target_pj,          // 目标输出模数
    uint32_t target_idx,                                       // QHatModp 列索引
    const std::vector<NativeInteger>& QHatInvModq,
    const std::vector<NativeInteger>& QHatInvModqPrecon,
    const std::vector<std::vector<NativeInteger>>& QHatModp,
    const DoubleNativeInt& modpBarrettMu_target
) const;
// 内部把 sum[] 退化为单值，删除 for j 循环
```

| 项 | 评估 |
|----|------|
| 改动范围 | 单文件，~80 行新代码 |
| 向后兼容性 | ✅ 纯新增 API，不破坏老接口 |
| OpenMP 并行性 | ✅ 仍按 ringDim 切，无退化 |
| 性能预期 | CPU 路径 BConv 时间从 ~469 μs → ~94 μs（5× 缩减，sizeP=5 时） |

### 5.2 FPGA 路径

读了 [FpgaManager.h:610-718](../../src/core/include/FpgaManager.h) 和 [dcrtpoly-impl.h:1167-1227](../../src/core/include/lattice/hal/default/dcrtpoly-impl.h)：

```cpp
// dcrtpoly-impl.h:1175
if (FpgaManager::GetInstance().IsReady() && sizeQ <= KERNEL_LIMB_Q && sizeP <= KERNEL_MAX_OUT_COLS) {
    // → 走 FPGA
}
```

**`sizeP=1` 已经在合法范围内**，FPGA path **不用改 HLS** 就能跑 single-tower BConv。

| 项 | 当前 (sizeP=5) | single-tower (sizeP=1) | 变化 |
|----|----------------|------------------------|------|
| H2D 输入 | 3 × 4096 × 8B = 96 KB | 96 KB | **不变**（输入仍是全 sizeQ） |
| D2H 输出 | 5 × 4096 × 8B = 160 KB | 32 KB | **5× 缩减** |
| Kernel cycles | 8,232 | ~3,600 | **2.3× 缩减** |
| DMA 固定开销 | ~20 μs | ~20 μs | **不变** ⚠️ |
| Systolic PE 利用率 | 100% | 20% (4/5 PE 浪费) | 浪费但不影响功能 |

**FPGA 侧 single-tower kernel 是 L4 可选优化**（专门的 1-tower BConv 微架构能把 PE 利用率拉回 100%，但当前不做也能跑通）。

### 5.3 综合改造成本

| Gap | 修复点 | 成本 | 预期增益 |
|-----|--------|------|---------|
| #2 (BConv 粒度) | 加 CPU `ApproxSwitchCRTBasisSingleTower` + 改 OC 分支调用 | 中（OpenFHE 改动） | 计算 5× / D2H 5× |
| #3 (Section1 bypass) | 纯 Host 逻辑，OC 外层判断 `if (p < sizeQl) skip digit p/alpha` | **低** | 计算 1/dnum 节省 |
| #4 (Section1/Section2 区分) | 改 OC 外层循环范围 `for p in [0, sizeQl + sizeP)`，配合 #3 | 低 | 与 #3 联动 |
| #5 (on-chip 累加) | OC 内层加 `tower_accumulator` 局部变量 | 低 | PCIe 流量缩减（次要） |

**推荐顺序**：#3 → #4 → #2 → #5。`#3 + #4` 是纯 Host 逻辑改动，**不动 OpenFHE 也不动 FPGA**，应该最先验证。#2 是最大收益项，但要碰 OpenFHE 库。

## 6. 待办（下一步）

- [ ] **Gap #3 最小补丁**：在 `keyswitch-hybrid.cpp:530` 的外层循环里加 bypass 判断（`p < sizeQl` 时第 `p/alpha` 个 digit 跳过 BConv 复用原始 EVAL 系数）。验收标准：实测 OC 延迟应从 5.99 ms 降到 ~4 ms 级。
- [ ] **Gap #4 改循环范围**：`for p in [0, sizeQl + sizeP)` 并区分 Section1/Section2，配合 #3 一起验证。
- [ ] **Gap #2 OpenFHE 重载**：在 `dcrtpoly.h` / `dcrtpoly-impl.h` 加 `ApproxSwitchCRTBasisSingleTower`，先验证 CPU 路径数值正确性（vs 现有全 sizeP 的对应输出 tower），再切到 OC 分支。预期最大收益。
- [ ] **修完后回归测试**：重跑 `hks-benchmark --strategy OC --repeat 100`，更新 [实验规划.md 表 5-4](../papers/实验规划.md) 实测数据，并核对答辩 PPT Slide 19 的"OC 延迟略优"论述。

## 7. 注意事项

- **当前论文 [content.md §5.3 调度策略](../papers/content.md) 的论述**与实测数据冲突。修复 OC 后需要回头同步章节文字与表 5-4。
- 答辩 PPT Slide 19 写"OC 延迟略优"是预期，不是当前实测——这是论文/答辩材料里的隐性 bug。
- CiFlow 原文测的 OC 加速比 1.30× ~ 4.16× 是在 RPU 向量处理器架构 + 大 on-chip SRAM (32MB) 场景下，**与我们的小参数 (N=4096) + PCIe-bound FPGA 场景不可直接套用**。我们能拿到的收益主要来自 #2 的计算量节省，而非 CiFlow 强调的"片外带宽节省"。
