# Archived Notes

这里存放项目早期、已被取代或与当前架构脱节的 md 笔记。**保留以备追溯，不再维护**。需要了解当前实现请回到 [../PROJECT_STATUS.md](../PROJECT_STATUS.md)。

| 文件 | 归档原因 | 当前应看哪里 |
|------|----------|-------------|
| `NTT_COMPREHENSIVE_REPORT.md` | 描述的是**标准变几何 NTT**（`ntt_kernel.cpp` 的 `Compute_NTT`），已被 CG-NTT 替代为主力路径 | [../CG-NTT-Migration-Report.md](../CG-NTT-Migration-Report.md) |
| `EXPLORATION_FINDINGS.md` | 早期标准 NTT 探索总结，206% BRAM 超标的现状已通过 CG-NTT 重构解决 | [../CG-NTT-Migration-Report.md](../CG-NTT-Migration-Report.md) |
| `ntt_exploration_summary.md` | 同上来源的展开版，与 `EXPLORATION_FINDINGS.md` 内容重复 | [../CG-NTT-Migration-Report.md](../CG-NTT-Migration-Report.md) |
| `CG_NTT_IMPLEMENTATION_PLAN.md` | CG-NTT **实施前**的规划稿，已被实施完成报告取代 | [../CG-NTT-Migration-Report.md](../CG-NTT-Migration-Report.md) |
| `cg_ntt.md` | CG-NTT 端到端重构思路稿，与 Migration-Report Section 3/4 重复 | [../CG-NTT-Migration-Report.md](../CG-NTT-Migration-Report.md) |
| `综合后修改方案.md` | 针对**旧标准 NTT** 综合后 BRAM 超标/时序违例的修复链路，13 个 P0/P1 问题大多在 CG-NTT 重构时一并解决 | [../CG-NTT-Migration-Report.md](../CG-NTT-Migration-Report.md) + [../cg_ntt_optimization.md](../cg_ntt_optimization.md) |
| `HW_XCLBIN_ANALYSIS.md` | 早期 Vitis HLS 综合失败分析（2026-04-13），列出的 BRAM 206% / `OP_MUL` 拼写 / `u55c.cfg part=` 空等问题已修复 | [../cg_ntt_optimization.md](../cg_ntt_optimization.md)（看现在的 L1/L2/L3 优化轨道） |
| `fix.md` | `MultMod` 动态移位+128bit 加法时序违例 / AXI MUX / `addr++` 突发优化的修复指南，建议在解决新的时序问题前作参考资料 | 同上 |

## 归档时间

2026-06-29 — CG-NTT 迁移已稳定、HKS 三策略实测已完成后整理。
