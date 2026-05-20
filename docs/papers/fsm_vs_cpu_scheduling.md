# FPGA FSM vs CPU 调度：决策框架与学习笔记

## 一、核心判据

FPGA 侧 FSM 和 CPU 调度的选择，本质上只取决于一个比值：

$$\text{开销比} = \frac{\text{CPU 调度开销}}{\text{单次计算周期}}$$

- 开销比 < 0.5%：CPU 调度完全无感。
- 开销比 1%~10%：边界区域，视调用频率决定。
- 开销比 > 10%：必须 FPGA 内部 FSM 或合并为大 op。

## 二、CPU 调度开销的构成

CPU 向 FPGA 发一次操作的实际代价（以 Zynq/PCIe 平台为例）：

| 开销项 | 量级（cycles） |
|--------|:---:|
| AXI-lite 写寄存器（opcode + 参数） | 10–20 |
| ap_start 握手 | 5–10 |
| 中断等待 + ISR 上下文切换 | 500–2,000 |
| 用户态↔内核态切换（UIO） | 1,000+ |

一条指令合计约 **2,000–3,000 cycles** 的纯开销。

## 三、什么任务适合 CPU 调度

### 3.1 粗粒度、确定性操作

当前项目各 op 的计算量：

| Op | 计算量（cycles） | 开销比 | 结论 |
|----|:---:|:---:|------|
| NTT (N=65,536) | ~300K | ~1% | CPU 调度 |
| BCONV | ~500K+ | <1% | CPU 调度 |
| MULT (逐 limb 模乘) | ~10K–50K | ~6–15% | 边界但可接受 |
| 单次模乘 | ~5 | 60,000% | 必须合并/FSM |

只要计算远大于调度开销，CPU 调度就是最优解——你的场景完全符合。

### 3.2 复杂决策逻辑

```cpp
if (security_level == 128)
    chain = chain_128;
else if (security_level == 192)
    chain = chain_192;
```

这类判断用 CPU C++ 表达比硬件 FSM 可维护得多。

### 3.3 多用户 / 多任务仲裁

多个 CPU 核共享一块 FPGA 时，CPU 侧调度天然拥有锁、队列、优先级机制，无需硬件重实现。

### 3.4 开发灵活性与调试

- 改操作顺序：CPU 改一行源码，几秒编译
- 加新 op：写一个 case，不用重综合
- 单步调试：GDB + cosim，周期级可观测

FPGA 内 FSM 每次改序列 → 重综合（几小时）。

## 四、什么任务必须用 FPGA 内部 FSM

### 4.1 细粒度迭代

```cpp
// 反例：1024 次调用，每次 50 cycles 计算 + 2,000 cycles 调度
for (int i = 0; i < 1024; i++)
    Top(..., OP_SINGLE_MODMULT, ...);
// 总时间 = 1024 × 2050 = 2.1M cycles，仅 2.4% 在做计算
```

解法：要么合并成批量 op（一次 `OP_BATCH_MODMULT`），要么在 FPGA 内部用 FSM 循环。核心思想：**把循环边界从 CPU/FPGA 通信边界拉到 FPGA 内部**。

### 4.2 流水线级联（无气泡背靠背）

```
Load(a) → NTT → MULT → INTT → Store
```

| 方案 | 阶段间延迟 |
|------|:---:|
| CPU 逐个调度 | ~2,000 cycle 气泡 |
| FPGA 内部 FSM | 0 cycle（Done 直连下一个 Start） |

当前项目已将各阶段提取为独立子函数（`Load_XXX`、`Execute_NTT`、`Compute_Mult`），如日后需要 FSM 流水线，在子函数间插 `state++` 即可——接入成本低。

### 4.3 数据依赖分支

```cpp
res = compute_norm(poly);
if (res < threshold)
    skip_rescale();
else
    do_rescale();
```

分支由 FPGA 计算结果决定。CPU 调度需：计算完 → 中断 → CPU 读回 → 判断 → 写新 opcode → 再启动，往返 5,000+ cycles。FSM 0 cycle。

### 4.4 多操作异步并发

```cpp
// 同时启动三个子模块
NTT_on_bank0_start && MULT_on_bank1_start && Store_on_bank2_start;
```

FSM 可同时发出多个 `ap_start`，CPU 调度只能串行发指令。

## 五、决策流程图

```
单次计算量 < 10,000 cycles?
  ├── YES → 合并成大 op，或 FPGA 内部 FSM
  └── NO
         ↓
  单次 10K–100K cycles + 高频调用?
    ├── YES → 考虑 FPGA 内部 FSM
    └── NO  → CPU 调度即可
            ↓
    复合重复调用的计算过程?
    ├── YES → 整合 opcode（阶段 2）
    └── NO  → 回到开销比判断：开销比 > 阈值 → FSM
                                 开销比 ≤ 阈值 → 接受阶段 1
```

## 六、架构演进路径

```
阶段 1（现状）:  CPU 逐个发 opcode
                 最灵活，每个 op 计算量 >> 调度开销
                 ↓
阶段 2（过渡）:  CPU 发 "复合 opcode"
                 如 OP_NTT_MULT_INTT_CHAIN，内部顺序调用已有子函数
                 不改硬件架构，消除热点路径往返开销，性价比最高
                 ↓
阶段 3（极致）:  FPGA 内部完整 FSM + 流水线
                 CPU 只发高层语义命令
                 需频繁重综合，灵活性最低
```

对本项目而言，**阶段 2 性价比最高**：在 `switch(opcode)` 中新增几个复合 case，顺序调用已有的 `Load()` / `Execute_NTT()` / `Compute_Mult()` / `Store()` 即可，无需改变 `top.cpp` 的顶层结构。

## 七、一句话总结

> FSM 解决的是"CPU 调度开销相对计算量不可忽略"的问题。
> 粗粒度 FHE 算子（NTT/BCONV）的计算量天然压倒通信开销，
> 因此 opcode-switch RPC 模型就是最优设计。

## 八、代码注释：当前 `top.cpp` 的 FSM 结构

```cpp
// 本质：opcode-driven RPC dispatcher，而非自主跳转 FSM
void Top(const uint8_t opcode, ...) {
    switch(opcode) {
        case OP_INIT: Load_Init_Params(...); break;     // 只初始化，不计算
        case OP_ADD:  Load→Compute→Store;   break;     // 经典三段式
        case OP_NTT:  Load→Execute_NTT→Store; break;
        // ... 每个 case 独立，无跨 case 状态保持
    }
}
```

设计要点：
1. **子函数提取 + `INLINE off`**：将 AXI 访问隔离到子模块内部，顶层只保留 ap_start/ap_done 握手，消除 gmem 顶层巨型 MUX（缓解布线拥塞）。
2. **静态 BRAM 缓冲区**：`poly_buffer_1/2`、`result_buffer` 作为片上 scratchpad，Load→Compute→Store 三段之间不经过 AXI。
3. **无跨调用状态**：每次 `Top()` 执行一个 opcode 即返回，编排逻辑完全在 CPU 侧。