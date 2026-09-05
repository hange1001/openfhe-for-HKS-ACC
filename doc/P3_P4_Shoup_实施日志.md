# P3 / P4 / Shoup 实施日志

## 目标与提交节点

用户授权：依次完成 P3、P4、Shoup；每个节点验证后独立 Git 提交，不推送。
保持当前 PE_PARALLEL=4；不修改、不停止另一个任务的 PE 扫参。

| 节点 | 目标 | 状态 | Git |
|---|---|---|---|
| P3 | 工作区直接作为变换 A bank；一条 scratch；预乘原位；消除变换装卸拷贝 | 完成 | `53cf757`（源码/功能/RTL），后续提交补物理证据 |
| P4 | 预乘复用现有蝶形模乘 lane，顺序切换并排空流水线 | 实施中 | — |
| Shoup | BConv 定常数模乘改为 Shoup，更新两个入口及 OpenFHE 元数据/验证 | 待开始 | — |

证据分层：native/OpenFHE 精确数值 → HLS 端口/II/资源 → RTL 精确结果/周期 → OOC 布线。
没有板卡，不将仿真换算时间称作实测端到端速度，也不把综合 slack 当作布线结果。

## 2026-09-04：开工与预测

- 基线提交 `32c04b5`（P2），工作区硬件源码初始干净。保留已有未提交决策和分析文档。
- 当前另一个任务正在对 PE32 做物理实现；本轮采用独立项目标签和实现工作目录。
- P3：`poly_buffer_1` 作为 A_work，输入物理行为 digit-local，模数编号保持 global；
  BConv 输出物理行 `LIMB_Q+p`。输出选择在 AXI 写回边界完成，不增加片上整塔复制。
- 保留通用双操作数所需的 `poly_buffer_2` 和原始 EVAL / 通用结果所需的 `result_buffer`；
  不新增整套工作数组。增加一条 scratch 的同时删除旧 bank_a、bank_b。
- 测量前预测（PE4）：P3 暖态约 95,000 周期，±15% 为 80,750～109,250；
  DSP 仍 1,160；BRAM_18K 预计约 424（440 减一个 16-block tower），允许 ±15%。
  实际 RAM 端口绑定与 HLS 分区传播可能使 BRAM 预测偏离，届时逐项解释。
- P4 预测：在 P3 上再减少约 9,216 个预乘主体周期；Barrett 模乘实例 20→19，DSP 1,160→1,102。
- Shoup 预测：只替换 BConv 的 15 路；每路 DSP 58→约36，合计再减约330 DSP。
  这些均为预测，不能作为验收结果；BConv 原有 II=1，不预言其主体循环等比例加速。

## 验收门槛

每个节点都要保持 Top/HKS/独立 BConv 回归、真实 OpenFHE 精确余数对比、非法请求与哨兵检查；
保存 HLS 资源、II、可达 RTL 实例审计与 RTL smoke/perf。结构优化必须以生成 RTL 为证据。
P3 检查直接工作区、单 scratch、无 TRANSFORM_LOAD/STORE；P4 检查无独立预乘模乘实例；
Shoup 检查 x≥目标模数、边界/随机向量、元数据版本保护、两个 BConv 入口均使用新模乘。

## P3 第一次综合：发现并拒绝调度退化

- native Top 18/18、HKS 22、OpenFHE Release 与 ASan 各 472 checks / 1,523,712 residues 通过。
- r1 BRAM_18K=424，符合存储预测；但蝶形 II=5，单变换 31,045 周期，**不予验收**。
- 原因定位：Vitis 把运行时方向控制的写地址判为跨迭代 WAW，串行化后仅保留一路蝶形模乘；
  不是 BRAM 实际带宽不足。两个地址映射都是 [0,N) 的双射，独立检查器已穷举全部槽位/方向/级。
- 修正只解除上述 inter-WAW，不复制旧版整套 inter/intra false；保留 RAW 和 stage 排空约束。
  使用新标签 r2 重新综合，保留 r1 作为失败证据，不能用功能通过替代性能验收。
- r2 仅解除 WAW 后 II=3，定位到读端口调度：先选择地址再访问 RAM 使编译器看不清互斥 bank。
  r3 把方向判断保留在 RAM 访问分支，让每条访问的 bank 低位可静态推导；并行度不变。
- r3 结构审计通过：4 路、II=1、单 scratch、AXI256；BRAM424 / DSP1160 / FF81549 / LUT180505 / URAM96。
  核心 6,481 周期，含 wrapper 6,483 周期。LUT 比 P2 增 25,653，来自多塔工作区的运行时选通；待布线检验。
- 尝试在重新打开的 csynth 工程中更换 TB 复测，Vitis 出现 wrapper 链接缺失；复制部分工程做 export 也失败。
  已撤下这个不可靠的辅助入口，恢复既有完整 cosim 工程流程；失败不计入功能结果。
- 最终r3：C-sim通过；RTL smoke 35次调用/4种digit通过；OpenFHE RTL 40,960余数通过。
  两个digit为52,779+47,643=100,422周期，相对P2减少20,400（16.88%）；INIT295,063不变。
  smoke/perf生成Verilog逐文件相同。物理实现正在独立目录运行，等待严格时序结果再签核P3。

## 等待P3布线时的独立Shoup预验证（尚未接入Top）

- 新原语仅增加独立文件，没有进入P3的KERNEL_SRC，不改变P3生成RTL。
- native及RTL共20,704个算术向量精确通过；18,013个输入≥目标模数，4,702个需要末次减法，
  另检验p=0禁用返回零、零权重、最大uint64输入、接近2^63的模数及逐调用参数切换。
- 实际HLS原语：36 DSP（16+10+10）、898 FF、948 LUT、II=1、latency12周期；
  Barrett原语为58 DSP/latency19周期。与资源预测吻合，但整机BConv还未接入，不能据此宣称整机收益。

## Shoup 协议准备与物理任务协调

- 采用显式ABI v2：HKS为2字头+15权重+15precon+3逆元=35字（280B）；
  独立BCONV为2字头+15权重+15precon+5目标模数=37字（296B）。原先33字估计仅包含HKS有效载荷，未含头。
- 共享协议常量使用不依赖ap_int/OpenFHE类型的头文件。当前仅预备此文件，尚未接入P3。
- 规划新增只读能力查询opcode9，保留0..8编号与退役7；新Host先确认ABI/固定参数形状再运行，
  防止旧bitstream误解释新元数据。新kernel先检查高位magic和版本/长度，拒绝旧元数据而不越界读取新增字段。
- 预计算只在Host的参数准备阶段发生；按目标列的模数计算，缓存键覆盖有序上下文/分片/权重/逆元。
  不改OpenFHE CPU数学参考，不将QHatInv预乘搬回CPU，不引入Montgomery域转换。
- P3物理综合通过，布局布线前与用户PE16物理扫参共用16GiB内存，交换区已接近满载。
  仅暂停本轮P3的Vivado进程PID376659（项目hks_digit_cosim-perf_mem_p3_r3）；用户扫参保持运行。
  待内存释放后必须恢复该进程并继续物理验收，当前不能宣称P3物理通过。
- 随后确认暂停进程仍持有约4GB，系统交换区约7.3/8GB；为释放内存，已正常终止本轮PID376659。
  综合DCP和Vivado工程完整保留；新增 `hks_resume_impl.tcl` 仅从已完成synth_1继续impl_1，禁止自动reset/resynth。
  后续应通过该入口继续P3，不再向已结束的PID发CONT。用户PE16任务未被干预。

## 2026-09-05：P3 最终签核

- P3 源码、功能与 RTL 证据已保存在 `53cf757`。该提交同时包含尚未接入 Top 的 Shoup
  原语 WIP；不改写已有历史，后续仍按 P4、Shoup 整机接入分别提交。
- perf fixture 在 `config_cosim -random_stall` 下重放通过，说明三个 AXI master 出现随机
  背压时仍能保持结果精确；周期比较继续采用无 stall fixture，避免把随机等待混入架构周期。
- 恢复路径完成 OOC implementation：235,225/235,225 条可路由 net 全部布通，routing
  error 为 0。默认 6 ns 的 WNS/TNS 为 +0.029/0 ns；6 ns +0.75 ns uncertainty
  为 -0.721/-595.871 ns，不通过；7 ns +0.75 ns 为 +0.279/0 ns，通过。
  三种情景 hold 均为 WHS +0.007 ns、THS 0、0 个失败端点。
- 布线后资源：CLB LUT 108,392、register 64,058、BRAM tile 272、DSP 1,160、
  URAM 96。HLS LUT 180,505 不能与布线后 LUT 直接相减，两套口径分别报告。
- 最差 setup 路径转到 BConv `Load_W`，5.951 ns 中 route 占 5.436 ns（91.346%）；
  P3 删除的旧整塔装卸路径没有重新出现。至此 P3 按保守 7 ns 情景完成验收，进入 P4。
