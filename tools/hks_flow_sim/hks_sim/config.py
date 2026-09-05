"""配置 schema、校验与 config hash。

三种 dataflow 必须共享同一份 HardwareConfig。fairness 由 hardware_config_hash
保证：三次运行的 hash 不同即判定对比无效。

数据来源标记（贯穿全仿真器）：
  measured    -- RTL co-sim / HLS 报告里直接读到的数
  calibrated  -- 由 measured 解出的模型参数
  projected   -- 外推到未标定的参数点，或未实现的模块
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace
from typing import Any, Dict, Optional

CONFIG_SCHEMA_VERSION = "0.1.0"

# 证据基线：P4 = commit e948a69，报告 docs/reports/hls/hks_mem_p4_20260905/
CALIBRATION_BASELINE = "P4/e948a69"


class ConfigError(ValueError):
    """非法配置。仿真器拒绝运行，不静默修正。"""


@dataclass(frozen=True)
class CostParams:
    """算子周期模型参数。

    measured（直接来自 hks_mem_p4_20260905/）：
      CG_Transform_Work  NTT/INTT  = 6481   transform-work-csynth.rpt
      CG_Transform_Work  SCALE     = 1047   同上，latency min
      Compute_BConv_Systolic       = 4145   top-csynth.rpt
      Load_Init_Params             = 196674 top-csynth.rpt
      OP_INIT 事务                 = 295063 perf-transactions.rpt

    calibrated（由 alpha=2 / alpha=1 两个事务解出，零残差）：
      axi_tower_extra    = 29    每塔 AXI = ceil(N*limb_bits/axi_bits) + 29
      control_per_digit  = 656   每 digit 固定控制开销

    hold-out：把这两个参数原样套到结构不同的 P3（预乘是独立硬件、1 系数/拍），
    两 digit 合计预测 100431 vs 实测 100422，误差 +0.009%。
    """

    # -- 变换引擎（measured @ P4 / PE=4 / N=4096） --
    butterfly_overhead: int = 26      # 538 = 2048/4 + 26
    stage_pair_overhead: int = 4      # WORK_STAGE_LOOP 单次迭代 = 2*538 + 4 = 1080
    transform_call_overhead: int = 1  # 6*1080 + 1 = 6481
    scale_loop_overhead: int = 22     # 1046 = 4096/4 + 22
    scale_call_overhead: int = 1      # 引擎口径 1047

    # -- BConv 阵列（measured @ P4） --
    bconv_overhead: int = 49          # 4145 = 4096 * ceil(2/3) * ceil(3/5) + 49

    # -- AXI / 控制（calibrated） --
    axi_tower_extra: int = 29
    control_per_digit: int = 656

    # -- 模乘原语（measured，见 doc/HKS_digit间并行与BConv资源拆解.md §2） --
    modmul_latency: int = 19
    modmul_ii: int = 1

    # -- KeyMult / Accumulate（projected：硬件未实现） --
    # T_KM(d,t) = passes * ceil(N/lanes) * II + keymult_pipe + keymult_mem
    # T_ACC(t)  = ceil(N/lanes) + accumulate_overhead
    # 方案 §5.2 把 T_acc 写进 T_KM；这里拆成两个 Op，二者之和等价，
    # 便于 M2 单独观察 accumulator 端口冲突。不得重复计数。
    keymult_mul_passes: int = 2       # c_j*a_j 与 c_j*b_j
    keymult_pipe: int = 0
    keymult_mem: int = 0
    accumulate_overhead: int = 0

    # -- 上下文级固定段（measured） --
    init_cycles: int = 295063

    def evidence_level(self) -> str:
        """该参数组的证据等级。KeyMult 未设定时不得输出单点结论。"""
        if (self.keymult_pipe, self.keymult_mem, self.accumulate_overhead) == (0, 0, 0):
            return "calibrated;keymult_unset"
        return "calibrated;keymult_projected"


@dataclass(frozen=True)
class WorkloadConfig:
    """CKKS/HKS workload 形状。

    手写大参数只表示 workload shape，不代表通过了 CKKS 安全性或精度验证。
    正式 industrial case 应由 OpenFHE CryptoContext 导出（source=openfhe_context）。
    """

    ring_dimension: int = 4096          # N
    q_towers: int = 3                   # L
    p_towers: int = 2                   # K
    num_part_q: int = 2                 # D
    alpha: Optional[int] = None         # None -> ceil(L / D)
    limb_bits: int = 64
    batch_count: int = 1
    level: Optional[int] = None
    source: str = "handwritten"         # handwritten | openfhe_context
    label: str = "p4_reproduction"

    def validate(self) -> None:
        n = self.ring_dimension
        if n < 2 or (n & (n - 1)) != 0:
            raise ConfigError(f"ring_dimension 必须是 >=2 的 2 的幂，得到 {n}")
        if self.q_towers < 1:
            raise ConfigError(f"q_towers 必须 >=1，得到 {self.q_towers}")
        if self.p_towers < 1:
            raise ConfigError(f"p_towers 必须 >=1，得到 {self.p_towers}")
        if not 1 <= self.num_part_q <= self.q_towers:
            raise ConfigError(
                f"num_part_q 必须在 [1, q_towers={self.q_towers}]，得到 {self.num_part_q}"
            )
        if self.alpha is not None and not 1 <= self.alpha <= self.q_towers:
            raise ConfigError(f"alpha 必须在 [1, q_towers]，得到 {self.alpha}")
        if not 1 <= self.limb_bits <= 64:
            raise ConfigError(f"limb_bits 必须在 [1, 64]，得到 {self.limb_bits}")
        if self.batch_count < 1:
            raise ConfigError(f"batch_count 必须 >=1，得到 {self.batch_count}")
        if self.source not in ("handwritten", "openfhe_context"):
            raise ConfigError(f"未知 workload source: {self.source}")

    @property
    def stages(self) -> int:
        """CG-NTT 级数 = log2(N)。"""
        return self.ring_dimension.bit_length() - 1

    @property
    def bytes_per_tower(self) -> int:
        """B_tower = ceil(limb_bits/8) * N。64-bit limb 下为 8N。"""
        return ((self.limb_bits + 7) // 8) * self.ring_dimension


@dataclass(frozen=True)
class HardwareConfig:
    """固定硬件配置 H。三种 dataflow 必须共享同一份。"""

    clock_period_ns: float = 7.0        # P4 保守签核情景（6ns+0.75ns 并未闭合）
    transform_engines: int = 1
    transform_lanes: int = 4            # PE_PARALLEL
    bconv_arrays: int = 1
    bconv_rows: int = 3                 # LIMB_Q
    bconv_cols: int = 5                 # MAX_OUT_COLS
    keymult_binding: str = "shared_transform_mul"   # 默认不新增 DSP
    axi_bits: int = 256
    axi_channels: int = 1               # P4 的 AXI 事务边界是串行的

    sram_capacity_bytes: int = 0        # 0 = 不限容量（M1 只统计 peak，不 spill）
    sram_banks: int = 8

    ddr_bandwidth_GBps: float = 0.0     # 0 = 不建模 DDR
    pcie_fixed_us: float = 0.0          # 无板卡：扫描维，不是标定量
    pcie_bandwidth_GBps: float = 0.0
    dma_compute_overlap: bool = False   # P4 的 XRT 路径完全阻塞（task.yaml q8）
    # 不同计算引擎能否并发。P4 的 Top 是 opcode-driven RPC dispatcher（ADR-008），
    # BConv 与变换引擎由顶层顺序调用，没有 DATAFLOW，因此默认 False。
    # 置 True 等于假设已经做了 phase4 step 4.1 的 DATAFLOW 拆分——那是另一个
    # 硬件配置点，hash 会变。
    #
    # 这一条对结论影响极大且**标定点分辨不出来**：DC 的单次 BConv 必然先于它的
    # 全部 NTT，永远不会有两个引擎同时 ready，所以 P4 的 91242 在 True/False 下
    # 都能对上；但 OC 的 single-tower 粒度会让 BConv/NTT 重叠。若默认给 True，
    # OC 会仅因为模型白送了一个硬件并不具备的并发能力而胜出。
    allow_engine_overlap: bool = False

    allow_tiling: bool = False
    strict_p4: bool = False             # 只允许 P4 已 RTL 验证过的形状

    cost: CostParams = field(default_factory=CostParams)

    def validate(self, wl: WorkloadConfig) -> None:
        if self.clock_period_ns <= 0:
            raise ConfigError("clock_period_ns 必须 > 0")
        for name in ("transform_engines", "transform_lanes",
                     "bconv_arrays", "bconv_rows", "bconv_cols", "axi_channels"):
            if getattr(self, name) < 1:
                raise ConfigError(f"{name} 必须 >=1，得到 {getattr(self, name)}")
        if self.axi_bits % wl.limb_bits != 0:
            raise ConfigError(
                f"axi_bits={self.axi_bits} 必须是 limb_bits={wl.limb_bits} 的整数倍"
            )
        if self.keymult_binding not in ("shared_transform_mul", "dedicated"):
            raise ConfigError(f"未知 keymult_binding: {self.keymult_binding}")
        if self.sram_capacity_bytes < 0:
            raise ConfigError("sram_capacity_bytes 不能为负")
        for name in ("ddr_bandwidth_GBps", "pcie_fixed_us", "pcie_bandwidth_GBps"):
            if getattr(self, name) < 0:
                raise ConfigError(f"{name} 不能为负")
        if self.strict_p4:
            self._assert_p4_shape(wl)

    def _assert_p4_shape(self, wl: WorkloadConfig) -> None:
        """--strict-p4：只允许 P4 已经过 RTL 精确比对的形状。"""
        expect = {
            "ring_dimension": 4096, "q_towers": 3, "p_towers": 2,
            "num_part_q": 2, "limb_bits": 64,
        }
        bad = {k: getattr(wl, k) for k, v in expect.items() if getattr(wl, k) != v}
        if bad:
            raise ConfigError(f"strict_p4 下 workload 必须为 {expect}，越界项：{bad}")
        if (self.transform_lanes, self.bconv_rows, self.bconv_cols) != (4, 3, 5):
            raise ConfigError("strict_p4 下 lanes/rows/cols 必须为 4/3/5")
        if self.allow_tiling:
            raise ConfigError("strict_p4 与 allow_tiling 互斥")

    def cycles_to_us(self, cycles: float) -> float:
        return cycles * self.clock_period_ns / 1000.0


@dataclass(frozen=True)
class SimConfig:
    """一次运行的完整配置。"""

    workload: WorkloadConfig = field(default_factory=WorkloadConfig)
    hardware: HardwareConfig = field(default_factory=HardwareConfig)
    # ready 队列的 tie-break 策略。方案 §5.2 未规定，但它会系统性偏袒某种
    # dataflow（LIFO 偏袒 DC 的深度优先，FIFO 偏袒 MP 的广度优先），
    # 因此必须固定并计入 hardware hash。trace_order = 严格按 trace 发射序，
    # 是唯一不引入额外自由度的选择。
    tie_break: str = "trace_order"
    # 方案 §3 把「kernel invocation boundary 和输入输出格式」列进了公平性条件，
    # 所以这两项与 tie_break 一样计入 hardware_config_hash，三种策略必须相同。
    boundary: str = "modup_only"        # modup_only | full_hks
    invocation: str = "per_digit"       # per_hks | per_digit | per_output_tower
    include_init: bool = False          # 是否把 OP_INIT 计入总时间
    note: str = ""

    def validate(self) -> None:
        self.workload.validate()
        self.hardware.validate(self.workload)
        if self.tie_break not in ("trace_order", "fifo", "lifo"):
            raise ConfigError(f"未知 tie_break: {self.tie_break}")
        if self.boundary not in ("modup_only", "full_hks"):
            raise ConfigError(f"未知 boundary: {self.boundary}")
        if self.invocation not in ("per_hks", "per_digit", "per_output_tower"):
            raise ConfigError(f"未知 invocation: {self.invocation}")

    def fairness_key(self) -> Dict[str, Any]:
        """构成 hardware_config_hash 的全部字段（方案 §3 的等同条件清单）。"""
        return {
            "hardware": asdict(self.hardware),
            "tie_break": self.tie_break,
            "boundary": self.boundary,
            "invocation": self.invocation,
            "include_init": self.include_init,
        }


def _canonical(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def config_hash(obj: Any, prefix: str = "") -> str:
    """稳定的配置指纹。用于证明三种策略跑在同一配置上。"""
    body = asdict(obj) if hasattr(obj, "__dataclass_fields__") else obj
    payload = _canonical({"schema": CONFIG_SCHEMA_VERSION, "prefix": prefix, "body": body})
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def hardware_config_hash(cfg: "SimConfig") -> str:
    """硬件 + 调度自由度 + 调用边界的指纹。

    三种策略必须给出同一个 hash，否则对比无效（方案 §3）。tie_break 之所以
    在里面，是因为它是调度器的自由度而非 dataflow 的性质；boundary/invocation
    在里面，是因为方案 §3 明确把调用边界列进了等同条件。
    """
    return config_hash(cfg.fairness_key(), "hw")


def workload_config_hash(wl: WorkloadConfig) -> str:
    return config_hash(wl, "wl")


def strategy_config_hash(strategy: str, tile_width: int = 1) -> str:
    """**调度**自由度的指纹，与硬件分开。

    output tile width 改变的是 dataflow 的调度粒度，不是硬件——同一套 3x5 阵列、
    同样的引擎和端口。因此它绝不能进 hardware_config_hash（那会让 OC-w1..w5
    之间的对比显示成「配置不同」而失去可比性），但必须单独记录，否则报告里就分不清
    某个数字是哪个 tile width 跑出来的。
    """
    return config_hash({"strategy": strategy, "oc_output_tile_width": tile_width},
                       "strat")


def _filter_known(cls, data: Dict[str, Any]) -> Dict[str, Any]:
    known = set(cls.__dataclass_fields__)
    unknown = set(data) - known
    if unknown:
        raise ConfigError(f"{cls.__name__} 不认识这些字段：{sorted(unknown)}")
    return {k: v for k, v in data.items() if k in known}


def sim_config_from_dict(data: Dict[str, Any]) -> SimConfig:
    """从 YAML/JSON 载入的 dict 构造 SimConfig。未知字段直接报错，不静默忽略。"""
    data = dict(data or {})
    wl = WorkloadConfig(**_filter_known(WorkloadConfig, data.pop("workload", {}) or {}))
    hw_raw = dict(data.pop("hardware", {}) or {})
    cost = CostParams(**_filter_known(CostParams, hw_raw.pop("cost", {}) or {}))
    hw = HardwareConfig(cost=cost, **_filter_known(HardwareConfig, hw_raw))
    cfg = SimConfig(workload=wl, hardware=hw, **_filter_known(SimConfig, data))
    cfg.validate()
    return cfg


def with_overrides(cfg: SimConfig, **kv: Any) -> SimConfig:
    """按 `workload.x` / `hardware.y` / `hardware.cost.z` 覆盖字段，供 --sweep 使用。"""
    wl_kv, hw_kv, cost_kv, top_kv = {}, {}, {}, {}
    for key, val in kv.items():
        if key.startswith("workload."):
            wl_kv[key[len("workload."):]] = val
        elif key.startswith("hardware.cost."):
            cost_kv[key[len("hardware.cost."):]] = val
        elif key.startswith("hardware."):
            hw_kv[key[len("hardware."):]] = val
        else:
            top_kv[key] = val
    wl = replace(cfg.workload, **_filter_known(WorkloadConfig, wl_kv)) if wl_kv else cfg.workload
    hw = cfg.hardware
    if cost_kv:
        hw = replace(hw, cost=replace(hw.cost, **_filter_known(CostParams, cost_kv)))
    if hw_kv:
        hw = replace(hw, **_filter_known(HardwareConfig, hw_kv))
    out = replace(cfg, workload=wl, hardware=hw, **_filter_known(SimConfig, top_kv))
    out.validate()
    return out
