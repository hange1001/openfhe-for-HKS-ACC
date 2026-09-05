"""算子周期模型。

标定基线 P4 / e948a69，证据见 docs/reports/hls/hks_mem_p4_20260905/。
每个函数都标注它是 measured 复现、calibrated 还是 projected。

零残差校验（tests/test_cost_model.py 会跑）：
    alpha=2 digit: 5*6481 + 2*1047 + 4145 + 7*1053 + 656 = 46671   实测 46671
    alpha=1 digit: 5*6481 + 1*1047 + 4145 + 6*1053 + 656 = 44571   实测 44571
hold-out（P3，预乘是独立硬件、1 系数/拍）：
    沿用同一组 axi/control 参数，两 digit 预测 100431 vs 实测 100422，+0.009%
"""

from __future__ import annotations

from dataclasses import dataclass

from .config import CostParams, HardwareConfig, WorkloadConfig


def _ceil_div(a: int, b: int) -> int:
    return -(-a // b)


@dataclass(frozen=True)
class CostModel:
    """把 workload/hardware 形状翻译成周期数。"""

    wl: WorkloadConfig
    hw: HardwareConfig

    @property
    def c(self) -> CostParams:
        return self.hw.cost

    # ---------- 变换引擎 ----------
    def butterfly_cycles(self) -> int:
        """单级蝶形。measured @P4：2048/4 + 26 = 538。"""
        half_n = self.wl.ring_dimension // 2
        return _ceil_div(half_n, self.hw.transform_lanes) + self.c.butterfly_overhead

    def transform_cycles(self) -> int:
        """一次 NTT 或 INTT（单塔）。measured @P4：6481。

        WORK_STAGE_LOOP 每次迭代处理两级：2*538 + 4 = 1080，共 6 次 = 6480，
        加一拍调用开销 = 6481。写成 stages*butterfly + (stages//2)*pair_overhead
        以便推广到 log2(N) 为奇数的参数点（此时为 projected）。
        """
        stages = self.wl.stages
        return (stages * self.butterfly_cycles()
                + (stages // 2) * self.c.stage_pair_overhead
                + self.c.transform_call_overhead)

    def scale_cycles(self) -> int:
        """QHatInv 预乘（单塔）。measured @P4：4096/4 + 22 + 1 = 1047。

        P4 之后 SCALE 复用变换引擎的 4 路模乘，因此与 NTT/INTT 互斥。
        P3 时它是独立硬件、1 系数/拍（hold-out 解出约 4110 拍/塔）。
        """
        return (_ceil_div(self.wl.ring_dimension, self.hw.transform_lanes)
                + self.c.scale_loop_overhead + self.c.scale_call_overhead)

    # ---------- BConv ----------
    def bconv_cycles(self, alpha: int, beta: int) -> int:
        """一次 BConv 调用。measured @P4（alpha=2,beta=3 与 alpha=1,beta=4 均为 4145）。

        T = N * ceil(alpha/rows) * ceil(beta/cols) + overhead

        注意 ceil(beta/cols)：固定 rows x cols 阵列**并行产出全部 cols 列**，
        所以在 beta<=cols 时，产出 1 个输出塔和产出 5 个输出塔耗时相同。
        这正是 OC 的 single-tower BConv 在当前 3x5 固定阵列上不省时间的原因，
        与 ADR-010 / OC_strategy_gap_analysis 描述的「sizeP 倍冗余」是同一件事。
        要让 OC 真正省下 BConv 时间，必须缩窄阵列或做列分块（allow_tiling）。
        """
        if alpha < 1 or beta < 1:
            raise ValueError(f"bconv_cycles 需要 alpha,beta >= 1，得到 {alpha},{beta}")
        row_tiles = _ceil_div(alpha, self.hw.bconv_rows)
        col_tiles = _ceil_div(beta, self.hw.bconv_cols)
        return self.wl.ring_dimension * row_tiles * col_tiles + self.c.bconv_overhead

    # ---------- AXI / 搬运 ----------
    def axi_beats_per_tower(self) -> int:
        """单塔 AXI 拍数。N=4096/64-bit/256-bit 总线 -> 1024。"""
        bits = self.wl.ring_dimension * self.wl.limb_bits
        return _ceil_div(bits, self.hw.axi_bits)

    def axi_tower_cycles(self) -> int:
        """单塔 AXI 传输周期。calibrated @P4：1024 + 29 = 1053。"""
        return self.axi_beats_per_tower() + self.c.axi_tower_extra

    def transfer_cycles(self, towers: int) -> int:
        return towers * self.axi_tower_cycles()

    def spill_cycles(self, nbytes: int) -> int:
        """片外 spill。ddr_bandwidth 未给时退化为 AXI 口径（projected）。"""
        if self.hw.ddr_bandwidth_GBps > 0:
            secs = nbytes / (self.hw.ddr_bandwidth_GBps * 1e9)
            return max(1, int(secs * 1e9 / self.hw.clock_period_ns))
        towers = _ceil_div(nbytes, max(1, self.wl.bytes_per_tower))
        return self.transfer_cycles(towers)

    # ---------- 控制 ----------
    def invocation_cycles(self) -> int:
        """一次 kernel 调用的固定控制开销。calibrated @P4：656。

        标定时 1 个 digit = 1 次调用。OC 若按输出塔逐个发起调用，这一项会乘上
        调用次数——方案 §12 点名的「固定开销抵消 OC 收益」在模型里就体现于此。
        """
        return self.c.control_per_digit

    def init_cycles(self) -> int:
        """OP_INIT。measured @P4：295063（上下文级固定段，不随 digit 缩放）。"""
        return self.c.init_cycles

    def pcie_fixed_cycles(self) -> int:
        """PCIe 固定开销。无板卡，只能作为扫描维，默认 0（不建模）。"""
        if self.hw.pcie_fixed_us <= 0:
            return 0
        return int(self.hw.pcie_fixed_us * 1000.0 / self.hw.clock_period_ns)

    # ---------- KeyMult / Accumulate（projected：硬件未实现） ----------
    def keymult_cycles(self) -> int:
        """T_KM = passes * ceil(N/lanes) * II + pipe + mem。

        passes=2 对应 c_j*b_j 与 c_j*a_j。lanes 取变换引擎的模乘 lane 数
        （默认绑定 shared_transform_mul，即不新增 DSP）。
        """
        per_pass = _ceil_div(self.wl.ring_dimension, self.hw.transform_lanes)
        return (self.c.keymult_mul_passes * per_pass * self.c.modmul_ii
                + self.c.keymult_pipe + self.c.keymult_mem)

    def accumulate_cycles(self) -> int:
        """T_ACC = ceil(N/lanes) + overhead。与 keymult 之和等价于方案 §5.2 的 T_KM。"""
        return (_ceil_div(self.wl.ring_dimension, self.hw.transform_lanes)
                + self.c.accumulate_overhead)

    # ---------- 自检 ----------
    def p4_reference(self) -> dict:
        """P4 标定点的复现值，供测试与报告引用。"""
        return {
            "transform": self.transform_cycles(),
            "scale": self.scale_cycles(),
            "bconv_a2b3": self.bconv_cycles(2, 3),
            "bconv_a1b4": self.bconv_cycles(1, 4),
            "axi_tower": self.axi_tower_cycles(),
            "invocation": self.invocation_cycles(),
            "init": self.init_cycles(),
        }
