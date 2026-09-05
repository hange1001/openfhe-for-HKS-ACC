"""Workload 展开：digit 布局、tower 身份与算子计数。

塔的全局编号约定（贯穿全仿真器）：
    0 .. L-1        Q towers
    L .. L+K-1      P towers（special primes）
digit d 的 native 塔 = [d*alpha, min((d+1)*alpha, L))
digit d 的 complement 塔 = 全部 L+K 个塔去掉它的 native 塔

由此 beta_d = L + K - alpha_d，与 doc/HKS_大参数扩展与单digit并行分析.md §2 一致。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

from .config import ConfigError, WorkloadConfig


@dataclass(frozen=True)
class Digit:
    """一个 digit 的塔布局。"""

    index: int
    native_towers: Tuple[int, ...]      # 该 digit 直接持有的 Q 塔（EVAL 域输入）
    complement_towers: Tuple[int, ...]  # 需要 BConv 补出来的塔

    @property
    def alpha(self) -> int:
        return len(self.native_towers)

    @property
    def beta(self) -> int:
        return len(self.complement_towers)


@dataclass(frozen=True)
class Workload:
    """展开后的 workload。所有 trace 生成器共享这一份。"""

    cfg: WorkloadConfig
    alpha: int
    digits: Tuple[Digit, ...]

    # ---- 基本形状 ----
    @property
    def N(self) -> int:
        return self.cfg.ring_dimension

    @property
    def L(self) -> int:
        return self.cfg.q_towers

    @property
    def K(self) -> int:
        return self.cfg.p_towers

    @property
    def D(self) -> int:
        return len(self.digits)

    @property
    def total_towers(self) -> int:
        """ModUp 输出基的塔数 = L + K。"""
        return self.L + self.K

    @property
    def bytes_per_tower(self) -> int:
        return self.cfg.bytes_per_tower

    def is_native(self, digit_index: int, tower: int) -> bool:
        return tower in self.digits[digit_index].native_towers

    # ---- 算子计数（不依赖任何周期参数，M1 用） ----
    def operation_counts(self) -> Dict[str, int]:
        """三种 dataflow 都必须产生这一组计数。

        注意：DC / MP / OC 的算术量在数学上恒等——
            sum_d beta_d = D*(L+K) - L = L*(D-1) + K*D
        方案 §6.3 给的 N_BConv_OC 与 DC 的 sum_d beta_d 是同一个数。
        所以 operation-count 报告本身没有区分度，区分在 lifetime 与调度上。
        """
        n_intt = sum(d.alpha for d in self.digits)          # = L
        n_scale = n_intt
        n_bconv = sum(d.beta for d in self.digits)          # = L*(D-1) + K*D
        n_ntt = n_bconv
        n_keymult = self.D * self.total_towers              # 每 digit 每输出塔一次
        n_accumulate = n_keymult
        return {
            "intt": n_intt,
            "scale": n_scale,
            "bconv": n_bconv,
            "ntt": n_ntt,
            "keymult": n_keymult,
            "accumulate": n_accumulate,
            "bypass": sum(d.alpha for d in self.digits),    # native 塔走旁路，不做 BConv/NTT
        }

    def oc_bconv_closed_form(self) -> int:
        """方案 §6.3 的闭式：N_BConv_OC = L*(D-1) + K*D。用于交叉校验。"""
        return self.L * (self.D - 1) + self.K * self.D

    def describe(self) -> str:
        parts = ", ".join(
            f"d{d.index}(alpha={d.alpha},beta={d.beta})" for d in self.digits
        )
        return f"N={self.N} L={self.L} K={self.K} D={self.D} alpha={self.alpha} [{parts}]"


def build_workload(cfg: WorkloadConfig) -> Workload:
    """从 WorkloadConfig 推导 digit 布局。非法参数直接拒绝。"""
    cfg.validate()
    L, K, D = cfg.q_towers, cfg.p_towers, cfg.num_part_q

    alpha = cfg.alpha if cfg.alpha is not None else -(-L // D)  # ceil(L/D)
    if alpha < 1:
        raise ConfigError(f"推导出的 alpha={alpha} 非法")

    # digit 数由 alpha 与 L 唯一确定；与 num_part_q 不符说明配置自相矛盾。
    derived_d = -(-L // alpha)
    if derived_d != D:
        raise ConfigError(
            f"alpha={alpha} 推出 ceil(L/alpha)={derived_d} 个 digit，"
            f"与 num_part_q={D} 不一致。请只指定其中一个。"
        )

    all_towers = tuple(range(L + K))
    digits: List[Digit] = []
    for d in range(D):
        lo = d * alpha
        alpha_d = min(alpha, L - lo)
        if alpha_d < 1:
            raise ConfigError(f"digit {d} 的 alpha_d={alpha_d} 非法（L={L}, alpha={alpha}）")
        native = tuple(range(lo, lo + alpha_d))
        complement = tuple(t for t in all_towers if t not in native)
        assert len(complement) == L + K - alpha_d, "beta_d 推导与定义不符"
        digits.append(Digit(index=d, native_towers=native, complement_towers=complement))

    covered = sorted(t for dg in digits for t in dg.native_towers)
    if covered != list(range(L)):
        raise ConfigError(f"digit 划分没有恰好覆盖 Q 塔 0..{L-1}，得到 {covered}")

    wl = Workload(cfg=cfg, alpha=alpha, digits=tuple(digits))

    # 闭式交叉校验：方案 §6.3 与逐 digit 展开必须一致。
    counts = wl.operation_counts()
    if counts["bconv"] != wl.oc_bconv_closed_form():
        raise ConfigError(
            f"BConv 计数不自洽：展开 {counts['bconv']} vs 闭式 {wl.oc_bconv_closed_form()}"
        )
    return wl
