"""资源池定义与绑定解析。

资源名是 trace 与调度器之间的契约。trace 只写资源名，容量与绑定由
HardwareConfig 决定，因此三种 dataflow 天然共享同一组资源。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

from .config import HardwareConfig

TRANSFORM = "transform_engine"
BCONV = "bconv_array"
KEYMULT = "keymult_engine"
#: **所有**搬运共用一个 AXI 池，容量 = axi_channels。
#:
#: 早先按用途拆成 dma_h2d / dma_d2h / dma_spill 三个各自容量 1 的池，等于白送
#: 3 倍并发带宽，而且按「用途」划分通道本身就是臆造的——真实绑定是按 gmem 端口。
#: 实测能证伪那个模型：P4 单 digit 事务 46671 只有在 alpha 个 load 塔加 L+K 个
#: store 塔**全部串行**（7x1053）时才对得上；任意两笔重叠，实测值都会更小。
#: 因此 digit 路径上的 AXI 是串行的，默认 axi_channels=1。
#:
#: （INIT 的两张 twiddle 表确实在 gmem0/gmem1 上重叠，但那是另一条代码路径，
#: 且 295063 是直接实测的常数，整体计入，不由本池建模。）
DMA = "dma_axi"
CONTROL = "control"

# 兼容旧名字：三者现在都指向同一个共享池
DMA_H2D = DMA
DMA_D2H = DMA
DMA_SPILL = DMA


@dataclass(frozen=True)
class ResourceSet:
    """名字 -> 容量。容量为并发单元数。"""

    capacity: Dict[str, int]
    #: KeyMult/Accumulate 实际落到哪个池上（共享变换模乘时就是 transform_engine）
    keymult_resource: str
    #: 变换 lane 数，供 cost_model 使用
    transform_lanes: int
    #: BConv 阵列的输出列数。一次调用最多能同时产出这么多输出塔，
    #: 因此它也是 OC output tile width 的物理上限。
    bconv_cols: int = 5

    def cap(self, name: str) -> int:
        if name not in self.capacity:
            raise KeyError(f"未知资源 {name}；已定义 {sorted(self.capacity)}")
        return self.capacity[name]


def build_resources(hw: HardwareConfig) -> ResourceSet:
    """按 HardwareConfig 构造资源池。

    keymult_binding="shared_transform_mul"（默认）时，KeyMult/Accumulate 与
    NTT/INTT/SCALE 争同一个变换引擎——这正是方案 §2「默认不增加额外 DSP」的含义。
    换成 "dedicated" 会新增一个独立引擎，属于另一个硬件配置点，hash 会变。
    """
    keymult_resource = TRANSFORM if hw.keymult_binding == "shared_transform_mul" else KEYMULT
    capacity = {
        TRANSFORM: hw.transform_engines,
        BCONV: hw.bconv_arrays,
        DMA: hw.axi_channels,
        CONTROL: 1,
    }
    if keymult_resource == KEYMULT:
        capacity[KEYMULT] = hw.transform_engines
    return ResourceSet(
        capacity=capacity,
        keymult_resource=keymult_resource,
        transform_lanes=hw.transform_lanes,
        bconv_cols=hw.bconv_cols,
    )
