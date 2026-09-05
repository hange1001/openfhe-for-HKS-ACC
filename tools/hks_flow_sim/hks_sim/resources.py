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
DMA_H2D = "dma_h2d"
DMA_D2H = "dma_d2h"
DMA_SPILL = "dma_spill"
CONTROL = "control"


@dataclass(frozen=True)
class ResourceSet:
    """名字 -> 容量。容量为并发单元数。"""

    capacity: Dict[str, int]
    #: KeyMult/Accumulate 实际落到哪个池上（共享变换模乘时就是 transform_engine）
    keymult_resource: str
    #: 变换 lane 数，供 cost_model 使用
    transform_lanes: int

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
        DMA_H2D: hw.axi_channels,
        DMA_D2H: hw.axi_channels,
        DMA_SPILL: hw.axi_channels,
        CONTROL: 1,
    }
    if keymult_resource == KEYMULT:
        capacity[KEYMULT] = hw.transform_engines
    return ResourceSet(
        capacity=capacity,
        keymult_resource=keymult_resource,
        transform_lanes=hw.transform_lanes,
    )
