import matplotlib.pyplot as plt
import numpy as np
import os

# ── 真实的端到端与 Kernel 耗时 (μs) ──────────────────────────────────────
strategies = ['DC', 'MP', 'OC']

# 真实的 Kernel 硬件执行总耗时 (直接取自日志)
kernel_exec = np.array([5537, 5555, 5397])

# 真实的 PCIe 传输耗时 (H2D + D2H)
pcie_transfer = np.array([1829 + 863, 1975 + 927, 1775 + 865]) # [2692, 2902, 2640]

# 期望的总端到端耗时应严格等于: DC=8229, MP=8457, OC=8037

# 原始提取的各子模块累加耗时 (纯计算时间)
raw_sub_ops = {
    'INTT':   np.array([789,  1228,   717]),
    'BConv':  np.array([938,   925,  2656]),
    'NTT':    np.array([2512, 2715,  2306]),
    'ModMul': np.array([231,   302,   310]),
}

# 核心修正：计算归一化比例系数
# 将各子模块的时间按比例映射到真实的 kernel_exec 时间内，以消除流水线重叠或stall的影响
raw_sums = sum(raw_sub_ops.values())
scale_factors = kernel_exec / raw_sums

# 应用缩放
data = {}
for op, vals in raw_sub_ops.items():
    data[op] = vals * scale_factors

data['PCIe Transfer'] = pcie_transfer

# 颜色映射
colors = {
    'INTT':          '#5C85D6',
    'BConv':         '#E07B39',
    'NTT':           '#4CAF50',
    'ModMul':        '#9C6BB5',
    'PCIe Transfer': '#A9A9A9',
}

# ── 绘图 ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6.5, 5.5))

x = np.arange(len(strategies))
bar_width = 0.45
bottoms = np.zeros(len(strategies))

bars = {}
for op, vals in data.items():
    b = ax.bar(x, vals, bar_width, bottom=bottoms,
               color=colors[op], edgecolor='white', linewidth=0.8, 
               label=op, zorder=3)
    bars[op] = (b, bottoms.copy(), vals)
    bottoms += vals

# 在每段中间标注数值
for op in ['INTT', 'BConv', 'NTT', 'ModMul', 'PCIe Transfer']:
    b, bot, vals = bars[op]
    for i, (rect, bv, v) in enumerate(zip(b, bot, vals)):
        # 为了美观，太窄的块不标数字，显示的数字取整
        if v > 400: 
            cy = bv + v / 2
            # PCIe 的灰色背景用深色字更清晰
            text_color = '#333' if op == 'PCIe Transfer' else 'white'
            ax.text(rect.get_x() + rect.get_width() / 2, cy,
                    f'{int(round(v))}', ha='center', va='center',
                    fontsize=8.5, color=text_color, fontweight='bold')

# 总延迟标注在柱顶 (强行对齐真实的总耗时)
totals = bottoms
for i, (xi, total) in enumerate(zip(x, totals)):
    ax.text(xi, total + 150, f'{int(round(total)):,} μs\n({total/1000:.2f} ms)',
            ha='center', va='bottom', fontsize=9, color='#333', fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(strategies, fontsize=12, fontweight='bold')
ax.set_ylabel('End-to-End Latency (μs)', fontsize=11, fontweight='bold')

ax.set_title('Fig. 5-3 End-to-End Latency Breakdown by Scheduling Strategy\n'
             r'($N=4096$, $\mathrm{sizeP}=2$, 2 digits, FPGA Offload)',
             fontsize=11, fontweight='bold', pad=15)

ax.set_ylim(-max(totals) * 0.15, max(totals) * 1.25)

# 图例移出图表避免遮挡
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9, framealpha=0.9, edgecolor='#ccc')

ax.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)
ax.set_axisbelow(True)

# 峰值缓冲副标注
peak_buf = ['Peak buf:\n64 KB', 'Peak buf:\n128 KB', 'Peak buf:\n32 KB']
for i, (xi, pb) in enumerate(zip(x, peak_buf)):
    ax.text(xi, -max(totals) * 0.03, pb,
            ha='center', va='top', fontsize=8, color='#444',
            bbox=dict(facecolor='#f9f9f9', edgecolor='#ccc', boxstyle='round,pad=0.4'))

plt.tight_layout()

out_dir = os.path.dirname(os.path.abspath(__file__))
for ext in ('pdf', 'png'):
    path = os.path.join(out_dir, f'fig5-3_latency_breakdown_fixed.{ext}')
    plt.savefig(path, bbox_inches='tight', dpi=300 if ext == 'png' else None)
    print(f'Saved: {path}')

plt.show()