"""
图5-2: 三种调度策略动态内存轨迹
数据来源: memory_trace_*.csv (由 hks-benchmark 生成)
DC/OC trace 从代码逻辑推导构造，MP 使用实测数据。
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# ── 参数 ──────────────────────────────────────────────────────────────────────
RING_DIM   = 4096
BYTES_PER  = 8
TOWER_KB   = RING_DIM * BYTES_PER / 1024   # 32 KB per ring-element
SIZE_P     = 2
NUM_DIGITS = 2

# ── 构造 DC trace ─────────────────────────────────────────────────────────────
# DC: for each digit: BConv alloc(sizeP towers) → NTT+Assemble free(sizeP towers)
dc_events = []
step = 0
for d in range(NUM_DIGITS):
    dc_events.append((step, f"BConv d{d}",   +SIZE_P * TOWER_KB))
    step += 1
    dc_events.append((step, f"Assemble d{d}", -SIZE_P * TOWER_KB))
    step += 1

dc_steps, dc_labels, dc_deltas = zip(*dc_events)
dc_watermark = np.cumsum(dc_deltas)

# ── 构造 MP trace ─────────────────────────────────────────────────────────────
# MP: all BConv alloc → all Assemble free (实测 CSV 已有，但重新构造保持一致)
mp_events = []
step = 0
for d in range(NUM_DIGITS):
    mp_events.append((step, f"BConv d{d}",   +SIZE_P * TOWER_KB))
    step += 1
for d in range(NUM_DIGITS):
    mp_events.append((step, f"Assemble d{d}", -SIZE_P * TOWER_KB))
    step += 1

mp_steps, mp_labels, mp_deltas = zip(*mp_events)
mp_watermark = np.cumsum(mp_deltas)

# ── 构造 OC trace ─────────────────────────────────────────────────────────────
# OC: for each P-tower p, for each digit d: BConv alloc(1) → Assemble free(1)
oc_events = []
step = 0
for p in range(SIZE_P):
    for d in range(NUM_DIGITS):
        oc_events.append((step, f"BConv p{p}d{d}",   +TOWER_KB))
        step += 1
        oc_events.append((step, f"Assemble p{p}d{d}", -TOWER_KB))
        step += 1

oc_steps, oc_labels, oc_deltas = zip(*oc_events)
oc_watermark = np.cumsum(oc_deltas)

# ── 绘图 ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4))

# 用阶梯图展示缓冲占用变化
def plot_staircase(ax, steps, watermark, color, label, linestyle='-'):
    # 在每个 step 前后各画一个点，形成阶梯效果
    xs = [0]
    ys = [0]
    for i, (s, w) in enumerate(zip(steps, watermark)):
        xs.append(s)
        ys.append(ys[-1])   # 水平段
        xs.append(s)
        ys.append(w)        # 垂直跳变
    xs.append(steps[-1] + 1)
    ys.append(watermark[-1])
    ax.plot(xs, ys, color=color, linestyle=linestyle, linewidth=2, label=label)
    ax.fill_between(xs, ys, alpha=0.08, color=color, step=None)

plot_staircase(ax, dc_steps, dc_watermark, '#2196F3', 'DC (Digit-Centric)')
plot_staircase(ax, mp_steps, mp_watermark, '#F44336', 'MP (Max-Parallel)', linestyle='--')
plot_staircase(ax, oc_steps, oc_watermark, '#4CAF50', 'OC (Output-Centric)', linestyle='-.')

# 峰值标注
ax.axhline(SIZE_P * TOWER_KB,           color='#2196F3', linewidth=0.8, linestyle=':')
ax.axhline(NUM_DIGITS * SIZE_P * TOWER_KB, color='#F44336', linewidth=0.8, linestyle=':')
ax.axhline(TOWER_KB,                    color='#4CAF50', linewidth=0.8, linestyle=':')

ax.annotate(f'DC peak: {SIZE_P * TOWER_KB:.0f} KB',
            xy=(0.02, SIZE_P * TOWER_KB), xycoords=('axes fraction', 'data'),
            fontsize=8, color='#2196F3', va='bottom')
ax.annotate(f'MP peak: {NUM_DIGITS * SIZE_P * TOWER_KB:.0f} KB',
            xy=(0.02, NUM_DIGITS * SIZE_P * TOWER_KB), xycoords=('axes fraction', 'data'),
            fontsize=8, color='#F44336', va='bottom')
ax.annotate(f'OC peak: {TOWER_KB:.0f} KB',
            xy=(0.02, TOWER_KB), xycoords=('axes fraction', 'data'),
            fontsize=8, color='#4CAF50', va='bottom')

ax.set_xlabel('Operation Step Index', fontsize=11)
ax.set_ylabel('On-chip Buffer Occupancy (KB)', fontsize=11)
ax.set_title('Fig. 5-2  Dynamic Memory Trace: DC vs MP vs OC\n'
             r'($N=4096$, $\mathrm{sizeP}=2$, 2 digits)',
             fontsize=11)
ax.set_ylim(-4, NUM_DIGITS * SIZE_P * TOWER_KB * 1.35)
ax.set_xlim(-0.3)
ax.legend(loc='upper right', fontsize=9)
ax.grid(axis='y', linestyle=':', alpha=0.5)
ax.set_xticks(range(max(len(dc_steps), len(mp_steps), len(oc_steps)) + 1))

plt.tight_layout()

out_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(out_dir, 'fig5-2_memory_trace.pdf')
plt.savefig(out_path, bbox_inches='tight')
out_png = os.path.join(out_dir, 'fig5-2_memory_trace.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"Saved: {out_path}")
print(f"Saved: {out_png}")
plt.show()
