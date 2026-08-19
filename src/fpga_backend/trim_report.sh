#!/usr/bin/env sh
# =============================================================================
# trim_report.sh — 把 Vitis HLS 的 csynth.rpt 裁到可入库的体积，
#                  并把顶层那一行数字抽进 docs/reports/summary.csv
#
# 用法（由 Makefile 的 report target 调用，一般不用手敲）：
#   ./trim_report.sh <src csynth.rpt> <dst .rpt> <summary.csv> <module>
#
# 为什么要裁：
#   完整 csynth.rpt 里 "Bind Op Report" + "SW I/O Information" 占 ~70% 篇幅
#   （Compute_NTT: 614/866 行），内容是逐个算子的 DSP 绑定明细。
#   跨 commit diff 时这两段全是噪声，会把真正要看的资源/时序变化淹掉。
#   保留：Synthesis Summary（资源+slack）/ HW Interfaces（AXI 位宽）/
#         Storage Report（BRAM-URAM 绑定）/ Pragma Report（抓被忽略的 pragma）。
#
# 要改保留哪几段：设环境变量 HLS_RPT_DROP="段名1|段名2"
# =============================================================================
set -e

SRC="$1"
DST="$2"
CSV="$3"
MODULE="$4"

if [ -z "$SRC" ] || [ -z "$DST" ] || [ -z "$CSV" ] || [ -z "$MODULE" ]; then
    echo "usage: $0 <src.rpt> <dst.rpt> <summary.csv> <module>" >&2
    exit 2
fi
if [ ! -f "$SRC" ]; then
    echo "[trim_report] source not found: $SRC" >&2
    exit 1
fi

DROP="${HLS_RPT_DROP:-SW I/O Information|Bind Op Report}"

# --- 1. 裁剪：按 "=== / == 段名 / ===" 三行头切段，丢掉 DROP 里的段 ----------
awk -v drop="$DROP" '
    { line[NR] = $0 }
    END {
        n = NR
        # 定位所有段头（形如 ====== / == Name / ======）
        ns = 0
        for (i = 2; i < n; i++) {
            if (line[i] ~ /^== / && line[i-1] ~ /^=+$/ && line[i+1] ~ /^=+$/) {
                ns++
                start[ns] = i - 1
                name[ns]  = substr(line[i], 4)
                sub(/[ \t]+$/, "", name[ns])
            }
        }
        # 段头之前的引言（Date/Version/Target device 那块）永远保留
        head_end = (ns > 0) ? start[1] - 1 : n
        for (i = 1; i <= head_end; i++) print line[i]

        for (s = 1; s <= ns; s++) {
            stop = (s < ns) ? start[s+1] - 1 : n
            keep = 1
            m = split(drop, d, "|")
            for (k = 1; k <= m; k++)
                if (name[s] == d[k]) keep = 0
            if (keep) {
                for (i = start[s]; i <= stop; i++) print line[i]
            } else {
                printf "%s\n== %s  [trimmed by trim_report.sh]\n%s\n\n", \
                       line[start[s]], name[s], line[start[s]]
            }
        }
    }
' "$SRC" > "$DST"

# --- 2. 抽顶层那一行进 summary.csv ------------------------------------------
# 顶层模块行以 "|+ " 开头；子模块是 "|  + "、循环是 "| o "，正则天然区分。
ROW=$(awk -v mod="$MODULE" '
    /^[ \t]*\* Date:/           { d = $0; sub(/^[^:]*:[ \t]*/, "", d) }
    /^[ \t]*\* Version:/        { v = $0; sub(/^[^:]*:[ \t]*/, "", v); sub(/ .*/, "", v) }
    /^[ \t]*\* Target device:/  { p = $0; sub(/^[^:]*:[ \t]*/, "", p) }
    /^[ \t]*\|\+ / && !done {
        split($0, f, "|")
        for (i = 1; i <= 15; i++) { gsub(/^[ \t]+|[ \t]+$/, "", f[i]) }
        # f[2]="+ Name" f[3]=Issue f[4]=Slack f[5]=Lat(cyc) f[6]=Lat(ns)
        # f[11]=BRAM f[12]=DSP f[13]=FF f[14]=LUT f[15]=URAM
        issue = f[3]; slack = f[4]; cyc = f[5]; ns = f[6]
        # 只留绝对值，丢掉 "(4%)"——百分比在 .rpt 里看，CSV 要能直接画图
        for (i = 11; i <= 15; i++) { sub(/[ \t]*\(.*/, "", f[i]); if (f[i] == "") f[i] = "-" }
        # 时钟周期由 latency_ns / latency_cycles 反推：能抓到 5ns/6ns 口径漂移
        clk = "-"
        if (cyc ~ /^[0-9]+$/ && cyc + 0 > 0 && ns ~ /^[0-9.eE+]+$/)
            clk = sprintf("%.2f", (ns + 0) / (cyc + 0))
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", \
               d, mod, p, v, clk, issue, slack, cyc, f[11], f[12], f[13], f[14], f[15]
        done = 1
    }
' "$SRC")

HDR="synth_date,module,part,vitis,clk_ns,issue,slack_ns,latency_cycles,BRAM,DSP,FF,LUT,URAM"

if [ -n "$ROW" ]; then
    mkdir -p "$(dirname "$CSV")"
    [ -f "$CSV" ] || echo "$HDR" > "$CSV"
    # upsert：同名 module 的旧行删掉再追加，保持一模块一行
    TMP="$CSV.tmp$$"
    awk -F, -v mod="$MODULE" 'NR==1 || $2 != mod' "$CSV" > "$TMP"
    echo "$ROW" >> "$TMP"
    # 表头之外按模块名排序，diff 才稳定
    { head -1 "$TMP"; tail -n +2 "$TMP" | sort -t, -k2,2; } > "$CSV"
    rm -f "$TMP"
fi

SRC_LINES=$(wc -l < "$SRC")
DST_LINES=$(wc -l < "$DST")
echo "  [+] $(basename "$DST")  ($SRC_LINES -> $DST_LINES 行)"
[ -n "$ROW" ] && echo "  [+] $(basename "$CSV")  <- $MODULE"
exit 0
