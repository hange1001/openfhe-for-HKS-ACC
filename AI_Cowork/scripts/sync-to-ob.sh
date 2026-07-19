#!/bin/bash
# sync-to-ob.sh — 将项目文件夹的学习笔记同步到 Obsidian vault
# 用法: ./sync-to-ob.sh [--dry-run]

set -e

# ====== 配置 ======
PROJECT_DIR="/f/project/openfhe-for-HKS-ACC"
OB_VAULT="/f/Obsidian/Repos/IC_Research/IC_Research"
OB_TARGET_DIR="${OB_VAULT}/HKS-ACC"

# 需要同步的文件列表（项目路径 → OB 相对路径）
declare -A SYNC_FILES=(
    ["openfhe-for-HKS-ACC学习笔记.md"]="openfhe-for-HKS-ACC学习笔记.md"
    ["问题回答.md"]="问题回答.md"
)

# ====== 参数 ======
DRY_RUN=false
if [ "$1" = "--dry-run" ]; then
    DRY_RUN=true
    echo "[DRY RUN] 仅预览，不实际复制"
fi

# ====== 检查 ======
if [ ! -d "$OB_VAULT" ]; then
    echo "错误: OB vault 路径不存在: $OB_VAULT"
    echo "请确认 Obsidian vault 路径正确"
    exit 1
fi

if [ ! -d "$PROJECT_DIR" ]; then
    echo "错误: 项目路径不存在: $PROJECT_DIR"
    exit 1
fi

# ====== 同步 ======
mkdir -p "$OB_TARGET_DIR"

echo "项目: $PROJECT_DIR"
echo "OB:   $OB_TARGET_DIR"
echo "---"

SYNC_COUNT=0
SKIP_COUNT=0

for src in "${!SYNC_FILES[@]}"; do
    dst="${SYNC_FILES[$src]}"
    src_path="${PROJECT_DIR}/${src}"
    dst_path="${OB_TARGET_DIR}/${dst}"

    if [ ! -f "$src_path" ]; then
        echo "⚠ 跳过（源文件不存在）: $src"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    # 比较文件是否相同
    if [ -f "$dst_path" ] && cmp -s "$src_path" "$dst_path"; then
        echo "○ 未变化: $src"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    if [ "$DRY_RUN" = true ]; then
        echo "→ 将复制: $src → $dst"
        if [ -f "$dst_path" ]; then
            echo "  (覆盖已有文件)"
        fi
    else
        cp "$src_path" "$dst_path"
        if [ -f "$dst_path" ]; then
            echo "✓ 已更新: $src"
        else
            echo "✓ 已新增: $src"
        fi
    fi
    SYNC_COUNT=$((SYNC_COUNT + 1))
done

echo "---"
if [ "$DRY_RUN" = true ]; then
    echo "[DRY RUN] 将同步 $SYNC_COUNT 个文件，跳过 $SKIP_COUNT 个"
else
    echo "同步完成: 更新 $SYNC_COUNT 个文件，跳过 $SKIP_COUNT 个（未变化）"
fi

# ====== 提示 ======
echo ""
echo "下一步（在 OB 里手动完成）："
echo "  1. 将新模块章节拆成独立笔记（模块笔记/ or 概念笔记/）"
echo "  2. 加 [[]] 双向链接"
echo "  3. 更新 MOC 索引"
