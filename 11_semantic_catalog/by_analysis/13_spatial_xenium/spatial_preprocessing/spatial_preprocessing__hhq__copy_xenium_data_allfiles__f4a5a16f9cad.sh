#!/bin/bash
#===============================================================================
# copy_xenium_data.sh
#
# 功能:
#   从指定的 Xenium 项目路径中查找 Xenium 输出目录，并把输出目录下的
#   全部文件和子文件夹复制到目标目录。
#
# 用法:
#   bash /mnt/NexData_file/home/hhq/func/copy_xenium_data_allfiles.sh /data/xenium/XE2508254414-02
#   bash copy_xenium_data.sh /data/xenium/XE2508254414-02 --dry-run
#
# 目标路径:
#   /mnt/NexData_file/home/hhq/project/HE_FM_data/<项目号>/<片子号>/
#   如果输入路径本身就是 Xenium 输出目录，则复制到:
#   /mnt/NexData_file/home/hhq/project/HE_FM_data/<项目号>/
#===============================================================================

set -euo pipefail

# ===== 配置 =====
DEST_BASE="/mnt/NexData_file/home/hhq/project/HE_FM_data"

# 这些文件只用于识别 Xenium 输出目录，不再作为拷贝白名单。
MARKER_FILES=(
    "cell_feature_matrix.h5"
    "cells.parquet"
    "transcripts.parquet"
)

# 项目目录下常见的辅助目录，不作为片子目录处理。
SKIP_DIRS=(
    "aux_outputs"
    "analysis"
)

# ===== 颜色输出 =====
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

usage() {
    echo "用法: $0 <xenium_project_path> [--dry-run]"
    echo "示例: $0 /data/xenium/XE2508254414-02"
}

contains_marker_file() {
    local dir="$1"
    local marker

    for marker in "${MARKER_FILES[@]}"; do
        if [ -f "${dir%/}/${marker}" ]; then
            return 0
        fi
    done

    return 1
}

find_xenium_output_dir() {
    local dir="$1"

    find "$dir" -maxdepth 2 -type f \( \
        -name "cell_feature_matrix.h5" -o \
        -name "cells.parquet" -o \
        -name "transcripts.parquet" \
    \) -exec dirname {} \; -quit 2>/dev/null
}

is_skip_dir() {
    local name="$1"
    local skip

    for skip in "${SKIP_DIRS[@]}"; do
        if [ "$name" = "$skip" ]; then
            return 0
        fi
    done

    return 1
}

copy_all_contents() {
    local src_dir="${1%/}"
    local dest_dir="${2%/}"
    local label="$3"
    local dir_size
    local file_count

    dir_size=$(du -sh "$src_dir" 2>/dev/null | cut -f1 || echo "未知")
    file_count=$(find "$src_dir" -type f 2>/dev/null | wc -l | tr -d ' ')

    echo -e "${GREEN}------------------------------------------${NC}"
    echo -e "${GREEN}处理: ${label}${NC}"
    echo -e "${BLUE}  源路径: ${src_dir}${NC}"
    echo -e "${BLUE}  目标:   ${dest_dir}${NC}"
    echo -e "${BLUE}  内容:   ${file_count} 个文件, ${dir_size}${NC}"

    if [ "$DRY_RUN" = true ]; then
        echo -e "  ${YELLOW}[DRY-RUN] 将复制该目录下的全部内容${NC}"
        return 0
    fi

    mkdir -p "$dest_dir"
    cp -a "${src_dir}/." "${dest_dir}/"
    echo -e "  ${GREEN}完成: 已复制全部内容${NC}"
}

# ===== 参数解析 =====
if [ $# -lt 1 ]; then
    echo -e "${RED}错误: 请提供 Xenium 项目路径${NC}"
    usage
    exit 1
fi

SOURCE_DIR="${1%/}"
DRY_RUN=false

if [ "${2:-}" = "--dry-run" ]; then
    DRY_RUN=true
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW}  DRY-RUN 模式: 仅预览，不实际复制${NC}"
    echo -e "${YELLOW}========================================${NC}"
    echo ""
elif [ $# -gt 1 ]; then
    echo -e "${RED}错误: 未知参数 ${2}${NC}"
    usage
    exit 1
fi

if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}错误: 项目路径不存在: ${SOURCE_DIR}${NC}"
    exit 1
fi

PROJECT_ID=$(basename "$SOURCE_DIR")

echo -e "${CYAN}===============================================${NC}"
echo -e "${CYAN}  Xenium 数据拷贝脚本${NC}"
echo -e "${CYAN}===============================================${NC}"
echo -e "${BLUE}项目路径:${NC}     ${SOURCE_DIR}"
echo -e "${BLUE}项目号:${NC}       ${PROJECT_ID}"
echo -e "${BLUE}目标根路径:${NC}   ${DEST_BASE}"
echo -e "${BLUE}拷贝策略:${NC}     复制 Xenium 输出目录下的全部内容"
echo ""

# ===== 统计变量 =====
TOTAL_SLIDES=0
SUCCESS_COUNT=0
SKIP_COUNT=0

# 如果输入目录本身就是 Xenium 输出目录，直接复制它的全部内容。
if contains_marker_file "$SOURCE_DIR"; then
    echo -e "${YELLOW}检测到输入路径本身就是 Xenium 输出目录，将复制其全部内容。${NC}"
    copy_all_contents "$SOURCE_DIR" "${DEST_BASE}/${PROJECT_ID}" "$PROJECT_ID"
    SUCCESS_COUNT=1
else
    # 否则，把第一层子目录当作片子目录逐个处理。
    SLIDE_FOUND=false

    for SLIDE_DIR in "${SOURCE_DIR}"/*/; do
        if [ ! -d "$SLIDE_DIR" ]; then
            continue
        fi

        SLIDE_FOUND=true
        SLIDE_ID=$(basename "$SLIDE_DIR")
        TOTAL_SLIDES=$((TOTAL_SLIDES + 1))

        if is_skip_dir "$SLIDE_ID"; then
            echo -e "${YELLOW}跳过辅助目录: ${SLIDE_ID}${NC}"
            SKIP_COUNT=$((SKIP_COUNT + 1))
            continue
        fi

        COPY_SOURCE="${SLIDE_DIR%/}"
        if ! contains_marker_file "$COPY_SOURCE"; then
            XENIUM_SUBDIR=$(find_xenium_output_dir "$COPY_SOURCE")
            if [ -n "$XENIUM_SUBDIR" ]; then
                echo -e "${YELLOW}在子目录中找到 Xenium 输出: ${XENIUM_SUBDIR}${NC}"
                COPY_SOURCE="$XENIUM_SUBDIR"
            else
                echo -e "${YELLOW}警告: ${SLIDE_ID} 下未找到 Xenium 输出标识文件，跳过${NC}"
                SKIP_COUNT=$((SKIP_COUNT + 1))
                continue
            fi
        fi

        copy_all_contents "$COPY_SOURCE" "${DEST_BASE}/${PROJECT_ID}/${SLIDE_ID}" "$SLIDE_ID"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    done

    if [ "$SLIDE_FOUND" = false ]; then
        echo -e "${RED}错误: 未在 ${SOURCE_DIR} 下找到片子子目录，也未检测到根目录是 Xenium 输出目录${NC}"
        exit 1
    fi
fi

if [ "$SUCCESS_COUNT" -eq 0 ]; then
    echo -e "${RED}错误: 未复制任何 Xenium 输出目录${NC}"
    exit 1
fi

# ===== 汇总报告 =====
echo ""
echo -e "${CYAN}===============================================${NC}"
echo -e "${CYAN}  拷贝完成 - 汇总报告${NC}"
echo -e "${CYAN}===============================================${NC}"
echo -e "${BLUE}项目号:${NC}       ${PROJECT_ID}"
echo -e "${BLUE}成功处理数:${NC}   ${SUCCESS_COUNT}"
echo -e "${BLUE}跳过目录数:${NC}   ${SKIP_COUNT}"
echo -e "${GREEN}目标路径: ${DEST_BASE}/${PROJECT_ID}/${NC}"

if [ "$DRY_RUN" = true ]; then
    echo ""
    echo -e "${YELLOW}这是 DRY-RUN 模式，没有文件被实际复制。${NC}"
    echo -e "${YELLOW}去掉 --dry-run 参数即可执行实际复制。${NC}"
fi
