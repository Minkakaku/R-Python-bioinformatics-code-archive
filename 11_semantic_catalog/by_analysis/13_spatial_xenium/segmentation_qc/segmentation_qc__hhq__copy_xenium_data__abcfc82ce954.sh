#!/bin/bash
#===============================================================================
# copy_xenium_data.sh
#
# 功能: 从指定的 Xenium 项目路径中遍历所有片子（子文件夹），
#       将关键分析文件拷贝到目标目录，按 "项目号/片子号" 的结构组织。
#
# 用法:
#   bash /mnt/NexData_file/home/hhq/func/copy_xenium_data.sh /data/xenium/XE2508254414-02
#   bash copy_xenium_data.sh /data/xenium/XE2508254414-02 --dry-run   # 仅预览，不实际拷贝
#
# 需要拷贝的文件:
#   1. cell_feature_matrix.h5       — cell × gene 原始计数矩阵
#   2. cells.parquet                — 每个 cell 的坐标与 metadata
#   3. cell_boundaries.parquet      — 细胞分割边界 (boundary polygon)
#   4. nucleus_boundaries.parquet   — 细胞核分割边界
#   5. morphology_focus/            — DAPI / morphology 聚焦图像（文件夹）
#   6. transcripts.parquet          — transcript-level 分子坐标表
#
# 目标路径:
#   /mnt/NexData_file/home/hhq/project/HE_FM_data/<项目号>/<片子号>/
#===============================================================================

# 需要拷贝的目录
# /data/xenium/XE2508254414-02（20）
# /data/xenium/20250402__040500__FWSP20251137_1136 （12）
# /data/xenium/FWSP20251137 （3）
# /data/xenium/FWSP20251137-03 （16）
# /data/xenium/FWSP20251137-05 （31）
# /data/xenium/XE2411273372-01 （5）


set -euo pipefail

# ===== 配置 =====
DEST_BASE="/mnt/NexData_file/home/hhq/project/HE_FM_data"

# 需要拷贝的文件列表（相对于每个片子目录）
FILES_TO_COPY=(
    "cell_feature_matrix.h5"
    "cells.parquet"
    "cell_boundaries.parquet"
    "nucleus_boundaries.parquet"
    "transcripts.parquet"
)

# 需要拷贝的文件夹
DIRS_TO_COPY=(
    "morphology_focus"
)

# ===== 颜色输出 =====
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# ===== 参数解析 =====
if [ $# -lt 1 ]; then
    echo -e "${RED}错误: 请提供 Xenium 项目路径${NC}"
    echo "用法: $0 <xenium_project_path> [--dry-run]"
    echo "示例: $0 /data/xenium/XE2508254414-02"
    exit 1
fi

SOURCE_DIR="$1"
DRY_RUN=false

if [ "${2:-}" = "--dry-run" ]; then
    DRY_RUN=true
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW}  DRY-RUN 模式: 仅预览，不实际拷贝${NC}"
    echo -e "${YELLOW}========================================${NC}"
    echo ""
fi

# 检查源路径是否存在
if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}错误: 项目路径不存在: ${SOURCE_DIR}${NC}"
    exit 1
fi

# 提取项目号（取路径最后一段）
PROJECT_ID=$(basename "$SOURCE_DIR")

echo -e "${CYAN}===============================================${NC}"
echo -e "${CYAN}  Xenium 数据拷贝脚本${NC}"
echo -e "${CYAN}===============================================${NC}"
echo -e "${BLUE}项目路径:${NC} ${SOURCE_DIR}"
echo -e "${BLUE}项目号:${NC}   ${PROJECT_ID}"
echo -e "${BLUE}目标根路径:${NC} ${DEST_BASE}"
echo ""

# ===== 统计变量 =====
TOTAL_SLIDES=0
SUCCESS_COUNT=0
SKIP_COUNT=0
MISSING_FILES=()

# ===== 遍历项目下的子文件夹（每个子文件夹是一个片子） =====
# 这里只找第一层的子目录
SLIDE_FOUND=false

for SLIDE_DIR in "${SOURCE_DIR}"/*/; do
    # 检查是否有子目录
    if [ ! -d "$SLIDE_DIR" ]; then
        continue
    fi

    SLIDE_FOUND=true
    SLIDE_ID=$(basename "$SLIDE_DIR")
    TOTAL_SLIDES=$((TOTAL_SLIDES + 1))

    # 跳过 aux_outputs 等辅助目录
    if [[ "$SLIDE_ID" == "aux_outputs" || "$SLIDE_ID" == "analysis" ]]; then
        echo -e "${YELLOW}跳过辅助目录: ${SLIDE_ID}${NC}"
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi

    DEST_DIR="${DEST_BASE}/${PROJECT_ID}/${SLIDE_ID}"

    echo -e "${GREEN}------------------------------------------${NC}"
    echo -e "${GREEN}处理片子: ${SLIDE_ID}${NC}"
    echo -e "${BLUE}  源路径: ${SLIDE_DIR}${NC}"
    echo -e "${BLUE}  目标:   ${DEST_DIR}${NC}"

    # 检查该子目录是否包含 Xenium 输出文件
    # 通过检查是否至少存在一个关键文件来判断
    HAS_XENIUM_FILES=false
    for f in "${FILES_TO_COPY[@]}"; do
        if [ -f "${SLIDE_DIR}${f}" ]; then
            HAS_XENIUM_FILES=true
            break
        fi
    done

    if [ "$HAS_XENIUM_FILES" = false ]; then
        # 尝试递归查找（可能 Xenium 输出在子目录的 output-xxx 文件夹中）
        XENIUM_SUBDIR=$(find "$SLIDE_DIR" -maxdepth 2 -name "cell_feature_matrix.h5" -printf "%h\n" 2>/dev/null | head -1)
        if [ -n "$XENIUM_SUBDIR" ]; then
            echo -e "${YELLOW}  在子目录中找到 Xenium 输出: ${XENIUM_SUBDIR}${NC}"
            SLIDE_DIR="${XENIUM_SUBDIR}/"
        else
            echo -e "${YELLOW}  警告: 该目录下未找到 Xenium 输出文件，跳过${NC}"
            SKIP_COUNT=$((SKIP_COUNT + 1))
            continue
        fi
    fi

    # 创建目标目录
    if [ "$DRY_RUN" = false ]; then
        mkdir -p "$DEST_DIR"
    else
        echo -e "  ${YELLOW}[DRY-RUN] 将创建目录: ${DEST_DIR}${NC}"
    fi

    # 拷贝文件
    for FILE in "${FILES_TO_COPY[@]}"; do
        SRC_FILE="${SLIDE_DIR}${FILE}"
        if [ -f "$SRC_FILE" ]; then
            FILE_SIZE=$(du -sh "$SRC_FILE" 2>/dev/null | cut -f1)
            if [ "$DRY_RUN" = false ]; then
                echo -e "  ${GREEN}✓ 拷贝 ${FILE} (${FILE_SIZE})${NC}"
                cp "$SRC_FILE" "$DEST_DIR/"
            else
                echo -e "  ${GREEN}✓ [DRY-RUN] 将拷贝 ${FILE} (${FILE_SIZE})${NC}"
            fi
        else
            echo -e "  ${RED}✗ 缺失: ${FILE}${NC}"
            MISSING_FILES+=("${SLIDE_ID}/${FILE}")
        fi
    done

    # 拷贝文件夹（morphology_focus）
    for DIR in "${DIRS_TO_COPY[@]}"; do
        SRC_SUBDIR="${SLIDE_DIR}${DIR}"
        if [ -d "$SRC_SUBDIR" ]; then
            DIR_SIZE=$(du -sh "$SRC_SUBDIR" 2>/dev/null | cut -f1)
            FILE_COUNT=$(find "$SRC_SUBDIR" -type f | wc -l)
            if [ "$DRY_RUN" = false ]; then
                echo -e "  ${GREEN}✓ 拷贝 ${DIR}/ (${DIR_SIZE}, ${FILE_COUNT} 个文件)${NC}"
                cp -r "$SRC_SUBDIR" "$DEST_DIR/"
            else
                echo -e "  ${GREEN}✓ [DRY-RUN] 将拷贝 ${DIR}/ (${DIR_SIZE}, ${FILE_COUNT} 个文件)${NC}"
            fi
        else
            # morphology_focus 也可能是单个 ome.tif 文件（某些旧版本）
            ALT_FILE="${SLIDE_DIR}morphology_focus.ome.tif"
            if [ -f "$ALT_FILE" ]; then
                FILE_SIZE=$(du -sh "$ALT_FILE" 2>/dev/null | cut -f1)
                if [ "$DRY_RUN" = false ]; then
                    echo -e "  ${GREEN}✓ 拷贝 morphology_focus.ome.tif (${FILE_SIZE}) [替代 morphology_focus 文件夹]${NC}"
                    mkdir -p "$DEST_DIR/morphology_focus"
                    cp "$ALT_FILE" "$DEST_DIR/morphology_focus/"
                else
                    echo -e "  ${GREEN}✓ [DRY-RUN] 将拷贝 morphology_focus.ome.tif (${FILE_SIZE})${NC}"
                fi
            else
                echo -e "  ${RED}✗ 缺失: ${DIR}/ (文件夹和 .ome.tif 均未找到)${NC}"
                MISSING_FILES+=("${SLIDE_ID}/${DIR}/")
            fi
        fi
    done

    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
done

# 如果源目录本身就是一个 Xenium 输出目录（没有子文件夹包含 Xenium 数据）
if [ "$SLIDE_FOUND" = false ] || ( [ "$SUCCESS_COUNT" -eq 0 ] && [ "$SKIP_COUNT" -eq "$TOTAL_SLIDES" ] ); then
    # 检查源目录本身是否包含 Xenium 文件
    HAS_ROOT_FILES=false
    for f in "${FILES_TO_COPY[@]}"; do
        if [ -f "${SOURCE_DIR}/${f}" ]; then
            HAS_ROOT_FILES=true
            break
        fi
    done

    if [ "$HAS_ROOT_FILES" = true ]; then
        echo -e "${YELLOW}注意: 源目录本身就是一个 Xenium 输出目录（无子片子目录）${NC}"
        DEST_DIR="${DEST_BASE}/${PROJECT_ID}"

        if [ "$DRY_RUN" = false ]; then
            mkdir -p "$DEST_DIR"
        fi

        for FILE in "${FILES_TO_COPY[@]}"; do
            SRC_FILE="${SOURCE_DIR}/${FILE}"
            if [ -f "$SRC_FILE" ]; then
                FILE_SIZE=$(du -sh "$SRC_FILE" 2>/dev/null | cut -f1)
                if [ "$DRY_RUN" = false ]; then
                    echo -e "  ${GREEN}✓ 拷贝 ${FILE} (${FILE_SIZE})${NC}"
                    cp "$SRC_FILE" "$DEST_DIR/"
                else
                    echo -e "  ${GREEN}✓ [DRY-RUN] 将拷贝 ${FILE} (${FILE_SIZE})${NC}"
                fi
            else
                echo -e "  ${RED}✗ 缺失: ${FILE}${NC}"
                MISSING_FILES+=("${PROJECT_ID}/${FILE}")
            fi
        done

        for DIR in "${DIRS_TO_COPY[@]}"; do
            SRC_SUBDIR="${SOURCE_DIR}/${DIR}"
            if [ -d "$SRC_SUBDIR" ]; then
                DIR_SIZE=$(du -sh "$SRC_SUBDIR" 2>/dev/null | cut -f1)
                if [ "$DRY_RUN" = false ]; then
                    echo -e "  ${GREEN}✓ 拷贝 ${DIR}/ (${DIR_SIZE})${NC}"
                    cp -r "$SRC_SUBDIR" "$DEST_DIR/"
                else
                    echo -e "  ${GREEN}✓ [DRY-RUN] 将拷贝 ${DIR}/ (${DIR_SIZE})${NC}"
                fi
            fi
        done

        SUCCESS_COUNT=1
    else
        echo -e "${RED}错误: 未在 ${SOURCE_DIR} 下找到任何 Xenium 输出数据${NC}"
        exit 1
    fi
fi

# ===== 汇总报告 =====
echo ""
echo -e "${CYAN}===============================================${NC}"
echo -e "${CYAN}  拷贝完成 - 汇总报告${NC}"
echo -e "${CYAN}===============================================${NC}"
echo -e "${BLUE}项目号:${NC}       ${PROJECT_ID}"
echo -e "${BLUE}处理片子数:${NC}   ${SUCCESS_COUNT}"
echo -e "${BLUE}跳过目录数:${NC}   ${SKIP_COUNT}"

if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo ""
    echo -e "${RED}缺失文件列表:${NC}"
    for mf in "${MISSING_FILES[@]}"; do
        echo -e "  ${RED}• ${mf}${NC}"
    done
fi

echo ""
echo -e "${GREEN}目标路径: ${DEST_BASE}/${PROJECT_ID}/${NC}"

if [ "$DRY_RUN" = true ]; then
    echo ""
    echo -e "${YELLOW}这是 DRY-RUN 模式，没有文件被实际拷贝。${NC}"
    echo -e "${YELLOW}去掉 --dry-run 参数以执行实际拷贝。${NC}"
fi
