#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制免疫细胞类型的空间分布图（每个免疫细胞类型独立颜色，其他细胞浅灰色背景，白色边框线宽0.05）
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import scanpy as sc
from pathlib import Path

# ==================== 配置参数 ====================
H5AD_PATH = "TMA2.h5ad"               # 你的 h5ad 文件
BOUNDARIES_PATH = "cell_boundaries.csv"  # 细胞边界文件
OUTPUT_DIR = "spatial_immu"          # 输出目录

# 为每个免疫细胞类型独立分配一种颜色（共16种，互不重复）
SELECTED_TYPE_COLORS = {
    "CD8+ T": "#E6194B",                # 红色
    "CXCL9+ CD8+ T": "#3CB44B",         # 绿色
    "MKI67+ T": "#FFE119",              # 黄色
    "CD4+T": "#4363D8",                 # 蓝色
    "Treg": "#F58231",                  # 橙色
    "Mature B": "#911EB4",              # 紫色
    "B naive": "#46F0F0",               # 青色
    "plasma": "#F032E6",                # 品红
    "Mast": "#BCF60C",                  # 亮绿
    "CXCL9+ M1": "#FABEBE",             # 浅粉
    "M2": "#E6BEFF",                    # 浅紫
    "Monocyte-derived Macrophages": "#FFFAC8",  # 浅黄
    "M1": "#008080",                    # 深青
    "LAMP3+DC": "#E6DB74",              # 橄榄黄
    "cDC": "#800000",                   # 栗色
    "Neutrophils": "#AAFFC3"            # 薄荷绿
}

# ==================== 数据加载 ====================
def load_and_preprocess_adata(h5ad_path):
    """读取 h5ad，确保 obs 中有 celltype 列，并将 cell_id 设为索引（如果存在）"""
    print(f"读取数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

    if 'celltype' not in adata.obs.columns:
        raise KeyError("obs 中缺少 'celltype' 列，无法进行细胞类型筛选")

    # 如果有 cell_id 列，处理重复并设为索引
    if 'cell_id' in adata.obs.columns:
        dup_mask = adata.obs.duplicated(subset='cell_id', keep=False)
        if dup_mask.sum() > 0:
            print(f"发现 {dup_mask.sum()} 个重复 cell_id，保留第一个出现")
            adata = adata[~adata.obs.duplicated(subset='cell_id', keep='first'), :].copy()
        adata.obs.set_index('cell_id', inplace=True)
        print("已将 cell_id 设为索引")
    else:
        print("警告: 未找到 'cell_id' 列，将使用当前索引作为细胞 ID，请确保其与边界文件中的 cell_id 一致")

    return adata

def load_boundaries(boundaries_path):
    """读取细胞边界文件，自动识别坐标列"""
    print(f"读取边界文件: {boundaries_path}")
    boundaries = pd.read_csv(boundaries_path)

    # 自动检测坐标列名
    x_col, y_col = None, None
    for col in ['vertex_x', 'x', 'X']:
        if col in boundaries.columns:
            x_col = col
            break
    for col in ['vertex_y', 'y', 'Y']:
        if col in boundaries.columns:
            y_col = col
            break

    if x_col is None or y_col is None:
        num_cols = boundaries.select_dtypes(include=np.number).columns
        if len(num_cols) >= 2:
            x_col, y_col = num_cols[0], num_cols[1]
            print(f"自动选择坐标列: {x_col}, {y_col}")
        else:
            raise ValueError("无法确定坐标列，请检查列名是否为 vertex_x/vertex_y 或 x/y")

    print(f"边界文件包含 {boundaries['cell_id'].nunique()} 个独特细胞")
    return boundaries, x_col, y_col

# ==================== 构建多边形（包含所有细胞） ====================
def build_all_polygons(boundaries, adata, sample_name, x_col, y_col):
    """
    构建指定样本中所有细胞的多边形。
    若细胞类型在 SELECTED_TYPE_COLORS 中则使用独立颜色，否则统一使用浅灰色。
    """
    sample_cells = adata[adata.obs['sample'] == sample_name].obs.index.tolist()
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cells)]

    if len(sample_boundaries) == 0:
        return [], [], []

    grouped = sample_boundaries.groupby('cell_id')
    polygons, celltypes, fill_colors = [], [], []
    
    for cell_id, group in grouped:
        points = group[[x_col, y_col]].values
        if len(points) < 3:
            continue
        if cell_id not in adata.obs.index:
            continue
        ct = adata.obs.loc[cell_id, 'celltype']
        polygons.append(patches.Polygon(points, closed=True))
        celltypes.append(ct)
        if ct in SELECTED_TYPE_COLORS:
            fill_colors.append(SELECTED_TYPE_COLORS[ct])
        else:
            fill_colors.append("#D3D3D3")  # 浅灰色

    print(f"  样本 {sample_name}: 构建 {len(polygons)} 个多边形 (包含所有细胞)")
    return polygons, celltypes, fill_colors

# ==================== 绘图 ====================
def plot_immune_celltype_map(adata, boundaries, x_col, y_col, output_dir):
    """为每个样本绘制免疫细胞空间图，其他细胞为浅灰色背景"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制免疫细胞空间图 ({len(samples)} 个样本) =====")

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes, fill_colors = build_all_polygons(boundaries, adata, sample, x_col, y_col)
        if not polygons:
            print(f"  样本 {sample} 无细胞数据，跳过")
            continue

        fig, ax = plt.subplots(figsize=(16, 12))
        # 白色边框，线宽 0.05
        pc = PatchCollection(polygons, facecolor=fill_colors,
                             edgecolor='white', linewidth=0.05, alpha=0.9)
        ax.add_collection(pc)

        # 坐标范围
        all_verts = np.vstack([p.get_xy() for p in polygons])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        margin = 0.06
        ax.set_xlim(x_min - margin * (x_max - x_min), x_max + margin * (x_max - x_min))
        ax.set_ylim(y_min - margin * (y_max - y_min), y_max + margin * (y_max - y_min))
        ax.set_aspect('equal')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'Immune Celltypes Map - Sample: {sample}', pad=20)

        # 图例：显示实际出现的免疫细胞类型，以及 "Other"
        present_types = set(celltypes)
        handles = []
        has_other = False
        # 按字典顺序显示图例，便于查找
        for ct_name, color in SELECTED_TYPE_COLORS.items():
            if ct_name in present_types:
                handles.append(patches.Patch(color=color, label=ct_name))
        # 检查是否存在非免疫细胞
        for ct in present_types:
            if ct not in SELECTED_TYPE_COLORS:
                has_other = True
                break
        if has_other:
            handles.append(patches.Patch(color="#D3D3D3", label="Other"))

        if handles:
            ax.legend(handles=handles, title='Cell Types',
                      bbox_to_anchor=(1.02, 1), loc='upper left',
                      framealpha=0.9, edgecolor='black')

        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(output_dir / f"immune_spatial_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"immune_spatial_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 immune_spatial_{safe_sample}.png/pdf")

# ==================== 主流程 ====================
def main():
    try:
        adata = load_and_preprocess_adata(H5AD_PATH)
        boundaries, x_col, y_col = load_boundaries(BOUNDARIES_PATH)

        adata_ids = set(adata.obs.index)
        boundary_ids = set(boundaries['cell_id'].unique())
        matched = adata_ids & boundary_ids
        print(f"细胞ID匹配: {len(matched)}/{len(adata_ids)} ({len(matched)/len(adata_ids)*100:.1f}%)")
        if len(matched) == 0:
            raise ValueError("无匹配细胞ID，无法绘图")

        plot_immune_celltype_map(adata, boundaries, x_col, y_col, OUTPUT_DIR)

        print("\n所有任务完成！")

    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
