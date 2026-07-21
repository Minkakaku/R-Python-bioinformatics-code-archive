#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
空间转录组可视化 - 仅绘制选定细胞类型的空间分布
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
H5AD_PATH = "TMA.h5ad"
BOUNDARIES_PATH = "cell_boundaries.csv"
OUTPUT_DIR = "select_spatial_maps"

# 选定的细胞类型及其颜色（顺序决定图例顺序）
SELECTED_CELLTYPE_COLORS = {
    "MKI67+ Endo": "#A5D6A7",
    "iCAF": "#388E3C",
    "CD8+ T": "#FF5733",
    "CXCL9+ CD8+ T": "#FF7043",
    "MKI67+ T": "#FF8A65",
    "Proliferative Epithelials": "#0D47A1",
    "Hypoxic Epithelials": "#1565C0",
    "IFN-Epithelials": "#2196F3",
    "TA-like Epithelials": "#63B8FF",
    "Wnt-active stem": "#93C5FD",
    "Stem-like Tumor": "#C7E0FC",
}

# 其他所有细胞类型的颜色
OTHER_COLOR = "#D3D3D3"

# ==================== 数据加载与预处理 ====================
def load_and_preprocess_adata(h5ad_path):
    print(f"读取数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

    if 'celltype' not in adata.obs.columns:
        raise KeyError("obs 中缺少 'celltype' 列")
    if 'sample' not in adata.obs.columns:
        raise KeyError("obs 中缺少 'sample' 列")

    if 'cell_id' in adata.obs.columns:
        dup_mask = adata.obs.duplicated(subset='cell_id', keep=False)
        if dup_mask.sum() > 0:
            print(f"发现 {dup_mask.sum()} 个重复 cell_id，保留第一个出现")
            adata = adata[~adata.obs.duplicated(subset='cell_id', keep='first'), :].copy()
        adata.obs.set_index('cell_id', inplace=True)
        print("已将 cell_id 设为索引")
    else:
        print("警告: 未找到 'cell_id' 列，将使用当前索引作为细胞 ID")
    return adata

def load_boundaries(boundaries_path):
    print(f"读取边界文件: {boundaries_path}")
    boundaries = pd.read_csv(boundaries_path)
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
            raise ValueError("无法确定坐标列")
    print(f"边界文件包含 {boundaries['cell_id'].nunique()} 个独特细胞")
    return boundaries, x_col, y_col

def build_cell_polygons(boundaries, adata, sample_name, x_col, y_col):
    """返回 polygons 和 celltypes 列表"""
    sample_cells = adata[adata.obs['sample'] == sample_name].obs.index.tolist()
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cells)]
    if len(sample_boundaries) == 0:
        return [], []

    grouped = sample_boundaries.groupby('cell_id')
    polygons, celltypes = [], []
    for cell_id, group in grouped:
        if cell_id not in adata.obs.index:
            continue
        points = group[[x_col, y_col]].values
        if len(points) < 3:
            continue
        polygons.append(patches.Polygon(points, closed=True))
        celltypes.append(adata.obs.loc[cell_id, 'celltype'])
    print(f"  样本 {sample_name}: 构建 {len(polygons)} 个多边形")
    return polygons, celltypes

# ==================== 绘制细胞类型图 ====================
def plot_selected_celltype_map(adata, boundaries, x_col, y_col, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制选定细胞类型空间图 ({len(samples)} 个样本) =====")

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes = build_cell_polygons(boundaries, adata, sample, x_col, y_col)
        if not polygons:
            continue

        # 确定每个细胞的颜色和显示标签
        fill_colors = []
        present_selected = set()    # 记录本次样本中出现的选定类型
        has_other = False
        for ct in celltypes:
            if ct in SELECTED_CELLTYPE_COLORS:
                fill_colors.append(SELECTED_CELLTYPE_COLORS[ct])
                present_selected.add(ct)
            else:
                fill_colors.append(OTHER_COLOR)
                has_other = True

        fig, ax = plt.subplots(figsize=(16, 12))
        pc = PatchCollection(polygons, facecolor=fill_colors, edgecolor='none', alpha=0.9)
        ax.add_collection(pc)

        # 坐标范围
        all_verts = np.vstack([p.get_xy() for p in polygons])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        margin = 0.06
        ax.set_xlim(x_min - margin*(x_max-x_min), x_max + margin*(x_max-x_min))
        ax.set_ylim(y_min - margin*(y_max-y_min), y_max + margin*(y_max-y_min))
        ax.set_aspect('equal')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'Selected Cell Types - Sample: {sample}', pad=20)

        # 图例：严格按 SELECTED_CELLTYPE_COLORS 的顺序添加存在的选定类型，最后加 Others
        handles = []
        for ct_name, color in SELECTED_CELLTYPE_COLORS.items():
            if ct_name in present_selected:
                handles.append(patches.Patch(color=color, label=ct_name))
        if has_other:
            handles.append(patches.Patch(color=OTHER_COLOR, label='Others'))

        if handles:
            ax.legend(handles=handles, title='Cell Types', bbox_to_anchor=(1.02, 1),
                      loc='upper left', framealpha=0.9, edgecolor='black')

        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(output_dir / f"selected_celltype_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"selected_celltype_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 selected_celltype_{safe_sample}.png/pdf")

# ==================== 主流程 ====================
def main():
    try:
        adata = load_and_preprocess_adata(H5AD_PATH)
        boundaries, x_col, y_col = load_boundaries(BOUNDARIES_PATH)

        # 匹配统计
        adata_ids = set(adata.obs.index)
        boundary_ids = set(boundaries['cell_id'].unique())
        matched = adata_ids & boundary_ids
        print(f"细胞ID匹配: {len(matched)}/{len(adata_ids)} ({len(matched)/len(adata_ids)*100:.1f}%)")
        if len(matched) == 0:
            raise ValueError("无匹配细胞ID，无法绘图")

        plot_selected_celltype_map(adata, boundaries, x_col, y_col, OUTPUT_DIR)
        print("\n所有任务完成！")

    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
