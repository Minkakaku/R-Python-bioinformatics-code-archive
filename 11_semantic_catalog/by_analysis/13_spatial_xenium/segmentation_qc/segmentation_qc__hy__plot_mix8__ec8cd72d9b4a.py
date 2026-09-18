#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制指定细胞类型的空间分布图（上皮红系、CAF绿系、CD8+T蓝系，其余灰色）
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
H5AD_PATH = "TMA2.h5ad"
BOUNDARIES_PATH = "cell_boundaries.csv"
OUTPUT_DIR = "spatial_selected_mix8"          # 输出目录

# 新色彩方案：上皮（6种红色系）、CAF（2种绿色系）、CD8+ T（2种蓝色系）
SELECTED_TYPE_COLORS = {
    # 红色系 – 上皮细胞
    "Proliferative Epithelials": "#67000d",   # 深红
    "Hypoxic Epithelials":        "#a50f15",   # 暗红
    "IFN-Epithelials":            "#cb181d",   # 正红
    "TA-like Epithelials":        "#ef3b2c",   # 亮红
    "Wnt-active stem":            "#fb6a4a",   # 橙红
    "Stem-like Tumor":            "#fcbba1",   # 浅粉红
    # 绿色系 – CAF
    "iCAF":                       "#238b45",   # 深绿
    "Matrix Fibroblasts":         "#74c476",   # 浅绿
    # 蓝色系 – CD8+ T 细胞
    "CD8+ T":                     "#08519c",   # 深蓝
    "CXCL9+ CD8+ T":              "#6baed6"    # 浅蓝
}

# ==================== 数据加载（不变） ====================
def load_and_preprocess_adata(h5ad_path):
    """读取 h5ad，确保 obs 中有 celltype 列，并将 cell_id 设为索引"""
    print(f"读取数据: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

    if 'celltype' not in adata.obs.columns:
        raise KeyError("obs 中缺少 'celltype' 列")

    if 'cell_id' in adata.obs.columns:
        dup_mask = adata.obs.duplicated(subset='cell_id', keep=False)
        if dup_mask.sum() > 0:
            print(f"发现 {dup_mask.sum()} 个重复 cell_id，保留第一个出现")
            adata = adata[~adata.obs.duplicated(subset='cell_id', keep='first'), :].copy()
        adata.obs.set_index('cell_id', inplace=True)
        print("已将 cell_id 设为索引")
    else:
        print("警告: 未找到 'cell_id' 列，将使用当前索引，请确保其与边界文件中的 cell_id 一致")

    return adata

def load_boundaries(boundaries_path):
    """读取细胞边界文件，自动识别坐标列"""
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

# ==================== 构建多边形（包含所有细胞） ====================
def build_all_polygons(boundaries, adata, sample_name, x_col, y_col):
    """
    构建指定样本中所有细胞的多边形。
    若细胞类型在 SELECTED_TYPE_COLORS 中则使用独立颜色，否则浅灰色。
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
            fill_colors.append("#D3D3D3")  # 浅灰

    print(f"  样本 {sample_name}: 构建 {len(polygons)} 个多边形")
    return polygons, celltypes, fill_colors

# ==================== 绘图（标题与图例更新） ====================
def plot_selected_celltype_map(adata, boundaries, x_col, y_col, output_dir):
    """为每个样本绘制空间图（目标细胞类型彩色，其余灰色）"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制选定细胞类型空间图 ({len(samples)} 个样本) =====")

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes, fill_colors = build_all_polygons(boundaries, adata, sample, x_col, y_col)
        if not polygons:
            print(f"  样本 {sample} 无细胞数据，跳过")
            continue

        fig, ax = plt.subplots(figsize=(16, 12))
        pc = PatchCollection(polygons, facecolor=fill_colors,
                             edgecolor='white', linewidth=0.05, alpha=0.9)
        ax.add_collection(pc)

        all_verts = np.vstack([p.get_xy() for p in polygons])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        margin = 0.06
        ax.set_xlim(x_min - margin * (x_max - x_min), x_max + margin * (x_max - x_min))
        ax.set_ylim(y_min - margin * (y_max - y_min), y_max + margin * (y_max - y_min))
        ax.set_aspect('equal')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'Selected Celltypes Map - Sample: {sample}', pad=20)

        # 图例：实际存在的目标细胞类型 + 可能的 Other
        present_types = set(celltypes)
        handles = []
        has_other = False
        for ct_name, color in SELECTED_TYPE_COLORS.items():
            if ct_name in present_types:
                handles.append(patches.Patch(color=color, label=ct_name))
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
        plt.savefig(output_dir / f"selected_spatial_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"selected_spatial_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 selected_spatial_{safe_sample}.png/pdf")

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

        plot_selected_celltype_map(adata, boundaries, x_col, y_col, OUTPUT_DIR)

        print("\n所有任务完成！")

    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
