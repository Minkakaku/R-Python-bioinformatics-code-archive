#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
空间转录组可视化
- 图1：细胞类型空间图（所有细胞，others为白色，其余按类型着色，无边框）
- 图2：干性评分空间图（所有细胞，others=0为白色，其余按评分四等分着色：白→浅红→红→深红）
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import scanpy as sc
from pathlib import Path

# ==================== 配置参数 ====================
H5AD_PATH = "subset_data.h5ad"
BOUNDARIES_PATH = "cell_boundaries.csv"
OUTPUT_DIR = "spatial_maps"

STEMNESS_GENES = [
    "LGR5", "ASCL2", "PROM1", "CD44", "SMOC2",
    "SOX9", "KLF4", "MYC", "HES1", "KLF5",
    "AXIN2", "RNF43", "EPCAM", "MSI2", "HMGA1"
]

# 细胞类型 → 颜色映射
CELLTYPE_COLORS = {
    "CXCL5+ Inf stem tumor Epi": "#FF5733",
    "IGF2+ prolif stem tumor Epi": "#FF7043",
    "Hypo Stem tumor Epi": "#FF8A65",
    "WNT+ ISC": "#FFAB91",
    "ESC-like Epi": "#FFCCBC",
    "Prol Metab Epi": "#1976D2",
    "EMT-like Inv Ent": "#42A5F5",
    "Inf Hyp Epi": "#64B5F6",
    "Early stress Epi": "#90CAF9",
    "Metab stress Epi": "#BBDEFB",
    "Mat Ent": "#388E3C",
    "MUC5B+ CGC": "#4CAF50",
    "MUC5AC+ MGC": "#66BB6A",
    "PCSK1+ AGC": "#81C784",
    "Hypo Metab Epi": "#A5D6A7",
    "others": "#FFFFFF"           # 图1中others用白色
}

# 干性评分四分组颜色（由低到高：白、浅红、红、深红）
STEMNESS_GROUP_COLORS = ["#FFFFFF", "#FFA0A0", "#FF0000", "#8B0000"]

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

    adata.layers["counts"] = adata.X.copy()
    print("正在标准化计数矩阵...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    available_genes = [g for g in STEMNESS_GENES if g in adata.var_names]
    missing_genes = set(STEMNESS_GENES) - set(available_genes)
    if missing_genes:
        print(f"警告: 以下基因不在数据集中，将被忽略: {missing_genes}")
    print(f"使用 {len(available_genes)}/{len(STEMNESS_GENES)} 个干性基因计算评分")

    sc.tl.score_genes(adata, gene_list=available_genes, score_name='stemness_score')
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

def build_cell_polygons(boundaries, adata, sample_name, x_col, y_col, include_others=True):
    """构建多边形，可选择是否包含 others 细胞（此处都包含）"""
    cell_subset = adata  # 始终使用全部细胞

    sample_cells = cell_subset[cell_subset.obs['sample'] == sample_name].obs.index.tolist()
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cells)]
    if len(sample_boundaries) == 0:
        return [], [], [], []

    grouped = sample_boundaries.groupby('cell_id')
    polygons, celltypes, scores, valid_ids = [], [], [], []
    for cell_id, group in grouped:
        if cell_id not in cell_subset.obs.index:
            continue
        points = group[[x_col, y_col]].values
        if len(points) < 3:
            continue
        polygon = patches.Polygon(points, closed=True)
        polygons.append(polygon)
        celltypes.append(cell_subset.obs.loc[cell_id, 'celltype'])
        scores.append(cell_subset.obs.loc[cell_id, 'stemness_score'])
        valid_ids.append(cell_id)

    print(f"  样本 {sample_name}: 构建 {len(polygons)} 个多边形 "
          f"(匹配率 {len(polygons)/len(sample_cells)*100:.1f}%)")
    return polygons, celltypes, scores, valid_ids

def plot_celltype_map(adata, boundaries, x_col, y_col, output_dir):
    """绘制细胞类型空间图（所有细胞，others为白色，其余按类型着色）"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制细胞类型空间图 ({len(samples)} 个样本) =====")

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes, _, _ = build_cell_polygons(
            boundaries, adata, sample, x_col, y_col, include_others=True
        )
        if not polygons:
            continue

        # 填充颜色：others → 白色，其余使用定义的颜色
        fill_colors = []
        for ct in celltypes:
            if ct == 'others':
                fill_colors.append('#FFFFFF')
            else:
                fill_colors.append(CELLTYPE_COLORS.get(ct, '#D3D3D3'))

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
        ax.set_title(f'Cell Type Map - Sample: {sample}', pad=20)

        # 图例（包含 others 为白色）
        handles = []
        visited = set()
        for ct in celltypes:
            if ct not in visited:
                visited.add(ct)
                if ct == 'others':
                    color = '#FFFFFF'
                else:
                    color = CELLTYPE_COLORS.get(ct, '#D3D3D3')
                handles.append(patches.Patch(color=color, label=ct))
        if handles:
            ax.legend(handles=handles, title='Cell Types', bbox_to_anchor=(1.02, 1),
                      loc='upper left', framealpha=0.9, edgecolor='black')

        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(output_dir / f"celltype_spatial_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"celltype_spatial_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 celltype_spatial_{safe_sample}.png/pdf")

def plot_stemness_map(adata, boundaries, x_col, y_col, output_dir):
    """
    绘制干性评分空间图（所有细胞，others 固定白色；
    上皮细胞按干性评分四等分，由低到高赋予白、浅红、红、深红）
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制干性评分空间图 ({len(samples)} 个样本) =====")

    # 仅基于上皮细胞（非others）计算四分位边界
    epi_mask = adata.obs['celltype'] != 'others'
    epi_scores = adata.obs.loc[epi_mask, 'stemness_score'].values
    if len(epi_scores) == 0:
        print("警告: 未找到任何上皮细胞，无法绘制干性图")
        return

    # 计算四个等分组的边界（依据分位数：0, 25%, 50%, 75%, 100%）
    quantiles = np.percentile(epi_scores, [0, 25, 50, 75, 100])
    # 边界做微调，确保四个区间：(q[0], q[1]]], (q[1], q[2]]], (q[2], q[3]]], (q[3], q[4]]]
    # 这里直接使用分位数作为区间边界，并将最低边界设为 -inf 以包含所有最小值
    bins = [-np.inf] + list(quantiles[1:-1]) + [np.inf]
    print(f"上皮细胞干性评分四分位边界: {quantiles}")
    print(f"分组区间: (<{quantiles[1]}] 白色, ({quantiles[1]}~{quantiles[2]}] 浅红, "
          f"({quantiles[2]}~{quantiles[3]}] 红, >{quantiles[3]} 深红)")

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes, scores, _ = build_cell_polygons(
            boundaries, adata, sample, x_col, y_col, include_others=True
        )
        if not polygons:
            continue

        # 根据评分分组确定颜色
        fill_colors = []
        for ct, s in zip(celltypes, scores):
            if ct == 'others':
                fill_colors.append('#FFFFFF')   # 白色
            else:
                # 找到所属分组 index (0,1,2,3)
                group_idx = np.digitize(s, bins) - 1
                # 防止越界（理论上不会，因为bins覆盖全实数）
                group_idx = min(max(group_idx, 0), 3)
                fill_colors.append(STEMNESS_GROUP_COLORS[group_idx])

        fig, ax = plt.subplots(figsize=(16, 12))
        pc = PatchCollection(polygons, facecolor=fill_colors, edgecolor='none', alpha=0.9)
        ax.add_collection(pc)

        all_verts = np.vstack([p.get_xy() for p in polygons])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        margin = 0.06
        ax.set_xlim(x_min - margin*(x_max-x_min), x_max + margin*(x_max-x_min))
        ax.set_ylim(y_min - margin*(y_max-y_min), y_max + margin*(y_max-y_min))
        ax.set_aspect('equal')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'Stemness Score (Quartile) Map - Sample: {sample}', pad=20)

        # 手动添加图例表示四个分组
        legend_labels = [
            f'Low (≤{quantiles[1]:.2f})',
            f'Mid-Low ({quantiles[1]:.2f}–{quantiles[2]:.2f})',
            f'Mid-High ({quantiles[2]:.2f}–{quantiles[3]:.2f})',
            f'High (>{quantiles[3]:.2f})'
        ]
        legend_handles = [patches.Patch(color=c, label=lab)
                          for c, lab in zip(STEMNESS_GROUP_COLORS, legend_labels)]
        # 再添加一个others条目
        legend_handles.append(patches.Patch(color='#FFFFFF', label='Others'))
        ax.legend(handles=legend_handles, title='Stemness Group',
                  bbox_to_anchor=(1.02, 1), loc='upper left',
                  framealpha=0.9, edgecolor='black')

        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(output_dir / f"stemscore_spatial_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"stemscore_spatial_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 stemscore_spatial_{safe_sample}.png/pdf")

def main():
    try:
        adata = load_and_preprocess_adata(H5AD_PATH)
        boundaries, x_col, y_col = load_boundaries(BOUNDARIES_PATH)

        epi_count = (adata.obs['celltype'] != 'others').sum()
        other_count = (adata.obs['celltype'] == 'others').sum()
        print(f"\n细胞统计: 上皮 {epi_count} / 其他 {other_count} / 总计 {adata.n_obs}")

        adata_ids = set(adata.obs.index)
        boundary_ids = set(boundaries['cell_id'].unique())
        matched = adata_ids & boundary_ids
        print(f"细胞ID匹配: {len(matched)}/{len(adata_ids)} ({len(matched)/len(adata_ids)*100:.1f}%)")
        if len(matched) == 0:
            print("错误: 无匹配细胞ID")
            return

        plot_celltype_map(adata, boundaries, x_col, y_col, OUTPUT_DIR)
        plot_stemness_map(adata, boundaries, x_col, y_col, OUTPUT_DIR)
        print("\n所有任务完成！")

    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
