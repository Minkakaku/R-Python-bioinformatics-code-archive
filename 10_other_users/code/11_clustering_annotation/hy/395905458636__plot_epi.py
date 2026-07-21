#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
空间转录组干性评分可视化
- 边框颜色：细胞类型（线宽加粗）
- 填充颜色：干性评分（白→红渐变）
- 图例位于图形右侧，与颜色条不重叠
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
OUTPUT_DIR = "spatial_stemness_maps"

STEMNESS_GENES = [
    "LGR5", "ASCL2", "PROM1", "CD44", "SMOC2",
    "SOX9", "KLF4", "MYC", "HES1", "KLF5",
    "AXIN2", "RNF43", "EPCAM", "MSI2", "HMGA1"
]

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
    "others": "#D3D3D3"
}

STEMNESS_CMAP = LinearSegmentedColormap.from_list("stemness", ["#FFFFFF", "#FF0000"])

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

def build_cell_polygons(boundaries, adata, sample_name, x_col, y_col):
    sample_cells = adata[adata.obs['sample'] == sample_name].obs.index.tolist()
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cells)]
    if len(sample_boundaries) == 0:
        return [], [], [], []
    grouped = sample_boundaries.groupby('cell_id')
    polygons, celltypes, scores, valid_ids = [], [], [], []
    for cell_id, group in grouped:
        if cell_id not in adata.obs.index:
            continue
        points = group[[x_col, y_col]].values
        if len(points) < 3:
            continue
        polygon = patches.Polygon(points, closed=True)
        polygons.append(polygon)
        celltypes.append(adata.obs.loc[cell_id, 'celltype'])
        scores.append(adata.obs.loc[cell_id, 'stemness_score'])
        valid_ids.append(cell_id)
    print(f"  样本 {sample_name}: 构建 {len(polygons)} 个多边形 "
          f"(匹配率 {len(polygons)/len(sample_cells)*100:.1f}%)")
    return polygons, celltypes, scores, valid_ids

def plot_sample_spatial(adata, boundaries, x_col, y_col, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = adata.obs['sample'].unique()
    print(f"\n处理 {len(samples)} 个样本...")

    all_scores = adata.obs['stemness_score'].values
    vmin, vmax = np.percentile(all_scores, [1, 99])

    for sample in samples:
        print(f"\n=== 样本: {sample} ===")
        polygons, celltypes, scores, _ = build_cell_polygons(
            boundaries, adata, sample, x_col, y_col
        )
        if not polygons:
            continue

        border_colors = []
        fill_colors = []
        for ct, score in zip(celltypes, scores):
            border_colors.append(CELLTYPE_COLORS.get(ct, CELLTYPE_COLORS["others"]))
            norm_score = np.clip((score - vmin) / (vmax - vmin), 0, 1)
            fill_colors.append(STEMNESS_CMAP(norm_score))

        # 增加图形宽度，为右侧元素留足空间
        fig, ax = plt.subplots(figsize=(16, 12))

        # 填充多边形
        pc_fill = PatchCollection(polygons, facecolor=fill_colors, edgecolor='none', alpha=0.9)
        ax.add_collection(pc_fill)

        # 绘制边框（线宽加粗）
        BORDER_LINEWIDTH = 0.05
        for poly, bcolor in zip(polygons, border_colors):
            ax.add_patch(patches.PathPatch(poly.get_path(),
                                           facecolor='none',
                                           edgecolor=bcolor,
                                           linewidth=BORDER_LINEWIDTH))

        # 坐标范围
        all_verts = np.vstack([p.get_xy() for p in polygons])
        x_min, x_max = all_verts[:, 0].min(), all_verts[:, 0].max()
        y_min, y_max = all_verts[:, 1].min(), all_verts[:, 1].max()
        margin = 0.05
        ax.set_xlim(x_min - margin*(x_max-x_min), x_max + margin*(x_max-x_min))
        ax.set_ylim(y_min - margin*(y_max-y_min), y_max + margin*(y_max-y_min))
        ax.set_aspect('equal')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'Spatial Stemness Map - Sample: {sample}', pad=20)

        # ---- 构建图例（置于更右侧，避免与颜色条重叠）----
        legend_types = []
        for ct in celltypes:
            if ct in CELLTYPE_COLORS and ct != "others":
                legend_types.append(ct)
            else:
                legend_types.append("others")
        unique_types = sorted(set(legend_types) - {"others"})
        if "others" in set(legend_types):
            unique_types.append("others")

        legend_handles = []
        for ct in unique_types:
            color = CELLTYPE_COLORS.get(ct, CELLTYPE_COLORS["others"])
            legend_handles.append(patches.Patch(color=color, label=ct))

        # 将图例放置在更靠右的位置 (bbox_to_anchor x 增大)
        ax.legend(handles=legend_handles,
                  loc='center left',
                  bbox_to_anchor=(1.15, 0.5),   # 向右移动
                  title='Cell Type (border)',
                  framealpha=0.9)

        # 颜色条（保持原位置）
        sm = ScalarMappable(cmap=STEMNESS_CMAP, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Stemness Score', rotation=270, labelpad=15)

        # 调整布局：压缩主图区域，给右侧更多空间
        plt.subplots_adjust(right=0.72)   # 减小 right 值，让主图区域变窄

        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        out_png = output_dir / f"stemness_spatial_{safe_sample}.png"
        out_pdf = output_dir / f"stemness_spatial_{safe_sample}.pdf"
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_pdf, bbox_inches='tight')
        plt.close()
        print(f"  图片已保存: {out_png}")

    print("\n所有样本处理完成！")

def main():
    try:
        adata = load_and_preprocess_adata(H5AD_PATH)
        boundaries, x_col, y_col = load_boundaries(BOUNDARIES_PATH)
        adata_ids = set(adata.obs.index)
        boundary_ids = set(boundaries['cell_id'].unique())
        matched = adata_ids & boundary_ids
        print(f"\n细胞ID匹配: {len(matched)}/{len(adata_ids)} ({len(matched)/len(adata_ids)*100:.1f}%)")
        if len(matched) == 0:
            print("错误: 无匹配细胞ID")
            return
        plot_sample_spatial(adata, boundaries, x_col, y_col, OUTPUT_DIR)
    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
