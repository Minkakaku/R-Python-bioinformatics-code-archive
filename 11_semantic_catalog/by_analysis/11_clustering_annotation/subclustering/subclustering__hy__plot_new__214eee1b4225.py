#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
空间转录组可视化（多通路版本，评分仅限指定上皮类型，统一边框，图例有序）
- 图1：细胞类型空间图（图例按CELLTYPE_COLORS顺序，其余归为others）
- 图2：各通路功能评分空间图（评分四等分，白/浅红/红/深红；others白色；
       所有细胞带浅灰边框，线宽0.01）
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
H5AD_PATH = "subset_data.h5ad"
BOUNDARIES_PATH = "cell_boundaries.csv"
GENE_FUNC_PATH = "gene_func.csv"
OUTPUT_DIR = "spatial_maps"

# 细胞类型 → 颜色映射（请注意保留字典顺序，Python 3.7+ 自动维护插入顺序）
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
    "others": "#FFFFFF"       # 所有未定义类型统一为 others，白色
}

# 所有被选择的上皮类型（用于评分限制）
VALID_EPI_TYPES = set(CELLTYPE_COLORS.keys()) - {"others"}

# 评分四分组颜色（由低到高：白、浅红、红、深红）
SCORE_GROUP_COLORS = ["#FFFFFF", "#FFA0A0", "#FF0000", "#8B0000"]

# ==================== 通路-基因读取 ====================
def load_gene_func(file_path):
    print(f"读取通路-基因文件: {file_path}")
    with open(file_path, 'r') as f:
        first_line = f.readline()
    sep = ',' if first_line.count(',') >= 2 else '\t'

    df = pd.read_csv(file_path, sep=sep)
    if not {'Function', 'Gene'}.issubset(df.columns):
        raise ValueError("gene_func文件必须包含 'Function' 和 'Gene' 列")

    func_gene_dict = df.groupby('Function')['Gene'].apply(list).to_dict()
    print(f"发现 {len(func_gene_dict)} 个功能通路:")
    for func, genes in func_gene_dict.items():
        print(f"  - {func}: {len(genes)} 个基因")
    return func_gene_dict

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
    return adata

def compute_function_scores(adata, func_gene_dict):
    for func_name, gene_list in func_gene_dict.items():
        available = [g for g in gene_list if g in adata.var_names]
        missing = set(gene_list) - set(available)
        if missing:
            print(f"  功能 '{func_name}': 缺失基因 {missing}，将使用 {len(available)}/{len(gene_list)} 个基因")
        if len(available) == 0:
            print(f"  功能 '{func_name}': 无有效基因，跳过评分计算")
            continue

        score_name = f"{func_name}_score"
        sc.tl.score_genes(adata, gene_list=available, score_name=score_name)
        print(f"  功能 '{func_name}': 评分列 '{score_name}' 已生成")
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
    """返回 polygons, celltypes, scores_df"""
    sample_cells = adata[adata.obs['sample'] == sample_name].obs.index.tolist()
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cells)]
    if len(sample_boundaries) == 0:
        return [], [], pd.DataFrame()

    grouped = sample_boundaries.groupby('cell_id')
    polygons, celltypes, valid_ids = [], [], []
    for cell_id, group in grouped:
        if cell_id not in adata.obs.index:
            continue
        points = group[[x_col, y_col]].values
        if len(points) < 3:
            continue
        polygons.append(patches.Polygon(points, closed=True))
        celltypes.append(adata.obs.loc[cell_id, 'celltype'])
        valid_ids.append(cell_id)

    score_cols = [c for c in adata.obs.columns if c.endswith('_score')]
    scores_df = adata.obs.loc[valid_ids, score_cols].copy()
    print(f"  样本 {sample_name}: 构建 {len(polygons)} 个多边形 "
          f"(匹配率 {len(polygons)/len(sample_cells)*100:.1f}%)")
    return polygons, celltypes, scores_df

# ==================== 细胞类型图 ====================
def plot_celltype_map(adata, boundaries, x_col, y_col, output_dir):
    """细胞类型空间图：图例按 CELLTYPE_COLORS 顺序，其余归为 others"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制细胞类型空间图 ({len(samples)} 个样本) =====")

    # 辅助函数：将 celltype 映射为颜色及图例标签
    def map_color_and_label(ct):
        if ct in CELLTYPE_COLORS:   # 包括 "others"
            return CELLTYPE_COLORS[ct], ct
        else:
            return "#FFFFFF", "others"

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes, _ = build_cell_polygons(boundaries, adata, sample, x_col, y_col)
        if not polygons:
            continue

        # 颜色渲染
        fill_colors = []
        mapped_labels = []
        for ct in celltypes:
            color, label = map_color_and_label(ct)
            fill_colors.append(color)
            mapped_labels.append(label)

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

        # 图例：按照 CELLTYPE_COLORS 顺序，只显示实际出现的类型
        present_types = set(mapped_labels)
        handles = []
        for ct_name in CELLTYPE_COLORS.keys():
            if ct_name in present_types:
                color = CELLTYPE_COLORS[ct_name]
                label = ct_name
                handles.append(patches.Patch(color=color, label=label))
        # 如果出现未定义的类型（此时它们的 label 已经是 'others'），
        # 则 'others' 已经在 CELLTYPE_COLORS 的键中，会被上面的循环捕获。
        # 所以不再需要额外添加。

        if handles:
            ax.legend(handles=handles, title='Cell Types', bbox_to_anchor=(1.02, 1),
                      loc='upper left', framealpha=0.9, edgecolor='black')

        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(output_dir / f"celltype_spatial_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"celltype_spatial_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 celltype_spatial_{safe_sample}.png/pdf")

# ==================== 通路评分图（四分组 + 浅灰边框） ====================
def plot_function_score_map(adata, boundaries, x_col, y_col, output_dir, function_name):
    score_col = f"{function_name}_score"
    if score_col not in adata.obs:
        print(f"  列 '{score_col}' 不存在，跳过")
        return

    samples = adata.obs['sample'].unique()
    print(f"\n===== 绘制 {function_name} 评分空间图 ({len(samples)} 个样本) =====")

    # 仅基于有效上皮细胞（字典中除 others 以外的类型）计算分位数
    epi_mask = adata.obs['celltype'].isin(VALID_EPI_TYPES)
    epi_scores = adata.obs.loc[epi_mask, score_col].dropna()
    if len(epi_scores) == 0:
        print(f"  功能 '{function_name}': 无有效上皮细胞评分，跳过")
        return

    quantiles = np.percentile(epi_scores, [0, 25, 50, 75, 100])
    bins = [-np.inf] + list(quantiles[1:-1]) + [np.inf]
    group_labels = [
        f'Low (≤{quantiles[1]:.2f})',
        f'Mid-Low ({quantiles[1]:.2f}–{quantiles[2]:.2f})',
        f'Mid-High ({quantiles[2]:.2f}–{quantiles[3]:.2f})',
        f'High (>{quantiles[3]:.2f})'
    ]
    print(f"  上皮细胞评分分位点: {quantiles}")

    for sample in samples:
        print(f"\n样本: {sample}")
        polygons, celltypes, scores_df = build_cell_polygons(boundaries, adata, sample, x_col, y_col)
        if polygons is None or len(polygons) == 0:
            continue

        sample_scores = scores_df[score_col].values
        fill_colors = []
        for ct, s in zip(celltypes, sample_scores):
            # others 或者未定义类型统一白色
            if ct not in VALID_EPI_TYPES:   # 包括 'others' 或任何不在字典中的类型
                fill_colors.append('#FFFFFF')
            else:
                group_idx = np.digitize(s, bins) - 1
                group_idx = min(max(group_idx, 0), 3)
                fill_colors.append(SCORE_GROUP_COLORS[group_idx])

        fig, ax = plt.subplots(figsize=(16, 12))
        # 关键：添加浅灰色边框，线宽 0.01
        pc = PatchCollection(polygons, facecolor=fill_colors,
                             edgecolor='lightgrey', linewidth=0.01, alpha=0.9)
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
        ax.set_title(f'{function_name} Score Map - Sample: {sample}', pad=20)

        # 图例：四分组 + Others
        legend_handles = [patches.Patch(color=c, label=lab)
                          for c, lab in zip(SCORE_GROUP_COLORS, group_labels)]
        legend_handles.append(patches.Patch(color='#FFFFFF', label='Others'))
        ax.legend(handles=legend_handles, title=f'{function_name} Group',
                  bbox_to_anchor=(1.02, 1), loc='upper left',
                  framealpha=0.9, edgecolor='black')

        safe_function = str(function_name).replace('/', '_').replace(' ', '_')
        safe_sample = str(sample).replace('/', '_').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(output_dir / f"score_{safe_function}_spatial_{safe_sample}.png", dpi=300)
        plt.savefig(output_dir / f"score_{safe_function}_spatial_{safe_sample}.pdf")
        plt.close()
        print(f"  已保存 score_{safe_function}_spatial_{safe_sample}.png/pdf")

# ==================== 主流程 ====================
def main():
    try:
        # 1. 数据加载
        adata = load_and_preprocess_adata(H5AD_PATH)
        boundaries, x_col, y_col = load_boundaries(BOUNDARIES_PATH)

        # 2. 通路‑基因列表与评分计算
        func_gene_dict = load_gene_func(GENE_FUNC_PATH)
        adata = compute_function_scores(adata, func_gene_dict)

        # 3. 统计信息
        epi_count = adata.obs['celltype'].isin(VALID_EPI_TYPES).sum()
        other_count = (~adata.obs['celltype'].isin(VALID_EPI_TYPES)).sum()
        print(f"\n细胞统计: 选定上皮类型 {epi_count} / 其他 (含others) {other_count} / 总计 {adata.n_obs}")

        adata_ids = set(adata.obs.index)
        boundary_ids = set(boundaries['cell_id'].unique())
        matched = adata_ids & boundary_ids
        print(f"细胞ID匹配: {len(matched)}/{len(adata_ids)} ({len(matched)/len(adata_ids)*100:.1f}%)")
        if len(matched) == 0:
            raise ValueError("无匹配细胞ID，无法绘图")

        # 4. 仅绘制一次细胞类型图（按用户需求）
        plot_celltype_map(adata, boundaries, x_col, y_col, OUTPUT_DIR)

        # 5. 为每个通路绘制评分图（所有细胞带浅灰边框）
        for func_name in func_gene_dict.keys():
            plot_function_score_map(adata, boundaries, x_col, y_col, OUTPUT_DIR, func_name)

        print("\n所有任务完成！")

    except Exception as e:
        print(f"运行出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
