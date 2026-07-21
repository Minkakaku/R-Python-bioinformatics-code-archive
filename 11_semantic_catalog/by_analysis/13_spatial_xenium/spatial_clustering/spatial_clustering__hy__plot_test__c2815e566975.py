"""
Cell Type Enrichment Heatmap (without clustering)
===================================================
Input:  AnnData with columns 'celltype' and 'CNs_10' (or custom names)
Output: PDF/PNG heatmap using the style from the new script,
        but with no hierarchical clustering (rows/columns kept in order).
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List

# =============================================================
# 1. 计算富集矩阵（与两个原脚本算法一致）
# =============================================================
def compute_enrichment(adata, celltype_col='celltype', CN_col='CNs_10',
                       excluded_cells=None):
    """
    计算每个 CN 中细胞类型的 log2 富集分数（带平滑项）。
    返回 (fc_df, counts_df, global_props)
    """
    if excluded_cells is None:
        excluded_cells = []

    # 提取所有细胞类型并过滤
    all_cts = sorted(adata.obs[celltype_col].dropna().unique())
    ct_order = [ct for ct in all_cts if ct not in excluded_cells]

    # 提取所有 CN 并按原始顺序排序（数值排序）
    all_cns = sorted(adata.obs[CN_col].dropna().unique(),
                     key=lambda x: int(x.replace('CN', '')) if x.startswith('CN') else str(x))

    # 构建计数矩阵 (CN × celltype)
    obs = adata.obs[[CN_col, celltype_col]].dropna()
    counts = pd.crosstab(obs[CN_col], obs[celltype_col])
    counts = counts.reindex(index=all_cns, columns=ct_order, fill_value=0).values

    # 全局比例（每个细胞类型在所有细胞中的比例）
    global_counts = counts.sum(axis=0)
    global_props = global_counts / global_counts.sum()

    # 计算富集分数：log2( ((counts + global_props) / row_sum) / global_props )
    # 注意：这里将 global_props 作为平滑项，避免除零
    fc = np.log2(
        ((counts + global_props) /
         (counts + global_props).sum(axis=1, keepdims=True)) /
        global_props
    )

    fc_df = pd.DataFrame(fc, index=all_cns, columns=ct_order)
    counts_df = pd.DataFrame(counts, index=all_cns, columns=ct_order)
    return fc_df, counts_df, global_props


# =============================================================
# 2. 绘制无聚类热图（采用新脚本样式）
# =============================================================
def plot_heatmap_no_cluster(fc_df, vmin=-2, vmax=2, cmap='coolwarm',
                            linewidths=0, figsize=(12, 6),
                            xlabel_rotation=45, fontsize=6,
                            save_dir=None, filename='ct_enrichment'):
    """
    绘制热图，行列顺序保持原始，不进行任何聚类。
    参数与 new script 风格一致。
    """
    plt.figure(figsize=figsize)
    sns.set_style("white")

    # 使用 sns.heatmap（无聚类，直接绘制）
    ax = sns.heatmap(fc_df,
                     vmin=vmin, vmax=vmax,
                     cmap=cmap,
                     linewidths=linewidths,
                     cbar_kws={"shrink": 0.6, "aspect": 10})

    # 设置 x 轴标签：旋转，右对齐，小字体
    ax.set_xticks(np.arange(len(fc_df.columns)) + 0.5)
    ax.set_xticklabels(fc_df.columns, rotation=xlabel_rotation,
                       ha='right', fontsize=fontsize)

    # 设置 y 轴标签（右对齐，常规）
    ax.set_yticks(np.arange(len(fc_df.index)) + 0.5)
    ax.set_yticklabels(fc_df.index, rotation=0, fontsize=fontsize)

    # 调整布局，防止标签被裁剪
    plt.tight_layout()

    # 保存文件
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        for ext in ('pdf', 'png'):
            fpath = os.path.join(save_dir, f"{filename}.{ext}")
            plt.savefig(fpath, dpi=300, bbox_inches='tight')
        print(f"[INFO] Figures saved to: {save_dir}")

    plt.show()
    return ax


# =============================================================
# 3. 包装函数，一站式执行
# =============================================================
def ct_enrichment_heatmap_no_cluster(
    adata,
    celltype_col='celltype',
    CN_col='CNs_10',
    excluded_cells=None,
    vmin=-2, vmax=2,
    cmap='coolwarm',
    linewidths=0,
    figsize=(12, 6),
    xlabel_rotation=45,
    fontsize=6,
    save_dir=None,
    filename='ct_enrichment_no_cluster',
    return_data=False
):
    """
    主函数：计算富集并绘制无聚类热图（新脚本风格）。
    """
    fc_df, counts_df, global_props = compute_enrichment(
        adata, celltype_col, CN_col, excluded_cells
    )

    plot_heatmap_no_cluster(
        fc_df,
        vmin=vmin, vmax=vmax,
        cmap=cmap,
        linewidths=linewidths,
        figsize=figsize,
        xlabel_rotation=xlabel_rotation,
        fontsize=fontsize,
        save_dir=save_dir,
        filename=filename
    )

    if return_data:
        return fc_df, counts_df, global_props


# =============================================================
# 示例调用
# =============================================================
if __name__ == "__main__":
    # 读取数据（假定 adata 已存在）
    adata = sc.read('adata_with_CNs_10.h5ad')

    # 直接运行，输出到当前目录
    fc, counts, props = ct_enrichment_heatmap_no_cluster(
        adata,
        celltype_col='celltype',
        CN_col='CNs_10',
        excluded_cells=['Unknown'],   # 可选的排除细胞类型
        save_dir='./CT_heatmap_no_cluster',
        filename='ct_enrichment',
        return_data=True
    )
