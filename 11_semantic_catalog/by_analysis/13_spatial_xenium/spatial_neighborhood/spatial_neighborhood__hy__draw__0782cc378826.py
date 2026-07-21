import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Tuple, List

adata = sc.read('adata_with_CNs_10.h5ad')


def _build_label_counts(adata, celltype_col, neighborhood_col,
                        ct_order, cn_order):
    """
    Build the (n_cns × n_cts) count matrix from AnnData.
    Uses pd.crosstab for fast vectorized counting.
    """
    obs = adata.obs[[celltype_col, neighborhood_col]].dropna()
    ct = pd.crosstab(obs[neighborhood_col], obs[celltype_col])

    # Align to requested order, fill missing with 0
    label_counts = pd.DataFrame(
        0, index=cn_order, columns=ct_order, dtype=float
    )
    for cn in cn_order:
        for c in ct_order:
            if cn in ct.index and c in ct.columns:
                label_counts.loc[cn, c] = ct.loc[cn, c]

    return label_counts.values, cn_order


def cal_heatmap(adata, celltype_col='singler', neighborhood_col='CNs_9',
                excluded_cells=None, cn_index=None):
    """
    Calculate cell type enrichment heatmap data from an AnnData object.
    [保持原有文档字符串不变]
    """
    if excluded_cells is None:
        excluded_cells = []

    # [保持原有函数体不变...]
    # ---------------------------------------------------------------
    # 1. Derive cell types and CN labels from AnnData
    # ---------------------------------------------------------------
    all_cts = sorted(adata.obs[celltype_col].dropna().unique())
    ct_order = [ct for ct in all_cts if ct not in excluded_cells]

    import re

    def _sort_key(x):
        """Sort by trailing numeric part (e.g. 'CN1'->1, 'CN10'->10)."""
        m = re.search(r'(\d+)$', str(x))
        if m:
            return (0, int(m.group(1)))
        try:
            return (0, int(x))
        except (ValueError, TypeError):
            return (1, str(x))

    all_cns = sorted(adata.obs[neighborhood_col].dropna().unique(),
                     key=_sort_key)

    n_cns = len(all_cns)
    n_cts = len(ct_order)

    # ---------------------------------------------------------------
    # 2. Build count matrix (n_cns × n_cts)
    # ---------------------------------------------------------------
    label_counts, cn_labels = _build_label_counts(
        adata, celltype_col, neighborhood_col, ct_order, all_cns
    )

    # ---------------------------------------------------------------
    # 3. Tissue-wide average proportions
    # ---------------------------------------------------------------
    tissue_avgs = label_counts.sum(axis=0)
    tissue_avgs = tissue_avgs / tissue_avgs.sum()

    # ---------------------------------------------------------------
    # 4. Log2 fold-change enrichment
    # ---------------------------------------------------------------
    fc = np.log2(
        ((label_counts + tissue_avgs)
         / (label_counts + tissue_avgs).sum(axis=1, keepdims=True))
        / tissue_avgs
    )

    if cn_index is not None:
        fc = pd.DataFrame(fc, columns=ct_order, index=cn_index)
    else:
        fc = pd.DataFrame(
            fc, columns=ct_order,
            index=[f'CN{i+1}' for i in range(n_cns)],
        )

    # ---------------------------------------------------------------
    # 5. CN frequency & cell frequency
    # ---------------------------------------------------------------
    cn_freq = label_counts / label_counts.sum(axis=1, keepdims=True)
    cn_cell_freq = label_counts.sum(axis=1) / label_counts.sum()

    return fc, cn_freq, cn_cell_freq


def draw_heatmap(
    fc,
    cn_freq,
    cn_cell_freq,
    n_cts,
    n_cns,
    filename=None,
    title=None,
    ax=None,
    cbar=True,
    cbar_ax=None,
    vert=False,
    standalone=False,
    is_top=False,
    is_bottom=False,
    color_cn=False,
    figsize=None,
):
    """
    Draw a cell type enrichment heatmap.
    [保持原有文档字符串不变]
    """
    # [保持原有函数体不变...]
    if ax is None:
        if figsize is None:
            figsize = (5, 10) if vert else (10, 5)
        fig, ax = plt.subplots(figsize=figsize)

    if not vert:
        sns.heatmap(
            fc, vmin=-4, vmax=4, cmap='bwr', linewidths=0.5,
            ax=ax, cbar=cbar, cbar_ax=cbar_ax,
            cbar_kws={'orientation': 'horizontal'},
        )
    else:
        sns.heatmap(
            fc.T, vmin=-4, vmax=4, cmap='bwr', linewidths=0.5,
            ax=ax, cbar=cbar, cbar_ax=cbar_ax,
            cbar_kws={'orientation': 'vertical'},
        )

    matx, maty = np.zeros((n_cns, n_cts)), np.zeros((n_cns, n_cts))
    for i in range(len(matx)):
        matx[i] += i
    for j in range(len(matx[0])):
        maty[:, j] += j

    if not vert:
        scatter = ax.scatter(
            0.5 + maty, 0.5 + matx,
            marker='o', s=cn_freq * 60, color='black',
        )
        legend = ax.legend(
            *scatter.legend_elements(
                **dict(prop='sizes', num=4, func=lambda s: s / 60)
            ),
            title='CT frequency per CN',
            loc='upper left', frameon=False, edgecolor='k',
            bbox_to_anchor=(1, 0),
        )
        if standalone or is_top:
            ax.annotate('Cell%', (-1.9, -0.4), annotation_clip=False)
            ax.annotate('CN', (n_cts + .5, -0.4), annotation_clip=False)
            if cbar_ax is not None:
                cbar_ax.tick_params(
                    top=True, bottom=False,
                    labeltop=True, labelbottom=False, labelsize=8,
                )
                cbar_ax.set_title('CN Enrichment Score')
        if standalone or is_bottom:
            ax.annotate('CT:', (-3, n_cns + 1), annotation_clip=False)
        else:
            ax.get_xaxis().set_visible(False)
            ax.get_legend().remove()

        if title:
            ax.set_title(
                title, rotation='vertical',
                x=-0.08, y=0.5, ha='right', va='center',
            )

        freq_range = max(cn_cell_freq) - min(cn_cell_freq)
        if freq_range > 0:
            alpha_list = [
                (i - min(cn_cell_freq)) / freq_range * 0.75 + 0.25
                for i in cn_cell_freq
            ]
        else:
            alpha_list = [1.0] * len(cn_cell_freq)
        for i, (x, y) in enumerate(
            zip(
                np.ones(len(cn_cell_freq)) * 0 - 1.9,
                0.8 + np.arange(len(cn_cell_freq)),
            )
        ):
            ax.annotate(
                round(cn_cell_freq[i], 3), (x, y),
                annotation_clip=False, alpha=alpha_list[i],
            )
    else:
        scatter = ax.scatter(
            0.5 + matx, 0.5 + maty,
            marker='o', s=cn_freq * 60, color='black',
        )
        legend = ax.legend(
            *scatter.legend_elements(
                **dict(prop='sizes', num=4, func=lambda s: s / 60)
            ),
            title='CT frequency per CN',
            loc='upper left', frameon=False, edgecolor='k',
            bbox_to_anchor=(1, 0),
        )

    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=30 if vert else 45, ha='right',
    )
    ax.yaxis.set_ticks_position('right')
    ax.tick_params(axis='y', pad=25)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    palette = list(sns.color_palette('bright'))
    for i in range(n_cns):
        c_idx = i % len(palette)
        c_alpha = 0.8 if color_cn else 0
        ax.annotate(
            f'{i + 1}', (n_cts + 0.8, i + 0.5),
            bbox={
                'boxstyle': 'circle',
                'color': palette[c_idx],
                'pad': 0.05, 'alpha': c_alpha,
            },
            ha='center', va='center',
            annotation_clip=False, size=10,
        )

    if filename:
        plt.savefig(filename, bbox_inches='tight')

    return ax


def ct_enrichment(
    adata,
    celltype_col: str = 'singler',
    neighborhood_col: str = 'CNs_9',
    excluded_cells: Optional[List[str]] = None,
    cn_index: Optional[List[str]] = None,
    title: Optional[str] = None,
    color_cn: bool = True,
    vert: bool = False,
    cluster_cols: bool = False,
    ax=None,
    cbar_ax=None,
    standalone: bool = True,
    is_top: bool = False,
    is_bottom: bool = False,
    out_dir: Optional[str] = None,
    return_data: bool = False,
    figsize: Optional[Tuple[float, float]] = None,
):
    """
    Perform cell type enrichment analysis across cellular neighborhoods
    and visualize results as a heatmap.
    [保持原有文档字符串不变]
    """
    if excluded_cells is None:
        excluded_cells = []

    # ---------------------------------------------------------------
    # 1. Input validation
    # ---------------------------------------------------------------
    required_cols = {
        'celltype_col': celltype_col,
        'neighborhood_col': neighborhood_col,
    }
    missing = [
        name for name, col in required_cols.items()
        if col not in adata.obs.columns
    ]
    if missing:
        raise ValueError(
            "Missing columns in adata.obs: "
            + ", ".join(f"{name}='{required_cols[name]}'" for name in missing)
        )

    # ---------------------------------------------------------------
    # 2. Calculate enrichment
    # ---------------------------------------------------------------
    fc, cn_freq, cn_cell_freq = cal_heatmap(
        adata,
        celltype_col=celltype_col,
        neighborhood_col=neighborhood_col,
        excluded_cells=excluded_cells,
        cn_index=cn_index,
    )

    if cluster_cols:
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, leaves_list
        dist = pdist(fc.T.values, metric='euclidean')
        link = linkage(dist, method='ward')
        col_order = leaves_list(link)
        
        fc = fc.iloc[:, col_order]
        cn_freq = cn_freq[:, col_order]

    n_cts = fc.shape[1]
    n_cns = fc.shape[0]

    # ---------------------------------------------------------------
    # 3. Plot heatmap
    # ---------------------------------------------------------------
    fig = None
    if ax is None:
        if figsize is None:
            figsize = (5, 10) if vert else (10, 5)
        fig, ax = plt.subplots(figsize=figsize)
        if standalone and not vert:
            cbar_ax = ax.inset_axes([0.0, 1.05, 1.0, 0.04])
    else:
        fig = ax.get_figure()

    draw_heatmap(
        fc=fc,
        cn_freq=cn_freq,
        cn_cell_freq=cn_cell_freq,
        n_cts=n_cts,
        n_cns=n_cns,
        title=title,
        ax=ax,
        cbar=True,
        cbar_ax=cbar_ax,
        vert=vert,
        standalone=standalone,
        is_top=is_top,
        is_bottom=is_bottom,
        color_cn=color_cn,
    )

    plt.tight_layout()

    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        fname = "ct_enrichment"
        for ext in ('pdf', 'png'):
            fpath = os.path.join(out_dir, f"{fname}.{ext}")
            fig.savefig(fpath, bbox_inches='tight', dpi=300)
        print(f"[INFO] Figures saved to: {out_dir}")

    plt.show()

    if return_data:
        return fc, cn_freq, cn_cell_freq


# ============================================================================
# 使用示例
# ============================================================================
if __name__ == "__main__":
    # 假设 adata 已经被加载（需要根据你的实际情况修改）
    # import scanpy as sc
    # adata = sc.read_h5ad('your_data.h5ad')
    
    # 示例1: 基础用法
    print("Running basic enrichment analysis...")
    fc, cn_freq, cn_cell_freq = ct_enrichment(
        adata,
        celltype_col='celltype',
        neighborhood_col='CNs_10',
        out_dir='./CT_enrichment',
        return_data=True,
        standalone=True,
        cluster_cols=True
    )
    
    # 打印结果摘要
    print(f"Enrichment matrix shape: {fc.shape}")
    print(f"Top enriched cell types in CN1:")
    print(fc.iloc[0].sort_values(ascending=False).head())
