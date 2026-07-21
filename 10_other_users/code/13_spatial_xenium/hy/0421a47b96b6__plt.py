"""
CT Enrichment Analysis
======================
Usage:
    # Basic usage
    fc, cn_freq, cn_cell_freq = ct_enrichment(adata, return_data=True)

    # With custom columns and saving
    fc, cn_freq, cn_cell_freq = ct_enrichment(
        adata,
        celltype_col='celltype',
        neighborhood_col='CNs_10',
        excluded_cells=['Unknown'],
        out_dir='./results/',
        return_data=True,
    )
"""

import os
import re
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List

adata = sc.read('adata_with_CNs_10.h5ad')


def _build_label_counts(adata, celltype_col, neighborhood_col,
                        ct_order, cn_order):
    """
    Build the (n_cns × n_cts) count matrix from AnnData.
    Uses pd.crosstab for fast vectorized counting.
    """
    obs = adata.obs[[celltype_col, neighborhood_col]].dropna()
    ct = pd.crosstab(obs[neighborhood_col], obs[celltype_col])

    label_counts = pd.DataFrame(
        0, index=cn_order, columns=ct_order, dtype=float
    )
    for cn in cn_order:
        for c in ct_order:
            if cn in ct.index and c in ct.columns:
                label_counts.loc[cn, c] = ct.loc[cn, c]

    return label_counts.values


def cal_heatmap(adata, celltype_col='celltype', neighborhood_col='CNs_10',
                excluded_cells=None, cn_index=None):
    """
    Calculate cell type enrichment heatmap data from an AnnData object.

    The enrichment is computed as:
        fc = log2( ((counts + tissue_avg) / row_sum) / tissue_avg )
    where tissue_avg is the global cell-type proportion.

    Parameters
    ----------
    adata : AnnData
    celltype_col : str, default 'celltype'
    neighborhood_col : str, default 'CNs_10'
    excluded_cells : list of str, optional
    cn_index : list of str, optional
        If None, CNs are sorted naturally (CN0, CN1, ..., CN10, ...).

    Returns
    -------
    fc : pd.DataFrame
    cn_freq : np.ndarray
    cn_cell_freq : np.ndarray
    """
    if excluded_cells is None:
        excluded_cells = []

    all_cts = sorted(adata.obs[celltype_col].dropna().unique())
    ct_order = [ct for ct in all_cts if ct not in excluded_cells]

    # Natural sort CN labels by the numeric part
    raw_cns = adata.obs[neighborhood_col].dropna().unique().tolist()
    def _cn_sort_key(label):
        match = re.search(r'\d+', str(label))
        return int(match.group()) if match else float('inf')
    all_cns = sorted(raw_cns, key=_cn_sort_key)

    n_cns = len(all_cns)
    n_cts = len(ct_order)

    label_counts = _build_label_counts(
        adata, celltype_col, neighborhood_col, ct_order, all_cns
    )

    tissue_avgs = label_counts.sum(axis=0)
    tissue_avgs = tissue_avgs / tissue_avgs.sum()

    fc = np.log2(
        ((label_counts + tissue_avgs)
         / (label_counts + tissue_avgs).sum(axis=1, keepdims=True))
        / tissue_avgs
    )
    print(fc.max(), fc.min())

    if cn_index is not None:
        fc = pd.DataFrame(fc, columns=ct_order, index=cn_index)
    else:
        fc = pd.DataFrame(fc, columns=ct_order, index=all_cns)

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
):
    """
    Draw a cell type enrichment heatmap.

    Parameters
    ----------
    fc : pd.DataFrame
    cn_freq : np.ndarray
    cn_cell_freq : np.ndarray
    n_cts : int
    n_cns : int
    filename : str, optional
    title : str, optional
    ax : matplotlib.axes.Axes, optional
    cbar : bool, default True
    cbar_ax : matplotlib.axes.Axes, optional
    vert : bool, default False
    standalone : bool, default False
    is_top : bool, default False
    is_bottom : bool, default False
    color_cn : bool, default False

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        if not vert:
            fig, ax = plt.subplots(figsize=(10, 5))
        else:
            fig, ax = plt.subplots(figsize=(5, 10))

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
            ax.annotate('Cell%', (-2.5, -0.4), annotation_clip=False)
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

        # Left cell frequency column
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
                np.ones(len(cn_cell_freq)) * -2.5,
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

    # Right-side colored circles with CN labels — made even smaller
    palette = list(sns.color_palette('bright'))
    for i in range(1, n_cns + 1):
        c_idx = (i - 1) % len(palette)
        c_alpha = 0.8 if color_cn else 0
        ax.annotate(
            f'{fc.index[i-1]}',
            (n_cts + 0.8, i - 0.5),
            bbox={
                'boxstyle': 'circle',
                'color': palette[c_idx],
                'pad': 0.01,          # tighter padding → smaller circle
                'alpha': c_alpha,
            },
            ha='center', va='center',
            annotation_clip=False, size=4,   # smaller font
        )

    if filename:
        plt.savefig(filename, bbox_inches='tight')

    return ax


def ct_enrichment(
    adata,
    celltype_col: str = 'celltype',
    neighborhood_col: str = 'CNs_10',
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
):
    """
    Perform cell type enrichment analysis across cellular neighborhoods
    and visualize results as a heatmap.

    Parameters
    ----------
    adata : AnnData
    celltype_col : str, default 'celltype'
    neighborhood_col : str, default 'CNs_10'
    excluded_cells : list of str, optional
    cn_index : list of str, optional
    title : str, optional
    color_cn : bool, default True
    vert : bool, default False
    cluster_cols : bool, default False
    ax : matplotlib.axes.Axes, optional
    cbar_ax : matplotlib.axes.Axes, optional
    standalone : bool, default True
    is_top : bool, default False
    is_bottom : bool, default False
    out_dir : str or None, default None
    return_data : bool, default False

    Returns
    -------
    fc : pd.DataFrame (if return_data=True)
    cn_freq : np.ndarray (if return_data=True)
    cn_cell_freq : np.ndarray (if return_data=True)
    """
    if excluded_cells is None:
        excluded_cells = []

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

    fig = None
    if ax is None:
        if standalone and not vert:
            fig, ax = plt.subplots(figsize=(10, 5))
            cbar_ax = ax.inset_axes([0.0, 1.05, 1.0, 0.04])
        elif not vert:
            fig, ax = plt.subplots(figsize=(10, 5))
        else:
            fig, ax = plt.subplots(figsize=(5, 10))
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


fc, cn_freq, cn_cell_freq = ct_enrichment(
    adata,
    celltype_col='celltype',
    neighborhood_col='CNs_10',
    out_dir='./CT_enrichment',
    return_data=True,
    standalone=True,
    cluster_cols=True
)
