from pathlib import Path
from urllib.request import urlretrieve
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import omicverse as ov
from omicverse.single import Monocle

ov.plot_set(font_path='Arial')
np.random.seed(42)

# %reload_ext autoreload
# %autoreload 2
import scanpy as sc
adata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/B_Plasma/B_Plasma.h5ad")
print(f'B cells: {adata.n_obs}')
adata.obs['cell_type'].value_counts()

mono = Monocle(adata)
mono.preprocess()
mono.select_ordering_genes(max_genes=500)
print(f"Ordering genes: {mono.adata.var['use_for_ordering'].sum()}")

import os, sys
out_dir = "/mnt/NexData_file/home/sh/Projects/SC2604012626-02/02.Mission/Bcell_monocle2/"
mono.plot_ordering_genes(figsize=(5, 4))
plt.savefig(os.path.join(out_dir, 'plot_ordering_genes.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'plot_ordering_genes.pdf'),bbox_inches='tight', dpi=300)

mono.reduce_dimension(max_components=2, verbose=False)
available_cols = ['cell_type','group']
for col in available_cols:
    print(f"Plotting trajectory: {col}")
    ov.pl.trajectory(
        mono.adata,
        method='monocle',
        basis='X_DDRTree',
        color=col,
        figsize=(7,5),
    )
    ov.plt.savefig(os.path.join(out_dir, f'plot_trajectory_{col}.pdf'),bbox_inches='tight', dpi=300)
    ov.plt.savefig(os.path.join(out_dir, f'plot_trajectory_{col}.png'),bbox_inches='tight', dpi=300)

mono.order_cells(root_by_column='cell_type', root_by_value='Naive B')
print(mono)
available_cols = ['cell_type','group','Pseudotime','State']
for col in available_cols:
    print(f"Plotting trajectory: {col}")
    ov.pl.trajectory(
        mono.adata,
        method='monocle',
        basis='X_DDRTree',
        color=col,
        figsize=(7,5),
    )
    ov.plt.savefig(os.path.join(out_dir, f'plot_trajectory_{col}.pdf'),bbox_inches='tight', dpi=300)
    ov.plt.savefig(os.path.join(out_dir, f'plot_trajectory_{col}.png'),bbox_inches='tight', dpi=300)

for col in available_cols:
    print(f"Plotting trajectory: {col}")
    ov.pl.trajectory_tree(
        mono.adata,
        method='monocle',
        color=col,
        figsize=(7,5),
    )
    ov.plt.savefig(os.path.join(out_dir, f'plot_trajectory_tree_{col}.pdf'),bbox_inches='tight', dpi=300)
    ov.plt.savefig(os.path.join(out_dir, f'plot_trajectory_tree_{col}.png'),bbox_inches='tight', dpi=300)


fig, ax = ov.pl.branch_streamplot(
    mono.adata,
    group_key='cell_type',
    pseudotime_key='Pseudotime',
    show=False,
)
ov.plt.savefig(os.path.join(out_dir, 'plot_branch_streamplot.pdf'),bbox_inches='tight', dpi=300)
ov.plt.savefig(os.path.join(out_dir, 'plot_branch_streamplot.png'),bbox_inches='tight', dpi=300)

ordering_genes = mono.adata.var_names[mono.adata.var['use_for_ordering']].tolist()
mono_ord = Monocle(mono.adata[:, ordering_genes].copy())

de = mono_ord.differential_gene_test(cores=-1)
sig = de[(de['qval'] < 0.01) & (de['status'] == 'OK')]

top30 = sig.sort_values('pval').head(30).index.tolist()

g = ov.pl.dynamic_heatmap(
    mono.adata,
    pseudotime='Pseudotime',
    var_names=top30,
    cell_annotation='State',
    use_cell_columns=False,
    use_fitted=True,
    cell_bins=200,
    figsize=(7, 7),
    show_row_names=True,
    standard_scale='var',
    cmap='RdBu_r',
    order_by='peak',
    show=False,
)
ov.plt.savefig(os.path.join(out_dir, 'plot_dynamic_heatmap.pdf'),bbox_inches='tight', dpi=300)
ov.plt.savefig(os.path.join(out_dir, 'plot_dynamic_heatmap.png'),bbox_inches='tight', dpi=300)


marker_genes = [g for g in ["Fao","S1pr2","Bcl6","Prdm1"] if g in mono.adata.var_names]
res = ov.single.dynamic_features(
    mono.adata,
    genes=marker_genes,
    pseudotime='Pseudotime',
    store_raw=True,
    raw_obs_keys=['cell_type'],
)

ov.pl.dynamic_trends(
    res,
    genes=marker_genes,
    add_point=True,
    point_color_by='cell_type',
    figsize=(7, 5),
    legend_loc='right margin',
    legend_fontsize=8,
)
ov.plt.savefig(os.path.join(out_dir, 'plot_gene_dynamic_trends.pdf'),bbox_inches='tight', dpi=300)
ov.plt.savefig(os.path.join(out_dir, 'plot_gene_dynamic_trends.png'),bbox_inches='tight', dpi=300)


olsson_global_dyn = ov.single.dynamic_features(
    mono.adata,
    genes=marker_genes,
    pseudotime='Pseudotime',
    use_raw=False,
    distribution='normal',
    link='identity',
    n_splines=8,
    store_raw=True,
    raw_obs_keys=['cell_type', 'group'],
)

ov.pl.dynamic_trends(
    olsson_global_dyn,
    genes=marker_genes,
    compare_features=True,
    add_point=True,
    point_color_by='cell_type',
    line_style_by='features',
    figsize=(6, 3.2),
    linewidth=2.2,
    legend_loc='right margin',
    legend_fontsize=8,
)
ov.plt.savefig(os.path.join(out_dir, 'plot_gene_dynamic_trends_all.pdf'),bbox_inches='tight', dpi=300)
ov.plt.savefig(os.path.join(out_dir, 'plot_gene_dynamic_trends_all.png'),bbox_inches='tight', dpi=300)


branch_subtypes = ['A1_', 'A2_']
branch_split_mask = mono.adata.obs['cell_type'].astype(str).isin(["Naive B","Memory B","Activate B"])

olsson_branch_dyn = ov.single.dynamic_features(
    mono.adata,
    genes=marker_genes,
    pseudotime='Pseudotime',
    groupby='group',
    groups=branch_subtypes,
    use_raw=False,
    distribution='normal',
    link='identity',
    n_splines=8,
    store_raw=True,
)
branch_split_time = float(np.nanmedian(mono.adata.obs.loc[branch_split_mask, 'Pseudotime'])) if branch_split_mask.any() else float(np.nanmedian(mono.adata.obs['Pseudotime']))

ov.pl.dynamic_trends(
    olsson_branch_dyn,
    genes=marker_genes,
    compare_groups=True,
    split_time=branch_split_time,
    shared_trunk=True,
    add_point=True,
    point_color_by='group',
    figsize=(4.6, 3),
    linewidth=2.2,
    ncols=3,
    legend_loc='right margin',
    legend_fontsize=8
    # title='Hematopoietic branch-aware marker trends',
)
ov.plt.savefig(os.path.join(out_dir, 'plot_gene_dynamic_trends_group.pdf'),bbox_inches='tight', dpi=300)
ov.plt.savefig(os.path.join(out_dir, 'plot_gene_dynamic_trends_group.png'),bbox_inches='tight', dpi=300)