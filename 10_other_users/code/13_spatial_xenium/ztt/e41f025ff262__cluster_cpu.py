import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.set_figure_params(dpi=100, dpi_save=200, figsize=(5, 5), facecolor='white')


def dim_reduction_cluster(adata, n_hvg=2000, n_neighbors=20, n_pcs=30, res=0.8, batch_remove='harmony'):
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir, exist_ok=True)

    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, flavor='seurat', batch_key='sample', n_top_genes=n_hvg)
    sc.pl.highly_variable_genes(adata, show=False)
    adata = adata[:, adata.var['highly_variable']]

    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True, show=False)

    # 去批次，使用harmony或者bbknn
    if batch_remove == 'harmony':
        sc.external.pp.harmony_integrate(adata, 'sample', adjusted_basis='X_pca')
        sc.pp.neighbors(adata, metric='cosine', n_neighbors=n_neighbors, n_pcs=n_pcs) 
    elif batch_remove == 'bbknn':
        sc.external.pp.bbknn(adata, batch_key="sample", n_pcs=30)

    sc.tl.umap(adata, min_dist=0.3, spread=1.5)
    sc.tl.leiden(adata, resolution=res, key_added=f'res_{res}')

    if 120000/adata.obs.shape[0] < 1:
        size = (120000/adata.obs.shape[0]) * 5
    else:
        size = 120000/adata.obs.shape[0] * 1.5

    fig_umap, axs = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    sc.pl.umap(
        adata, 
        color=[f'res_{res}'],
        title=['Cluster'],
        # legend_loc='right',
        size=size,
        legend_fontsize=10,
        legend_fontoutline=2, 
        frameon=False,
        show=False,
        ax=axs
    )
    # plt.savefig(os.path.join(output_dir, 'dimplot_umap_cluster.png'), bbox_inches = 'tight', pad_inches = 0.1, facecolor='white')
    # plt.savefig(os.path.join(output_dir, 'dimplot_umap_cluster.pdf'), bbox_inches = 'tight', pad_inches = 0.1, facecolor='white')

    groups = adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig_umap_group, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            sc.pl.umap(
                adata[adata.obs["group"]==groups[i]], 
                color=f'res_{res}',
                title=groups[i], 
                ax=axs[i], 
                wspace=0.1, 
                size=size,
                legend_loc='on data', 
                legend_fontsize=12, 
                legend_fontweight='normal', 
                legend_fontoutline=2, 
                frameon=False,
                show=False
            )
        # plt.savefig(os.path.join(output_dir, 'dimplot_umap_cluster_group.png'), bbox_inches = 'tight', pad_inches = 0.1, facecolor='white')
        # plt.savefig(os.path.join(output_dir, 'dimplot_umap_cluster_group.pdf'), bbox_inches = 'tight', pad_inches = 0.1, facecolor='white')
        
    return fig_umap, fig_umap_group, adata