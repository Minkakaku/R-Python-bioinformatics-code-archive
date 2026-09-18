import os
import random
import argparse
import cupy as cp
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import rapids_singlecell as rsc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import warnings
warnings.filterwarnings("ignore")
plt.style.use('default')  # 使用默认样式
plt.rcParams['figure.facecolor'] = 'white'  # 确保图像背景为白色
sc.settings.set_figure_params(dpi=100, dpi_save=350, figsize=(6, 6), facecolor='white', vector_friendly=True)

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False, # default is False
    devices=0,            # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)
random.seed(0)


def cluster_gpu(adata_h5ad, n_hvg=2000, n_neighbors=20, n_pcs=15, res=1, out_dir=None):
    cluster_col_name = f'res_{res:g}'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    adata = sc.read_h5ad(adata_h5ad)
    rsc.get.anndata_to_GPU(adata, convert_all=True)
    rsc.pp.normalize_total(adata, target_sum=10000, inplace=True)
    rsc.pp.log1p(adata)
    adata.raw = adata.copy()
    # HVG
    rsc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg,
        flavor="seurat_v3",
        layer='counts',
        batch_key="sample",
    )
    adata = adata[:, adata.var['highly_variable']]
    rsc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    rsc.pp.scale(adata, max_value=10)
    rsc.pp.pca(adata)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=30)
    rsc.pp.harmony_integrate(adata, key="sample")
    # Cluster
    rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric='cosine', use_rep='X_pca_harmony', random_state=10)
    rsc.tl.umap(adata, min_dist=0.5, spread=1.5, random_state=10)
    rsc.tl.leiden(adata, resolution=res, key_added=cluster_col_name, random_state=10)
    rsc.get.anndata_to_CPU(adata, convert_all=True)

    if 120000/adata.obs.shape[0] < 1:
        size = (120000/adata.obs.shape[0]) * 5
    else:
        size = 120000/adata.obs.shape[0]
        
    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.umap(
            adata, 
            color=['sample', cluster_col_name],
            legend_loc='on data', 
            legend_fontsize=9, 
            legend_fontoutline=3,
            frameon=True,
            size=size,
            show=False
        )
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_cluster.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

    library_names = adata.obs['sample'].drop_duplicates().to_list()
    for library in library_names:
        with rc_context({'figure.figsize': (8, 8)}):
            sc.pl.spatial(
                adata[adata.obs['sample'] == library,:],
                color=cluster_col_name,
                library_id=library,
                title=f'{library}',
                size=1.5,
                alpha_img=0.1,
                marker='s',
                show=False
            )
            plt.savefig(os.path.join(out_dir, f'spatial_{library}_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)
            plt.savefig(os.path.join(out_dir, f'spatial_{library}_cluster.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

    groups = adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(
                adata[adata.obs["group"]==groups[i]], 
                color=cluster_col_name,
                title=groups[i], 
                ax=axs[i], 
                wspace=0.1, 
                size=8,
                legend_loc='on data', 
                legend_fontsize=12, 
                legend_fontweight='normal', 
                legend_fontoutline=2, 
                show=False
            )
        plt.savefig(os.path.join(out_dir, 'umap_split_group_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, 'umap_split_group_cluster.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

    result = adata.obs.groupby([cluster_col_name, 'sample']).size().reset_index(name='number').pivot(index=cluster_col_name, columns='sample', values='number').reset_index()
    result.columns.name = None
    # result['Total'] = result.sum(axis=1)
    # result['Percent'] = round((result['Total'] / sum(result['Total'])) * 100, 1)
    result.to_csv(os.path.join(out_dir, 'cluster_statistics.xls'), sep='\t', index=None)
    adata.obs.to_csv(os.path.join(out_dir, 'all_cluster_cells_features.xls'), sep='\t')
    pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP_1', 'UMAP_2']).to_csv(os.path.join(out_dir, 'data_umap.xls'), sep='\t')

    ct_by_sample=pd.crosstab(adata.obs['sample'], adata.obs[cluster_col_name])
    cts=ct_by_sample.div(ct_by_sample.sum(axis=1),axis=0)
    cls_dic=dict(zip(np.unique(adata.obs[cluster_col_name]), adata.uns[f'{cluster_col_name}_colors']))

    cts.plot(
        kind='bar',
        stacked=True,
        color=cts.columns.map(cls_dic),
        legend=False,
        figsize=(len(cts.index)*0.75, 4.5),
        width=0.98
    )
    plt.grid(False)
    plt.title('Cluster cell proportions')
    plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0, ncol=2)
    plt.savefig(os.path.join(out_dir, 'cluster_cell_proportion.png'), facecolor='white', bbox_inches="tight")
    plt.savefig(os.path.join(out_dir, 'cluster_cell_proportion.pdf'), facecolor='white', bbox_inches="tight")
    adata.write(os.path.join(Path(out_dir).parent, 'adata_cluster.h5ad'))
    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform reduction and cluster analysis")
    parser.add_argument("--adata_h5ad", type=str, help="Adata to the reduction and cluster")
    parser.add_argument("--n_hvg", type=int, default=3000, required=True, help="Number of higher variance gene")
    parser.add_argument("--n_neighbors", type=int, default=20, required=True, help="Number of neighbors")
    parser.add_argument("--n_pcs", type=int, default=15, required=True, help="Number of PCs")
    parser.add_argument("--resolution", type=float, default=1, help="Leiden cluster resolution (e.g., 0.3 0.5 0.8, 1)")  
    parser.add_argument("--output_dir", type=str, required=False, help="Path to the output directory")
    args = parser.parse_args()

    adata = cluster_gpu(
        adata_h5ad=args.adata_h5ad, 
        n_hvg=args.n_hvg, 
        n_neighbors=args.n_neighbors, 
        n_pcs=args.n_pcs, 
        res=args.resolution,
        out_dir=args.output_dir
    )
