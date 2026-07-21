import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")


def data_qc_and_filter_by_quantile(adata, qt_min=0.1, qt_max=0.99, pct_mt_max=20, log10GenesPerUMI_min=0.75, log10GenesPerUMI_max=0.99):
    # Data QC
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-') 
    adata.var['ribo'] = adata.var_names.str.upper().str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes_by_counts']) / np.log10(adata.obs['total_counts'])
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'log10GenesPerUMI'], groupby='sample', jitter=False, rotation=90, multi_panel=True, show=False)

    # Cell Filtering
    adata_list = []
    for sample in adata.obs['sample'].drop_duplicates():
        tmp = adata[adata.obs['sample']==sample]  
        min_counts = int(np.quantile(tmp.obs['total_counts'], qt_min))
        if min_counts < 30:
            min_counts = 30
        elif min_counts > 50:
            min_counts = 50
        max_counts = int(np.quantile(tmp.obs['total_counts'], qt_max))
        tmp = tmp[(tmp.obs['total_counts'] > min_counts) & (tmp.obs['total_counts'] < max_counts)]
        tmp = tmp[(tmp.obs['log10GenesPerUMI'] > log10GenesPerUMI_min) & (tmp.obs['log10GenesPerUMI'] < log10GenesPerUMI_max)]
        tmp = tmp[tmp.obs['pct_counts_mt'] < pct_mt_max]
        sc.pp.filter_genes(tmp, min_cells=10)
        adata_list.append(tmp)
        print(f"{sample} -- min_counts: {min_counts}  max_counts: {max_counts}; min_genes: {tmp.obs['n_genes_by_counts'].min()}; max_genes: {tmp.obs['n_genes_by_counts'].max()}")
    adata_filter = ad.concat(adata_list, join='inner', merge='same', uns_merge='same', index_unique=None)
    sc.pl.violin(adata_filter, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'log10GenesPerUMI'], groupby='sample', jitter=False, rotation=90, multi_panel=True, show=False)

    return adata_filter


def dim_reduction_cluster(adata, n_hvg=2000, n_neighbors=20, n_pcs=30, res=0.8, batch_remove='harmony'):
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
        sc.external.pp.bbknn(adata, batch_key="sample", n_pcs=n_pcs)
    else:
        sc.pp.neighbors(adata, metric='cosine', n_neighbors=n_neighbors, n_pcs=n_pcs) 

    sc.tl.umap(adata, min_dist=0.3, spread=1.5)
    sc.tl.leiden(adata, resolution=res, key_added=f'res_{res}')

    if 120000/adata.obs.shape[0] < 1:
        size = (120000/adata.obs.shape[0]) * 5
    else:
        size = 120000/adata.obs.shape[0] * 1.5

    sc.pl.umap(
        adata, 
        color=['sample', f'res_{res}'],
        title=['Sample', 'Cluster'],
        # legend_loc='right',
        size=size,
        legend_fontsize=10,
        legend_fontoutline=2, 
        frameon=False,
        show=False
    )

    groups = adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
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

    return adata


def find_all_markers(adata, groupby='cluster', cmap='Spectral_r', fc=0.25, output_dir=None):       
    sc.tl.rank_genes_groups(
        adata, 
        groupby=groupby, 
        use_raw=True, 
        method="t-test", # method='t-test','wilcoxon',key_added='rank_genes_groups' 
        pts=True
    )
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=5,
        values_to_plot="logfoldchanges",
        # min_logfoldchange=1,
        vmax=1,
        vmin=-1,
        cmap=cmap,
    )
    # result = adata.uns['rank_genes_groups']
    # groups = result['names'].dtype.names
    # all_diff = []
    # for c in groups:
    #     out = os.path.join(output_dir, f'cluster{c}')
    #     os.makedirs(out, exist_ok=True)
    #     tmp = sc.get.rank_genes_groups_df(adata, group=c)
    #     tmp = tmp[(np.abs(tmp['logfoldchanges'])>fc) & (tmp['pvals_adj']<0.05)]
    #     tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
    #     tmp['cluster'] = c
    #     all_diff.append(tmp)

    # all_diff = pd.concat(all_diff, axis=0)
    # all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    # all_diff.to_csv(os.path.join(output_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    # all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(output_dir, 'cluster_diff_statistics.xls'), sep='\t')

    return adata



















