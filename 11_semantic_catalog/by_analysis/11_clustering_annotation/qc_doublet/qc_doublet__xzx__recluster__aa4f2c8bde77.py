import os
import shutil
import random
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import streamlit as st
from PIL import Image
import celltypist
from celltypist import models
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

# import pandas as pd
import scipy.io
import anndata
from scipy.sparse import csr_matrix


import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA_platform/scVision5/tmp/.scanpy_cached"
sc.settings.set_figure_params(figsize=(5,5), facecolor='white')

def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

adata_counts = sc.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/SC2411133007-07/03-Cell_filter/adata_counts.h5ad')


def dim_reduction_cluster(_adata, method, n_pcs, n_neighbors, resolution, batch_effect, out_dir):
    _adata.layers["counts"] = _adata.X.copy()
    # data normalization and log10 transform
    sc.pp.normalize_total(_adata)
    sc.pp.log1p(_adata)
    _adata.raw = _adata
    # select highly variable gene and data scale 
    sc.pp.highly_variable_genes(_adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    sc.pl.highly_variable_genes(_adata, show=False)
    save_fig(plt, out_dir, file_name='highly_variable_genes', fig_size={'width':16, 'height':8})
    _adata = _adata[:, _adata.var['highly_variable']]
    sc.pp.regress_out(_adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(_adata, max_value=10)
    # PCA
    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pl.pca(_adata, color='sample', show=False)
    save_fig(plt, out_dir, file_name='pca_sample')
    sc.pl.pca(_adata, color='group', show=False)
    save_fig(plt, out_dir, file_name='pca_group')
    sc.pl.pca_variance_ratio(_adata, log=True, show=False)
    save_fig(plt, out_dir, file_name='pca_variance_ratio')

    if batch_effect == 'harmony':
        sce.pp.harmony_integrate(_adata, 'sample', adjusted_basis='X_pca')
    # cluster
    sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)    

    # Remove batch effect
    if batch_effect == 'bbknn':
        sce.pp.bbknn(_adata, batch_key='sample', n_pcs=n_pcs)  # running bbknn 1.3.6


    sc.tl.umap(_adata)
    if method == 'leiden':
        sc.tl.leiden(_adata, resolution=resolution, key_added='cluster')
    if method == 'louvain':
        sc.tl.louvain(_adata, resolution=resolution, key_added='cluster')



    # UMAP
    sc.pl.umap(_adata, 
               color='cluster',
               title='Cluster',
               legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=2, 
               show=False)
    save_fig(plt, out_dir, file_name='umap_dimplot')
    sc.pl.umap(_adata, color='sample', title='', show=False)
    save_fig(plt, out_dir, file_name='umap_sample_dimplot')
    sc.pl.umap(_adata, color='group', title='', show=False)
    save_fig(plt, out_dir, file_name='umap_group_dimplot')
      
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, 
                   color=['cluster', "total_counts", "n_genes_by_counts", "total_counts_mt"],
                   title=['Cluster', 'total_counts', 'n_genes_by_counts', "total_counts_mt"],
                   legend_loc='on data',
                   frameon=True,
                   legend_fontsize=9, 
                   legend_fontoutline=3, 
                   wspace=0.2, 
                   ncols=2, 
                   cmap='coolwarm', 
                   use_raw=True, 
                   show=False)
               
        plt.savefig(os.path.join(out_dir, 'cluster_total_counts_n_genes_by_counts.png'), facecolor='white', bbox_inches="tight", dpi=300)
    
    # write all cells features
    _adata.obs.to_csv(os.path.join(out_dir,'all_cells_features.xls'), sep='\t')
    
    # write cluster cell statistics result
    result = _adata.obs.groupby(['cluster', 'sample']).size().reset_index(name='number').pivot(index='cluster', columns='sample', values='number').reset_index()
    result.columns.name = None    
    result.to_csv(os.path.join(out_dir,'cluster_cells_statistics.xls'), sep='\t', index=None)
    
    # cluster correlation
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.correlation_matrix(_adata, 'cluster', show=False)
        plt.savefig(os.path.join(out_dir, 'cluster_correlation.png'), facecolor='white', bbox_inches="tight", dpi=300)

    groups = _adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["group"]==groups[i]], color="cluster", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_group_cluster', fig_size={'width':len(groups)*5, 'height':5})

    samples = _adata.obs['sample'].drop_duplicates().tolist()
    if len(samples) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
        fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
        for i in range(len(samples)):
            if i == len(samples) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["sample"]==samples[i]], color="cluster", title=samples[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_sample_cluster', fig_size={'width':len(samples)*5, 'height':5})
    _adata.write(os.path.join(out_dir, 'result.h5ad'))

    return _adata

outdir = '/mnt/NexData_file/home/xzx/project/SC2411133007-07/report'
out = os.path.join(outdir,'04-Dimension_reduction_cluster')
os.makedirs(out,exist_ok=True)

adata = dim_reduction_cluster(_adata=adata_counts,method='leiden',n_pcs=30,n_neighbors=20,resolution=0.5,batch_effect='harmony',out_dir=out)



def find_marker_gene(_adata, method, out_dir, logfc=1, fdr=0.05, pct=0.25):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
   
    sc.tl.rank_genes_groups(_adata, 'cluster', method=method, pts=True)
    sc.pl.rank_genes_groups(_adata, n_genes=25, sharey=False,show=False)
    save_fig(plt, out_dir, file_name='rank_genes_groups', fig_size={'width':20, 'height':12})
    
    result = _adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, f'cluster{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(_adata, group=c)
        tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)

        # enr_go = gp.enrichr(gene_list=list(tmp['names']),
        #         gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
        #         background='hsapiens_gene_ensembl',
        #         organism='Human', 
        #         cutoff=0.05,
        #         no_plot=True,
        #         outdir=out,
        #         verbose=False
        #         )

        # for f in os.listdir(out):
        #     if f.startswith('gseapy'):
        #         os.remove(os.path.join(out, f))

        sc.pl.rank_genes_groups_violin(_adata, groups=c, gene_names=list(tmp['names'].head(5)), jitter=0, size=0, show=False)
        save_fig(plt, out, file_name='top_5_genes_violin')

    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    # if st.session_state['pipe']['organism'] == 'human':
    #     gene_annot = os.path.join(annotDir, '9606_human_go_kegg_annotation.xls')
    # elif st.session_state['pipe']['organism'] == 'mouse':
    #     gene_annot = os.path.join(annotDir, '10090_mouse_go_kegg_annotation.xls')
    # gene_annot = pd.read_csv(gene_annot, sep='\t')
    # all_diff = pd.merge(all_diff, gene_annot, on='Symbol', how='left')
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'cluster_diff_statistics.xls'), sep='\t')

    sc.pl.rank_genes_groups_dotplot(
        _adata,
        n_genes=3,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=1,
        cmap='bwr',
        colorbar_title='log fold change',
        show=False
    )
    save_fig(plt, out_dir, file_name='dotplot_top3_marker_gene', fig_size={'width':20, 'height':8})

    sc.pl.rank_genes_groups_heatmap(
        _adata,
        n_genes=10,
        values_to_plot="logfoldchanges",
        min_logfoldchange=1,
        swap_axes=True,
        cmap='bwr',
        standard_scale='var',
        show_gene_labels=True,
        dendrogram=True,
        show=False
    )
    save_fig(plt, out_dir, file_name='heatmap_top10_marker_gene', fig_size={'width':26, 'height':20})


outdir = '/mnt/NexData_file/home/xzx/project/SC2411133007-07/report'
out = os.path.join(outdir,'05-Marker_gene_analysis')
os.makedirs(out,exist_ok=True)

find_marker_gene(_adata=adata,method='wilcoxon',out_dir=out)


def celltype_annotation(_adata, database, organism, out_dir, cell_type=None):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)

    if organism == 'mouse':
        _adata.var['gene'] = _adata.var.index
        _adata.var.index = [i.upper() for i in _adata.var.index]
            
    if database == 'PanglaoDB':
        markers = dc.get_resource('PanglaoDB')        
        if cell_type:
            annotation_dict = cell_type
        else:
            # markers = markers[(markers[organism]=='True')&(markers['canonical_marker']=='True')&(markers['organ']!='Pancreas')]
            markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]
            dc.run_ora(mat=_adata, net=markers, source='cell_type', target='genesymbol', min_n=3, verbose=True, use_raw=False)
            acts = dc.get_acts(_adata, obsm_key='ora_estimate')
            mean_enr = dc.summarize_acts(acts, groupby="cluster", min_std=1)
            annotation_dict = dc.assign_groups(mean_enr)
        _adata.obs['cell_type'] = [annotation_dict[clust] for clust in _adata.obs["cluster"]]
    if database == 'CellTypist':
        # models.download_models(force_update = True)
        # models.models_path
        # models.models_description()
        if cell_type:
            annotation_dict = cell_type
            _adata.obs['cell_type'] = [annotation_dict[clust] for clust in _adata.obs["cluster"]]
        else:
            if organism == 'human':
                predictions = celltypist.annotate(_adata, model = 'Immune_All_High.pkl', majority_voting = True)
            if organism == 'mouse':   
                predictions = celltypist.annotate(_adata, model = 'Adult_Mouse_Gut.pkl', majority_voting = True)
            
            _adata = predictions.to_adata()
            _adata.obs['cell_type'] = _adata.obs['majority_voting']
            df = _adata.obs[["cluster", "cell_type"]]
            grouped = df.groupby('cluster')['cell_type'].apply(list).to_dict()
            annotation_dict = {k:'|'.join(set(v)) for k, v in grouped.items()}

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(_adata, color='cell_type', show=False, title='', legend_loc='on data', legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2)
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    # with rc_context({'figure.figsize': (5, 5)}):
    #     sc.pl.tsne(_adata, color='cell_type', show=False, title='', legend_loc='on data', legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2)
    #     plt.savefig(os.path.join(out_dir, 'tsne_dimplot_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    groups = _adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["group"]==groups[i]], color="cell_type", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_group_celltype', fig_size={'width':len(groups)*5, 'height':5})

        # axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        # fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 6))
        # for i in range(len(groups)):
        #     if i == len(groups) - 1:
        #         legend = 'right margin'
        #     else:
        #         legend = None
        #     sc.pl.tsne(_adata[_adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
        # save_fig(plt, out_dir, file_name='tsne_split_group_celltype', fig_size={'width':len(groups)*5, 'height':6})

    result = _adata.obs.groupby(['cell_type', 'sample']).size().reset_index(name='number').pivot(index='cell_type', columns='sample', values='number').reset_index()
    result.columns.name = None
    # result['Total'] = result.sum(axis=1)
    # result['Percent'] = round((result['Total'] / sum(result['Total'])) * 100, 1)

    result.to_csv(os.path.join(out_dir, 'celltype_statistics.xls'), sep='\t', index=None)

    _adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')
    #plot_celltype_proportion(_adata, out_dir=out_dir,type='group')
    #plot_celltype_proportion(_adata, out_dir=out_dir,type='sample')
    # pd.DataFrame(_adata.obsm['X_tsne'], index=_adata.obs_names, columns=['tSNE_1', 'tSNE_2']).to_csv(os.path.join(out_dir, 'data_tsne.xls'), sep='\t')
    pd.DataFrame(_adata.obsm['X_umap'], index=_adata.obs_names, columns=['UMAP_1', 'UMAP_2']).to_csv(os.path.join(out_dir, 'data_umap.xls'), sep='\t')

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(_adata, color=["cluster","cell_type"], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    
    # with rc_context({'figure.figsize': (5, 5)}):
    #     sc.pl.tsne(_adata, color=["cluster","cell_type"], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
    #     plt.savefig(os.path.join(out_dir, 'tsne_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    if organism == 'mouse':
        _adata.var.index = _adata.var['gene']

    _adata.write(os.path.join(out_dir, 'result.h5ad'))

    return _adata, annotation_dict


outdir = '/mnt/NexData_file/home/xzx/project/SC2411133007-07/report'
out = os.path.join(outdir,'06-Cell_type_annotation')
os.makedirs(out,exist_ok=True)

celltype_annotation(_adata=adata,database='PanglaoDB',organism='human',out_dir=out)
