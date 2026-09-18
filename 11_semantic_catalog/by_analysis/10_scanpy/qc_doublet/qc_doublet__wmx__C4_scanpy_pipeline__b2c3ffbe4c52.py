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

@st.fragment
def make_report_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for i in ['01-Dnbc4tools', '02-Data_QC', '03-Cell_filter', '04-Dimension_reduction_cluster', '05-Marker_gene_analysis', '06-Cell_type_annotation']:
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)
    return out_dir



def read_anndata_C4(path):
    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat ,dtype="float32")
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[0] .values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    adata.obs_names = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.var_names_make_unique()
    return adata



@st.fragment
def get_C4_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    for i in os.listdir(data_path):
        name = i.split('.')[0]
        for k,v in sample_dict.items():
            if name in v:
                path = os.path.join(data_path, name, 'output', 'filter_matrix')
                data = read_anndata_C4(path)
                # data = sc.read_10x_mtx(os.path.join(data_path, name, 'outs', 'filtered_feature_bc_matrix'), var_names='gene_symbols', cache=True)
                data.var_names_make_unique()
                data.obs['sample'] = name
                data.obs['group'] = k
                adatas_list.append(data)
                os.makedirs(os.path.join(out_dir, "01-Dnbc4tools", name), exist_ok=True)
                os.system(f'cp {os.path.join(data_path, name, "output", "metrics_summary.xls")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                # os.system(f'cp {os.path.join(data_path, name, "output", "cloupe.cloupe")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                os.system(f'cp {os.path.join(data_path, name, "output", name+"_scRNA_report.html")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                os.system(f'cp -r {os.path.join(data_path, name, "output", "filter_matrix")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
    adata = adatas_list[0].concatenate(adatas_list[1:], index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    del adata.obs['batch']

    return adata

def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

@st.fragment
def data_qc(_adata, out_dir):   
    _adata.var['mt'] = _adata.var_names.str.upper().str.startswith('MT-') 
    _adata.var['ribo'] = _adata.var_names.str.upper().str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(_adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
    axes = plt.gcf().axes
    for ax in axes:
        ax.set_xticklabels(ax.get_xticklabels(), horizontalalignment='right')
    save_fig(plt, out_dir, file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(_adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
    sc.pl.scatter(_adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
    save_fig(fig, out_dir, file_name='orig_nCount_nFeature_scatter')
    
    pd.DataFrame(_adata.obs).to_csv(os.path.join(out_dir, 'all_cell_features_raw.xls'), sep='\t')

    for s in _adata.obs['sample'].drop_duplicates().tolist():
        os.makedirs(os.path.join(out_dir, s), exist_ok=True)
        sc.pl.violin(_adata[_adata.obs['sample']==s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
        save_fig(plt, os.path.join(out_dir, s), file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, os.path.join(out_dir, s), file_name='orig_nCount_nFeature_scatter')
        pd.DataFrame(_adata[_adata.obs['sample']==s].obs).to_csv(os.path.join(out_dir, s, 'cell_features_raw.xls'), sep='\t')
    
    return _adata

@st.fragment
def cell_filtering(_adata, min_genes, min_cells, pct_counts_mt, pct_counts_ribo, upper_lim, total_counts, expected_doublet_rate, threshold, out_dir):    
    sc.pp.filter_cells(_adata, min_genes=min_genes)
    sc.pp.filter_genes(_adata, min_cells=min_cells)
    _adata = _adata[(_adata.obs.n_genes_by_counts < upper_lim)]
    _adata = _adata[_adata.obs['pct_counts_mt'] < pct_counts_mt, :]
    _adata = _adata[_adata.obs['pct_counts_ribo'] < pct_counts_ribo, :]
    _adata = _adata[_adata.obs['total_counts'] < total_counts, :]
    
    sc.external.pp.scrublet(_adata, expected_doublet_rate=expected_doublet_rate / 100, threshold=threshold, batch_key='sample', random_state=1234)
    
    samples = _adata.obs['sample'].drop_duplicates().tolist()
    if len(samples) > 1:
        sc.external.pl.scrublet_score_distribution(_adata)
        save_fig(plt, out_dir, file_name='scrublet_score_distribution', fig_size={'width':8, 'height':6})  
    _adata = _adata[_adata.obs['predicted_doublet']==False, :]

    pd.DataFrame(_adata.obs).to_csv(os.path.join(out_dir, 'all_cell_features_filtered.xls'), sep='\t')
    sc.pl.violin(_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby = 'sample', rotation= 45, multi_panel=True, show=False)
    save_fig(plt, out_dir, file_name='filtered_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})  
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(_adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
    sc.pl.scatter(_adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
    save_fig(fig, out_dir, file_name='filtered_nCount_nFeature_scatter')      
    
    for s in _adata.obs['sample'].drop_duplicates().tolist():
        os.makedirs(os.path.join(out_dir, s), exist_ok=True)
        sc.pl.violin(_adata[_adata.obs['sample']==s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby = 'sample', rotation= 45, multi_panel=True, show=False)
        save_fig(plt, os.path.join(out_dir, s), file_name='filtered_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(_adata[_adata.obs['sample']==s], x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, os.path.join(out_dir, s), file_name='filtered_nCount_nFeature_scatter')       
        pd.DataFrame(_adata[_adata.obs['sample']==s].obs).to_csv(os.path.join(out_dir, s, 'cell_features_filtered.xls'), sep='\t')
    _adata.write(os.path.join(out_dir, 'adata_counts.h5ad'))

    return _adata

@st.fragment
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
    
    # cluster
    if batch_effect == '**harmony**':
        sce.pp.harmony_integrate(_adata, 'sample', adjusted_basis='X_pca')
    sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # Remove batch effect
    if batch_effect == '**bbknn**':
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
    
    # TSNE
    # sc.tl.tsne(_adata)
    # sc.pl.tsne(_adata, 
    #            color='cluster', 
    #            title='Cluster', 
    #            legend_loc='on data', 
    #            legend_fontsize=12, 
    #            legend_fontoutline=2, 
    #            show=False)
    # save_fig(plt, out_dir, file_name='tsne_dimplot')   
    # sc.pl.tsne(_adata, color='sample', title='', show=False)
    # save_fig(plt, out_dir, file_name='tsne_sample_dimplot')
    # sc.pl.tsne(_adata, color='group', title='', show=False)
    # save_fig(plt, out_dir, file_name='tsne_group_dimplot')
    
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



    return _adata

@st.fragment
def find_marker_gene(_adata, method, out_dir, logfc=1, fdr=0.05, pct=0.25):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
   
    sc.tl.rank_genes_groups(_adata, 'cluster', method=method, pts=True)
    sc.pl.rank_genes_groups(_adata, n_genes=25, sharey=False)
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

    #_adata.write(os.path.join(out_dir, 'adata_cluster.h5ad'))

    return _adata

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
    os.system(f"Rscript /mnt/NexData_file/share/scRNA_platform/sc_C4/moudle/Rscript/cell_proportion.R {os.path.join(out_dir, 'all_cells_features.xls')} cell_type {out_dir}")

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

    #_adata.write(os.path.join(out_dir, 'result.h5ad'))

    return _adata, annotation_dict

def visua_gene_express(_adata, genes, out_dir):
    st.set_option('deprecation.showPyplotGlobalUse', False)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color=genes, s=15, ncols=3, cmap='coolwarm', use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, '_'.join(genes)+'.png'), facecolor='white', bbox_inches="tight", dpi=200)
    gene_expression = Image.open(os.path.join(out_dir, '_'.join(genes)+'.png'))

    if len(genes) > 1:
        fig, axs = plt.subplots(len(genes), 1, figsize=(10, len(genes)*4), sharey=True)
        for i in range(len(genes)):
            sc.pl.violin(_adata, genes[i], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=axs[i], show=False)
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.violin(_adata, genes[0], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=ax, show=False)

    return gene_expression, fig


def visua_gene_express(_adata, genes, out_dir):
    #st.set_option('deprecation.showPyplotGlobalUse', False)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color=genes, s=15, ncols=3, cmap='coolwarm', use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, '_'.join(genes)+'.png'), facecolor='white', bbox_inches="tight", dpi=200)
    gene_expression = Image.open(os.path.join(out_dir, '_'.join(genes)+'.png'))

    if len(genes) > 1:
        fig, axs = plt.subplots(len(genes), 1, figsize=(10, len(genes)*4), sharey=True)
        for i in range(len(genes)):
            sc.pl.violin(_adata, genes[i], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=axs[i], show=False)
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.violin(_adata, genes[0], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=ax, show=False)

    return gene_expression, fig