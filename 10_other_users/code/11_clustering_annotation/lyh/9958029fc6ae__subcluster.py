import os
import math
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import cosg as cosg
from PIL import Image
import seaborn as sns
import streamlit as st
import gseapy as gp
from gseapy.plot import barplot, dotplot
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import warnings
warnings.filterwarnings("ignore")

sc.settings.verbosity = 1
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, dpi_save=300, figsize=(5, 5), facecolor='white')

plt.style.use('seaborn-white')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=12)
plt.rcParams['pdf.fonttype'] = 42


def make_report_dir(out_dir):   
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

def image_view(_container, png_file, caption='', resize=None):
    image = Image.open(png_file)
    if resize:
        new_image = image.resize(resize)
        _container.image(new_image, caption=caption)
    else:
        _container.image(image, caption=caption, use_column_width=True)

def sub_cluster(_adata, n_pcs, n_neighbors, resolution, organism, dot_size, out_dir='./'):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # sc.pp.normalize_total(_adata)
    # sc.pp.log1p(_adata)

    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # UMAP
    sc.tl.umap(_adata)
    sc.tl.leiden(_adata, resolution=resolution)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color='leiden', title='Cluster', size=dot_size, show=False)
        plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)

    # TSNE
    sc.tl.tsne(_adata)
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.tsne(_adata, color='leiden', title='Cluster', size=dot_size, show=False)
        plt.savefig(os.path.join(out_dir, 'dimplot_tsne_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)

    # Cluster correlation
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.correlation_matrix(_adata, 'leiden', show=False)
        plt.savefig(os.path.join(out_dir, 'cluster_correlation.png'), facecolor='white', bbox_inches="tight", dpi=300)


    method='wilcoxon'
    logfc=1
    fdr=0.05
    pct=0.25
    # Find marker
    sc.tl.rank_genes_groups(_adata, "leiden", method='wilcoxon', pts=True)
    result = _adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, 'cluster_marker_gene', f'cluster{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(_adata, group=c)
        tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)

        try:
            if organism != 'No':
                gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023']
                if organism == 'mouse':
                    gene_sets=['KEGG_2019_Mouse', 'GO_Biological_Process_2023']

                enr = gp.enrichr(gene_list=tmp['names'],
                                gene_sets=gene_sets,
                                organism=organism,
                                #  background=_adata.var_names.tolist(),
                                outdir=None, # don't write to disk
                                )

                enr_kegg = enr.results.loc[enr.results["Gene_set"] == gene_sets[0],]
                enr_kegg[enr_kegg['Adjusted P-value']<0.05].to_csv(os.path.join(out, 'KEGG_Enrichment.xls'), sep='\t', index=None)
                barplot(enr_kegg,
                        title='KEGG Enrichment', 
                        cutoff=0.05,
                        color='darkred',
                        top_term=20, 
                        figsize=(6, 8),
                        ofname=os.path.join(out, 'KEGG_Enrichment.png'))

                enr_bp = enr.results.loc[enr.results["Gene_set"] == gene_sets[1],]
                enr_bp[enr_bp['Adjusted P-value']<0.05].to_csv(os.path.join(out, 'GO_Biological_Process_Enrichment.xls'), sep='\t', index=None)
                barplot(enr_bp,
                        title='GO Biological Process Enrichment', 
                        cutoff=0.05,
                        color='darkred',
                        top_term=20, 
                        figsize=(6, 8),
                        ofname=os.path.join(out, 'GO_Biological_Process_Enrichment.png'))
        except:
            pass

    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    if organism != 'No':
        if organism == 'human':
            gene_annot = os.path.join('config', 'annotation', '9606_human_go_kegg_annotation.xls')
        else:
            gene_annot = os.path.join('config', 'annotation', '10090_mouse_go_kegg_annotation.xls')
        gene_annot = pd.read_csv(gene_annot, sep='\t')
        all_diff = pd.merge(all_diff, gene_annot, on='Symbol', how='left')
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'cluster_diff_statistics.xls'), sep='\t')

    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.rank_genes_groups_dotplot(
            _adata,
            n_genes=5,
            values_to_plot="logfoldchanges",
            vmin=-3,
            vmax=3,
            min_logfoldchange=1,
            cmap='bwr',
            colorbar_title='log fold change',
            show=False
        )
        plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight", dpi=300)

    with rc_context({'figure.figsize': (6, 6)}):
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
        plt.savefig(os.path.join(out_dir, 'heatmap_top10_marker_gene.png'), facecolor='white', bbox_inches="tight", dpi=300)


    # # Marker (COSG)
    # cosg.cosg(_adata,
    #         groups='all',
    #         key_added='cosg',
    #         mu=1,
    #         n_genes_user=200,
    #         groupby='leiden')
    # colnames = ['names', 'scores']
    # test = [pd.DataFrame(_adata.uns["cosg"][c]) for c in colnames]
    # test = pd.concat(test, axis=1, names=[None, 'group'], keys=colnames)

    # markers = {}
    # all_cluster_marker = []
    # cats = _adata.obs.leiden.cat.categories
    # for i, c in enumerate(cats):
    #     tmp_out_dir = os.path.join(out_dir, 'cluster_marker_gene', f'cluster{c}')
    #     os.makedirs(tmp_out_dir)
    #     cell_type_df = test.loc[:, 'names'][c]
    #     scores_df = test.loc[:, 'scores'][c]
    #     markers[c] = cell_type_df.values.tolist()[:4]
    #     all_markers = pd.DataFrame({'names': list(cell_type_df), 'scores': list(scores_df)})
    #     all_markers['cluster'] = c
    #     all_markers.to_csv(os.path.join(tmp_out_dir, 'differentially_expressed_genes.xls'), sep='\t', index=None)

    #     try:
    #         if organism != 'No':
    #             gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023']
    #             if organism == 'mouse':
    #                 gene_sets=['KEGG_2019_Mouse', 'GO_Biological_Process_2023']

    #             enr = gp.enrichr(gene_list=all_markers['names'],
    #                             gene_sets=gene_sets,
    #                             organism=organism,
    #                             #  background=_adata.var_names.tolist(),
    #                             outdir=None, # don't write to disk
    #                             )

    #             enr_kegg = enr.results.loc[enr.results["Gene_set"] == gene_sets[0],]
    #             enr_kegg[enr_kegg['Adjusted P-value']<0.05].to_csv(os.path.join(tmp_out_dir, 'KEGG_Enrichment.xls'), sep='\t', index=None)
    #             barplot(enr_kegg,
    #                     title='KEGG Enrichment', 
    #                     cutoff=0.05,
    #                     color='darkred',
    #                     top_term=20, 
    #                     figsize=(6, 8),
    #                     ofname=os.path.join(tmp_out_dir, 'KEGG_Enrichment.png'))

    #             enr_bp = enr.results.loc[enr.results["Gene_set"] == gene_sets[1],]
    #             enr_bp[enr_bp['Adjusted P-value']<0.05].to_csv(os.path.join(tmp_out_dir, 'GO_Biological_Process_Enrichment.xls'), sep='\t', index=None)
    #             barplot(enr_bp,
    #                     title='GO Biological Process Enrichment', 
    #                     cutoff=0.05,
    #                     color='darkred',
    #                     top_term=20, 
    #                     figsize=(6, 8),
    #                     ofname=os.path.join(tmp_out_dir, 'GO_Biological_Process_Enrichment.png'))
    #     except:
    #         pass

    #     all_cluster_marker.append(all_markers)
    
    # tmp = pd.concat(all_cluster_marker, axis=0)
    # tmp.rename(columns={'names': 'Symbol'}, inplace=True)
    # tmp.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    
    # sc.pl.dotplot(_adata,
    #             var_names = markers, 
    #             groupby='leiden',
    #             cmap='Spectral_r',
    #             standard_scale='var',
    #             show=False)
    # plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight", dpi=300)
    # _adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')

    samples = _adata.obs['group'].drop_duplicates().tolist()
    if len(samples) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
        fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 6))
        for i in range(len(samples)):
            if i == len(samples) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["group"]==samples[i]], color="leiden", size=dot_size, legend_loc=legend, title=samples[i], ax=axs[i], show=False)
        plt.savefig(os.path.join(out_dir, 'umap_group_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)

    sc.set_figure_params(dpi=300, color_map = 'viridis_r')
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color=["n_genes_by_counts", "total_counts", "total_counts_mt"], wspace=0.2, size=dot_size, ncols=3, use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_n_genes_by_counts.png'), facecolor='white', bbox_inches="tight", dpi=300)

    return _adata

def visua_gene_express(_adata, genes, dot_size, out_dir):
    st.set_option('deprecation.showPyplotGlobalUse', False)
    sc.set_figure_params(dpi=100, color_map = 'viridis_r')
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color=genes, ncols=3, cmap='coolwarm', size=dot_size, use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, '_'.join(genes)+'.png'), facecolor='white', bbox_inches="tight", dpi=300)
    
    # with rc_context({'figure.figsize': (8, 8)}):
    #     sc.pl.dotplot(_adata, genes, 'leiden', dendrogram=False)
    #     plt.savefig(os.path.join(out_dir, '_'.join(genes)+'_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)

    return _adata

def celltype_rename(_adata, annotation_dict, out_dir):
    # rename
    _adata.obs['cell_type'] = [annotation_dict[clust] for clust in _adata.obs["leiden"]]

    sc.pl.umap(_adata, color='cell_type', size=25, show=False, title='')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    sc.pl.tsne(_adata, color='cell_type', size=25, show=False, title='')
    plt.savefig(os.path.join(out_dir, 'dimplot_tsne_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    from matplotlib.pyplot import rc_context
    sc.set_figure_params(dpi=80, color_map = 'viridis_r')
    groups = _adata.obs['group'].drop_duplicates().tolist()
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 6))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
        plt.savefig(os.path.join(out_dir, 'umap_groups_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.umap(_adata, color=["leiden","cell_type"], wspace=0.2, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    
    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.tsne(_adata, color=["leiden","cell_type"], wspace=0.2, show=False)
        plt.savefig(os.path.join(out_dir, 'tsne_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    pd.DataFrame(_adata.obsm['X_tsne'], index=_adata.obs_names, columns=['tSNE_1', 'tSNE_2']).to_csv(os.path.join(out_dir, 'data_tsne.xls'), sep='\t')
    pd.DataFrame(_adata.obsm['X_umap'], index=_adata.obs_names, columns=['UMAP_1', 'UMAP_2']).to_csv(os.path.join(out_dir, 'data_umap.xls'), sep='\t')
    _adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')
    os.system(f"Rscript ./utils/Rscript/cell_proportion.R {os.path.join(out_dir, 'all_cells_features.xls')} cell_type {out_dir}")

    _adata.write(os.path.join(out_dir, 'sub_adata.h5ad'))

    return _adata
