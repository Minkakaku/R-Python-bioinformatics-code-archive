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
from matplotlib import font_manager as fm
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

def make_report_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for i in ['01-seeksoultools', '02-Data_QC', '03-Cell_filter', '04-Dimension_reduction_cluster', '05-Marker_gene_analysis', '06-Cell_type_annotation']:
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)
    return out_dir

def set_scanpy_color_palette(_adata, color_key, color_palette=None):
    """
    为 Scanpy 对象设置统一配色方案。
    
    参数:
    - adata: AnnData 对象
    - color_key: str，obs 中的列名，用于分组着色（如 'cell_type'）
    - color_palette: list，可选，自定义的颜色列表（最多 50 个）

    默认色板支持最多 50 种颜色，超出会报错。
    """
    default_palette = [
        "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
        "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
        "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
        "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
        "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
        "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
        "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
        "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
        "#f17d80", "#737495"
    ]
    palette = color_palette if color_palette else default_palette

    if color_key not in _adata.obs:
        raise KeyError(f"'{color_key}' 不存在于 adata.obs 中")

    categories = _adata.obs[color_key].astype(str).unique().tolist()
    n_categories = len(categories)

    if n_categories > len(palette):
        raise ValueError(f"类别数为 {n_categories}，超过了可用颜色数 {len(palette)}，请扩展色板。")

    # 设置为 Categorical 类型并指定顺序
    _adata.obs[color_key] = _adata.obs[color_key].astype("category")
    _adata.obs[color_key] = _adata.obs[color_key].cat.set_categories(categories)

    # 设置颜色
    _adata.uns[f"{color_key}_colors"] = palette[:n_categories]

def plot_celltype_proportion(_adata, out_dir,type='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = _adata.obs['cell_type'].cat.categories.to_list()
    values = list(_adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 创建堆积条形图
    ax = pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )

    # 设置X轴和Y轴标签
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')  # 设置X轴标签
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')  # 设置Y轴标签

    # 设置图例
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )

    # 设置坐标轴标签的字体大小和颜色
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    # 设置背景颜色
    ax.set_facecolor('white')

    # 设置图表边框加粗
    for spine in ax.spines.values():
        spine.set_edgecolor('black')  # 设置边框颜色
        spine.set_linewidth(1.5)

    # 设置图像大小，如果没有传入fig_size，则使用默认大小
    if fig_size is None:
        fig_size = {'width': len(set(_adata.obs[type])) * 2, 'height': 6}

    # 调整布局，确保图例不遮挡图形
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    # 保存图像
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)


## step 5.reduction_cluster

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

    _adata.write(os.path.join(out_dir, 'adata_cluster.h5ad'))

    return _adata


adata_filtered=sc.read_h5ad('/mnt/NexData_file/home/xzx/project/AE2604082792-01/report/03-Cell_filter/adata_filtered.h5ad')
out_dir = '/mnt/NexData_file/home/xzx/project/AE2604082792-01/report/04-Dimension_reduction_cluster'
cluster_method='leiden'
n_pcs = 20
n_neighbors = 20
resolution = 0.5
remove_batch='harmony'
adata_cluster= dim_reduction_cluster(adata_filtered, 
                                    cluster_method,
                                    n_pcs,
                                    n_neighbors,
                                    resolution,
                                    remove_batch,
                                    out_dir)


rank_genes_groups_method='wilcoxon'
rank_genes_groups_fc=0.5
rank_genes_groups_fdr=0.05
rank_genes_groups_pct=0.25
out_dir='/mnt/NexData_file/home/xzx/project/AE2604082792-01/report/05-Marker_gene_analysis'
adata_rank_genes_groups= find_marker_gene(adata_cluster,
                                        method=rank_genes_groups_method,
                                        logfc=rank_genes_groups_fc,
                                        fdr=rank_genes_groups_fdr,
                                        pct=rank_genes_groups_pct,
                                        out_dir=out_dir)  