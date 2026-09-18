import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from matplotlib import font_manager as fm
import matplotlib.colors as mcolors

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    '''
    out_dir : 输出文件夹地址
    file_name : 文件名称
    fig_size : 输出图片的大小
    '''
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None):
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    if fig_size is None:
        fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20252866/01.Mission/Monocytes/Macrophage.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20252866/01.Mission/Monocytes/"
sc.tl.rank_genes_groups(adata_raw, "cell_type", method="wilcoxon", pts=True)
rank_genes = adata_raw.uns['rank_genes_groups']
# 筛选满足条件的基因
# 初始化一个空的列表来存储符合条件的基因
selected_genes = []
# 遍历每个细胞类型组
for group in rank_genes['names'].dtype.names:
    # 获取当前组的基因信息
    rank_genes_sub = sc.get.rank_genes_groups_df(adata_raw, group=group)
    mask = (rank_genes_sub['logfoldchanges'] > 2) & (rank_genes_sub['pct_nz_group'] > 0.85) # 筛选满足条件的基因
    selected_genes.extend(rank_genes_sub['names'][mask])
# 去重
selected_genes = list(set(selected_genes))
# 绘制热图
sc.pl.heatmap(adata_raw, var_names=selected_genes, groupby='cell_type', 
              dendrogram=False, show_gene_labels=True, 
              figsize=(6, 15), cmap='viridis', swap_axes=True)
plt.savefig(os.path.join(out_dir, 'rank_genes_groups_heatmap.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'rank_genes_groups_heatmap.pdf'),bbox_inches='tight')

top_genes = []
for group in rank_genes['names'].dtype.names:
    # 获取当前组的基因信息
    rank_genes_sub = sc.get.rank_genes_groups_df(adata_raw, group=group)
    top_genes.extend(rank_genes_sub['names'].head(10))
    sc.pl.umap(adata_raw, color=rank_genes_sub['names'].head(10), show=False,cmap='Oranges',ncols=5)
    plt.savefig(os.path.join(f'{out_dir}umap_{group}_top10genes.png'))
    plt.savefig(os.path.join(f'{out_dir}umap_{group}_top10genes.pdf'))
# 绘制热图
sc.pl.heatmap(adata_raw, var_names=top_genes, groupby='cell_type', 
              dendrogram=False, show_gene_labels=True, 
              figsize=(12, 8), cmap='viridis')
plt.savefig(os.path.join(out_dir, 'top10_rank_genes_groups_heatmap.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_rank_genes_groups_heatmap.pdf'),bbox_inches='tight')

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20252866/01.Mission/M1_vs_M2/"
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(['M1','M2']), :]
adata_sub.uns['log1p']['base'] = 2
sc.tl.rank_genes_groups(adata_sub,groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata_sub,species='human',output=out_dir, logfc=0.25, fdr = 0.05, pct = 0.1)

samples = adata_raw.obs['sample'].drop_duplicates().tolist()
samples = ['AS1','AS2','AS3','C1','C2','C3']
if len(samples) > 1:
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
    axs = axs.flatten()
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["sample"]==samples[i]], color="MRC1", legend_loc=legend, title=samples[i], ax=axs[i], show=False,cmap='Oranges')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_MRC1.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_MRC1.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.umap(adata_raw, color=['CD86','IL1B','NOS2','CD163','MRC1'], show=False,cmap='Oranges',ncols=3)
plt.savefig(os.path.join(f'{out_dir}umap_{group}_top10genes.png'))
plt.savefig(os.path.join(f'{out_dir}umap_{group}_top10genes.pdf'))



adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20252866/01.Mission/Fibroblasts/Fibroblast.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20252866/01.Mission/Fibroblasts-AS_vs_control/"
adata_raw.uns['log1p']['base'] = 2
sc.tl.rank_genes_groups(adata_raw,groupby='group',method='wilcoxon',pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir, logfc=0.25, fdr = 0.05, pct = 0.1)

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20252866/01.Mission/Fibroblasts/"
sc.tl.rank_genes_groups(adata_raw, "cell_type", method="wilcoxon", pts=True)
rank_genes = adata_raw.uns['rank_genes_groups']
# 筛选满足条件的基因
# 初始化一个空的列表来存储符合条件的基因
selected_genes = []
# 遍历每个细胞类型组
for group in rank_genes['names'].dtype.names:
    # 获取当前组的基因信息
    rank_genes_sub = sc.get.rank_genes_groups_df(adata_raw, group=group)
    mask = (rank_genes_sub['logfoldchanges'] > 2) & (rank_genes_sub['pct_nz_group'] > 0.85) # 筛选满足条件的基因
    selected_genes.extend(rank_genes_sub['names'][mask])
# 去重
selected_genes = list(set(selected_genes))
# 绘制热图
sc.pl.heatmap(adata_raw, var_names=selected_genes, groupby='cell_type', 
              dendrogram=False, show_gene_labels=True, 
              figsize=(6, 24), cmap='viridis', swap_axes=True)
plt.savefig(os.path.join(out_dir, 'rank_genes_groups_heatmap.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'rank_genes_groups_heatmap.pdf'),bbox_inches='tight')

top_genes = []
for group in rank_genes['names'].dtype.names:
    # 获取当前组的基因信息
    rank_genes_sub = sc.get.rank_genes_groups_df(adata_raw, group=group)
    top_genes.extend(rank_genes_sub['names'].head(10))
    sc.pl.umap(adata_raw, color=rank_genes_sub['names'].head(10), show=False,cmap='Oranges',ncols=5)
    plt.savefig(os.path.join(f'{out_dir}umap_{group}_top10genes.png'))
    plt.savefig(os.path.join(f'{out_dir}umap_{group}_top10genes.pdf'))
# 绘制热图
sc.pl.heatmap(adata_raw, var_names=top_genes, groupby='cell_type', 
              dendrogram=False, show_gene_labels=True, 
              figsize=(12, 8), cmap='viridis')
plt.savefig(os.path.join(out_dir, 'top10_rank_genes_groups_heatmap.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_rank_genes_groups_heatmap.pdf'),bbox_inches='tight')
