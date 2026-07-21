import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
from scanpy_fun_subcluster import add_umap_arrows
import pandas as pd
from matplotlib import font_manager as fm

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

# 调整大类注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Bcells', 'Distal tubular cells', 'Endothelial cells', 'Fibroblasts', 'Intercalated cells', 'Macrophages', 'Podocytes', 'Principal cells', 'Proximal tubule cells', 'Tcells', 'Mesangial cell'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/08.Mission/Main/"

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace('Bcells', 'B cells')
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace('Tcells', 'T cells')
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace('Mesangial cell', 'Mesangial cells')

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data',frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cluster', title='Cluster', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='Sample', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='group', title='group', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

gene_dict = {
    'B cells': ['Cd79a','Ly6d','Cd22'],
    'Distal tubular cells': ['Kng2','Umod','Egf'],
    'Endothelial cells': ['Egfl7','Eng','Pecam1'],
    'Fibroblasts': ['Bgn','Cfh','Col1a2'],
    'Intercalated cells': ['Atp6v1g3','Car2','Atp6v0d2'],
    'Macrophages': ['Ctss','C1qb','C1qa'],
    'Podocytes': ['Podxl','Nphs2','Ptpro'],
    'Principal cells': ['Hsd11b2','Cdh16','Ivns1abp'],
    'Proximal tubule cells': ['Spink1','Pdzk1','Acsm2'],
    'T cells': ['Cd3d','Cd3e','Cd3g'],
    'Mesangial cells': ['Cald1', 'Mylk', 'Acta2']
}
sc.pl.dotplot(adata_raw, var_names=gene_dict, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_scale.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=gene_dict, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')

groups = adata_raw.obs['group'].cat.categories.to_list()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False,frameon=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

samples = adata_raw.obs['sample'].cat.categories.to_list()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["sample"]==samples[i]], color="cell_type", legend_loc=None, title=samples[i], ax=axs[i], show=False,frameon=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

count_table = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
percent_table = count_table.div(count_table.sum(axis=1), axis=0)
percent_table.to_csv(os.path.join(out_dir, 'cell_type_percentage_percents.csv'))
adata_raw.write(os.path.join(out_dir, 'recluster.h5ad'))

out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/09.Mission/UMAP_plot/"
# groups = adata_raw.obs['group'].cat.categories.to_list()
# if len(groups) > 1:
#     axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
#     fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
#     for i in range(len(groups)):
#         if i == len(groups) - 1:
#             legend = 'right margin'
#         else:
#             legend = None
#         sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False,frameon=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)



# 调整巨噬注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_recluster_Macrophages/recluster_Macrophages.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'DC', 'Immature Mac C1', 'Immature Mac C2', 'Monocyte', 'Pro-inflammatory Mac', 'Pro-repair Mac', 'Proliferative Mac', 'Remodeling Mac-Svep1'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/08.Mission/Macrophages/"

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data',frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='leiden', title='Cluster', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='Sample', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='group', title='group', show=False,frameon=False)
    add_umap_arrows(arrow_length=0.15)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

groups = adata_raw.obs['group'].cat.categories.to_list()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False,frameon=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

samples = adata_raw.obs['sample'].cat.categories.to_list()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["sample"]==samples[i]], color="cell_type", legend_loc=None, title=samples[i], ax=axs[i], show=False,frameon=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

count_table = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
percent_table = count_table.div(count_table.sum(axis=1), axis=0)
percent_table.to_csv(os.path.join(out_dir, 'cell_type_percentage_percents.csv'))

sc.tl.rank_genes_groups(adata_raw, groupby="cell_type", method='wilcoxon', pts=True)
# find_group_marker(adata=adata_raw,species='mouse',output=out_dir)
sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=20, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top20genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top20genes.pdf'),bbox_inches='tight')
sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/09.Mission/heatmap_mac/"