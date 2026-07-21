import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
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

# 对巨噬进行细分
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/subcluster_Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_subcluster_Macrophages/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Pro-inflammatory Mac',
    '1':'Immature Mac C1',
    '2':'Immature Mac C1',
    '3':'Remodeling Mac-Itga9',
    '4':'Pro-repair Mac',
    '5':'Immature Mac C1',
    '6':'Immature Mac C2',
    '7':'DC',
    '8':'Monocyte',
    '9':'Proliferative Mac',
    '10':'Pro-inflammatory Mac'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

adata.write(os.path.join(out_dir, 'subcluster_Macrophages.h5ad'))

from scanpy_fun_subcluster import find_group_marker
sc.tl.rank_genes_groups(adata,groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)



# 对巨噬进行重新注释
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_subcluster_Macrophages/subcluster_Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_recluster_Macrophages/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Pro-inflammatory Mac',
    '1':'Immature Mac C1',
    '2':'Immature Mac C1',
    '3':'Remodeling Mac-Svep1',
    '4':'Pro-repair Mac',
    '5':'Immature Mac C1',
    '6':'Immature Mac C2',
    '7':'DC',
    '8':'Monocyte',
    '9':'Proliferative Mac',
    '10':'Pro-inflammatory Mac'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

new_categories_order = ['control3', 'model4', 'PF11-H3', 'PF11-H4']  # 定义新的顺序
adata.obs['sample'] = adata.obs['sample'].cat.set_categories(new_categories_order)

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

adata.write(os.path.join(out_dir, 'recluster_Macrophages.h5ad'))
out_dir = "/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/Markerplot/"
marker_genes=['Svep1','Itga9']
sc.pl.umap(adata, color=marker_genes, show=False, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'Macrophages_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Macrophages_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'Macrophages_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Macrophages_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', color_map='Purples', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Macrophages_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Macrophages_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*9, 5))
    for i in range(len(groups)):
         sc.pl.dotplot(adata[adata.obs["group"]==groups[i]], var_names=marker_genes, groupby='cell_type', color_map='Purples', dendrogram=False, show=False, ax=axs[i], title=groups[i], vmax=0.25, dot_max=0.5)
    plt.savefig(os.path.join(out_dir, 'Macrophages_dotplot_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'Macrophages_dotplot_group.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.close()

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.h5ad")
adata = adata[adata.obs['cell_type'].isin(['Fibroblasts','Mesangial cell','Macrophages']), :]
sc.pl.umap(adata, color=marker_genes, show=False, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'Main_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Main_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'Main_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Main_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', color_map='Purples', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Main_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Main_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*9, 5))
    for i in range(len(groups)):
         sc.pl.dotplot(adata[adata.obs["group"]==groups[i]], var_names=marker_genes, groupby='cell_type', color_map='Purples', dendrogram=False, show=False, ax=axs[i], title=groups[i], vmax=0.75, dot_max=0.6)
    plt.savefig(os.path.join(out_dir, 'Main_dotplot_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'Main_dotplot_group.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.close()
    
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/subcluster_Fibroblasts.h5ad")
sc.pl.umap(adata, color=marker_genes, show=False, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'Fibroblasts_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Fibroblasts_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'Fibroblasts_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Fibroblasts_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', color_map='Purples', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Fibroblasts_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Fibroblasts_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*8, 4))
    for i in range(len(groups)):
         sc.pl.dotplot(adata[adata.obs["group"]==groups[i]], var_names=marker_genes, groupby='cell_type', color_map='Purples', dendrogram=False, show=False, ax=axs[i], title=groups[i], vmax=0.5, dot_max=0.7)
    plt.savefig(os.path.join(out_dir, 'Fibroblasts_dotplot_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'Fibroblasts_dotplot_group.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.close()