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

adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWAE20251249/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts','Leydig cells', 'Macrophages', 'Neutrophils', 'T cells'

# 对Epithelial cells注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/Epithelial cells/"
os.makedirs(out_dir, exist_ok=True)
adata_sub = adata_raw[adata_raw.obs.cell_type=="Epithelial cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.2,n_pcs = 9,output=out_dir)
adata = adata_recluster.sub_cluster()

adata.obs['cell_type'] = adata.obs['leiden']
adata.uns['cell_type_colors'] = adata.uns['leiden_colors']
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['leiden']].to_dict()['leiden']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '0':'Epithelial cells',
    '1':'Epithelial cells',
    '2':'Epithelial cells',
    '3':'Epithelial cells',
    '4':'Macrophages',
    '5':'Epithelial cells',
    '6':'Epithelial cells'
}
mask = adata_raw.obs['sub_cell_type'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
# 更新 categories
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'sub_cell_type'].map(anno_dict)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


adata_raw.write("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/result.h5ad")
adata_sub = adata_raw[adata_raw.obs.cell_type=="Epithelial cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.2,n_pcs = 7,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

adata.obs['cell_type'] = adata.obs['leiden']
adata.uns['cell_type_colors'] = adata.uns['leiden_colors']
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['leiden']].to_dict()['leiden']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.write("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/result.h5ad")

anno_dict={
    '0':'Tumor-like cells',
    '1':'Tumor cells C1',
    '2':'Tumor cells C2',
    '3':'Tumor cells C3',
    '4':'Tumor cells C1',
    '5':'Tumor cells C2'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Epithelialcells.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')



adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/result.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/Epithelial cells/Epithelialcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata.obs['sub_cell_type'] = adata.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata.obs['sub_cell_type'] = adata.obs['sub_cell_type'].astype('category')
adata.obs['sub_cell_type'].cat.categories
adata.write("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/result.h5ad")
sc.pl.umap(adata, color='sub_cell_type', legend_loc='right margin')
plt.savefig("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/Epithelial cells/umap_mapping.png", facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig("/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/Epithelial cells/umap_mapping.pdf", facecolor='white', bbox_inches="tight", dpi=300)

adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251249/01.Mission/Fibroblasts/Fibroblasts.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251249/02.Mission/Fibroblasts/"
os.makedirs(out_dir, exist_ok=True)
anno_dict={
    '0':'Fibroblasts C1',
    '1':'Fibroblasts C2',
    '2':'Fibroblasts C3',
    '3':'Fibroblasts C1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')