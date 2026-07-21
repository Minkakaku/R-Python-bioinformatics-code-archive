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

# Fibroblasts细分 重新注释
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/04.Mission/subcluster_Fibroblasts/"
os.makedirs(out_dir, exist_ok=True)
adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250326/01.Mission/subcluster_Fibroblasts/Fibroblasts.h5ad")
adata.obs['cell_type'].cat.categories
adata.obs['sample'].cat.categories
adata.obs['group'].cat.categories
group_dict={
    'control1':'Control',
    'control2':'Control',
    'BLM1':'Case_control',
    'BLM2':'Case_control',
    'BLM3':'Case_control',
    'BLM4':'Case_control',
    'Treatment1':'Treatment1',
    'Treatment2':'Treatment1',
    'Treatment3':'Treatment2',
    'Treatment4':'Treatment2',
    'Treatment5':'Treatment2',
    'Treatment6':'Treatment2'
}
adata.obs['group'] = [group_dict[samples] for samples in adata.obs['sample']]
adata.obs['group'] = adata.obs['group'].astype("category")

anno_dict={
    '0':'Serpina3n Fibroblasts',
    '1':'Spp1 Fibroblasts',
    '2':'Fth1 Fibroblasts',
    '3':'Inmt Fibroblasts',
    '4':'Dcn Fibroblasts',
    '5':'Cd74 Fibroblasts',
    '6':'Msln Fibroblasts'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cat_type = pd.CategoricalDtype(categories= ['Control', 'Case_control', 'Treatment1', 'Treatment2'], ordered=True)
adata.obs['group'] = adata.obs['group'].astype(cat_type)
groups = adata.obs['group'].cat.categories.tolist()
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

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))


# Macrophages细分 重新注释
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/04.Mission/subcluster_Macrophages/"
os.makedirs(out_dir, exist_ok=True)
adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250326/01.Mission/subcluster_Macrophages/Macrophages.h5ad")
adata.obs['cell_type'].cat.categories
adata.obs['sample'].cat.categories
adata.obs['group'].cat.categories
group_dict={
    'control1':'Control',
    'control2':'Control',
    'BLM1':'Case_control',
    'BLM2':'Case_control',
    'BLM3':'Case_control',
    'BLM4':'Case_control',
    'Treatment1':'Treatment1',
    'Treatment2':'Treatment1',
    'Treatment3':'Treatment2',
    'Treatment4':'Treatment2',
    'Treatment5':'Treatment2',
    'Treatment6':'Treatment2'
}
adata.obs['group'] = [group_dict[samples] for samples in adata.obs['sample']]
adata.obs['group'] = adata.obs['group'].astype("category")

anno_dict={
    '0':'C1qc Macrophages',
    '1':'Ctsb Macrophages',
    '2':'Chil3 Macrophages',
    '3':'Ccr2 Macrophages',
    '4':'Napsa Macrophages',
    '5':'DC',
    '6':'Ear1 Macrophages',
    '7':'Samsn1 Macrophages',
    '9':'DC',
    '10':'Sparc Macrophages'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cat_type = pd.CategoricalDtype(categories= ['Control', 'Case_control', 'Treatment1', 'Treatment2'], ordered=True)
adata.obs['group'] = adata.obs['group'].astype(cat_type)
groups = adata.obs['group'].cat.categories.tolist()
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

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250326/02.Mission/mapping_raw/mapping_raw.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/04.Mission/mapping_raw/"
os.makedirs(out_dir, exist_ok=True)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250326/04.Mission/subcluster_Fibroblasts/Fibroblasts.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

group_dict={
    'control1':'Control',
    'control2':'Control',
    'BLM1':'Case_control',
    'BLM2':'Case_control',
    'BLM3':'Case_control',
    'BLM4':'Case_control',
    'Treatment1':'Treatment1',
    'Treatment2':'Treatment1',
    'Treatment3':'Treatment2',
    'Treatment4':'Treatment2',
    'Treatment5':'Treatment2',
    'Treatment6':'Treatment2'
}
adata_raw.obs['group'] = [group_dict[samples] for samples in adata_raw.obs['sample']]
adata_raw.obs['group'] = adata_raw.obs['group'].astype("category")

cat_type = pd.CategoricalDtype(categories= ['Control', 'Case_control', 'Treatment1', 'Treatment2'], ordered=True)
adata_raw.obs['group'] = adata_raw.obs['group'].astype(cat_type)
groups = adata_raw.obs['group'].cat.categories.tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="sub_cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'umap_mapping_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_mapping_group.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'mapping_raw.h5ad'))