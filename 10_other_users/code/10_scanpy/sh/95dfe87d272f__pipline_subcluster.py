import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from harmonypy import run_harmony
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

# 对B cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/00.data/Main/Main.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Erythroid cells', 'HSPCs', 'Megakaryocytes', 'Myeloid cells', 'NK cells', 'Neutrophils', 'Plasma cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/B cells/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="B cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

# 运行 Harmony
# adata_harmony = run_harmony(adata, adata.obs, ['platform'])
# sc.pp.neighbors(adata_harmony)
# sc.tl.umap(adata_harmony)
# sc.pl.umap(adata, color='leiden')
# plt.savefig(os.path.join(out_dir, 'harmony_umap_celltype.png'))
# plt.savefig(os.path.join(out_dir, 'harmony_umap_celltype.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['CD27','FCER2','CD19','CD80','NT5E','PDCD1LG2'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_B_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_B_diffgenes.pdf'))

anno_dict={
    '0':'Naive B',
    '1':'Naive B',
    '2':'Pre-B',
    '3':'Pre-B',
		'4':'Memory B',
    '5':'Memory B',
    '6':'Naive B',
    '7':'Activated B'
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
adata.write(os.path.join(out_dir, 'Bcells.h5ad'))



# 对Plasma cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/00.data/Main/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/Plasma cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 10,output=out_dir)
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

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHM','IGHD','IGHE','CD3D','CD3E'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'))

cell_type_mapping = adata.obs[['leiden']].to_dict()['leiden']
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '1':'T cells',
    '5':'T cells',
    '10':'T cells'
}
mask = adata_raw.obs['sub_cell_type'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
# 更新 categories
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'sub_cell_type'].map(anno_dict)
adata_raw.write("/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/Main.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.2,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHM','IGHD','IGHE'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'))

anno_dict={
    '0':'IgA Plasma',
    '1':'IgG Plasma',
    '2':'IgA Plasma',
    '3':'IgA Plasma',
		'4':'IgG Plasma',
    '5':'IgG Plasma',
    '6':'lowQ',
    '7':'IgM Plasma',
    '8':'IgG Plasma',
    '9':'IgG Plasma',
    '10':'IgG Plasma',
    '11':'IgA Plasma',
		'12':'IgM Plasma'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_map = adata[adata.obs['leiden'].isin(['6']), :]
sc.tl.leiden(adata_map, resolution=0.1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata_map, color=['IGHA1','IGHG1','IGHM','IGHD','IGHE'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'))

anno_dict={
    '0':'IgA Plasma',
    '1':'IgG Plasma',
    '2':'IgA Plasma',
    '3':'IgG Plasma'
}
adata_map.obs['cell_type'] = [anno_dict[clust] for clust in adata_map.obs['leiden']]
adata_map.obs['cell_type'] = adata_map.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_mapping = adata_map.obs[['cell_type']].to_dict()['cell_type']
adata.obs['cell_type'] = adata.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
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
adata.write(os.path.join(out_dir, 'Plasma.h5ad'))

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/mapping_raw/"
adata_raw.obs['test'] = adata_raw.obs['cell_type']

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/B cells/Bcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['test'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['test']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/01.Mission/Plasma cells/Plasma.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['test'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['test']), axis=1)

adata_raw.obs['test'] = adata_raw.obs['test'].astype('category')
adata_raw.obs['test'].cat.categories

sc.pl.umap(adata_raw, color='test', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)