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

# 对Myeloid注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252616/07.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252616/08.Mission/mapping_raw/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata_raw, color=['ACTA2','TAGLN','COL1A1','COL3A1','FN1','PDGFRA','VIM'], cmap='coolwarm',show=False,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'), bbox_inches = 'tight')
sc.pl.dotplot(adata_raw, var_names=['ACTA2','TAGLN','COL1A1','COL3A1','FN1','PDGFRA','VIM'], groupby='cluster', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.pdf'), bbox_inches = 'tight')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cluster', title='cluster', show=False, legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.pdf'))


adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252616/07.Mission/Myeloid/Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252616/08.Mission/Myeloid/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata_raw, color=['ACTA2','TAGLN','COL1A1','COL3A1','FN1','PDGFRA','VIM'], cmap='coolwarm',show=False,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'), bbox_inches = 'tight')
sc.pl.dotplot(adata_raw, var_names=['ACTA2','TAGLN','COL1A1','COL3A1','FN1','PDGFRA','VIM'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.pdf'), bbox_inches = 'tight')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='leiden', title='leiden', show=False, legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.pdf'))

adata=adata_raw.copy()
anno_dict={
    '0':'SPAAR- Macrophages',
    '1':'SPAAR- Macrophages',
    '2':'SPAAR- Macrophages',
    '3':'SPAAR- Macrophages',
    '4':'SPAAR- Macrophages',
    '5':'SPAAR- Macrophages',
    '6':'SPAAR- Macrophages',
    '7':'SPAAR- Macrophages',
    '8':'SPAAR+ Macrophages',
    '9':'Myofibroblast',
    '10':'SPAAR- Macrophages',
    '11':'SPAAR+ Macrophages'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot_cell_type.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot_cell_type.pdf'),bbox_inches='tight')
plt.close()

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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))

# 映射所有细分 后面没做
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252616/07.Mission/mapping_raw/"
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252616/02.Mission/mapping_raw/Main_mapping.h5ad")

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')