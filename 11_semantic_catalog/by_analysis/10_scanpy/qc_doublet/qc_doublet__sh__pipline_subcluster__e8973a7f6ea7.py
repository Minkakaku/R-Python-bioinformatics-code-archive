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
import seaborn as sns

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

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20251028/02.Mission/FWSC20251028/06-Cell_type_annotation/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'AT1', 'AT2', 'B cells', 'Basal cells', 'Ciliated cells', 'Doublets', 'Endothelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophil', 'Plasma cells', 'SMC', 'T/NK', 'lowQ'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20251028/03.Mission/Main/"

mainmarkers=['ITGB4',#LNEP
                 'TP63','KRT5','KRT17',#BASAL
                 'CALCA','CHGA',#NEUROENDOCRINE 16？ √
                 'SCGB3A1','MUC5AC','MUC5B',#GOBLET ×
                 'CCDC153','FOXJ1','CDHR3',#CILIATED √
                 'SCGB3A2','SCGB1A1',#CLUB
                 'KRT8',#KRT8+ 17？ √
                 'SFTPD','SFTPC','SFTPA1','CHI3L1',#AT2 √
                 'PDPN','AGER','AQP5']#AT1 √
sc.pl.dotplot(adata_raw, var_names=mainmarkers, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_mainmarkers.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_mainmarkers.pdf'),bbox_inches='tight')
sc.pl.umap(adata_raw, color=mainmarkers, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_mainmarkers.png'))
plt.savefig(os.path.join(out_dir, 'umap_mainmarkers.pdf'))

adata_raw = adata_raw[~adata_raw.obs['cell_type'].isin(['Doublets','lowQ']), :]

sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.tl.umap(adata_raw)
sc.tl.leiden(adata_raw, resolution=0.8)
# plot cluster
sc.pl.umap(adata_raw, color='leiden', show=False,legend_loc='on data')
plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '11':'Club cells',
    '20':'Neuroendocrine'
}
mask = adata_raw.obs['leiden'].isin(anno_dict.keys())
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'leiden'].map(anno_dict)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_raw.write(os.path.join(out_dir, 'Main_recluster.h5ad'))

groups =  adata_raw.obs['group'].drop_duplicates().tolist()
groups =  ['NS', 'S', 'EP', 'EPCA','CA']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':9, 'h':10})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':9, 'h':10})