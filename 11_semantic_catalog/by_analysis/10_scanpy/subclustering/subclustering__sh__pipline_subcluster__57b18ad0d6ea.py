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

# Macrophage注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2512157481-01-wmx/06-Cell_type_annotation/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages', 'Mast cells', 'NK cells', 'Neutrophils', 'Pericyte', 'Proliferating cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512157481-01/01.Mission/Macrophages/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['LYZ','CD14','VCAN'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))

adata_sub = adata[~adata.obs.leiden.isin(['13','14']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.3,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

sc.pl.umap(adata, color=['NOS2','CD86','CCL5','IL1B','CD274','MRC1','CD163','CCND1','CKB','GSN'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['NOS2','CD86','CCL5','IL1B','CD274','MRC1','CD163','CCND1','CKB','GSN'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'M2',
    '1':'M2',
    '2':'M1',
    '3':'M1',
    '4':'M2',
    '5':'bulkM',
    '6':'M2',
    '7':'M2'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.dotplot(adata, var_names=['CD86','IL1B','MRC1','CD163','CCND1'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'),bbox_inches='tight')

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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Macrophage.h5ad'))

out_dir="/mnt/NexData_file/home/sh/Projects/SC2512157481-01/01.Mission/Macrophages/EPOR_plot/"
sc.pl.umap(adata, color=['EPOR'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_EPOR.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_EPOR.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['EPOR'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type_EPOR.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type_EPOR.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['EPOR'], groupby='group', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_group_EPOR.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_EPOR.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['EPOR'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_EPOR.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_EPOR.pdf'),bbox_inches='tight')