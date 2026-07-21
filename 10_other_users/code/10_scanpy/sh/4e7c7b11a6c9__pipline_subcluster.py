import os, sys
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
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

# 每个leiden cluster0-16汇总其cell_type构成比例
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2512227685-01/07.Mission/05_cluster_fine_annotation/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512227685-01/09.Mission/result_1/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='leiden', title='leiden', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_leiden.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_leiden.pdf'))

plot_celltype_proportion(adata, out_dir=out_dir,type='leiden',fig_size={'w':10, 'h':5})
type='leiden'
crosstab = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index')
crosstab.to_csv(os.path.join(out_dir, 'leiden_celltype_proportions.csv'))

type='sample'
keys = adata.obs['leiden'].cat.categories.to_list()
values = list(adata.uns['leiden_colors'])
leiden_colors = dict(zip(keys, values))
ax = pd.crosstab(adata.obs[type], adata.obs['leiden'], normalize='index').plot(
    kind='bar', xlabel='', stacked=True,
    color=[leiden_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['leiden'], normalize='index').columns],
    edgecolor='black', width=0.8,
    figsize=(12,8)
)
ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
legend_font = fm.FontProperties(family='sans-serif', size=10)
ax.legend(
    bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title="Cell Types", frameon=False, labelspacing=1.0, prop=legend_font, fontsize=10, title_fontsize=12, handletextpad=1.5, handleheight=2, 
    columnspacing=1.0, markerscale=1, alignment='left', handlelength=2)
ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
ax.tick_params(axis='y', labelsize=12, labelcolor='black')
ax.grid(False)
ax.set_facecolor('white')
for spine in ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(1.5)
    fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}

plt.tight_layout()
# plt.subplots_adjust(right=0.8)
file_name=f"leiden_proportion_{type}"
save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)
crosstab = pd.crosstab(adata.obs[type], adata.obs['leiden'], normalize='index')
crosstab.to_csv(os.path.join(out_dir, 'sample_leiden_proportions.csv'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color=['n_genes_by_counts','total_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo','doublet_score'], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_QC.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_QC.pdf'))

sc.pl.dotplot(adata, var_names=['n_genes_by_counts','total_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo','doublet_score'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_QC.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_QC.pdf'), bbox_inches="tight")
