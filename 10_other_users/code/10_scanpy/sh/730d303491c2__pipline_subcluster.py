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


adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511126506-03/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511126506-03/01.Mission/markerplot/"
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata_raw, var_names=['Cd276','Gzmb'], groupby='cell_type', dendrogram=False, show=False,swap_axes=True)
    plt.savefig(os.path.join(out_dir, 'cell_type_dotplot_genes.png'))
    plt.savefig(os.path.join(out_dir, 'cell_type_dotplot_genes.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata_raw, var_names=['Cd276','Gzmb'], groupby='group', dendrogram=False, show=False,swap_axes=True)
    plt.savefig(os.path.join(out_dir, 'group_dotplot_genes.png'))
    plt.savefig(os.path.join(out_dir, 'group_dotplot_genes.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata_raw, var_names=['Cd276','Gzmb'], groupby=['cell_type','group'], dendrogram=False, show=False,swap_axes=True)
    plt.savefig(os.path.join(out_dir, 'cell_type_group_dotplot_genes.png'))
    plt.savefig(os.path.join(out_dir, 'cell_type_group_dotplot_genes.pdf'))

sc.pl.umap(adata_raw, color=['Cd276'], show=False, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_Cd276.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Cd276.pdf'),bbox_inches='tight')
sc.pl.umap(adata_raw, color=['Gzmb'], show=False, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_Gzmb.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Gzmb.pdf'),bbox_inches='tight')

genes = ['Cd276', 'Gzmb']
# 直接提取表达矩阵并转换为 DataFrame
expr_df = pd.DataFrame(
    adata_raw.raw[:, genes].X.toarray() if hasattr(adata_raw.raw[:, genes].X, 'toarray') else adata_raw.raw[:, genes].X,
    columns=genes,
    index=adata_raw.obs_names
)
# 合并到 obs
adata_raw.obs = adata_raw.obs.join(expr_df)
# 输出表格
adata_raw.obs[['cell_type', 'group'] + genes].to_csv('gene_expression.csv')
