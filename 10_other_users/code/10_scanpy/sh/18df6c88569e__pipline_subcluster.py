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

# 对Macrophage注释细分
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/01.Mission/Macrophages/Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/04.Mission/Macrophages/"
sc.pl.umap(adata, color=['Trem2'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Trem2'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Trem2'], groupby=['group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_group_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Trem2'], groupby=['cell_type','sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-sample_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-sample_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Trem2'], groupby=['sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_sample_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_sample_genes.pdf'),bbox_inches='tight')

out_dir="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/04.Mission/M1_M2/"
sc.pl.umap(adata, color=['Il1b','Mrc1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Il1b','Mrc1'], groupby=['leiden','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Il1b','Mrc1'], groupby=['group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_group_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Il1b','Mrc1'], groupby=['leiden','sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-sample_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-sample_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Il1b','Mrc1'], groupby=['sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_sample_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_sample_genes.pdf'),bbox_inches='tight')


# 对Endothelial注释细分
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/01.Mission/Endothelial/Endothelial.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/04.Mission/Endothelial/"
sc.pl.umap(adata, color=['Dusp1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Dusp1'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-group_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Dusp1'], groupby=['group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_group_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Dusp1'], groupby=['cell_type','sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-sample_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cell_type-sample_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Dusp1'], groupby=['sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_sample_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_sample_genes.pdf'),bbox_inches='tight')