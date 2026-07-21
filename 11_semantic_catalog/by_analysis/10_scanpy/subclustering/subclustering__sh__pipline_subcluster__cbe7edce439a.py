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

out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/17.Mission/CD8T/"
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_C/Tsub_cells.h5ad")
sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

# marker_UMAP图
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# 'Themis', 'Lck', 'Klrc1', 'Ptpn22', 'Thy1', 'Klrk1', 'Klrd1', 'Cd226', 'Ifng', 'Pik3ap1'
colors = ['white', 'red', 'black']
cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)

fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.umap(adata, color=['Themis'], cmap=cmap, frameon=False,ax=ax,colorbar_loc=None,show=False,size=25)
plt.savefig(os.path.join(out_dir, 'umap_Themis.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Themis.pdf'),bbox_inches='tight')

groups=['asham','aday1','aday3','aday7','aday14']
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots(figsize=(6, 4))
    sc.pl.violin(adata, keys=['Klrc1'], groupby='group', show=False, ax=ax, 
             stripplot=False, density_norm='count', xlabel=None, 
             ylabel='Expression level', jitter=False)
# 计算每个组的均值并添加点
    mean_values = []
    for i, group in enumerate(groups):
        expr = adata.raw[adata.obs['group'] == group, 'Klrc1'].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        mean_val = np.mean(expr)
        mean_values.append(mean_val)
        ax.scatter(i, mean_val, color='black', s=50, zorder=10)
    ax.plot(range(len(groups)), mean_values, color='black', linestyle='-', linewidth=2, zorder=5)
    ax.grid(False)
    plt.savefig(os.path.join(out_dir, 'violin_Klrc1.png'))
    plt.savefig(os.path.join(out_dir, 'violin_Klrc1.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['Themis'], groupby='cell_type', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_Themis.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_Themis.pdf'))


with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color=['Themis'], cmap='coolwarm',show=False,vmax='p99')
    plt.savefig(os.path.join(out_dir, 'umap_Themis_default.png'), bbox_inches = 'tight')
    plt.savefig(os.path.join(out_dir, 'umap_Themis_default.pdf'), bbox_inches = 'tight')


adata_raw=sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/07.Mission/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/17.Mission/"
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype_legend.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype_legend.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)