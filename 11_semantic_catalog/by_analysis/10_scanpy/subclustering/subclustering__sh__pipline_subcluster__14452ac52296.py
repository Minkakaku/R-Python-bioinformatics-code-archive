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

# 对CD8T cells注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/12.Mission/CD8T/"
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/CD8T/CD8T.h5ad")

anno_dict={
    '0':'Cmss1+',
    '1':'Pf4+',
    '2':'Atp2b1+',
    '3':'Isg15+',
    '4':'Ccr7+',
    '5':'Top2a+',
    '6':'Gzma+'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
adata.uns['cell_type_colors']=['#ECA8A9', '#74AED4', '#D3E2B7', '#CFAFD4', '#F7C97E','#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#b5bd61', '#17becf', '#aec7e8']
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'),bbox_inches='tight')

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata.obs['group'] = (
    pd.Categorical(adata.obs['group'],
                   categories=groups,
                   ordered=True)
)
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
adata.obs['sample'] = (
    pd.Categorical(adata.obs['sample'],
                   categories=['asham1', 'asham2','asham3','aday11', 'aday12', 'aday13','aday31', 'aday32', 'aday33', 'aday71', 'aday72', 'aday73', 'aday141', 'aday142', 'aday143', ],
                   ordered=True)
)
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'CD8T.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

cross_tab = pd.crosstab(adata.obs['group'], adata.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/12.Mission/CD8T/celltype_proportions_by_group.csv")


# marker_UMAP图
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/12.Mission/gene_plot/"
colors = ['white', 'red', 'black']
cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)
sc.pl.umap(adata, color=['Prf1','Gzma','Gzmb','Gzmk'], cmap=cmap, frameon=False,colorbar_loc=None,show=False)
# plt.savefig(os.path.join(out_dir, 'umap_genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'),bbox_inches='tight')
save_fig(plt, out_dir = out_dir, file_name= 'umap_genes', fig_size={'w':20, 'h':5})

sc.pl.dotplot(adata, var_names=['Prf1','Gzma','Gzmb','Gzmk'], groupby='group', dendrogram=False, show=False, swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight', dpi=300, facecolor='white')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight', dpi=300, facecolor='white')


fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.umap(adata, color=['Ccr2'], cmap=cmap, frameon=False,ax=ax,colorbar_loc=None,show=False,size=25)
plt.savefig(os.path.join(out_dir, 'umap_Ccr2.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Ccr2.pdf'),bbox_inches='tight')

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots(figsize=(6, 4))
    sc.pl.violin(adata, keys=['Ccr2'], groupby='group', show=False, ax=ax, 
             stripplot=False, density_norm='count', xlabel=None, 
             ylabel='Expression level', jitter=False)
# 计算每个组的均值并添加点
    mean_values = []
    for i, group in enumerate(groups):
        expr = adata.raw[adata.obs['group'] == group, 'Ccr2'].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        mean_val = np.mean(expr)
        mean_values.append(mean_val)
        ax.scatter(i, mean_val, color='black', s=50, zorder=10)
    ax.plot(range(len(groups)), mean_values, color='black', linestyle='-', linewidth=2, zorder=5)
    ax.grid(False)
    plt.savefig(os.path.join(out_dir, 'violin_Ccr2.png'))
    plt.savefig(os.path.join(out_dir, 'violin_Ccr2.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['Ccr2'], groupby='cell_type', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_Ccr2.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_Ccr2.pdf'))


# marker_UMAP图
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/12.Mission/CD8T/"
colors = ['lightgrey', 'red', 'black']
cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)

fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.umap(adata_sub, color=['Il2ra'], cmap=cmap, frameon=False,ax=ax,show=False)
plt.savefig(os.path.join(out_dir, 'umap_Il2ra.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Il2ra.pdf'),bbox_inches='tight')