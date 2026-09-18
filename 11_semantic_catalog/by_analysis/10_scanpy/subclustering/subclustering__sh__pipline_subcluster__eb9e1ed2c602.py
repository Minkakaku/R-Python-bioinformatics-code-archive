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


#### Macrophages 注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Erythrocytes', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Melanocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Macrophages/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

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
    '0':'Lgals3 Macrophages',
    '1':'Apoe Macrophages',
    '2':'Pf4 Macrophages',
    '3':'Lgals3 Macrophages',
    '4':'Pf4 Macrophages',
    '5':'Cd52 Macrophages',
    '6':'Fnip2 Macrophages',
    '7':'Fnip2 Macrophages'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='mouse',output=out_dir)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

cross_tab = pd.crosstab(adata.obs['group'], adata.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Macrophages_ROE/celltype_proportions_by_group.csv")


# marker_UMAP图
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Cxcl16_plot/"
colors = ['white', 'red', 'black']
cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)
fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.umap(adata, color=['Cxcl16'], cmap=cmap, frameon=False,ax=ax,colorbar_loc=None,show=False)
plt.savefig(os.path.join(out_dir, 'umap_Cxcl16.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Cxcl16.pdf'),bbox_inches='tight')

fig = plt.figure(figsize=(9, 7))
main_ax = fig.add_axes([0.05, 0.05, 0.7, 0.9])  # [left, bottom, width, height]
sc.pl.umap( adata,  color=['Cxcl16'],  cmap=cmap,  frameon=False, ax=main_ax, colorbar_loc=None, show=False)
# 颜色条位置（右侧）
cbar_ax = fig.add_axes([0.78, 0.2, 0.03, 0.3])  # [left, bottom, width, height]
# 提取 raw 中 Cxcl16 的表达值
raw_expr = adata.raw[:, 'Cxcl16'].X
if hasattr(raw_expr, 'toarray'):
    raw_expr = raw_expr.toarray().flatten()
else:
    raw_expr = raw_expr.flatten()
    # 使用 raw 数据的 min/max
norm = mpl.colors.Normalize(vmin=raw_expr.min(), vmax=raw_expr.max())
# norm = mpl.colors.Normalize(vmin=adata.obs['Cxcl16'].min(), vmax=adata.obs['Cxcl16'].max())
cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
cb.ax.text(1.2, 0, 'Low', transform=cb.ax.transAxes, fontsize=10, va='center')
cb.ax.text(1.2, 1, 'High', transform=cb.ax.transAxes, fontsize=10, va='center')
cb.ax.set_yticks([])
quantiles = [0.75, 1.5, 2.25]
for q in quantiles:
    cb.ax.axhline(y=q, color='white', linewidth=1, linestyle='dashed', 
                  xmin=0, xmax=1.11, clip_on=False)

plt.savefig(os.path.join(out_dir, 'umap_Cxcl16_with_colorbar.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Cxcl16_with_colorbar.pdf'),bbox_inches='tight')


#### marker 气泡图
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Cxcl16_plot/"
with rc_context({'figure.figsize': (9, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['Cxcl16'], groupby='cell_type', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.violin(adata, ['Cxcl16'], jitter=0, groupby='group', rotation=45, multi_panel=True, show=False)
    plt.savefig(os.path.join(out_dir, 'violin_genes.png'))
    plt.savefig(os.path.join(out_dir, 'violin_genes.pdf'))


#### Fibroblasts 注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Erythrocytes', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Melanocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Fibroblasts/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

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
    '0':'S100a9 Fibroblasts',
    '1':'Cygb Fibroblasts',
    '2':'Igfbp5 Fibroblasts',
    '3':'Igfbp5 Fibroblasts',
    '4':'Prg4 Fibroblasts',
    '5':'Prg4 Fibroblasts',
    '6':'Cygb Fibroblasts',
    '7':'Tnc Fibroblasts',
    '8':'Prg4 Fibroblasts'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='mouse',output=out_dir)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

cross_tab = pd.crosstab(adata.obs['group'], adata.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Fibroblasts_ROE/celltype_proportions_by_group.csv")


#### Neutrophils 注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Erythrocytes', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Melanocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Neutrophils/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Neutrophils", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

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
    '0':'Cstdc4 Neutrophils',
    '1':'Fabp5 Neutrophils',
    '2':'Slfn4 Neutrophils',
    '3':'Cebpb Neutrophils',
    '4':'Tnf Neutrophils',
    '5':'Tnf Neutrophils',
    '6':'Fabp5 Neutrophils',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Neutrophils.h5ad'))

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='mouse',output=out_dir)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

cross_tab = pd.crosstab(adata.obs['group'], adata.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Neutrophils_ROE/celltype_proportions_by_group.csv")


# # T cells细胞数柱状图
# adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/mapping_T/Tcells.h5ad")
# adata_raw.obs['cell_type'].cat.categories
# # 'CD4T exhaust', 'CD4T memory', 'CD4T naive', 'CD8T cytotoxic', 'CD8T memory', 'CD8T naive', 'NK', 'Proliferative Tcells', 'Th2', 'Treg', 'γδT'
# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_B/"

# group_counts = adata_raw.obs['group'].value_counts()
# group_counts = group_counts.sort_index()
# colors = ['gray'] + ['#E95137'] * (len(group_counts) - 1)
# plt.figure(figsize=(7, 7))
# bars = plt.bar(group_counts.index, group_counts.values, color=colors, edgecolor='black', width=0.6)
# for bar in bars:
#     height = bar.get_height()
#     plt.text(bar.get_x() + bar.get_width()/2., height + 5,
#              f'{int(height)}', ha='center', va='bottom', fontsize=10)

# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.ylabel('Total T/NK cell number', fontsize=12)
# plt.xticks(rotation=45, ha='right')
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, 'TandNK_cell_number_barplot.png'))
# plt.savefig(os.path.join(out_dir, 'TandNK_cell_number_barplot.pdf'))

# adata_sub=adata_raw[~adata_raw.obs['cell_type'].isin(["NK"]), :]
# group_counts = adata_sub.obs['group'].value_counts()
# group_counts = group_counts.sort_index()
# colors = ['gray'] + ['#E95137'] * (len(group_counts) - 1)
# plt.figure(figsize=(7, 7))
# bars = plt.bar(group_counts.index, group_counts.values, color=colors, edgecolor='black', width=0.6)
# for bar in bars:
#     height = bar.get_height()
#     plt.text(bar.get_x() + bar.get_width()/2., height + 5,
#              f'{int(height)}', ha='center', va='bottom', fontsize=10)

# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.ylabel('Total T cell number', fontsize=12)
# plt.xticks(rotation=45, ha='right')
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, 'T_cell_number_barplot.png'))
# plt.savefig(os.path.join(out_dir, 'T_cell_number_barplot.pdf'))

# adata_sub=adata_sub[~adata_sub.obs['cell_type'].isin(['γδT']), :]
# group_counts = adata_sub.obs['group'].value_counts()
# group_counts = group_counts.sort_index()
# colors = ['gray'] + ['#E95137'] * (len(group_counts) - 1)
# plt.figure(figsize=(7, 7))
# bars = plt.bar(group_counts.index, group_counts.values, color=colors, edgecolor='black', width=0.6)
# for bar in bars:
#     height = bar.get_height()
#     plt.text(bar.get_x() + bar.get_width()/2., height + 5,
#              f'{int(height)}', ha='center', va='bottom', fontsize=10)

# ax = plt.gca()
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.ylabel('Total T cell number', fontsize=12)
# plt.xticks(rotation=45, ha='right')
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, 'Tsub_cell_number_barplot.png'))
# plt.savefig(os.path.join(out_dir, 'Tsub_cell_number_barplot.pdf'))


# #### umap
# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_C/"
# adata_sub.obs['raw_cell_type']=adata_sub.obs['cell_type']
# anno_dict={
#     'CD4T exhaust':'CD4T cell',
#     'CD4T memory':'CD4T cell',
#     'CD4T naive':'CD4T cell',
#     'CD8T cytotoxic':'CD8T cell',
#     'CD8T memory':'CD8T cell',
#     'CD8T naive':'CD8T cell',
#     'Proliferative Tcells':'CD8T cell',
#     'Th2':'Th2 cell',
#     'Treg':'Treg cell'
# }
# adata_sub.obs['cell_type'] = [anno_dict[clust] for clust in adata_sub.obs['raw_cell_type']]
# adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype('category')
# adata_sub.uns['cell_type_colors']=['#ECA8A9', '#74AED4', '#D3E2B7', '#CFAFD4', '#F7C97E','#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#b5bd61', '#17becf', '#aec7e8']
# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False,legend_loc='on data')
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

# adata_sub.write(os.path.join(out_dir, 'Tsub_cells.h5ad'))


# # 饼图
# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_D/"
# # 计算每个细胞类型的细胞数
# celltype_counts = adata_sub.obs['cell_type'].value_counts()
# cell_types = celltype_counts.index.tolist()
# counts = celltype_counts.values
# colors = adata_sub.uns['cell_type_colors']
# fig, ax = plt.subplots(figsize=(7, 8))
# wedges, texts, autotexts = ax.pie(
#     counts, 
#     colors=colors,
#     autopct=lambda pct: f'{int(pct/100.*sum(counts))}',  # 显示细胞数而不是百分比
#     startangle=90,
#     textprops={'fontsize': 15, 'fontweight': 'bold'},
#     pctdistance=0.5  # 将数字放在饼图中心
# )
# ax.axis('equal')
# ax.legend(wedges, cell_types, 
#           loc='upper center', bbox_to_anchor=(0.5, -0.1), 
#           ncol=1, frameon=False, fontsize=25)
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, 'pieplot_celltype.png'))
# plt.savefig(os.path.join(out_dir, 'pieplot_celltype.pdf'))


# # marker_UMAP图
# from matplotlib.colors import LinearSegmentedColormap
# import matplotlib as mpl
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_E/"
# colors = ['white', 'red', 'black']
# cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)
# fig, ax = plt.subplots(figsize=(7, 7))
# sc.pl.umap(adata_sub, color=['Cd44'], cmap=cmap, frameon=False,ax=ax,colorbar_loc=None,show=False)
# plt.savefig(os.path.join(out_dir, 'umap_Cd44.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'umap_Cd44.pdf'),bbox_inches='tight')

# fig = plt.figure(figsize=(9, 7))
# main_ax = fig.add_axes([0.05, 0.05, 0.7, 0.9])  # [left, bottom, width, height]
# sc.pl.umap( adata_sub,  color=['Cd44'],  cmap=cmap,  frameon=False, ax=main_ax, colorbar_loc=None, show=False)
# # 颜色条位置（右侧）
# cbar_ax = fig.add_axes([0.78, 0.2, 0.03, 0.3])  # [left, bottom, width, height]
# norm = mpl.colors.Normalize(vmin=adata_sub.obs['Cd44'].min(), vmax=adata_sub.obs['Cd44'].max())
# cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
# cb.ax.text(1.2, 0, 'Low', transform=cb.ax.transAxes, fontsize=10, va='center')
# cb.ax.text(1.2, 1, 'High', transform=cb.ax.transAxes, fontsize=10, va='center')
# cb.ax.set_yticks([])
# quantiles = [0.25, 0.5, 0.75]
# for q in quantiles:
#     cb.ax.axhline(y=q, color='white', linewidth=1, linestyle='dashed', 
#                   xmin=0, xmax=1.11, clip_on=False)

# plt.savefig(os.path.join(out_dir, 'umap_Cd44_with_colorbar.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'umap_Cd44_with_colorbar.pdf'),bbox_inches='tight')


# #### marker 气泡图
# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_F/"
# marker_list = {'Naive T cell': ['Sell', 'Ccr7'], 'Activated T cell': ['Cd44', 'Prf1', 'Gzmb'], 'Exhausted T cell': ['Pdcd1']}
# with rc_context({'figure.figsize': (9, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.dotplot(adata_sub, var_names=marker_list, groupby='group', dendrogram=False, show=False,swap_axes=True)
#     plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'))
#     plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'))


#### 分组染色umap+饼图
adata_total = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
adata_total = adata_total[~adata_total.obs["cell_type"].isin(['T cells']),]
adata_sub=sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_C/Tsub_cells.h5ad")
cell_type_mapping = adata_sub.obs[['cell_type']].to_dict()['cell_type']
adata_total.obs['cell_type'] = adata_total.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_total.obs['cell_type'] = adata_total.obs['cell_type'].astype('category')
adata_total.obs['cell_type'].cat.categories

out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_G/"
cross_tab = pd.crosstab(adata_total.obs['group'], adata_total.obs['cell_type'])
proportions = cross_tab.div(cross_tab.sum(axis=1), axis=0) * 100
groups = ['asham', 'aday1', 'aday3', 'aday7', 'aday14']
group_mapping = {}
for group in groups:
    cd4t_prop = proportions.loc[group, 'CD4T cell'] if 'CD4T cell' in proportions.columns else 0
    cd8t_prop = proportions.loc[group, 'CD8T cell'] if 'CD8T cell' in proportions.columns else 0
    Th2_prop = proportions.loc[group, 'Th2'] if 'Th2' in proportions.columns else 0
    Treg_prop = proportions.loc[group, 'Treg'] if 'Treg' in proportions.columns else 0
    total_prop = cd4t_prop + cd8t_prop + Th2_prop + Treg_prop
    group_mapping[group] = f"{group} ({total_prop:.2f}%)"

adata_sub.obs['group_display'] = adata_sub.obs['group'].map(group_mapping)
adata_sub.obs['raw_cell_type']=adata_sub.obs['cell_type']
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype(str)
groups_display= ['asham (0.51%)' , 'aday1 (0.49%)' , 'aday3 (0.54%)' , 'aday7 (0.09%)' , 'aday14 (1.00%)']
keep_cell_types = ['CD4T cell', 'CD8T cell']
color_map = {'CD4T cell': '#ECA8A9', 'CD8T cell': '#74AED4', 'other': '#8080801A'}
if len(groups) > 1:
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(21, 14))
    axs = axs.flatten()
    for i, group in enumerate(groups):
        adata_sub.obs.loc[adata_sub.obs["group"] != group, 'cell_type'] = 'other'
        group_mask = adata_sub.obs["group"] == group
        adata_sub.obs.loc[group_mask, 'cell_type'] = adata_sub.obs.loc[group_mask, 'raw_cell_type'].apply(
            lambda x: x if x in keep_cell_types else 'other'
        )
        adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype('category')
        categories = adata_sub.obs['cell_type'].cat.categories
        adata_sub.uns['cell_type_colors'] = [color_map[ct] for ct in categories]
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_sub, color="cell_type", legend_loc=legend, title=groups_display[i], ax=axs[i], show=False, frameon=False,legend_fontsize=15)
        axs[i].set_title(groups_display[i], fontsize=20)
    for i in range(len(groups), len(axs)):
        axs[i].axis('off')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


adata_sub=sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_C/Tsub_cells.h5ad")
adata_sub= adata_sub[adata_sub.obs["cell_type"].isin(['CD4T cell', 'CD8T cell']),]
colors = adata_sub.uns['cell_type_colors']
for group in groups:
    adata_group = adata_sub[adata_sub.obs['group'] == group]
    celltype_counts = adata_group.obs['cell_type'].value_counts()
    cell_types = celltype_counts.index.tolist()
    counts = celltype_counts.values
    
    fig, ax = plt.subplots(figsize=(7, 7))
    wedges, texts, autotexts = ax.pie(
        counts, 
        colors=colors,
        autopct=lambda pct: f'{int(pct/100.*sum(counts))}',  # 显示细胞数
        startangle=90,
        textprops={'fontsize': 15, 'fontweight': 'bold'},
        pctdistance=0.5  # 数字放在饼图中心
    )
    ax.axis('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'pieplot_celltype_{group}.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, f'pieplot_celltype_{group}.pdf'), bbox_inches='tight')
    plt.close()  # 关闭图形，释放内存
    print(f"Saved pie plot for group: {group}")



# 对CD8T cells注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/CD8T/"
adata_sub = adata_sub[adata_sub.obs['cell_type'].isin(['CD8T cell']), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

anno_dict={
    '0':'Cmss1+',
    '1':'Cxcl2+',
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

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='group', title='Group', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

sc.tl.dendrogram(adata,groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

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
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/CD8T/celltype_proportions_by_group.csv")



#### 改大类注释顺序
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
# 'B cells', 'Endothelial cells', 'Erythrocytes', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Melanocytes', 'Neutrophils','Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Main/"
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'on data'
        else:
            legend = 'on data'
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, legend_fontsize=5, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.umap(adata_raw, color=['Cd79a','Ms4a1','Jchain','Cdh5','Pecam1','Flt1','Hbb-bs','Hba-a1','Hba-a2','Dcn','Col1a1','Pdgfra','Krt14','Lgals7','Krt5','Vcan','Cd14','Adgre1','Dct','Mlana','Pmel','Ly6g','S100a8','Csf3r','Notch3', 'Pdgfrb', 'Rgs5','Mpz','S100b','Mag','Acta2','Myh11','Tagln','Cd3d','Cd3e','Trac'], cmap='coolwarm',ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genelist.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_genelist.pdf'),bbox_inches='tight')

sc.pl.stacked_violin(adata_raw, var_names=['Cd79a','Ms4a1','Jchain','Cdh5','Pecam1','Flt1','Hbb-bs','Hba-a1','Hba-a2','Dcn','Col1a1','Pdgfra','Krt14','Lgals7','Krt5','Vcan','Cd14','Adgre1','Dct','Mlana','Pmel','Ly6g','S100a8','Csf3r','Notch3', 'Pdgfrb', 'Rgs5','Mpz','S100b','Mag','Acta2','Myh11','Tagln','Cd3d','Cd3e','Trac'], groupby='cell_type', show=False, swap_axes=True, dendrogram=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'genelist_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'genelist_violin.png'),bbox_inches='tight')