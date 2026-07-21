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
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Macrophages/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.4,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

anno_dict={
    '0':'Macrophage 1',
    '1':'Macrophage 2',
    '2':'Macrophage 3',
    '3':'Macrophage 4',
    '4':'Macrophage 5',
    '5':'Macrophage 6'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

#### 分组染色umap
groups=['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata.obs['group'] = ( pd.Categorical(adata.obs['group'], categories=groups, ordered=True))
group_counts = adata.obs['group'].value_counts()
group_mapping = {group: f"{group} ({count/len(adata)*100:.1f}%)" for group, count in group_counts.items()}

groups_display= ['asham (13.5%)','aday1 (46.8%)','aday3 (27.7%)','aday7 (6.3%)','aday14 (5.7%)']
if len(groups) > 1:
    fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(35, 7))
    axs = axs.flatten()
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups_display[i], ax=axs[i], show=False, frameon=False,legend_fontsize=15)
        axs[i].set_title(groups_display[i], fontsize=20)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))

from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
colors = ['white', 'red', 'black']
cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)
fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.umap(adata, color=['Cxcl16'], cmap=cmap, frameon=False,ax=ax,colorbar_loc=None,show=False,size=12,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_Cxcl16.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Cxcl16.pdf'),bbox_inches='tight')

fig = plt.figure(figsize=(9, 7))
main_ax = fig.add_axes([0.05, 0.05, 0.7, 0.9])  # [left, bottom, width, height]
sc.pl.umap( adata,  color=['Cxcl16'],  cmap=cmap,  frameon=False, ax=main_ax, colorbar_loc=None, show=False)
# 颜色条位置（右侧）
cbar_ax = fig.add_axes([0.78, 0.2, 0.03, 0.3])  # [left, bottom, width, height]
gene_expr = adata.raw[:, 'Cxcl16'].X
if hasattr(gene_expr, 'toarray'):
    gene_expr = gene_expr.toarray().flatten()
else:
    gene_expr = gene_expr.flatten()
# 从基因表达计算归一化范围

norm = mpl.colors.Normalize(vmin=np.min(gene_expr), vmax=np.max(gene_expr))
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
plt.close()

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots(figsize=(6, 4))
    sc.pl.violin(adata, keys=['Cxcl16'], groupby='group', show=False, ax=ax, 
             stripplot=False, density_norm='count', xlabel=None, 
             ylabel='Expression level', jitter=False)
# 计算每个组的均值并添加点
    mean_values = []
    for i, group in enumerate(groups):
        expr = adata.raw[adata.obs['group'] == group, 'Cxcl16'].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        mean_val = np.mean(expr)
        mean_values.append(mean_val)
        ax.scatter(i, mean_val, color='black', s=50, zorder=10)
    ax.plot(range(len(groups)), mean_values, color='black', linestyle='-', linewidth=2, zorder=5)
    ax.grid(False)
    plt.savefig(os.path.join(out_dir, 'violin_Cxcl16.png'))
    plt.savefig(os.path.join(out_dir, 'violin_Cxcl16.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['Cxcl16'], groupby='cell_type', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_Cxcl16.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_Cxcl16.pdf'))

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
proportion_df.to_csv(os.path.join(out_dir, 'celltype_proportions_by_group.csv'))


#### Fibroblasts 注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Erythrocytes', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Melanocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Fibroblasts/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.3,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

anno_dict={
    '0':'Fibroblast 1',
    '1':'Fibroblast 2',
    '2':'Fibroblast 3',
    '3':'Fibroblast 4',
    '4':'Fibroblast 5'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

#### 分组染色umap
groups=['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata.obs['group'] = ( pd.Categorical(adata.obs['group'], categories=groups, ordered=True))
group_counts = adata.obs['group'].value_counts()
group_mapping = {group: f"{group} ({count/len(adata)*100:.1f}%)" for group, count in group_counts.items()}

groups_display= ['asham (31.8%)', 'aday1 (20.6%)', 'aday3 (9.0%)', 'aday7 (18.3%)', 'aday14 (20.4%)']
if len(groups) > 1:
    fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(35, 7))
    axs = axs.flatten()
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups_display[i], ax=axs[i], show=False, frameon=False,legend_fontsize=15)
        axs[i].set_title(groups_display[i], fontsize=20)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

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
proportion_df.to_csv(os.path.join(out_dir, 'celltype_proportions_by_group.csv'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/07.Mission/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Macrophages/Macrophages.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].cat.add_categories('CD8T')
mask = adata_raw.obs['sub_cell_type'].isin(['CD8T cytotoxic', 'CD8T memory', 'CD8T naive'])
adata_raw.obs.loc[mask, 'sub_cell_type'] = 'CD8T'
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].cat.remove_unused_categories()
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))