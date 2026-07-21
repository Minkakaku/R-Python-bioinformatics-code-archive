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

# 对T cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Erythrocytes', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Melanocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/Tcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

anno_dict={
    '0':'T cells',
    '1':'T cells',
    '2':'T cells',
    '3':'T cells',
    '4':'T cells',
    '5':'T cells',
    '6':'T cells',
    '7':'T cells',
    '8':'Mast cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[adata.obs.cell_type=="T cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 2)
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

sc.pl.umap(adata, color=['Cd4','Cd8a','Cd8b1','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Trdc',"Trgc1", "Trgc2",'Mki67'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Cd8b1','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Trdc',"Trgc1", "Trgc2",'Mki67'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'γδT',
    '1':'γδT',
    '2':'γδT',
    '3':'CD4T',
    '4':'NK',
    '5':'CD8T',
    '6':'Doublet',
    '7':'Doublet',
    '8':'CD4T',
    '9':'CD4T',
    '10':'Doublet'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

adata_sub = adata[~adata.obs.cell_type.isin(['Doublet']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1.2,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Cd4','Cd8a','Cd8b1','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Trdc',"Trgc1", "Trgc2",'Mki67'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Cd8b1','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Trdc',"Trgc1", "Trgc2",'Mki67'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'γδT',
    '1':'γδT',
    '2':'γδT',
    '3':'CD4T',
    '4':'NK',
    '5':'γδT',
    '6':'γδT',
    '7':'CD8T',
    '8':'CD4T',
    '9':'γδT',
    '10':'γδT',
    '11':'CD4T',
    '12':'CD4T',
    '13':'γδT'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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
adata.write(os.path.join(out_dir, 'Tcells.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)
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
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/Tcells_ROE/celltype_proportions_by_group.csv")


# 对CD4T\CD8T cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/Tcells/Tcells.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'CD4T', 'CD8T', 'NK', 'γδT'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/CD4_8T/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["CD4T",'CD8T']), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.8,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

anno_dict={
    '0':'CD8T',
    '1':'CD8T',
    '2':'CD4T',
    '3':'CD4T',
    '4':'CD8T',
    '5':'Treg',
    '6':'Th2',
    '7':'Proliferative Tcells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_raw.write(os.path.join(out_dir, 'Tcells.h5ad'))
#### 做亚群再细分 ####
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/CD4_8T_sub/"
adata_sub = adata[adata.obs['cell_type'].isin(["CD4T",'CD8T']), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=1,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

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

sc.pl.umap(adata, color=['Cd4','Cd8a','Cd8b1','Sell','Ccr7','Gzmb','Prf1','Pdcd1','Cd69','Rgs1','Tox','Pdcd1','Havcr2'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Cd8b1','Sell','Ccr7','Gzmb','Prf1','Pdcd1','Cd69','Rgs1','Tox','Pdcd1','Havcr2'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'CD8T cytotoxic',
    '1':'CD8T memory',
    '2':'CD4T naive',
    '3':'CD4T memory',
    '4':'CD8T cytotoxic',
    '5':'CD8T memory',
    '6':'CD4T memory',
    '7':'CD4T exhaust',
    '8':'CD8T naive',
    '9':'CD4T memory',
    '10':'CD4T naive',
    '11':'CD4T exhaust'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

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
adata.write(os.path.join(out_dir, 'CD4_8T_sub.h5ad'))

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='mouse',output=out_dir)
# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/mapping_T/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_raw.obs['group'] = (
    pd.Categorical(adata_raw.obs['group'],
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
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata_raw.write(os.path.join(out_dir, 'Tcells.h5ad'))

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)
sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/mapping_T/celltype_proportions_by_group.csv")