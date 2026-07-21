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
    keys = adata.obs['cluster'].cat.categories.to_list()
    values = list(adata.uns['cluster_colors'])
    cluster_colors = dict(zip(keys, values))
    ax = pd.crosstab(adata.obs[type], adata.obs['cluster'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cluster_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cluster'], normalize='index').columns],
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
    file_name=f"cluster_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20242199/result.h5ad")
adata_raw.obs['cluster'].cat.categories
# 'B cells', 'Cycling cells', 'Endothelial cells', 'Myeloid cells', 'T cells', 'Unkown'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/03.Mission/T cells/"

adata=adata_raw[adata_raw.obs['cluster'].isin(['1','2','3','4','5','7','11']), :]

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cluster', title='cluster', show=False ,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['control', 'Test']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cluster", legend_loc='on data', title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs['cluster_group'] = adata.obs['cluster']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        adata_sub = adata.copy()
        adata_sub.obs.loc[adata_sub.obs['group'] != groups[i], 'cluster_group'] = np.nan
        sc.pl.umap(adata_sub, color="cluster_group", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# adata.obs['cluster_group'] = adata.obs['cluster']
orig_cats = adata.obs['cluster'].cat.categories
orig_colors = dict(zip(orig_cats,adata.uns.get('cluster_colors',sns.color_palette('tab20', len(orig_cats)).as_hex())))
new_cats = list(orig_cats) + [""]
palette = [orig_colors[c] for c in orig_cats] + ['lightgray']
new_order = ['1', '2', '3', '4', '5', '7', '11', '']

if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        adata_sub = adata.copy()
        cluster_str = adata.obs['cluster'].astype(str)
        cluster_str[adata.obs['group'] != groups[i]] = ""      # 空字符串
        adata_sub.obs['cluster_group'] = cluster_str.astype(str).astype('category')
        adata_sub.obs['cluster_group'] = adata_sub.obs['cluster_group'].cat.reorder_categories(new_order, ordered=True)
        sc.pl.umap(adata_sub, color="cluster_group", palette=palette,legend_loc='on data', title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'result.h5ad'))

method = "cluster"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cluster')



adata=adata_raw[adata_raw.obs['cell_type'].isin(['T cells']), :]
adata=adata.raw.to_adata()
marker_genes=['Isg15','Ifit1','Bst2','Stat1','Irf1','Cxcr4']
output_path="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/03.Mission/markerplot_group/"

for gene in marker_genes:
    adata.obs[f'{gene}_gene'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'Test'
    adata.obs.loc[mask, f'{gene}_gene'] = np.nan

genelist=['Isg15_gene','Ifit1_gene','Bst2_gene','Stat1_gene','Irf1_gene','Cxcr4_gene']
sc.pl.umap(adata, color=genelist, show=False, cmap='Purples',size=15,ncols=6,na_color='lightgray')
save_fig(plt,file_name= 'Test_umap',out_dir=output_path,fig_size={'w':36,'h':6})

for gene in marker_genes:
    adata.obs[f'{gene}_gene'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'control'
    adata.obs.loc[mask, f'{gene}_gene'] = np.nan

sc.pl.umap(adata, color=genelist, show=False, cmap='Purples',size=15,ncols=6,na_color='lightgray')
save_fig(plt,file_name= 'control_umap',out_dir=output_path,fig_size={'w':36,'h':6})



genelist=['Cd4','Tcf7','Lef1','S1pr1']
genelist=['Cd8b1','Cd8a','Sell','Cd44']
genelist=['Foxp3','Icos','Il7r']
genelist=['Ccr7','Il7r','Cd27']
genelist=['Isg15','Ifit1','Bst2','Stat1','Irf1','Cxcr4']
genelist=['Eomes','Il2rb','Gzmk','Prf1']
genelist=['Il17re','Rorc','Il22','Il23r']
sc.pl.umap(adata, color=genelist, show=False, size=15,ncols=4)
save_fig(plt,file_name= 'Cluster1_umap',out_dir=output_path,fig_size={'w':24,'h':6})
save_fig(plt,file_name= 'Cluster2_umap',out_dir=output_path,fig_size={'w':24,'h':6})
save_fig(plt,file_name= 'Cluster3_umap',out_dir=output_path,fig_size={'w':18,'h':6})
save_fig(plt,file_name= 'Cluster4_umap',out_dir=output_path,fig_size={'w':18,'h':6})
save_fig(plt,file_name= 'Cluster5_umap',out_dir=output_path,fig_size={'w':36,'h':6})
save_fig(plt,file_name= 'Cluster7_umap',out_dir=output_path,fig_size={'w':24,'h':6})
save_fig(plt,file_name= 'Cluster11_umap',out_dir=output_path,fig_size={'w':24,'h':6})

marker_genes={'Cluster1':['Cd4','Tcf7','Lef1','S1pr1'],'Cluster2':['Cd8b1','Cd8a','Sell','Cd44'],'Cluster3':['Foxp3','Icos','Il7r'],'Cluster4':['Ccr7','Il7r','Cd27'],'Cluster5':['Isg15','Ifit1','Bst2','Stat1','Irf1','Cxcr4'],'Cluster7':['Eomes','Il2rb','Gzmk','Prf1'],'Cluster11':['Il17re','Rorc','Il22','Il23r']}
sc.tl.dendrogram(adata,groupby='cluster')
sc.pl.dotplot(adata, var_names=marker_genes, groupby='cluster', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'T_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'T_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.matrixplot(adata, marker_genes, groupby='cluster', cmap='coolwarm', dendrogram=False, show=False, swap_axes=False)
plt.savefig(os.path.join(output_path, 'T_matrixplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'T_matrixplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '1':'C1',
    '2':'C2',
    '3':'C3',
    '4':'C4',
    '5':'C5',
    '7':'C7',
    '11':'C11'
}
adata.obs['celltype'] = [anno_dict[clust] for clust in adata.obs['cluster']]
adata.obs['celltype'] = adata.obs['celltype'].astype('category')
new_order = ['C1','C2','C3','C4','C5','C7','C11']
adata.obs['celltype'] = adata.obs['celltype'].cat.reorder_categories(
    new_order, ordered=True)
marker_genes=['Cd3d','Cd3e','Cd3g','Cd4','Tcf7','Lef1','S1pr1','Cd8b1','Cd8a','Sell','Cd44','Foxp3','Icos','Il7r','Ccr7','Il7r','Cd27','Isg15','Ifit1','Bst2','Stat1','Irf1','Cxcr4','Eomes','Il2rb','Gzmk','Prf1','Il17re','Rorc','Il22','Il23r']
sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='celltype', show=False, swap_axes=False, dendrogram=False)
plt.savefig(os.path.join(output_path, 'T_stacked_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'T_stacked_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/02.Mission/pieplot/"
adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
sc.pl.umap(adata, color='cluster', show=False)
save_fig(plt,file_name= 'cluster_umap',out_dir=out_dir,fig_size={'w':7,'h':7})
sc.pl.umap(adata, color='group', show=False)
save_fig(plt,file_name= 'group_umap',out_dir=out_dir,fig_size={'w':7,'h':7})

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')


marker_genes={'Cluster1':['Cd4','Lef1','S1pr1'],'Cluster2':['Cd8b1','Cd8a','Sell','Cd44'],'Cluster3':['Foxp3','Icos','Il7r'],'Cluster4':['Ccr7','Cd27','Cd69','Tnfrsf18'],'Cluster5':['Isg15','Ifit1','Bst2','Stat1','Irf1','Cxcr4'],'Cluster7':['Eomes','Il2rb','Gzmk','Prf1'],'Cluster11':['Il17re','Rorc','Il22','Il23r']}
sc.tl.dendrogram(adata,groupby='cluster')
sc.pl.dotplot(adata, var_names=marker_genes, groupby='cluster', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'T_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'T_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)