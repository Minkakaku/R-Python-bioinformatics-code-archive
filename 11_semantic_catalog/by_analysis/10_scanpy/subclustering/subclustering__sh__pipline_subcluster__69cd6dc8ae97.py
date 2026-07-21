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
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/02.Mission/Main/"

adata_raw.obs['raw_cluster']=adata_raw.obs['cluster']
anno_dict={
    '0':'0',
    '1':'1',
    '2':'2',
    '3':'3',
    '4':'1',
    '5':'5',
    '6':'6',
    '7':'7',
    '8':'8',
    '9':'9',
    '10':'10',
    '11':'4',
    '12':'11'
}
adata_raw.obs['cluster'] = [anno_dict[clust] for clust in adata_raw.obs['cluster']]
adata_raw.obs['cluster'] = adata_raw.obs['cluster'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cluster', title='cluster', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata_raw.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cluster", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample')
adata_raw.write(os.path.join(out_dir, 'result.h5ad'))

adata=adata_raw[adata_raw.obs['cell_type'].isin(['T cells']), :]
adata=adata.raw.to_adata()
marker_genes=['Isg15','Ifit1','Ifi27l2a','Bst2','Ifi203']
output_path="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/02.Mission/markerplot_group/"

for gene in marker_genes:
    adata.obs[f'{gene}_gene'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'Test'
    adata.obs.loc[mask, f'{gene}_gene'] = np.nan

genelist=['Isg15_gene','Ifit1_gene','Ifi27l2a_gene','Bst2_gene','Ifi203_gene']
sc.pl.umap(adata, color=genelist, show=False, cmap='Purples',size=15,ncols=5,na_color='lightgray')
save_fig(plt,file_name= 'Test_umap',out_dir=output_path,fig_size={'w':25,'h':5})

for gene in marker_genes:
    adata.obs[f'{gene}_gene'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'control'
    adata.obs.loc[mask, f'{gene}_gene'] = np.nan

sc.pl.umap(adata, color=genelist, show=False, cmap='Purples',size=15,ncols=5,na_color='lightgray')
save_fig(plt,file_name= 'control_umap',out_dir=output_path,fig_size={'w':25,'h':5})

for gene in marker_genes:
    adata.obs[f'{gene}_gene'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'Test'
    adata.obs.loc[mask, f'{gene}_gene'] = np.nan

sc.pl.umap(adata, color=genelist, show=False,size=15,ncols=5,na_color='lightgray')
save_fig(plt,file_name= 'rawcolor_Test_umap',out_dir=output_path,fig_size={'w':25,'h':5})

for gene in marker_genes:
    adata.obs[f'{gene}_gene'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'control'
    adata.obs.loc[mask, f'{gene}_gene'] = np.nan

sc.pl.umap(adata, color=genelist, show=False,size=15,ncols=5,na_color='lightgray')
save_fig(plt,file_name= 'rawcolor_control_umap',out_dir=output_path,fig_size={'w':25,'h':5})

# sc.pl.umap(adata, color='group', show=False, cmap='coolwarm',ncols=5,groups=(['Test']),na_color='lightgray')
# save_fig(plt,file_name= 'group_umap',out_dir=output_path,fig_size={'w':5,'h':5})

adata_raw=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/02.Mission/Main/result.h5ad")
marker_genes={'B cell':['Cd79a','Cd79b','Ms4a1','Bank1'],'T cell':['Cd3d','Cd3e','Cd3g','Cd247'],'Endothelial cell':['Cdh5','Plvap','Eng','Apoe'],'Myeloid cell':['Tyrobp','Fcer1g','Spi1','Ctsz']}
output_path="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/02.Mission/Main_marker/"
sc.tl.dendrogram(adata_raw,groupby='cluster')
sc.pl.dotplot(adata_raw, var_names=marker_genes, groupby='cluster', dendrogram=True, show=False)
save_fig(plt,file_name= 'Main_dotplot',out_dir=output_path,fig_size={'w':7,'h':7})

sc.pl.matrixplot(adata_raw, marker_genes, groupby='cluster', cmap='coolwarm', dendrogram=True, show=False, swap_axes=False)
save_fig(plt,file_name= 'Main_matrixplot',out_dir=output_path,fig_size={'w':7,'h':7})


adata=adata_raw[adata_raw.obs['cell_type'].isin(['T cells']), :]
output_path="/mnt/NexData_file/home/sh/Projects/FWSC20242199-02/02.Mission/T_marker/"

genelist=['Cd4','Lef1','Tcf7','S1pr1']
genelist=['Cd8a','Cxcr6','Cd69','Bcl2']
genelist=['Foxp3','Icos','Ccr7','Il7r']
genelist=['Isg15','Ifit1','Stat1','Usp18','Irf7']
genelist=['Eomes','Il2rb','Gzmk','Prf1']
genelist=['Il17re','Rorc','Il22','Il23r']
sc.pl.umap(adata, color=genelist, show=False, size=15,ncols=4)
save_fig(plt,file_name= 'Cluster1_umap',out_dir=output_path,fig_size={'w':20,'h':5})
save_fig(plt,file_name= 'Cluster2_umap',out_dir=output_path,fig_size={'w':20,'h':5})
save_fig(plt,file_name= 'Cluster3_umap',out_dir=output_path,fig_size={'w':20,'h':5})
save_fig(plt,file_name= 'Cluster5_umap',out_dir=output_path,fig_size={'w':25,'h':5})
save_fig(plt,file_name= 'Cluster7_umap',out_dir=output_path,fig_size={'w':20,'h':5})
save_fig(plt,file_name= 'Cluster4_umap',out_dir=output_path,fig_size={'w':20,'h':5})

marker_genes={'Cluster1':['Cd4','Lef1','Tcf7','S1pr1'],'Cluster2':['Cd8a','Cxcr6','Cd69','Bcl2'],'Cluster3':['Foxp3','Icos','Ccr7','Il7r'],'Cluster4':['Il17re','Rorc','Il22','Il23r'],'Cluster5':['Isg15','Ifit1','Stat1','Usp18','Irf7'],'Cluster7':['Eomes','Il2rb','Gzmk','Prf1']}
sc.tl.dendrogram(adata,groupby='cluster')
sc.pl.dotplot(adata, var_names=marker_genes, groupby='cluster', dendrogram=False, show=False)
save_fig(plt,file_name= 'T_dotplot',out_dir=output_path,fig_size={'w':7,'h':6})

sc.pl.matrixplot(adata, marker_genes, groupby='cluster', cmap='coolwarm', dendrogram=False, show=False, swap_axes=False)
save_fig(plt,file_name= 'T_matrixplot',out_dir=output_path,fig_size={'w':7,'h':6})

anno_dict={
    '1':'C1',
    '2':'C2',
    '3':'C3',
    '4':'C4',
    '5':'C5',
    '7':'C7'
}
adata.obs['celltype'] = [anno_dict[clust] for clust in adata.obs['cluster']]
adata.obs['celltype'] = adata.obs['celltype'].astype('category')
marker_genes=['Cd3d','Cd3e','Cd3g','Cd4','Cd8a','Sell','Cd44','Lef1','Tcf7','Cxcr6','Cd69','Foxp3','Icos','Isg15','Stat1','Ifit1','Eomes','Gzmk','Prf1','Il17re','Rorc']
sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='celltype', show=False, swap_axes=False, dendrogram=False)
save_fig(plt,file_name= 'T_stacked_violin',out_dir=output_path,fig_size={'w':8,'h':5})


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