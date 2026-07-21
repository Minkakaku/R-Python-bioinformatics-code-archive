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

#调整大类注释umap
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/06-Cell_type_annotation/adata_anno.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/02.Mission/Main/"
adata_recluster = recluster(adata=adata_raw,species='human',neighbors=10,resolution=0.6,n_pcs = 7,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')
sc.pl.umap(adata, color='sample', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_sample.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_sample.pdf"),bbox_inches='tight')

sc.pl.umap(adata,color=['CD4','CD8A','CD8B','NKG7'], show=False)
plt.savefig(os.path.join(out_dir,"umap_NK.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"umap_NK.pdf"),bbox_inches='tight')
anno_dict={
    '14':'NK'
}
mask = adata.obs['leiden'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata.obs['cell_type'].cat.categories]
# 更新 categories
adata.obs['cell_type'] = adata.obs['cell_type'].cat.set_categories(new_categories)
adata.obs.loc[mask, 'cell_type'] = adata.obs.loc[mask, 'leiden'].map(anno_dict)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# 对Tumor-like cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/06-Cell_type_annotation/adata_anno.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Endothelial cells', 'Fibroblast', 'Mast cells', 'Myeloid cells', 'T cells', 'Tumor-like cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/02.Mission/Tumor-like cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Tumor-like cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.8,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

# sc.pl.dotplot(adata_test, var_names='S100A8', groupby='group', dendrogram=False,show=False)
# plt.savefig(os.path.join(out_dir, 'S100A8_adata_test-group.png'), facecolor='white', bbox_inches="tight", dpi=300)
# sc.pl.dotplot(adata_test, var_names='S100A8', groupby=['cell_type','group'], dendrogram=False, show=False)
# plt.savefig(os.path.join(out_dir, 'S100A8_adata_test-cell_type-group.png'), facecolor='white', bbox_inches="tight", dpi=300)
# pd.DataFrame(
#     adata_test.raw[:,['S100A8','S100A9']].X.toarray(),
#     index=adata_test.obs_names,
#     columns=['S100A8','S100A9']
# ).to_csv(os.path.join(out_dir, 'X_raw.xls'), sep='\t', index=True)

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sample = adata.obs['sample'].drop_duplicates().tolist()
n_sample = len(sample)
if len(sample) > 1:
    nrows = 4            # 想要的行数
    ncols = int(np.ceil(n_sample / nrows))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*6, nrows*6))
    axes = axes.flatten() 
    keys = adata.obs['leiden'].cat.categories.to_list()
    values = adata.uns['leiden_colors']
    color_dict = dict(zip(keys,values))
    for i, smp in enumerate(sample):
        tmp = adata[adata.obs["sample"] == smp]
        sc.pl.umap(tmp,size=20, color='leiden', legend_loc='on data', palette=color_dict, title=smp, ax=axes[i], show=False)
    save_fig(plt,out_dir=out_dir,file_name='umap_split_samples', fig_size={'w':ncols*6, 'h':nrows*6})

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
    '0':'HCC_C1',
    '1':'HCC_C2',
    '2':'HCC_C2',
    '3':'HCC_C3',
    '4':'Hepatocytes_C1',
    '5':'HCC_C2',
    '6':'Hepatocytes_C2',
    '7':'HCC_C4',
    '8':'HCC_C5',
    '9':'HCC_C5',
    '10':'HCC_C5',
    '11':'HCC_C3',
    '12':'HCC_C4',
    '13':'HCC_C3',
    '14':'HCC_C2',
    '15':'Hepatocytes_C1',
    '16':'HCC_C1',
    '17':'HCC_C1',
    '18':'HCC_C3',
    '19':'HCC_C2',
    '20':'HCC_C2'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Tumor-like cells.h5ad'))

groups = adata.obs['group'].drop_duplicates().tolist()
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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':22, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

sc.tl.rank_genes_groups(adata,groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)


# 对T cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/06-Cell_type_annotation/adata_anno.h5ad")
# 'Endothelial cells', 'Fibroblast', 'Mast cells', 'Myeloid cells', 'T cells', 'Tumor-like cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/T cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.8,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['CD4','CD8A','CD8B','FOXP3','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD4','CD8A','CD8B','FOXP3','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'))

anno_dict={
    '0':'CD8 Tex',
    '1':'CD8 CTLs',
    '2':'CD4 Tem',
    '3':'CD56bright NK',
    '4':'CD4 Tem',
    '5':'CD8 Tex',
    '6':'CD4 naive cells',
    '7':'CD56bright NK',
    '8':'CD56dim NK',
    '9':'CD4 Tcm',
    '10':'CD4 Tcm',
    '11':'CD56dim NK',
    '12':'Treg',
    '13':'Treg'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':22, 'h':5})
adata.write(os.path.join(out_dir, 'T.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Endothelial cells', 'Fibroblast', 'Mast cells', 'Myeloid cells', 'T cells', 'Tumor-like cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/Myeloid cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['CD1A','CD1C','CLEC4C','HLA-DRA','IL3RA','ITGAM','ITGAX','LAMP3','LILRA4','LYZ','CCL22','CLEC9A','FSCN1','ITGAE','TNIP3','WDFY4','BST2','CCR9','CD14','FCGR3A'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD1A','CD1C','CLEC4C','HLA-DRA','IL3RA','ITGAM','ITGAX','LAMP3','LILRA4','LYZ','CCL22','CLEC9A','FSCN1','ITGAE','TNIP3','WDFY4','BST2','CCR9','CD14','FCGR3A'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'))
sc.pl.umap(adata, color=['CSF3R','FCGR3B','CXCR2','VCAN','TNFRSF17','JCHAIN','CD79A'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_B.png'))
plt.savefig(os.path.join(out_dir, 'umap_B.pdf'))
sc.pl.umap(adata, color=['C1QA','CD68','CD86','APOE'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Macrophage.png'))
plt.savefig(os.path.join(out_dir, 'umap_Macrophage.pdf'))
sc.pl.umap(adata, color=['IGHG1', 'IGKC','JCHAIN','TNFRSF17'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma.pdf'))


anno_dict={
    '0':'Monocyte',
    '1':'DCs',
    '2':'Macrophage',
    '3':'Macrophage',
		'4':'Neutrophil',
    '5':'Neutrophil',
    '6':'Monocyte',
    '7':'Macrophage',
    '8':'B cells',
    '9':'Neutrophil',
    '10':'Neutrophil',
    '11':'Tumor-like cells',
    '12':'Macrophage',
    '13':'Plasma cells',
    '14':'Plasma cells',
    '15':'Hepatocyte'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/06-Cell_type_annotation/adata_anno.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/mapping_raw/"

adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/Tumor-like cells/Tumor-like cells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/T cells/T.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.obs['T_only'] = adata_raw.obs['sub_cell_type']
mask = adata_raw.obs['sub_cell_type'] != 'Tumor-like cells'
adata_raw.obs.loc[mask, 'T_only'] = np.nan
sc.pl.umap(adata_raw, color='T_only', show=False, size=15,ncols=5,na_color='lightgray')
save_fig(plt,file_name= 'T_only_umap',out_dir=out_dir,fig_size={'w':5,'h':5})

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['Tumor-like cells']), :]
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')


# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Endothelial cells', 'Fibroblast', 'Mast cells', 'Myeloid cells', 'T cells', 'Tumor-like cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/Myeloid cells/"
adata_sub = adata_raw[adata_raw.obs.sub_cell_type.isin(['Macrophage','Monocyte','DCs']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='level1_y', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pl.umap(adata, color=['CD1A','CD1C','CLEC4C','HLA-DRA','IL3RA','ITGAM','ITGAX','LAMP3','LILRA4','LYZ','CCL22','CLEC9A','FSCN1','ITGAE','TNIP3','WDFY4','BST2','CCR9','CD14','FCGR3A'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD1A','CD1C','CLEC4C','HLA-DRA','IL3RA','ITGAM','ITGAX','LAMP3','LILRA4','LYZ','CCL22','CLEC9A','FSCN1','ITGAE','TNIP3','WDFY4','BST2','CCR9','CD14','FCGR3A'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'))
sc.pl.umap(adata, color=['ITGAM','CD33','ANPEP','FUT4','CEACAM8','LIN7A'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_mdsc.png'))
plt.savefig(os.path.join(out_dir, 'umap_mdsc.pdf'))
sc.pl.dotplot(adata, var_names=['ITGAM','CD33','ANPEP','FUT4','CEACAM8','LIN7A'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_mdsc.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_mdsc.pdf'))
sc.pl.umap(adata, color=['PDLIM3','PAM','PDLIM7','FSCN1','DPYSL2','ARID5B','LGALS3','KLF2'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_M2.png'))
plt.savefig(os.path.join(out_dir, 'umap_M2.pdf'))

anno_dict={
    '0':'MDSCs',
    '1':'TAM_1',
    '2':'DCs',
    '3':'TAM_2',
		'4':'TAM_3',
    '5':'Monocyte',
    '6':'Monocyte',
    '7':'TAM_4',
    '8':'TAM_5',
    '9':'DCs',
    '10':'DCs',
    '11':'MDSCs',
    '12':'MDSCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':22, 'h':5})
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')


# 运行 cytotrace2
adata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/Myeloid cells/Myeloid.h5ad")
adata.to_df().T.to_csv('/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace2/Myeloid_count.txt',sep='\t')
adata.obs['cell_type'].to_csv('/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace2/Myeloid_celltype.txt',sep='\t')
#运行shell
# cytotrace2 \
#    -f "/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace2/Myeloid_count.txt" \
#    -a "/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace2/Myeloid_celltype.txt" \
#    -sp human
import cytotrace2 as ct
from cytotrace2_py.cytotrace2_py import *
#adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/Myeloid cells/Myeloid.h5ad")
#if adata.raw is not None:
#     adata = adata.raw.to_adata()
# adata.var['gene'] = adata.var_names.astype(str).str.strip()
# mask = adata.var['gene'].replace('', np.nan).notna()
# adata = adata[:, mask].copy()
#adata = adata[:, ~adata.var['gene'].duplicated(keep='first')].copy()
#mat = adata.to_df().T.astype(int)
#mat.index = adata.var['gene']
#mat.index.name = ''
#mat.to_csv("/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace/Myeloid_count.tsv", sep='\t')
input_path = "/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace/Myeloid_count.tsv"
# 运行 CytoTRACE 2 主函数
results =  cytotrace2(input_path,annotation_path=example_annotation_path,species=example_species)
# 查看预测结果
print(results.head())
# 可视化结果
adata.obs['CytoTRACE_Score'] = results['CytoTRACE2_Score']
adata.obs['Potency']    = results['CytoTRACE2_Potency']
adata.obs['CytoTRACE2_Relative'] = results['CytoTRACE2_Relative']
adata.obs['preKNN_CytoTRACE2_Score']    = results['preKNN_CytoTRACE2_Score']
adata.obs['preKNN_CytoTRACE2_Potency'] = results['preKNN_CytoTRACE2_Potency']
adata.write('/mnt/NexData_file/home/sh/Projects/FWAE20254246/01.Mission/cytotrace/Myeloid_CytoTRACE.h5ad')