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

# 对Epithelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510215750/result.h5ad")
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.add_categories('Plasma cells')
mask = adata_raw.obs['cluster'] == '22'
adata_raw.obs.loc[mask, 'cell_type'] = 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('GAStis', 'Tumor')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('GASPtis', 'NonTumor')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('GASpbmc', 'PBMC')
adata_raw.obs['group'] = (
    pd.Categorical(adata_raw.obs['group'],
                   categories=['Tumor','NonTumor','PBMC'],
                   ordered=True)
)
groups = adata_raw.obs['group'].drop_duplicates().tolist()
groups = ['Tumor','NonTumor','PBMC']
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

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
adata_raw.write(os.path.join(out_dir, 'Main.h5ad'))

# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Epithelial/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Epithelial cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

adata = adata[adata.obs.leiden!="11", :]
adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))


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
    '0':'SCGB2A1 malignant epi',
    '1':'non-malignant epi-like',
    '2':'PGC malignant epi',
    '3':'UBE2D2 malignant epi',
    '4':'SCGB2A1 malignant epi',
    '5':'PGC malignant epi',
    '6':'non-malignant epi-like',
    '7':'PGC malignant epi',
    '8':'UBE2D2 malignant epi',
    '9':'HBB malignant epi',
    '10':'non-malignant epi-like'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Epithelial cells.h5ad'))

groups = ['Tumor','NonTumor','PBMC']
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

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

# method = "cell_type"
# keys = adata.obs[method].cat.categories.to_list()
# values = list(adata.uns[method+'_colors'])
# color_dict = dict(zip(keys,values))
# df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# # df.columns.names = ["cell_type"]
# df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对T cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/Main.h5ad")
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Tcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
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

sc.pl.umap(adata, color=['CD4','CD8A','FOXP3','CCR7','SELL','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA','TOX','MKI67','IL17A','CD69','TRDC','TRGV9'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['CD4','CD8A','FOXP3','CCR7','SELL','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA','TOX','MKI67','IL17A','CD69','TRDC','TRGV9'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'))

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'CD4T naive',
    '1':'NK',
    '2':'CD8T memory',
    '3':'CD8T naive',
    '4':'CD4T cytotoxic',
    '5':'CD4T naive',
    '6':'NK',
    '7':'CD4T naive',
    '8':'GDT',
    '9':'CD8T cytotoxic',
    '10':'lowQ',
    '11':'CD8T memory'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata = adata[~adata.obs['leiden'].isin(["10"]), :]

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['Tumor','NonTumor','PBMC']
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

sc.pl.dotplot(adata, var_names=['CD4','CD8A','FOXP3','CCR7','SELL','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA','TOX','MKI67','IL17A','CD69','TRDC','TRGV9'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")

# method = "cell_type"
# keys = adata.obs[method].cat.categories.to_list()
# values = list(adata.uns[method+'_colors'])
# color_dict = dict(zip(keys,values))
# df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# # df.columns.names = ["cell_type"]
# df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')



# 对Monocytes注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/Main.h5ad")
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Monocytes/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Monocytes", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
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

sc.pl.umap(adata, color=['CD86','MRC1','CLEC4A','IL3RA','ITGAM','ITGAX','LAMP3','LYZ','CLEC9A','ITGAE','WDFY4','BST2','VCAN','CD14','FCGR3A'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['CD86','MRC1','CLEC4A','IL3RA','ITGAM','ITGAX','LAMP3','LYZ','CLEC9A','ITGAE','WDFY4','BST2','VCAN','CD14','FCGR3A'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Classical Monocyte',
    '1':'Classical Monocyte',
    '2':'Non-classical Monocyte',
    '3':'Macrophage',
		'4':'DC',
    '5':'lowQ',
    '6':'lowQ'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata = adata[~adata.obs['leiden'].isin(["5","6"]), :]

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['Tumor','NonTumor','PBMC']
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
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))


# # 对B注释细分
# adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Bcells/"
# adata_sub = adata_raw[adata_raw.obs.cell_type=="B", :]
# adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
# adata = adata_recluster.sub_cluster()

# adata.obs['cell_type']=adata.obs['leiden']
# cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
# adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# anno_dict={
#     '0':'Naive B',
#     '1':'Naive B',
#     '2':'Memory B',
#     '3':'lowQ',
#     '4':'GCB'
# }
# adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
# adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# adata_sub = adata[~adata.obs.cell_type.isin(['lowQ']), :]
# adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
# adata = adata_recluster.sub_cluster()

# sc.tl.dendrogram(adata,groupby='leiden')
# sc.pl.rank_genes_groups_dotplot(
#     adata,
#     groupby='leiden',
#     key='rank_genes_groups',
#     n_genes=10,
#     show=False)
# plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
# plt.close()

# sc.pl.umap(adata, color=['Cd27','Cd38','Fcer2a','Tcl1','Il4ra'], cmap='coolwarm')
# plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.png'))
# plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.pdf'))
# sc.pl.dotplot(adata, var_names=['Cd27','Cd38','Fcer2a','Tcl1','Il4ra'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
# plt.savefig(os.path.join(out_dir, 'dotplot_NaiveB.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'dotplot_NaiveB.pdf'),bbox_inches='tight')
# sc.pl.umap(adata, color=['Cd27','Cd38','Aim2','Tnfrsf13b'], cmap='coolwarm')
# plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.png'))
# plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.pdf'))
# sc.pl.dotplot(adata, var_names=['Cd27','Cd38','Aim2','Tnfrsf13b'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
# plt.savefig(os.path.join(out_dir, 'dotplot_MemoryB.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'dotplot_MemoryB.pdf'),bbox_inches='tight')
# sc.pl.umap(adata, color=['Pcna','Cd19','Ms4a1','Cd27','Ptprc','Cd38','Cd40','Cd83','Tnfrsf13b','Tnfrsf13c','Cd69'], cmap='coolwarm')
# plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.png'))
# plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.pdf'))
# sc.pl.dotplot(adata, var_names=['Pcna','Cd19','Ms4a1','Cd27','Ptprc','Cd38','Cd40','Cd83','Tnfrsf13b','Tnfrsf13c','Cd69'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
# plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.pdf'),bbox_inches='tight')

# anno_dict={
#     '0':'MemoryB',
#     '1':'NaiveB',
#     '2':'NaiveB',
#     '3':'ActivatedB'
# }
# adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
# adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
# if len(groups) > 1:
#     axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
#     fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
#     for i in range(len(groups)):
#         if i == len(groups) - 1:
#             legend = 'right margin'
#         else:
#             legend = None
#         sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
# plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
# plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
# adata.write(os.path.join(out_dir, 'B.h5ad'))


# 对Fibroblast注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/Main.h5ad")
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Fibroblasts/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
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

sc.pl.umap(adata, color=['MMP11', 'CDH11','POSTN','COL1A1','CFD','C3', 'CD34','CXCL12', 'CD248','ACTA2', 'MCAM','RGS5','PDPN', 'MME', 'ENO1', 'NDRG1','HSPH1','IDO1','CD74', 'CCL19','MKI67'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['MMP11', 'CDH11','POSTN','COL1A1','CFD','C3', 'CD34','CXCL12', 'CD248','ACTA2', 'MCAM','RGS5','PDPN', 'MME', 'ENO1', 'NDRG1', 'HSPH1','IDO1','CD74', 'CCL19','MKI67'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_gene.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Matrix CAFs',
    '1':'Tumour-like CAFs CAFs',
    '2':'Matrix CAFs',
    '3':'Matrix CAFs',
    '4':'Inflammatory CAFs',
    '5':'Vascular CAFs',
    '6':'Matrix CAFs',
    '7':'Tumour-like CAFs CAFs',
    '8':'Antigen-presenting CAFs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['Tumor','NonTumor','PBMC']
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
adata.write(os.path.join(out_dir, 'Fibroblast.h5ad'))
sc.pl.dotplot(adata, var_names=['CD74','HLA-DRA','CFD','C3', 'CXCL12','MMP11', 'CDH11','COL1A1', 'HSPH1','ACTA2', 'MCAM'], groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'),bbox_inches='tight')


# 对Smooth muscle cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/Main.h5ad")
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/SMC/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Smooth muscle cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

adata = adata[~adata.obs['leiden'].isin(["8"]), :]

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
    '0':'FOS SMC',
    '1':'CARMN SMC',
    '2':'TMSB10 SMC',
    '3':'CFH SMC',
    '4':'FOS SMC',
    '5':'FOS SMC',
    '6':'CARMN SMC',
    '7':'TMSB10 SMC',
    '8':'lowQ'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['Tumor','NonTumor']
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
adata.write(os.path.join(out_dir, 'SMC.h5ad'))


# 对Endothelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/Main.h5ad")
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Neutrophils', 'Pericytes', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Endothelial/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
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

sc.pl.umap(adata, color=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'ACKR1 vein ECs',
    '1':'TAGLN vein ECs',
    '2':'MGP artery ECs',
    '3':'HSPG2 artery ECs',
    '4':'HSPG2 artery ECs',
    '5':'SLPI artery ECs',
    '6':'SLPI artery ECs',
    '7':'SLPI artery ECs',
    '8':'lymphatic ECs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['Tumor','NonTumor','PBMC']
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
adata.write(os.path.join(out_dir, 'Endothelial.h5ad'))


# # 对Plasma注释细分
# adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# # 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
# out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Plasma/"
# adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma", :]
# # adata_sub.uns['log1p']['base'] = 2
# adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.1,n_pcs = 10,output=out_dir)
# adata = adata_recluster.sub_cluster()

# sc.pl.rank_genes_groups_dotplot(
#     adata,
#     groupby='leiden',
#     key='rank_genes_groups',
#     n_genes=10,
#     show=False)
# plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
# plt.close()

# anno_dict={
#     '0':'Plasma C1',
#     '1':'Plasma C2',
#     '2':'Plasma C3'
# }
# adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
# adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
# if len(groups) > 1:
#     axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
#     fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
#     for i in range(len(groups)):
#         if i == len(groups) - 1:
#             legend = 'right margin'
#         else:
#             legend = None
#         sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
# plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
# plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
# adata.write(os.path.join(out_dir, 'Plasma.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Main/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Endothelial/Endothelial.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Epithelial/Epithelial cells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Fibroblasts/Fibroblast.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/SMC/SMC.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/Tcells/Tcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['Epithelial cells','Smooth muscle cells','T cells','Monocytes']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')