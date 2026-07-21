import os, sys
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
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
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/05.Mission/Epithelial cells/"
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/Epithelial cells/Epithelial cells.h5ad")
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_celltype', fig_size={'w':6, 'h':5})

anno_dict={
    '0':'Stem cells',
    '1':'Stem cells',
    '2':'Stem cells',
    '3':'Tumor-like cells',
    '4':'Tumor-like cells',
    '5':'Stem cells',
    '6':'Tumor-like cells',
    '7':'Goblet cells',
    '8':'Colonocyte',
    '9':'Tumor-like cells',
    '10':'Stem cells',
    '11':'Colonocyte',
    '12':'Tumor-like cells',
    '13':'Stem cells',
    '14':'Tumor-like cells',
    '15':'Tumor-like cells',
    '16':'Tumor-like cells',
    '17':'Tuft cells',
    '18':'TA cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_celltype', fig_size={'w':6, 'h':5})
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':15, 'h':5})

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='human',output=out_dir)

# sc.pl.dotplot(adata, var_names=['ALPI','SULT1B1','MEP1A','REG4','TFF3','SMIM14','SPDEF','B3GNT6','HMGCS2','MKI67','TOP2A','CENPF','COL1A1','COL1A2','COL3A1','RNF43','MLXIPL','AC004687.2'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
# plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Epithelial cells.h5ad'))

groups = ['LARC', 'N']
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


# 对Tumor-like cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/05.Mission/Epithelial cells/Epithelial cells.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Colonocyte', 'Goblet cells', 'Stem cells', 'TA cells', 'Tuft cells', 'Tumor-like cells'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/05.Mission/Tumor/"
adata_sub=adata_raw[adata_raw.obs['cell_type']=='Tumor-like cells']
# adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/Epithelial cells/Epithelial cells.h5ad")
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color='sample', title='Sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':6, 'h':5})

#### 后面没跑
anno_dict={
    '0':'Stem cells',
    '1':'Stem cells',
    '2':'Stem cells',
    '3':'Tumor-like cells',
    '4':'Tumor-like cells',
    '5':'Stem cells',
    '6':'Tumor-like cells',
    '7':'Goblet cells',
    '8':'Colonocyte',
    '9':'Tumor-like cells',
    '10':'Stem cells',
    '11':'Colonocyte',
    '12':'Tumor-like cells',
    '13':'Stem cells',
    '14':'Tumor-like cells',
    '15':'Tumor-like cells',
    '16':'Tumor-like cells',
    '17':'Tuft cells',
    '18':'TA cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_celltype', fig_size={'w':6, 'h':5})
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':15, 'h':5})

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='human',output=out_dir)

# sc.pl.dotplot(adata, var_names=['ALPI','SULT1B1','MEP1A','REG4','TFF3','SMIM14','SPDEF','B3GNT6','HMGCS2','MKI67','TOP2A','CENPF','COL1A1','COL1A2','COL3A1','RNF43','MLXIPL','AC004687.2'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
# plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Epithelial cells.h5ad'))

groups = ['LARC', 'N']
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


# 映射所有细分 #### 后面还没映射
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Bcells/Bcells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD4T/CD4T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD8T/CD8T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Endothelial/EC.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Epithelial cells/Epithelial cells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Fibroblasts/Fibroblasts.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Plasma/Plasma.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')