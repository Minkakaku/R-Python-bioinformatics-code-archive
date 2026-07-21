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

# 对DC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Monocytes/Macrophage.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/pDC/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub = adata_raw[~adata_raw.obs['cell_type'].isin(["mregDC"]), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_raw,species='human',neighbors=10,resolution=2,n_pcs = 20,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata, color=['LILRA4','TCF4','PLD4','BCL11A'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['LILRA4','TCF4','PLD4','BCL11A'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

adata=sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/pDC/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/Monocytes/"

pdc_barcodes = adata.obs[adata.obs['leiden'] == '11'].index.tolist()  # 注意leiden可能是字符串
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.add_categories(['pDC'])
adata_raw.obs.loc[adata_raw.obs.index.isin(pdc_barcodes), 'cell_type'] = 'pDC'

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/03.Mission/DC/DC.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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
#     '0':'IGHA1 cDC',
#     '1':'COL1A1 cDC',
#     '2':'IGHA1 cDC',
#     '3':'COL1A1 cDC',
#     '4':'COL1A1 cDC',
#     '5':'IGHA1 cDC'
# }
# adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
# adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata_raw.write(os.path.join(out_dir, 'Macrophage.h5ad'))

sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对Tumor注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/05.Mission/Tumor/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/Tumor/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))


# 对CD8T注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD8T/CD8T.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/04.Mission/CD8T/all_cells_features.xls", index_col=0, sep='\t')  # 调整路径和索引列
umap_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/04.Mission/CD8T/umapreductions.csv", index_col=0)  # 调整路径
adata_raw.obs['cell_type'] = cell_type_df['cell_type']
umap_coords = umap_df.loc[adata_raw.obs.index, ['UMAP_1', 'UMAP_2']].values  # 根据实际列名调整
adata_raw.obsm['X_umap'] = umap_coords

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['LARC', 'N']
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
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata_raw.write(os.path.join(out_dir, 'CD8T.h5ad'))

sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

method = "cell_type"
method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("celltype_proportions_by_group.csv")

# 对CD4T注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD4T/CD4T.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD4T/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/all_cells_features.xls", index_col=0, sep='\t')  # 调整路径和索引列
umap_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/umapreductions.csv", index_col=0)  # 调整路径
adata_raw.obs['cell_type'] = cell_type_df['cell_type']
umap_coords = umap_df.loc[adata_raw.obs.index, ['UMAP_1', 'UMAP_2']].values  # 根据实际列名调整
adata_raw.obsm['X_umap'] = umap_coords

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['LARC', 'N']
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
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata_raw.write(os.path.join(out_dir, 'CD4T.h5ad'))

sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


#### 画marker plot ####
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/markerplot/"
sc.pl.umap(adata_raw, color=['TCF7','PDCD1'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata_raw, var_names=['TCF7','PDCD1'], groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata_raw, var_names=['TCF7','PDCD1'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'), bbox_inches="tight")


# 映射所有细分 #
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/mapping_raw/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/Monocytes/Macrophage.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD4T/CD4T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/CD8T.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')