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

# 对Monocytes注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2408071595-01/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Erythrocytes', 'Macrophages', 'Megakaryocytes', 'Monocytes', 'NK cells', 'Neutrophils', 'Plasma cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2408071595-01/02.Mission/Monocytes/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Macrophages",'Monocytes']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1,n_pcs = 15,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[~adata.obs['leiden'].isin(["7"]), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=1,n_pcs = 20,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

monocyte_macrophage_marker = {
    "Classical Monocyte": ["CD14", "FCGR3A", "CCR2", "CX3CR1"],
    "Non-classical Monocyte": ["CD14", "FCGR3A", "CX3CR1", "CCR2"],
    "M1 Macrophage": ["CD80", "CD86", "HLA-DRA", "TNF", "IL1B", "CD38"],
    "M2 Macrophage": ["CD163", "IL10", "ARG1", "MRC1"],
    "cDC": ["IRF8", "HLA-DRA", "CD1C", "CLEC9A", "XCR1"],
    "pDC": ["SIGLEC1", "IL3RA", "GZMB", "IRF7", "TCF4"]
}
marker=["CD14", "FCGR3A", "CCR2", "CX3CR1","CD80", "CD86", "HLA-DRA", "TNF", "IL1B", "CD38","CD163", "IL10", "ARG1", "MRC1","IRF8", "CD1C", "CLEC9A", "XCR1","SIGLEC1", "IL3RA", "GZMB", "IRF7", "TCF4","LYZ"]
sc.pl.umap(adata, color=marker, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=monocyte_macrophage_marker, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'S100A12 Monocyte',
    '1':'GBP1 Monocyte',
    '2':'S100A12 Monocyte',
    '3':'ZFP36L1 Monocyte',
    '4':'RGS2 Monocyte',
    '5':'FCGR3A Nonclassical Monocyte',
    '6':'GBP1 Monocyte',
    '7':'C1QA Nonclassical Monocyte',
    '8':'IL32 Monocyte',
    '9':'S100A12 Monocyte',
    '10':'ZFP36L1 Monocyte',
    '11':'S100A12 Monocyte',
    '12':'DC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['NC', 'HLH', 'iMCD-TAFRO']
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
adata.write(os.path.join(out_dir, 'Monocyte.h5ad'))

# sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
# save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
# sc.pl.umap(adata, color='sample', title='sample', show=False)
# save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

out_dir="/mnt/NexData_file/home/sh/Projects/AE2408071595-01/02.Mission/Monocytes_cytokine_plot/"
cell_type_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",index_col=0, header=None)  # 调整路径和索引列
cell_type_df = cell_type_df.index
common_cells = list(set(cell_type_df) & set(adata.raw.var_names))

gene_means = np.array(adata.raw[:, common_cells].X.mean(axis=0)).flatten()
gene_order_df = pd.DataFrame({
    'gene': common_cells,
    'mean_expression': gene_means
})
genes_sorted = gene_order_df.sort_values('mean_expression', ascending=False)['gene'].tolist()

sc.pl.dotplot(adata, var_names=genes_sorted, groupby='group', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_group_scale.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_group_scale.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata, var_names=genes_sorted, groupby='sample', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_sample_scale.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_sample_scale.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata, var_names=genes_sorted, groupby='group', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_group.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_group.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata, var_names=genes_sorted, groupby='sample', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_sample.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_sample.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata, var_names=genes_sorted, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_cell_type.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_cell_type.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata, var_names=genes_sorted, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_cell_type_group.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_cytokine_genes_cell_type_group.pdf'), bbox_inches="tight")

adata_all = adata.raw.to_adata()
group_mean_df = (
    adata_all.to_df()
    .groupby(adata_all.obs['group'])
    .mean()
    .T
    .reindex(genes_sorted)  # 重新排序基因
)
group_mean_df.to_csv(os.path.join(out_dir, 'group_mean_expression.csv'))
sample_mean_df = (
    adata_all.to_df()
    .groupby(adata_all.obs['sample'])
    .mean()
    .T
    .reindex(genes_sorted)  # 重新排序基因
)
sample_mean_df.to_csv(os.path.join(out_dir, 'sample_mean_expression.csv'))
cell_type_mean_df = (
    adata_all.to_df()
    .groupby([adata_all.obs['cell_type'], adata_all.obs['group']])
    .mean()
    .T
    .reindex(genes_sorted)  # 重新排序基因
)
cell_type_mean_df.to_csv(os.path.join(out_dir, 'cell_type_mean_expression.csv'))

sc.tl.score_genes(adata, gene_list=common_cells, score_name='cytokine_score')
sc.pl.umap(adata, color='cytokine_score', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_cytokine_score', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='cytokine_score', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_cytokine_score', fig_size={'w': 10, 'h': 6})
sc.pl.dotplot(adata, var_names=['cytokine_score'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_cytokine_score', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['cytokine_score'], groupby=['cell_type','group'],  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_cytokine_score_group', fig_size={'w': 6, 'h': 12})


# 对CD8T注释细分 ####后面没做
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