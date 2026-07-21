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

adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20251028/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B/Plasma', 'doublets', 'Mast cells', 'Fibroblasts', 'T cells', 'AT2', 'Endothelial cells', 'Macrophages', 'AT1', 'pDCs', 'Ciliated cells', 'Neutrophils'
output_path="/mnt/NexData_file/home/sh/Projects/FWSC20251028/01.Mission/Markerplot/"

adata=adata_raw.raw.to_adata()
marker_genes=['HGF','MET']
sc.pl.umap(adata, color=marker_genes, show=False, cmap='viridis_r')
save_fig(plt,file_name= 'Marker_umap',out_dir=output_path,fig_size={'w':12,'h':6})

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="HGF", legend_loc=legend, title=groups[i], ax=axs[i], show=False, cmap='viridis_r')
    plt.savefig(os.path.join(output_path, 'HGF_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(output_path, 'HGF_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="MET", legend_loc=legend, title=groups[i], ax=axs[i], show=False, cmap='viridis_r')
    plt.savefig(os.path.join(output_path, 'MET_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(output_path, 'MET_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'Marker_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Marker_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names='HGF', groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'HGF_dotplot_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'HGF_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.matrixplot(adata, marker_genes, groupby='cell_type', cmap='coolwarm', dendrogram=False, show=False, swap_axes=False)
plt.savefig(os.path.join(output_path, 'Marker_matrixplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Marker_matrixplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=False, dendrogram=False)
plt.savefig(os.path.join(output_path, 'Marker_stacked_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Marker_stacked_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs['MET_gene'] = adata[:, 'MET'].X.toarray().ravel()
adata.obs['HGF_gene'] = adata[:, 'HGF'].X.toarray().ravel()
adata.obs.to_csv(os.path.join(output_path, "all_cells_features.xls"), sep='\t')

for gene in marker_genes:
    adata.obs[f'{gene}_CON'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'CON'
    adata.obs.loc[mask, f'{gene}_CON'] = None

for gene in marker_genes:
    adata.obs[f'{gene}_CONS'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'CONS'
    adata.obs.loc[mask, f'{gene}_CONS'] = None

for gene in marker_genes:
    adata.obs[f'{gene}_EMP'] = adata[:, gene].X.toarray().ravel()
    mask = adata.obs['group'] != 'EMP'
    adata.obs.loc[mask, f'{gene}_EMP'] = None

genelist=['HGF_CON', 'HGF_CONS', 'HGF_EMP']
sc.pl.dotplot(adata, var_names=genelist, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'HGF_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'HGF_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

genelist=['MET_CON', 'MET_CONS', 'MET_EMP']
sc.pl.dotplot(adata, var_names=genelist, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'MET_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'MET_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

gene = 'MET'                       # 换成你的基因名
dot_df = sc.get.obs_df(
            adata,
            keys=[gene,'cell_type', 'group'],
            use_raw=False                     # 如需 raw 改为 True
         )
# 1) 平均表达量
mean_df = (dot_df
           .groupby(['cell_type', 'group'])[gene]
           .mean()
           .rename('mean'))

# 2) 检测比例（>0 视为检测到）
pct_df = (dot_df
          .assign(detected=dot_df[gene] > 0)
          .groupby(['cell_type', 'group'])['detected']
          .mean()
          .mul(100)
          .rename('pct'))

# 3) 合并成最终表格（行索引 cell_type，列索引 group）
plot_df = pd.concat([mean_df, pct_df], axis=1).reset_index()

# 4. 画气泡图
plt.figure(figsize=(4,5))
sns.scatterplot(
    data=plot_df,
    x='group', y='cell_type',
    size='pct', hue='mean',
    sizes=(0, 300), palette='Reds', edgecolor='black', linewidth=0.5
)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=None, frameon=False)
plt.title(''); plt.xlabel(''); plt.ylabel('')
plt.tight_layout()
plt.margins(x=0.25)
# sc.pl.dotplot(adata, var_names=genelist, groupby='cell_type', dendrogram=False, mean_only_expressed=True, show=False)
plt.savefig(os.path.join(output_path, 'MET_dotplot_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'MET_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata, var_names='MET', groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'MET_dotplot_group_celltype_legend.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'MET_dotplot_group_celltype_legend.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


import scanpy as sc
import pandas as pd
import numpy as np

gene_use = ['Irf4']
outdir_tmp = '/lvdata/wzb/scRNA/FWSC20252335/2025.06.30/All/IRF'

df = sc.get.obs_df(adata, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata.obs
ax1 = sc.pl.dotplot(
 adata,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)

expr_matrix = pd.DataFrame(
    np.nan,  # 初始化为NaN
    index=df.index,
    columns=meta['group'].cat.categories.to_list()
)

for group in expr_matrix.columns:
    group_cells = meta[meta['group'] == group].index
    expr_matrix.loc[group_cells, group] = df.loc[group_cells, gene_use[0]]

new_adata = sc.AnnData(
    X=expr_matrix.values,  # 使用numpy数组
    obs=meta.drop(columns=['group']),  # 保留其他metadata
    var=pd.DataFrame(index=expr_matrix.columns)  # var是分组名
)

a = sc.pl.dotplot(
    adata,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)
b = a.dot_size_df
b['cell_type_group'] = b.index
b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1,expand=True)
b = b.drop(columns='cell_type_group')
b.columns = ['PCT', 'cell_type', 'group']
b = b.pivot(index='cell_type', columns='group', values='PCT')
b.columns.name  = ""
b = b[meta['group'].cat.categories.to_list()]
b = b.loc[meta['cell_type'].cat.categories.to_list()]
sc.pl.dotplot(
    new_adata,
    var_names=new_adata.var_names,  # 显示所有分组
    groupby='cell_type',  # 按细胞类型分组
    dot_size_df = b,
    title=f'{gene_use[0]} expression',
    show=False
)
save_fig(plt, out_dir = outdir_tmp, file_name=f'{gene_use[0]}_dotplot')