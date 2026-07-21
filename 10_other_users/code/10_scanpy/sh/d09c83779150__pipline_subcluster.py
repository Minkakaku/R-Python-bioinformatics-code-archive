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

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2505152218-01/01.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2505152218-01/11.Mission/"

sc.pl.umap(adata_raw, color=['Gls'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_Gls.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Gls.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.umap(adata_raw, color=['Gls'], cmap='coolwarm',show=False,vmin=1)
plt.savefig(os.path.join(out_dir, 'umap_Gls_fiter.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Gls_fiter.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata_raw, var_names=['Gls'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-cell_type_group.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-cell_type_group.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata_raw, var_names=['Gls'], groupby='sub_cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-sub_cell_type.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-sub_cell_type.pdf'), facecolor='white', bbox_inches="tight")

sc.pl.matrixplot(adata_raw, var_names=['Gls'], groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_cell_type.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata_raw, var_names=['Gls'], groupby='group', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_group.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata_raw, var_names=['Gls'], groupby=['cell_type','group'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_cell_type_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_cell_type_group.pdf'),bbox_inches='tight')
plt.close()


gene_use = ['Gls']
df = sc.get.obs_df(adata_raw, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata_raw.obs
ax1 = sc.pl.dotplot(
 adata_raw,
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
    adata_raw,
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
save_fig(plt, out_dir = out_dir, file_name=f'{gene_use[0]}_dotplot')


sc.pl.matrixplot(new_adata, var_names=new_adata.var_names, groupby=['cell_type'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, f'{gene_use[0]}_heatmap.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'{gene_use[0]}_heatmap.png'),dpi=300, bbox_inches='tight')
plt.close()


ax1 = sc.pl.dotplot(
 adata_raw,
    var_names=gene_use,  # 显示所有分组
    groupby=['sub_cell_type','group'],  # 按细胞类型分组
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
    adata_raw,
    var_names=gene_use,  # 显示所有分组
    groupby=['sub_cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)
b = a.dot_size_df
b['sub_cell_type_group'] = b.index
b[['sub_cell_type', 'group']] = b['sub_cell_type_group'].str.split('_', n=1,expand=True)
b = b.drop(columns='sub_cell_type_group')
b.columns = ['PCT', 'sub_cell_type', 'group']
b = b.pivot(index='sub_cell_type', columns='group', values='PCT')
b.columns.name  = ""
b = b[meta['group'].cat.categories.to_list()]
b = b.loc[meta['sub_cell_type'].cat.categories.to_list()]
sc.pl.dotplot(
    new_adata,
    var_names=new_adata.var_names,  # 显示所有分组
    groupby='sub_cell_type',  # 按细胞类型分组
    dot_size_df = b,
    title=f'{gene_use[0]} expression',
    show=False
)
save_fig(plt, out_dir = out_dir, file_name=f'{gene_use[0]}_dotplot')

sc.pl.matrixplot(new_adata, var_names=new_adata.var_names, groupby=['sub_cell_type'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, f'{gene_use[0]}_heatmap.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'{gene_use[0]}_heatmap.png'),dpi=300, bbox_inches='tight')
plt.close()

import scanpy as sc
import pandas as pd

# 提取数据
genes = new_adata.var_names.tolist()

df_expr = pd.DataFrame(
    new_adata.X.toarray() if hasattr(new_adata.X, 'toarray') else new_adata.X,
    index=new_adata.obs_names,
    columns=genes
)
df_expr['sub_cell_type'] = new_adata.obs['sub_cell_type'].values

# 转换为长格式
df_long = df_expr.melt(
    id_vars=['sub_cell_type'],
    var_name='group',
    value_name='gene expression'
)

# 导出
df_long.to_csv(os.path.join(out_dir, 'expression_long_format.csv'), index=False)
grouped_mean = df_long.groupby(['group', 'sub_cell_type'])['gene expression'].mean()
grouped_mean.to_csv(os.path.join(out_dir, 'group_sub_cell_type_mean_expression.csv'), index=True)