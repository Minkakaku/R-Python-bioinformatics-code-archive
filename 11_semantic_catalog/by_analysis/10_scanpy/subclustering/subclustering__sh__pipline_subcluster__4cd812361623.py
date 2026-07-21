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

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/08.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/17.Mission/Main/"
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*6, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="Arg2", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_Arg2.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_Arg2.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.obs['group'] = (
    pd.Categorical(adata_raw.obs['group'],
                   categories=groups,
                   ordered=True)
)
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots(figsize=(6, 4))
    sc.pl.violin(adata_raw, keys=['Arg2'], groupby='group', show=False, ax=ax, 
             stripplot=False, density_norm='count', xlabel=None, 
             ylabel='Expression level', jitter=False)
# 计算每个组的均值并添加点
    mean_values = []
    for i, group in enumerate(groups):
        expr = adata_raw.raw[adata_raw.obs['group'] == group, 'Arg2'].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = expr.flatten()
        mean_val = np.mean(expr)
        mean_values.append(mean_val)
        ax.scatter(i, mean_val, color='black', s=50, zorder=10)
    ax.plot(range(len(groups)), mean_values, color='black', linestyle='-', linewidth=2, zorder=5)
    ax.grid(False)
    plt.savefig(os.path.join(out_dir, 'violin_Arg2.png'))
    plt.savefig(os.path.join(out_dir, 'violin_Arg2.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata_raw, var_names=['Arg2'], groupby='group', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_Arg2_group.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_Arg2_group.pdf'))

gene_use = ['Arg2']
df = sc.get.obs_df(adata_raw, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata_raw.obs

expr_matrix = pd.DataFrame(
    np.NaN,  # 初始化为np.nan
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
# b = b.fillna(0)
sc.pl.dotplot(
    new_adata,
    var_names=new_adata.var_names,  # 显示所有分组
    groupby='cell_type',  # 按细胞类型分组
    dot_size_df = b,
    title=f'{gene_use[0]} expression',
    show=False
)
save_fig(plt, out_dir = out_dir, file_name=f'{gene_use[0]}_dotplot')



adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Macrophages/Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/17.Mission/Macrophages/"
marker_genes=['Arg2','Ccl2','Vim','Hk1','Hk2','Pfkp','Pkm','Il1b']
sc.pl.stacked_violin(adata_raw, var_names=marker_genes, groupby='group', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'marker_genes_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'marker_genes_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata_raw, var_names=marker_genes, groupby='group', dendrogram=False, swap_axes=True, show=False)
plt.savefig(os.path.join(out_dir, 'marker_genes_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'marker_genes_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.close()
sc.pl.dotplot(adata_raw, var_names=marker_genes, groupby='group', dendrogram=False, swap_axes=True, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'marker_genes_bubble_scale.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'marker_genes_bubble_scale.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.close()

sc.pl.umap(adata_raw, color=marker_genes, show=False, cmap='coolwarm')
save_fig(plt,file_name= 'marker_genes_umap',out_dir=out_dir,fig_size={'w':20,'h':10})
sc.pl.umap(adata_raw, color=marker_genes, show=False, cmap='coolwarm',vmax='p99')
save_fig(plt,file_name= 'marker_genes_umap_scale',out_dir=out_dir,fig_size={'w':20,'h':10})

for idx, gene in enumerate(marker_genes):
    print(f"Processing {idx+1}/{len(marker_genes)}: {gene}")
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*6, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color=gene, legend_loc=legend, title=groups[i], ax=axs[i], show=False)
        plt.savefig(os.path.join(out_dir, f'dimplot_umap_group_{gene}.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, f'dimplot_umap_group_{gene}.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.close()

adata_raw.obs['group'] = (
    pd.Categorical(adata_raw.obs['group'],
                   categories=groups,
                   ordered=True)
)
for idx, gene in enumerate(marker_genes):
    print(f"Processing {idx+1}/{len(marker_genes)}: {gene}")
    with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
        fig, ax = plt.subplots(figsize=(6, 4))
        sc.pl.violin(adata_raw, keys=gene, groupby='group', show=False, ax=ax, 
             stripplot=False, density_norm='count', xlabel=None, 
             ylabel='Expression level', jitter=False)
# 计算每个组的均值并添加点
        mean_values = []
        for i, group in enumerate(groups):
            expr = adata_raw.raw[adata_raw.obs['group'] == group, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()
            mean_val = np.mean(expr)
            mean_values.append(mean_val)
            ax.scatter(i, mean_val, color='black', s=50, zorder=10)
        ax.plot(range(len(groups)), mean_values, color='black', linestyle='-', linewidth=2, zorder=5)
        ax.grid(False)
        plt.savefig(os.path.join(out_dir, f'violin_group_{gene}.pdf'))
        plt.savefig(os.path.join(out_dir, f'violin_group_{gene}.png'))

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 定义要循环的基因列表
marker_genes = ['Arg2', 'Ccl2', 'Vim', 'Hk1', 'Hk2', 'Pfkp', 'Pkm', 'Il1b']
meta = adata_raw.obs

# 预先计算分组信息
group_categories = meta['group'].cat.categories.to_list()
cell_type_categories = meta['cell_type'].cat.categories.to_list()

for gene_use in marker_genes:
    print(f"Processing gene: {gene_use}")
    
    # 1. 提取基因表达
    df = sc.get.obs_df(adata_raw, keys=[gene_use], use_raw=True)
    
    # 2. 构建表达矩阵（分组 × 细胞）
    expr_matrix = pd.DataFrame(
        np.NaN,
        index=df.index,
        columns=group_categories
    )
    
    for group in expr_matrix.columns:
        group_cells = meta[meta['group'] == group].index
        expr_matrix.loc[group_cells, group] = df.loc[group_cells, gene_use]
    
    # 3. 创建新的 AnnData 对象
    new_adata = sc.AnnData(
        X=expr_matrix.values,
        obs=meta.drop(columns=['group']),
        var=pd.DataFrame(index=expr_matrix.columns)
    )
    
    # 4. 获取 dotplot 数据
    a = sc.pl.dotplot(
        adata_raw,
        var_names=[gene_use],
        groupby=['cell_type', 'group'],
        show=True,
        return_fig=True
    )
    
    b = a.dot_size_df
    b['cell_type_group'] = b.index
    b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1, expand=True)
    b = b.drop(columns='cell_type_group')
    b.columns = ['PCT', 'cell_type', 'group']
    b = b.pivot(index='cell_type', columns='group', values='PCT')
    b.columns.name = ""
    b = b[group_categories]
    b = b.loc[cell_type_categories]
    
    # 5. 绘制 dotplot
    sc.pl.dotplot(
        new_adata,
        var_names=new_adata.var_names,
        groupby='cell_type',
        dot_size_df=b,
        title=f'{gene_use} expression',
        show=False
    )
    
    # 6. 保存图片
    save_fig(plt, out_dir=out_dir, file_name=f'{gene_use}_dotplot')
    plt.close()  # 关闭图形释放内存

print("All dotplots completed!")