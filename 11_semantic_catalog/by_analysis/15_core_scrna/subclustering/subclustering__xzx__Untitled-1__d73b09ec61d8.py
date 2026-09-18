# %%

import diopy
# 保存图片
import os
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import font_manager as fm
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

# 加载亚群细分的函数
import sys
sys.path.append('/mnt/NexData_file/home/lst/scrna-code/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion
# os.chdir('/lvdata/wzb/scRNA/FW2023-656/2024.01.02/')
import pandas as pd

# %%
sys.path.append('/mnt/NexData_file/home/lst/scrna-code/diff')
from find_marker_1 import find_marker_gene_1

# %%
adata=sc.read_h5ad('/mnt/NexData_file/home/lst/project/FW2023-529/adata.h5ad')

# %%
outdir='/mnt/NexData_file/home/lst/project/FW2023-529/售后8'

# %%
sc.pl.umap(adata,color = ['Saa3'],size=40)

# %%
sc.pl.umap(adata,color = ['cell_type'],size=40)

# %%
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import cosg as cosg
import numpy as np
#import gseapy as gp
import scanpy as sc
import os
import re

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def find_marker_gene_1(adata, group,species=None,out_dir=os.getcwd(),filter = True, logfoldchanges=0.25,pvals_adj=0.05,pct_nz_group=0.1, method='wilcoxon', reference = 'rest',gene_annot=None):
    """
    查找标记基因。

    Args:
        adata: 输入的 AnnData 对象。
        group: 细胞类型群组（例如 leiden）。
        species (str, optional): 物种，如 'human' 或 'mouse'。默认为 None。
        out_dir (str, optional): 输出目录。默认为当前工作目录。
        filter (bool, optional): 是否过滤。默认为 True。
        logfoldchanges (float, optional): 对数折叠变化阈值。默认为 0.25。
        pvals_adj (float, optional): 调整后的 p 值阈值。默认为 0.05。
        pct_nz_group (float, optional): 组中非零百分比阈值。默认为 0.1。
        method (str, optional): 使用的方法，如 'wilcoxon'。默认为 'wilcoxon'。
        reference (str, optional): 参考组。默认为 'rest'。
        gene_annot (optional): 基因注释。默认为 None。

    Returns:
        写出group类的所有marker gene。
    """
    sample_counts = adata.obs[group].value_counts()
    problematic_groups = sample_counts[sample_counts == 1].index.to_list()
    adata = adata[~adata.obs[group].isin(problematic_groups)]
    adata.uns['log1p']['base'] = 2
    sc.tl.rank_genes_groups(adata, group, method=method,reference= reference,pts=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    all_diff = []
    for c in groups:
        
        # c = 'CD4+TN'
        d = re.sub(r'[+/\\-]',"_",c)
        out = os.path.join(out_dir, f'{d}_up')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        tmp[(tmp['pct_nz_group'] >=0.10)].to_csv(os.path.join(out, 'All_expressed_genes.xls'), sep='\t', index=None)
        if filter :
            tmp = tmp[(tmp['logfoldchanges']>logfoldchanges) & (tmp['pvals_adj']<pvals_adj) & (tmp['pct_nz_group'] >=pct_nz_group)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)
    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    all_diff.to_csv(os.path.join(out_dir, 'differentially_expressed_genes.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'diff_statistics.xls'), sep='\t')
    sc.pl.rank_genes_groups_heatmap(
    adata,
    n_genes=10,
    values_to_plot="logfoldchanges",
    min_logfoldchange=1,
    swap_axes=True,
    cmap='bwr',
    standard_scale='var',
    show_gene_labels=True,
    dendrogram=True,
    show=False)

    save_fig(plt, out_dir, file_name='heatmap_top10_marker_gene', fig_size={'width':20, 'height':20})

    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=10,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=1,
        cmap='bwr',
        colorbar_title='log fold change',
        show=False,
    )
    save_fig(plt, out_dir, file_name='dotplot_top10_marker_gene', fig_size={'width':20, 'height':10})

# %%
find_marker_gene_1(adata[adata.obs['cell_type'].isin(['Monocytes','Macrophages'])],group = 'group',species='mouse',out_dir = outdir)

# %%
outdir

# %%
adata.obs['cell_type']=adata.obs['cell_type'].replace('Monocytes','Mono_Macro')
adata.obs['cell_type']=adata.obs['cell_type'].replace('Macrophages','Mono_Macro')


# %%


# %%
set(adata.obs['cell_type'])

# %%
import scanpy as sc
import pandas as pd
import numpy as np

gene_use = ['Mmp13']
outdir_tmp = '/mnt/NexData_file/home/lst/project/FW2023-529/售后8/genes_1'

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
#b = list(set(b[meta['group']]))
b = b.loc[meta['cell_type'].cat.categories.to_list()]
#b = list(set(b.loc[meta['cell_type']]))
sc.pl.dotplot(
    new_adata,
    var_names=new_adata.var_names,  # 显示所有分组
    groupby='cell_type',  # 按细胞类型分组
    dot_size_df = b,
    title=f'{gene_use[0]} expression',
    show=False
)
save_fig(plt, out_dir = outdir_tmp, file_name=f'{gene_use[0]}_dotplot')

# %%
sc.pl.violin(adata, groupby='cell_type', keys='Saa3', rotation=90)
sc.pl.violin(adata, groupby='cell_type', keys='Il1b', rotation=90)
sc.pl.violin(adata, groupby='cell_type', keys='Arg2', rotation=90)


# %%
sc.pl.umap(adata,color = ['Il1b','Arg2'])

# %%
sc.pl.umap(adata,color = ['cell_type','group'])

# %%
find_marker_gene(adata,group = 'group',species='mouse',out_dir = outdir)

# %%
sample_counts = adata.obs['group'].value_counts()
problematic_groups = sample_counts[sample_counts == 1].index.to_list()

# %%
result = adata.uns['rank_genes_groups']

# %%
groups = result['names'].dtype.names
print(groups)


