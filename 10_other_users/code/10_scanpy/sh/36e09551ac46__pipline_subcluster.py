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
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), facecolor='white', bbox_inches="tight", dpi=300, pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), facecolor='white', bbox_inches="tight", dpi=300, pad_inches = 0.1)
    plt.clf()

def read_gmt_to_dict(gmt_file):
    """
    读取 GMT 文件，返回字典 {基因集名称: [基因列表]}
    """
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                set_name = parts[0]
                description = parts[1]
                genes = parts[2:]
                gene_sets[set_name] = genes
    return gene_sets

#### 
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2604163024/result.h5ad")
# adata_raw.obs.cell_type.cat.categories
# 'Oligodendrocyte', 'Microglia', 'Endothelial cells', 'T cells', 'Macrophages', 'Neurons', 'Astrocytes', 'Neutrophils', 'OPC', 'Erythrocyte', 'Fibroblasts', 'B cells'
# adata_raw.obs['sample'].cat.categories
# 'PD', 'PDNTMSC2', 'control'

out_dir="/mnt/NexData_file/home/sh/Projects/SC2604163024-03/02.Mission/S100a9_plot/"

sc.pl.umap(adata_raw, color=['S100a9'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata_raw, var_names=['S100a9'], groupby='group', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_gene_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_gene_group.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata_raw, var_names=['S100a9'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.violin(adata_raw, keys=['S100a9'], groupby='group', show=False, stripplot=True, xlabel=None, ylabel='Expression level', jitter=False)
plt.savefig(os.path.join(out_dir, 'violin_gene_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'violin_gene_group.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
with rc_context({'figure.figsize': (12, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.violin(adata_raw, keys=['S100a9'], groupby='cell_type', show=False, stripplot=True, xlabel=None, ylabel='Expression level', jitter=False)
    plt.savefig(os.path.join(out_dir, 'violin_gene_cell_type.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'violin_gene_cell_type.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

gene_use = ['S100a9']
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