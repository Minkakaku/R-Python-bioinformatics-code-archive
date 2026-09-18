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

adata=sc.read("/mnt/NexData_file/scRNA/scVision_output/FWAE20251216/result.h5ad")
adata.obs['sample'].cat.categories

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/geneplot/"
sc.pl.dotplot(adata, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'MIF-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'MIF-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


adata_sub = adata[adata.obs['sample'].isin(['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5']), :]
adata_sub.obs['sample'].cat.categories
adata_sub.obs['sample'] = pd.Categorical(
    adata_sub.obs['sample'],                       # 原值保留
    categories=['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5'],                     # 新顺序
    ordered=False                                 # 无需数学大小比较就 False
)

sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'sample_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sample_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
# sc.pl.violin(adata_sub, ['MIF'], groupby='sample',jitter=0, rotation=45, multi_panel=True, show=False)
# plt.savefig(os.path.join(out_dir, 'sample_violin.pdf'), facecolor='white',bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir, 'sample_violin.png'), facecolor='white',bbox_inches='tight', dpi=300)

dp=sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['sample'], dendrogram=False, show=False,return_fig=True)
mean_df = dp.dot_color_df
mean_df.to_csv("/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/geneplot/sample_mean_exp.csv", index=True)
frac_df = dp.dot_size_df
frac_df.to_csv("/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/geneplot/sample__mean_pct.csv", index=True)

adata_sub.obs['group']=adata_sub.obs['sample']
gene_use = ['MIF']
df = sc.get.obs_df(adata_sub, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata_sub.obs

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
    adata_sub,
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


adata_after=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/MIF_after/"
sc.pl.dotplot(adata_after, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata_after, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type_group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_after=adata_after[adata_after.obs.cell_type=="Plasma cells", :]
sc.pl.dotplot(adata_after, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Plasma-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Plasma-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


adata_after.obs['sample'] = pd.Categorical(
    adata_after.obs['sample'],                       # 原值保留
    categories=['NC_1', 'NC_2', 'NC_3','iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5','iMCD_20_after','iMCD_21_after','iMCD_5_after'], # 新顺序
    ordered=False # 无需数学大小比较就 False
)
sc.pl.dotplot(adata_after, var_names=['MIF'], groupby=['sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Plasma-sample_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'Plasma-sample_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)