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

adata=sc.read("/mnt/NexData_file/scRNA/scVision_output/FWAE20251216/result.h5ad")
adata.obs['sample'].cat.categories

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/geneplot/"
sc.pl.dotplot(adata, var_names=['IL6','MIF'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


adata_sub = adata[adata.obs['sample'].isin(['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5']), :]
adata_sub.obs['sample'].cat.categories
adata_sub.obs['sample'] = pd.Categorical(
    adata_sub.obs['sample'],                       # 原值保留
    categories=['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5'],                     # 新顺序
    ordered=False                                 # 无需数学大小比较就 False
)
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/MIF_plot/"
sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'sample_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sample_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.stacked_violin(adata_sub, var_names=['MIF'], groupby=['sample'], show=False, swap_axes=False, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'sample_violin.pdf'), facecolor='white',bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'sample_violin.png'), facecolor='white',bbox_inches='tight', dpi=300)

sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['cell_type','sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

plot_celltype_proportion(adata_sub, out_dir=out_dir,type='sample')
ax = pd.crosstab(adata_sub.obs['sample'], adata.obs['cell_type'], normalize='index')
ax.to_csv(os.path.join(out_dir,'celltype_proportion.xls'),sep="\t",index=True,index_label='cell_type')

group_dict={
    'iMCD_3':'severe',
    'iMCD_5':'non-severe',
    'iMCD_7':'severe',
    'iMCD_17':'non-severe',
    'iMCD_20':'severe',
    'iMCD_21':'severe'
}
adata_sub.obs['group'] = [group_dict[samples] for samples in adata_sub.obs['sample']]
adata_sub.obs['group'] = adata_sub.obs['group'].astype("category")

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/severe_plot/"
sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.stacked_violin(adata_sub, var_names=['MIF'], groupby=['group'], show=False, swap_axes=False, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'group_violin.pdf'), facecolor='white',bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'group_violin.png'), facecolor='white',bbox_inches='tight', dpi=300)

sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


outdir_tmp="/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/geneplot/MIF_sample/"
adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/mapping_raw/Main_mapping.h5ad")
adata_sub = adata[adata.obs['sample'].isin(['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5']), :]
adata_sub.obs['sample'].cat.categories
adata_sub.obs['sample'] = pd.Categorical(
    adata_sub.obs['sample'],                       # 原值保留
    categories=['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5'],                     # 新顺序
    ordered=False                                 # 无需数学大小比较就 False
)
sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['cell_type','sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# adata_sub.obs['raw_cell_type']=adata_sub.obs['cell_type']
# adata_sub.obs['cell_type']=adata_sub.obs['group']
# adata_sub.obs['group']=adata_sub.obs['raw_cell_type']
adata_sub.obs['group']=adata_sub.obs['sample']
gene_use = ['MIF']
df = sc.get.obs_df(adata_sub, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata_sub.obs
ax1 = sc.pl.dotplot(
 adata_sub,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)

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
save_fig(plt, out_dir = outdir_tmp, file_name=f'{gene_use[0]}_dotplot')