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

# 对CD8T注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/CD8T.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/14.Mission/CD8T/"
adata_raw.obs['sample'].cat.categories
adata_raw = adata_raw[adata_raw.obs['group'].isin(['LARC'])]

anno_dict={
    'BCR0807C': 'LOW',
    'DYD0703C': 'HIGH',
    'GC0812C': 'HIGH',
    'GYH0708C': 'HIGH',
    'HJM0808C': 'LOW',
    'JL0728C': 'HIGH',
    'LYF0725C': 'HIGH',
    'NX0708C': 'LOW',
    'QHX0729C': 'HIGH',
    'TXL0725C': 'LOW',
    'WHL0717C': 'HIGH',
    'WHT0811C': 'LOW',
    'WZQ0811C': 'LOW',
    'YBJ0710C': 'LOW',
    'ZBZ0808C': 'LOW',
    'ZYH0808C': 'LOW'
}
adata_raw.obs['group'] = [anno_dict[clust] for clust in adata_raw.obs['sample']]
adata_raw.obs['group'] = adata_raw.obs['group'].astype('category')

adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
adata_raw.obs['sample'] = (
    pd.Categorical(adata_raw.obs['sample'],
                   categories=['DYD0703C', 'GYH0708C', 'LYF0725C', 'JL0728C', 'QHX0729C', 'GC0812C', 'WHL0717C', 'YBJ0710C', 'TXL0725C', 'BCR0807C', 'HJM0808C', 'WZQ0811C', 'ZYH0808C', 'ZBZ0808C', 'NX0708C', 'WHT0811C'],
                   ordered=True)
)
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata_raw.write(os.path.join(out_dir, 'CD8T.h5ad'))

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv("celltype_proportions_by_group.csv")