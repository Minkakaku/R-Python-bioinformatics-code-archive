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
    plt.savefig(os.path.join(out_dir, f'{file_name}.tiff'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None,prefix='F'):
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
    file_name=f"{prefix}_celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result1/"
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig1_A_umap/Main_mapping.h5ad")
name_map = {
    'B cells': 'B cell',
    'Megakaryocytes': 'Megakaryocyte',
    'Myeloid cells': 'Myeloid cell',
    'NK cells': 'NK cell',
    'NKT cells': 'NKT cell',
    'Plasma cells': 'Plasma cell',
    'pDCs': 'pDC',
    'T cells': 'T cell'
}
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].map(name_map).astype('category')
# 1. 定义目标顺序
desired_order = ['B cell', 'Megakaryocyte', 'Myeloid cell', 'NK cell', 'NKT cell','Plasma cell', 'pDC', 'T cell']
# 2. 转换为 categorical 并指定顺序
adata_raw.obs['cell_type'] = pd.Categorical( adata_raw.obs['cell_type'], categories=desired_order, ordered=True)
# 3. 验证
print(adata_raw.obs['cell_type'].cat.categories) 
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',prefix='Main')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',prefix='Main')


adata=adata_raw[adata_raw.obs.cell_type.isin(['Myeloid cell','pDC']), :]
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result2/"
name_map = {
    'CCL2 Classical monocyte': 'CCL2 Classical monocyte',
    'DUSP1 Classical monocyte': 'DUSP1 Classical monocyte',
    'Non-classical monocyte': 'Non-classical monocyte',
    'VCAN Classical monocyte': 'VCAN Classical monocyte',
    'cDCs': 'cDC',
    'pDCs': 'pDC'
}
adata.obs['cell_type'] = adata.obs['sub_cell_type'].map(name_map).astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,wspace=0.1, legend_fontsize=12)
    plt.savefig(os.path.join(out_dir, 'Myeloid_dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'Myeloid_dimplot_umap_celltype.pdf'))
    plt.savefig(os.path.join(out_dir, 'Myeloid_dimplot_umap_celltype.tiff'))

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
sc.tl.dendrogram(adata,groupby="cell_type")
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()
# enrich
find_group_marker(adata=adata,species='human',output=out_dir)