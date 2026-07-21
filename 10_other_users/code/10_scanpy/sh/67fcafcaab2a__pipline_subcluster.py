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
from scipy.stats import spearmanr
from itertools import combinations
from matplotlib.patches import Circle
import matplotlib.colors as mcolors

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

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/06.Mission/Macrophages/Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/07.Mission/Macrophages_geneplot/"
sc.pl.umap(adata, color=['EGFR','TLR4'], size=25, show=False,cmap='RdBu_r',vmax='p99')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['EGFR','TLR4'], groupby='leiden', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_leiden', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['EGFR','TLR4'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['EGFR','TLR4'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_group', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['EGFR','TLR4'], groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_cell_type-group', fig_size={'w': 6, 'h': 6})


adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Muscle/Muscle.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/07.Mission/Muscle_geneplot/"
sc.pl.umap(adata, color=['TNC','MLKL'], size=25, show=False,cmap='RdBu_r',vmax='p99')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby='leiden', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_leiden', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_group', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_cell_type-group', fig_size={'w': 6, 'h': 6})

from scipy.stats import pearsonr
from scipy.stats import spearmanr
gene1 = 'TNC'
gene2 = 'MLKL'
cell_types = adata.obs['cell_type'].unique()
results = []
for ct in cell_types:
    adata_ct = adata[adata.obs['cell_type'] == ct]
    
    # 提取并转换
    expr1 = adata_ct.raw[:, gene1].X
    expr2 = adata_ct.raw[:, gene2].X
    
    # 转换为密集数组并展平
    if hasattr(expr1, 'toarray'):
        expr1 = expr1.toarray().flatten()
        expr2 = expr2.toarray().flatten()
    else:
        expr1 = np.ravel(expr1)
        expr2 = np.ravel(expr2)
    
    corr, pval = spearmanr(expr1, expr2)
    results.append({
        'cell_type': ct,
        'n_cells': adata_ct.n_obs,
        'gene1': gene1,
        'gene2': gene2,
        'spearman_r': corr,
        'p_value': pval
    })

results_df = pd.DataFrame(results)

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/Main_mapping/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/07.Mission/Main_plot/"
sc.pl.umap(adata, color=['TNC','MLKL'], size=25, show=False,cmap='RdBu_r',vmax='p99')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby='leiden', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_leiden', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_group', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['TNC','MLKL'], groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_gene_cell_type-group', fig_size={'w': 6, 'h': 6})