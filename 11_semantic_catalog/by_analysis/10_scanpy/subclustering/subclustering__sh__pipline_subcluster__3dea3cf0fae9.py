import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from harmonypy import run_harmony
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
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').iloc[:, ::-1].plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').iloc[:, ::-1].columns],
        edgecolor='black', width=0.8
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(reverse=True,
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
    file_name=f"Macrophages_celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.h5ad")
out_dir = "/mnt/NexData_file/home/sh/Projects/FWSC20241287/14.Mission/fre/"
# adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_recluster_Macrophages/recluster_Macrophages.h5ad")
adata.obs['cell_type'] = adata.obs['cell_type'].cat.add_categories('Remodeling Mac')
mask = adata.obs['cell_type'] == 'Remodeling Mac-Svep1'
adata.obs.loc[mask, 'cell_type'] = 'Remodeling Mac'
adata.obs['cell_type'] = adata.obs['cell_type'].cat.remove_unused_categories()
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

gene_dict = {
    'DC': ['Sept3','Gpr141b','Gcsam'],
    'Immature Mac C1': ['Slc5a2','Slc5a12','Prss8'],
    'Immature Mac C2': ['Cd209a','S100a4','Traf1'],
    'Monocyte': ['Cebpb','Tyrobp','Lyz2'],
    'Pro-inflammatory Mac': ['Oas2','Ms4a7','Lilra5',],
    'Proliferative Mac': ['Pbk','Hmmr','Nek2'],
    'Pro-repair Mac': ['Bcl6b','Myct1','Pimreg'],
    'Remodeling Mac': ['Fbn1','Sema3f','Col1a1']
}
sc.pl.dotplot(adata, var_names=gene_dict, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')
