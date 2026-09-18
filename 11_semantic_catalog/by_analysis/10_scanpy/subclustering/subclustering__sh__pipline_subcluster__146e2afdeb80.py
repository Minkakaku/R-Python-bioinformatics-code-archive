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

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250070/01.Mission/mapping_recluster/mapping_recluster.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Collecting duct_intercalated cells','Collecting duct_principal cells', 'Distal tuble cells', 'EC_1', 'EC_2','Loop of henle', 'Macrophages', 'Myofibroblast', 'NK cells','Neutrophils', 'Podocyte', 'Proximal tubule cells', 'T cells'

out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250070/02.Mission/mapping_recluster"
os.makedirs(out_dir, exist_ok=True)
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':6,'h':12})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':6,'h':12})
crosstab_table = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'], normalize='index')
crosstab_table.to_csv(os.path.join(out_dir,"celltype_proportion_group.xls"), sep='\t')
crosstab_table = pd.crosstab(adata_raw.obs['sample'], adata_raw.obs['cell_type'], normalize='index')
crosstab_table.to_csv(os.path.join(out_dir,"celltype_proportion_sample.xls"), sep='\t')

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250070/01.Mission/PTC cells/PTC cells.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Healthy', 'PTC_Bcam', 'PTC_Btg2', 'PTC_Cdkn1a', 'PTC_Miox', 'PTC_Mt1'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250070/02.Mission/PTC cells_cell_type/"
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250070/02.Mission/PTC cells/"
os.makedirs(out_dir, exist_ok=True)
adata_raw.obs['cell_type'] = adata_raw.obs['leiden']
adata_raw.uns['cell_type_colors'] = adata_raw.uns['leiden_colors']
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample')
crosstab_table = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'], normalize='index')
crosstab_table.to_csv(os.path.join(out_dir,"celltype_proportion_group.xls"), sep='\t')
crosstab_table = pd.crosstab(adata_raw.obs['sample'], adata_raw.obs['cell_type'], normalize='index')
crosstab_table.to_csv(os.path.join(out_dir,"celltype_proportion_sample.xls"), sep='\t')