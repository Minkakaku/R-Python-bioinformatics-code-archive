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


adata=sc.read('/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Fibroblasts/Fibroblasts.h5ad')
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Fibroblasts/"
sc.tl.rank_genes_groups(adata, groupby='group', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)
marker_list = ['PDGFRB','IGF1','ESR1']
with rc_context({'figure.figsize': (7, 7), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=marker_list, groupby='group', dendrogram=False, show=False, standard_scale='var')
    plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'))

marker_list={
    'Antiogenesis-related genes':['THY1','EDNRA','COL8A1','SRPX2','TMEM211','COL27A1','LEF1','VSTM4','GREM1','ECM1','LOXL2','HSPG2','PXDN','UNC5B','TGFBI','PLCD3','COL8A2','ANGPT2','PROK2','NOTCH3','CSPG4','SFRP2','FGF1','VASH1','NRXN1','NPR2'],
    'Growth factor':['RPS2','GPC1','HGF','WNT5A','BMP4','PDGFA','WNT4','IGF1','PDGFC'],
    'Ch_ILR':['CXCL14','CXCL16','IL1B'],
    'Extracellular matrix organization':['MMP11','COL1A1','COL1A2','COL3A1','SPARC','EMILIN1','COL5A1','COL5A2','COL10A1','CREB3L1','MMP7','COL12A1','ADAMTS16','SULF1','SULF2','COMP','BMP1','LOXL1','P3H4','MATN3','KAZALD1','KLK4','OLFML2B','LOXL4','OLFML2A','MMP28','COL11A1','SPINT2','SCARA3','MMP3','MMP10'],
}
with rc_context({'figure.figsize': (7, 7), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=marker_list, groupby='group', dendrogram=False, show=False, standard_scale='var')
    plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'))


adata=sc.read("/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Epithelial/Epithelial.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Epithelial/"
sc.tl.rank_genes_groups(adata, groupby='group', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)