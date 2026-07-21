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

# 巨噬-平滑肌细分 重新注释
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/05.Mission/subcluster_Muti/"
os.makedirs(out_dir, exist_ok=True)
adata_raw=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250326/04.Mission/mapping_raw/mapping_raw.h5ad")
adata_raw.obs['cell_type'].cat.categories
adata_raw.obs['sample'].cat.categories
adata_raw.obs['group'].cat.categories

adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['Macrophages','Smooth muscle cells'])]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pl.umap(adata_sub, color=['Acta2','Mcam','Myh11','Myl9','Pdgfa','Tagln','Adgre1','Apoe','C1qa','C1qb','C1qc','Cd68','Cd86','Irf8','Ccr2','Cx3cr1','H2-Aa','H2-Ab1','H2-Eb1','Mafb','Ms4a7'], ncols=5,cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'marker_plot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'marker_plot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.3,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Acta2','Mcam','Myh11','Myl9','Pdgfa','Tagln','Adgre1','Apoe','C1qa','C1qb','C1qc','Cd68','Cd86','Irf8','Ccr2','Cx3cr1','H2-Aa','H2-Ab1','H2-Eb1','Mafb','Ms4a7'], ncols=5,cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'subcluster_marker_plot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'subcluster_marker_plot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)