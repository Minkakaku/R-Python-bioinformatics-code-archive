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


# 对Main注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/04.Mission/Main/"
adata_sub = adata_raw[adata_raw.obs['group'].isin(['Tumor','NonTumor']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.8,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

cell_markers = {
    # === 上皮细胞 ===
    "Epithelial cells": ["EPCAM", "KRT14", "KRT15"],
    "T cells": ["CD3G", "CD3D", "CD3E"],
    "NK cells": ["NCAM1", "NKG7",  "GNLY"],
    "B cells": ["MS4A1", "CD19", "CD79A"],
    "Plasma cells": ["JCHAIN", "TNFRSF17", "SDC1", "MZB1"],
    "Macrophages": ["CD68", "CSF1R", "ITGAM"],
    "Neutrophils": ["CEACAM3", "FCGR3B", "CXCR2"],
    "Mast cells": ["TPSAB1", "CPA3", "MS4A2", "HDC"],
    "Fibroblasts": ["COL1A1", "COL3A1", "DCN", "FN1", "BGN", "FGF7"],
    "Smooth muscle cells": ["ACTA2", "CNN1", "MYH11", "TAGLN"],
    "Endothelial cells": ["PECAM1", "VWF", "CLDN5", "FLT1"],
    "Pericytes": ["PDGFRA", "CSPG4", "RGS5", "MCAM"]
}
markers = [gene for markers in cell_markers.values() for gene in markers]
sc.pl.umap(adata, color=markers, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Fibroblasts',
    '1':'Smooth muscle cells',
    '2':'Fibroblasts',
    '3':'T cells',
    '4':'Endothelial cells',
    '5':'Fibroblasts',
    '6':'Epithelial cells',
    '7':'Smooth muscle cells',
    '8':'Fibroblasts',
    '9':'Smooth muscle cells',
    '10':'Fibroblasts',
    '11':'Fibroblasts',
    '12':'Pericytes',
    '13':'Pericytes',
    '14':'Pericytes',
    '15':'Epithelial cells',
    '16':'Epithelial cells',
    '17':'Macrophages',
    '18':'Fibroblasts',
    '19':'Epithelial cells',
    '20':'Mast cells',
    '21':'Plasma cells',
    '22':'Smooth muscle cells',
    '23':'Neutrophils'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')
# adata_raw.obs['group'] = adata_raw.obs['group'].replace('GAStis', 'Tumor')
# adata_raw.obs['group'] = adata_raw.obs['group'].replace('GASPtis', 'NonTumor')
# adata_raw.obs['group'] = adata_raw.obs['group'].replace('GASpbmc', 'PBMC')
# adata_raw.obs['group'] = (
#     pd.Categorical(adata_raw.obs['group'],
#                    categories=['Tumor','NonTumor','PBMC'],
#                    ordered=True)
# )
# groups = adata_raw.obs['group'].drop_duplicates().tolist()
groups = ['Tumor','NonTumor']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Main.h5ad'))

cell_markers = {
    "Endothelial cells": ["PECAM1", "VWF", "CLDN5", "FLT1"],
    "Epithelial cells": ["EPCAM", "KRT19", "KRT8"],
    "Fibroblasts": ["COL1A1", "COL3A1", "DCN"],
    "Macrophages": ["CD68", "CSF1R", "ITGAM"],
    "Mast cells": ["CPA3", "MS4A2", "TPSAB1"],
    "Neutrophils": ["S100A8", "S100A9", "CSF3R"],
    "Pericytes": ["RGS5", "NOTCH3", "PDGFRB"],
    "Plasma cells": ["JCHAIN", "TNFRSF17", "MZB1"],
    "Smooth muscle cells": ["ACTA2", "CNN1", "MYH11"],
    "T cells": ["CD3G", "CD3D", "CD3E"]
}
sc.pl.dotplot(adata, var_names=cell_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


# 对Main_PBMC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510215750-01/01.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510215750-01/04.Mission/Main_PBMC/"
adata_sub = adata_raw[adata_raw.obs['group'].isin(['PBMC']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.8,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

cell_markers = {
    "T cells": ["CD3G", "CD3D", "CD3E"],
    "Treg": ["FOXP3", "CTLA4", "IL2RA"],
    "NK cells": ["NCAM1", "NKG7",  "GNLY"],
    "B cells": ["MS4A1", "CD19", "CD79A"],
    "Plasma cells": ["JCHAIN", "TNFRSF17", "SDC1", "MZB1"],
    "Macrophages": ["CD68", "CSF1R", "ITGAM"],
    "Monocytes": ["CD74", "VCAN", "CD14"],
    "Dendritic cells": ["HLA-DRA", "CD1C", "CLEC9A"],
    "Neutrophils": ["CEACAM3", "FCGR3B", "CXCR2"],
    "Mast cells": ["CPA3", "MS4A2", "TPSAB1"]
}
markers = [gene for markers in cell_markers.values() for gene in markers]
sc.pl.umap(adata, color=markers, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '0':'Monocytes',
    '1':'T cells',
    '2':'NK cells',
    '3':'T cells',
    '4':'Monocytes',
    '5':'T cells',
    '6':'T cells',
    '7':'T cells',
    '8':'Monocytes',
    '9':'NK cells',
    '10':'T cells',
    '11':'B cells',
    '12':'T cells',
    '13':'Monocytes',
    '14':'Dendritic cells',
    '15':'Monocytes'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

groups = ['PBMC']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})
adata.write(os.path.join(out_dir, 'Main.h5ad'))

cell_markers = {
    "B cells": ["MS4A1", "CD19", "CD79A"],
    "Dendritic cells": ["HLA-DRA", "CD1C", "CLEC9A"],
    "Monocytes": ["CD74", "VCAN", "CD14"],
    "NK cells": ["NCAM1", "NKG7",  "GNLY"],
    "T cells": ["CD3G", "CD3D", "CD3E"]
}
sc.pl.dotplot(adata, var_names=cell_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), facecolor='white', bbox_inches="tight", dpi=300)