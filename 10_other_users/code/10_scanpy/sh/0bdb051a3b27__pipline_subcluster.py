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

def read_gmt_to_dict(gmt_file):
    """
    读取 GMT 文件，返回字典 {基因集名称: [基因列表]}
    """
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                set_name = parts[0]
                description = parts[1]
                genes = parts[2:]
                gene_sets[set_name] = genes
    return gene_sets

# 基因集评分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2512157492-01/SC2512157492-01/require1/all_data.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512157492-01/02.Mission/geneset_score"
gene_sets = read_gmt_to_dict("/mnt/NexData_file/home/sh/Projects/SC2512157492-01/02.Mission/script/select.gmt")

adata=adata_raw.raw.to_adata()
for set_name, gene_list in gene_sets.items():
    available_genes = [g for g in gene_list if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=available_genes, score_name=set_name)

for set_name, gene_list in gene_sets.items():
    with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
        sc.pl.umap(adata, color=set_name, title=set_name, show=False,vmax='p99', cmap='coolwarm')
        plt.savefig(os.path.join(out_dir, f'dimplot_umap_{set_name}.png'))
        plt.savefig(os.path.join(out_dir, f'dimplot_umap_{set_name}.pdf'))

sc.pl.dotplot(adata, var_names=list(gene_sets.keys()), groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_sub_cell_type.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_sub_cell_type.pdf'), bbox_inches = 'tight')

sc.pl.dotplot(adata, var_names=list(gene_sets.keys()), groupby='celltype', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.pdf'), bbox_inches = 'tight')

sc.pl.umap(adata_raw, color=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], cmap='coolwarm',show=False,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'), bbox_inches = 'tight')
sc.pl.dotplot(adata, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_sub_cell_type.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_sub_cell_type.pdf'), bbox_inches = 'tight')

sc.pl.dotplot(adata, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='celltype', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.pdf'), bbox_inches = 'tight')


sc.pl.matrixplot(adata, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='sample', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_gene_sample.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_gene_sample.pdf'),bbox_inches='tight')

adata_sub=adata[adata.obs['celltype'].isin(['SPC']), :]
sc.pl.matrixplot(adata_sub, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='sample', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_gene_SPC_sample.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_gene_SPC_sample.pdf'),bbox_inches='tight')
adata_sub=adata[adata.obs['celltype'].isin(['SPG']), :]
sc.pl.matrixplot(adata_sub, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='sample', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_gene_SPG_sample.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_gene_SPG_sample.pdf'),bbox_inches='tight')
adata_sub=adata[adata.obs['celltype'].isin(['Spermatids']), :]
sc.pl.matrixplot(adata_sub, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='sample', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_gene_Spermatids_sample.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_gene_Spermatids_sample.pdf'),bbox_inches='tight')
adata_sub=adata[adata.obs['cell_type'].isin(['Undiff SPG Pro']), :]
sc.pl.matrixplot(adata_sub, var_names=['Smcp', 'Ppp3ca', 'Ppp3cb', 'Ppp3cc', 'Tppp', 'Tppp2', 'Tppp3'], groupby='sample', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_gene_Undiff SPG Pro_sample.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_gene_Undiff SPG Pro_sample.pdf'),bbox_inches='tight')

out_dir="/mnt/NexData_file/home/sh/Projects/SC2512157492-01/02.Mission/Main"
# cell_markers = {
#     "Fibroblasts": ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Fap", "Pdgfra", "Vim"],  # 成纤维细胞
#     "Leydig cells": ["Star", "Cyp11a1", "Hsd3b1", "Lhcgr", "Cyp17a1", "Insl3"],  # 间质细胞 
#     "Macrophages": ["Cd68", "Adgre1", "Cd163", "Lyve1", "Mrc1", "Csf1r"],  # 巨噬细胞（睾丸巨噬细胞标志物）
#     "SPC": ["Spo11", "Ovol2", "Nme8", "Id4"],
#     "SPG": ["Kit", "Gfra1", "Dazl", "Stra8", "Sycp3","Dmrt1","Ret"],
#     "Sertoli cells": ["Pebp1","Dynll1","Isyna1"],  # 支持细胞
#     "Spermatids": ["Prm1", "Prm2", "Tnp1", "Tnp2", "Acrv1", "Spaca1", "Odf1", "Tekt1"],  # 精子细胞（Spermatids）- 精子形成阶段
#     "T cells": ["Cd3d", "Cd3e", "Cd4", "Cd8a", "Trbc2"]  # T细胞
# }
cell_markers = {
    "Fibroblasts": ["Dcn", "Col1a2","Col1a1"],  # 成纤维细胞
    "Leydig cells": ["Cyp17a1", "Cyp11a1", "Insl3"],  # 间质细胞 
    "Macrophages": ["Csf1r", "Mrc1", "Adgre1"],  # 巨噬细胞（睾丸巨噬细胞标志物）
    "SPC": ["Ovol2", "Spo11", "Id4"],
    "SPG": ["Dazl", "Stra8", "Dmrt1"],
    "Sertoli cells": ["Pebp1","Dynll1","Isyna1"],  # 支持细胞
    "Spermatids": ["Prm1", "Prm2", "Tnp1",],  # 精子细胞（Spermatids）- 精子形成阶段
    "T cells": ["Cd3d", "Cd3e", "Cd8a"]  # T细胞
}
sc.pl.dotplot(adata_raw, var_names=cell_markers, groupby='celltype', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.pdf'), bbox_inches = 'tight')
sc.pl.dotplot(adata_raw, var_names=cell_markers, groupby='celltype', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type_scale.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type_scale.pdf'), bbox_inches = 'tight')

sc.pl.umap(adata_raw, color=['Dcn', 'Cyp17a1', 'Csf1r', 'Ovol2', 'Dazl', 'Pebp1', 'Prm1',"Cd3d"], cmap='coolwarm',show=False,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'), bbox_inches = 'tight')
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'), bbox_inches = 'tight')