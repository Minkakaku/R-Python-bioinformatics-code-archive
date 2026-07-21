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

# 调整大类注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Bcells', 'Distal tubular cells', 'Endothelial cells', 'Fibroblasts', 'Intercalated cells', 'Macrophages', 'Podocytes', 'Principal cells', 'Proximal tubule cells', 'Tcells', 'Mesangial cell'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/Main/"

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace('Bcells', 'B cells')
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace('Tcells', 'T cells')
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace('Mesangial cell', 'Mesangial cells')

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cluster', title='Cluster', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.pdf'))

sc.tl.rank_genes_groups(adata_raw, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)

sc.pl.dotplot(adata_raw, var_names=['Cd79a','Ly6d','Cd22','Kng2','Umod','Egf','Egfl7','Eng','Pecam1','Bgn','Cfh','Col1a2','Atp6v1g3','Car2','Atp6v0d2','Ctss','C1qb','C1qa','Podxl','Nphs2','Ptpro','Hsd11b2','Cdh16','Ivns1abp','Spink1','Pdzk1','Acsm2','Cd3d','Cd3e','Cd3g','Cald1', 'Mylk', 'Acta2'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')
gene_dict = {
    'B cells': ['Cd79a','Ly6d','Cd22'],
    'Distal tubular cells': ['Kng2','Umod','Egf'],
    'Endothelial cells': ['Egfl7','Eng','Pecam1'],
    'Fibroblasts': ['Bgn','Cfh','Col1a2'],
    'Intercalated cells': ['Atp6v1g3','Car2','Atp6v0d2'],
    'Macrophages': ['Ctss','C1qb','C1qa'],
    'Podocytes': ['Podxl','Nphs2','Ptpro'],
    'Principal cells': ['Hsd11b2','Cdh16','Ivns1abp'],
    'Proximal tubule cells': ['Spink1','Pdzk1','Acsm2'],
    'Tcells': ['Cd3d','Cd3e','Cd3g'],
    'Mesangial cell': ['Cald1', 'Mylk', 'Acta2']
}
sc.pl.dotplot(adata_raw, var_names=gene_dict, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')



# 调整巨噬注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_recluster_Macrophages/recluster_Macrophages.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'DC', 'Immature Mac C1', 'Immature Mac C2', 'Monocyte', 'Pro-inflammatory Mac', 'Pro-repair Mac', 'Proliferative Mac', 'Remodeling Mac-Svep1'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/Macrophages/"

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata_raw.obs['group'].cat.categories.to_list()
# if len(groups) > 1:
#     axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
#     fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
#     for i in range(len(groups)):
#         if i == len(groups) - 1:
#             legend = 'right margin'
#         else:
#             legend = None
#         sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.tl.rank_genes_groups(adata_raw, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)

sc.pl.dotplot(adata_raw, var_names=['Sept3','Gpr141b','Gcsam','Slc5a2','Slc5a12','Prss8','Cd209a','S100a4','Traf1','Cebpb','Tyrobp','Lyz2','Oas2','Ms4a7','Lilra5','Bcl6b','Myct1','Pimreg','Pbk','Hmmr','Nek2','Fbn1','Sema3f','Col1a1'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

with rc_context({'figure.figsize': (5, 1), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata_raw, var_names=['Sept3','Gpr141b','Gcsam','Slc5a2','Slc5a12','Prss8','Cd209a','S100a4','Traf1','Cebpb','Tyrobp','Lyz2','Oas2','Ms4a7','Lilra5','Pbk','Hmmr','Nek2','Bcl6b','Myct1','Pimreg','Fbn1','Sema3f','Col1a1'], groupby='cell_type', dendrogram=False, show=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'))

gene_dict = {
    'DC': ['Sept3','Gpr141b','Gcsam'],
    'Immature Mac C1': ['Slc5a2','Slc5a12','Prss8'],
    'Immature Mac C2': ['Cd209a','S100a4','Traf1'],
    'Monocyte': ['Cebpb','Tyrobp','Lyz2'],
    'Pro-inflammatory Mac': ['Oas2','Ms4a7','Lilra5',],
    'Proliferative Mac': ['Pbk','Hmmr','Nek2'],
    'Pro-repair Mac': ['Bcl6b','Myct1','Pimreg'],
    'Remodeling Mac-Svep1': ['Fbn1','Sema3f','Col1a1']
}
sc.pl.dotplot(adata_raw, var_names=gene_dict, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')


sc.pl.violin(adata_raw, ['Svep1'], jitter=0, groupby='cell_type', rotation=45, multi_panel=True, show=False)
save_fig(plt, out_dir, file_name='Svep1_voilon', fig_size={'width':25, 'height':10})
plt.close()
sc.pl.dotplot(adata_raw, var_names='Svep1', groupby='cell_type', dendrogram=False, show=False)
save_fig(plt, out_dir, file_name='Svep1_dotplot', fig_size={'width':5, 'height':5})
plt.close()