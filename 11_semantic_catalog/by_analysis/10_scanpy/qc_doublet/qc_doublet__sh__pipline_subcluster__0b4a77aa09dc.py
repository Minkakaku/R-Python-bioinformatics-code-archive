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


# 对11.Myofibroblasts注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/Main_mapping/Main_mapping.h5ad")
# adata_raw.obs.cell_type.cat.categories
# '1.B cells', '10.pDC', '11.Myofibroblasts', '2.Endothelial cells', '3.Fibroblasts', '4.Macrophages', '5.Muscle cells', '6.Pericytes', '7.Plasma cells', '8.Satelite cells', '9.T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/05.Mission/Myofibroblasts/"
sc.pl.umap(adata_raw, color=['ACTA2','COL1A1','COL3A1','POSTN','MYH11','CNN1','RGS5'], size=25, show=False,cmap='RdBu_r')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.pdf'))
sc.pl.umap(adata_raw, color=['n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo'], size=25, show=False,cmap='RdBu_r')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_QC.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_QC.pdf'))

# sc.pl.violin(adata_raw, keys=['ACTA2','COL1A1','COL3A1','POSTN','MYH11','CNN1','RGS5'], groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
# save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_marker', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata_raw, var_names=['ACTA2','COL1A1','COL3A1','POSTN','MYH11','CNN1','RGS5'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_marker', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata_raw, var_names=['ACTA2','COL1A1','COL3A1','POSTN','MYH11','CNN1','RGS5'], groupby=['cell_type','group'],  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_marker_group', fig_size={'w': 6, 'h': 12})

groups = adata_raw.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})


# 对4.Macrophages注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/Main_mapping/Main_mapping.h5ad")
# adata_raw.obs.cell_type.cat.categories
# '1.B cells', '10.pDC', '11.Myofibroblasts', '2.Endothelial cells', '3.Fibroblasts', '4.Macrophages', '5.Muscle cells', '6.Pericytes', '7.Plasma cells', '8.Satelite cells', '9.T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/05.Mission/Macrophages/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="4.Macrophages", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['SPP1','TGFB1'], size=25, show=False,cmap='RdBu_r',vmax=1)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.pdf'))
sc.pl.violin(adata, keys=['SPP1','TGFB1'], rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_marker', fig_size={'w': 12, 'h': 6})

anno_dict={
    '0':'SCN9A Macrophage',
    '1':'DIRC3 Macrophage',
    '2':'PVT1 Macrophage',
    '3':'Profibrotic Macrophage',
    '4':'TPRG1 Macrophage',
    '5':'TEX41 Macrophage',
    '6':'TRDN Macrophage'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

genes=['SPP1','TGFB1']
sc.tl.score_genes(adata, gene_list=genes, score_name='dp_gene_score')
sc.pl.umap(adata, color='dp_gene_score', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_dp_gene_score', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='dp_gene_score', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_dp_gene_score', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['dp_gene_score'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_dp_gene_score', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['dp_gene_score'], groupby=['cell_type','group'],  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_dp_gene_score_group', fig_size={'w': 6, 'h': 12})

# cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
# adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
# adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
# adata_raw.obs['cell_type'].cat.categories
# adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

groups = adata.obs['group'].drop_duplicates().tolist()
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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

# sc.pl.matrixplot(adata, var_names=genes, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
# plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type.png'), dpi=300, bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type.pdf'),bbox_inches='tight')
# plt.close()
# sc.pl.matrixplot(adata, var_names=genes, groupby='group', standard_scale='var', dendrogram=False, show=False)
# plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_group.png'), dpi=300, bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_group.pdf'),bbox_inches='tight')
# plt.close()
# sc.pl.matrixplot(adata, var_names=genes, groupby=['cell_type','group'], standard_scale='var', dendrogram=False, show=False)
# plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type-group.png'), dpi=300, bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type-group.pdf'),bbox_inches='tight')
# plt.close()


out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/05.Mission/TNC_gene/"
sc.pl.umap(adata_raw, color=['ITGAV','ITGB1','ITGB6','ITGA5','ITGB3','TLR4'], size=25, show=False,cmap='RdBu_r')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_TNC.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_TNC.pdf'))
sc.pl.dotplot(adata_raw, var_names=['ITGAV','ITGB1','ITGB6','ITGA5','ITGB3','TLR4'], swap_axes=True, groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_TNC', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata_raw, var_names=['ITGAV','ITGB1','ITGB6','ITGA5','ITGB3','TLR4'], swap_axes=True, groupby=['cell_type','group'],  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_TNC_group', fig_size={'w': 12, 'h': 6})

#### 没做映射需要调整结果
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
# adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
sc.pl.dotplot(adata_raw, var_names=['ITGAV','ITGB1','ITGB6','ITGA5','ITGB3','TLR4'], swap_axes=True, groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_TNC', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata_raw, var_names=['ITGAV','ITGB1','ITGB6','ITGA5','ITGB3','TLR4'], swap_axes=True, groupby=['cell_type','group'],  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_TNC_group', fig_size={'w': 12, 'h': 6})