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


#### 大类注释sample——UMAP,sample_barplot
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2605263926-01/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'ATC', 'B cells', 'Endothelial cells', 'Fibroblasts', 'Follicular', 'Macrophages', 'NK cells', 'PTC', 'Pericytes', 'Plasma cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2605263926-01/01.Mission/geneset_plot/"

TDS_geneset=['PAX8','NKX2-1','TTF1','FOXE1','TTF2','GLIS3','TG','TPO','TSHR','IYD','SLC5A5','SLC26A4','SLC26A7','DUOX1','DUOX2','DIO1','DIO2','THADA']
BRAF_like_geneset=['MET','CXCL10','TIMP1','FN1','CLDN10','DPP4','SERPINA1','PROM1','CD133','S100A4','LGALS3','KRT19']
RAS_like_geneset=['TPO','DIO1','DIO2','SLC26A4','DUOX1','DUOX2','ESRRG','IGFBP5','TG']

adata=adata_raw.raw.to_adata()
sc.tl.score_genes(adata, gene_list=TDS_geneset, score_name='TDS_geneset')
sc.tl.score_genes(adata, gene_list=BRAF_like_geneset, score_name='BRAF_like_geneset')
sc.tl.score_genes(adata, gene_list=RAS_like_geneset, score_name='RAS_like_geneset')
sc.pl.umap(adata, color=['TDS_geneset','BRAF_like_geneset','RAS_like_geneset'], size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_genesets', fig_size={'w': 21, 'h': 6})
sc.pl.umap(adata, color=['TDS_geneset','BRAF_like_geneset','RAS_like_geneset'], size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_genesets_filter', fig_size={'w': 21, 'h': 6})

import seaborn as sns
from scipy.stats import mannwhitneyu
# 提取数据
df = adata.obs[['group', 'TDS_geneset', 'BRAF_like_geneset','RAS_like_geneset']].copy()
df = df.dropna()
# 获取分组顺序
groups = df['group'].unique().tolist()
# 计算所有组间比较的 p-value
def get_significance(data1, data2):
    if len(data1) > 0 and len(data2) > 0:
        stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
        if p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        else:
            return 'ns'
    return 'ns'

# 绘制小提琴图
fig, axes = plt.subplots(1, 3, figsize=(24, 8), gridspec_kw={'wspace': 0.5})
for idx, feature in enumerate(['TDS_geneset', 'BRAF_like_geneset', 'RAS_like_geneset']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.violinplot(data=df, x='group', y=feature, order=groups, ax=ax, 
               color='black', cut=0, linewidth=1.5,fill=False,inner="quart",common_norm=True,density_norm='width',
               saturation=1, alpha=1)
    # sns.stripplot(data=df, x='group', y=feature, order=groups, ax=ax,
              # palette='RdBu_r', alpha=0.3, size=5, edgecolor='none', jitter=True)
    # 计算并添加显著性
    y_max = df[feature].max()
    y_min = df[feature].min()
    y_range = y_max - y_min
    # 组间比较
    comparisons = []
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            data1 = df[df['group'] == groups[i]][feature]
            data2 = df[df['group'] == groups[j]][feature]
            p_sig = get_significance(data1, data2)
            comparisons.append(((groups[i], groups[j]), p_sig))
    # 添加显著性标注
    for k, ((g1, g2), sig) in enumerate(comparisons):
        x1 = groups.index(g1)
        x2 = groups.index(g2)
        y = y_max + (k + 1) * y_range * 0.08
        ax.plot([x1, x1, x2, x2], [y - y_range*0.02, y, y, y - y_range*0.02], 
                'k-', lw=1.5)
        ax.text((x1 + x2) / 2, y + y_range*0.01, sig, 
                ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.set_title(feature, fontsize=14)
    ax.set_xlabel('Group', fontsize=12)
    ax.set_ylabel(feature, fontsize=12)

plt.tight_layout()
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_genesets_sig', fig_size={'w': 12, 'h': 8})

with rc_context({'figure.figsize': (5, 15), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['TDS_geneset', 'BRAF_like_geneset', 'RAS_like_geneset'], groupby='cell_type', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.pdf'))

with rc_context({'figure.figsize': (5,15), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['TDS_geneset', 'BRAF_like_geneset', 'RAS_like_geneset'], groupby='group', dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_genesets_group.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_genesets_group.pdf'))

with rc_context({'figure.figsize': (5, 15), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['TDS_geneset', 'BRAF_like_geneset', 'RAS_like_geneset'], groupby=['cell_type','group'], dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type_group.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type_group.pdf'))


sc.pl.matrixplot(adata, var_names=TDS_geneset, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_TDS_geneset_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_TDS_geneset_cell_type.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=TDS_geneset, groupby='group', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_TDS_geneset_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_TDS_geneset_group.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=TDS_geneset, groupby=['cell_type','group'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_TDS_geneset_cell_type_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_TDS_geneset_cell_type_group.pdf'),bbox_inches='tight')
plt.close()

BRAF_like_geneset=['MET','CXCL10','TIMP1','FN1','CLDN10','DPP4','SERPINA1','PROM1','S100A4','LGALS3','KRT19']
sc.pl.matrixplot(adata, var_names=BRAF_like_geneset, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_BRAF_like_geneset_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_BRAF_like_geneset_cell_type.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=BRAF_like_geneset, groupby='group', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_BRAF_like_geneset_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_BRAF_like_geneset_group.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=BRAF_like_geneset, groupby=['cell_type','group'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_BRAF_like_geneset_cell_type_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_BRAF_like_geneset_cell_type_group.pdf'),bbox_inches='tight')
plt.close()

sc.pl.matrixplot(adata, var_names=RAS_like_geneset, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_RAS_like_geneset_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_RAS_like_geneset_cell_type.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=RAS_like_geneset, groupby='group', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_RAS_like_geneset_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_RAS_like_geneset_group.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=RAS_like_geneset, groupby=['cell_type','group'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'heatmap_RAS_like_geneset_cell_type_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_RAS_like_geneset_cell_type_group.pdf'),bbox_inches='tight')
plt.close()

marker_genes=TDS_geneset+BRAF_like_geneset+RAS_like_geneset+['PPARG']
sc.pl.umap(adata, color=marker_genes, show=False, cmap='coolwarm',ncols=5)
save_fig(plt,file_name= f'genes_umap',out_dir=out_dir,fig_size={'w':25,'h':40})