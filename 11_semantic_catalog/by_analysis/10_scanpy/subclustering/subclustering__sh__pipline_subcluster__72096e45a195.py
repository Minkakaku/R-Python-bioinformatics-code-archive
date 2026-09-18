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
import gseapy as gp

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


# 对Macrophage注释细分
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/05.Mission/Macrophages/Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/08.Mission/Macrophages/"
adata.obs['cell_type'].cat.categories
# 'Cd74 Macrophages', 'Celf2 Macrophages', 'Fabp4 Macrophages', 'Fn1 Macrophages', 'Monocyte', 'Trem2 Macrophages', 'cDCs'
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cell_type.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cell_type.pdf'))

gmt_file = "/mnt/NexData_file/home/sh/Projects/AE2511106372-01/08.Mission/script/genelist.gmt"
gene_sets = gp.read_gmt(gmt_file)
# 对每个基因集计算评分
adata_sub = adata[adata.obs['group'].isin(['MI','MI_DZ']), :]
adata_sub = adata_sub.raw.to_adata()
for set_name, genes in gene_sets.items():
    # 检查哪些基因存在
    available = [g for g in genes if g in adata_sub.var_names]
    missing = set(genes) - set(available)
    if missing:
        print(f"Warning: {set_name} missing genes: {missing}")
    if len(available) >= 3:
        sc.tl.score_genes(
            adata_sub, 
            gene_list=available,
            score_name=set_name
        )
        print(f"✓ {set_name}: {len(available)}/{len(genes)} genes used")
    else:
        print(f"✗ {set_name}: insufficient genes ({len(available)}/{len(genes)})")

gene_sets_keys = list(gene_sets.keys())
gene_sets_keys = ['VEGF signaling pathway', 'HALLMARK_ANGIOGENESIS', 'GOBP_ANGIOGENESIS']
sc.pl.dotplot(adata_sub, var_names=gene_sets_keys, groupby=['group'], dendrogram=False, show=False,standard_scale="var",swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist_scale.pdf'),bbox_inches='tight')
# gene_sets_keys = ['VEGF signaling pathway','GOBP_ANGIOGENESIS_INVOLVED_IN_WOUND_HEALING']
sc.pl.dotplot(adata_sub, var_names=gene_sets_keys, groupby=['group'], dendrogram=False, show=False,swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist.pdf'),bbox_inches='tight')
sc.pl.matrixplot(adata_sub, var_names=gene_sets['VEGF signaling pathway'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_VEGF_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_VEGF_pathway.pdf'),bbox_inches='tight')
sc.pl.matrixplot(adata_sub, var_names=gene_sets['HALLMARK_ANGIOGENESIS'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_HALLMARK_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_HALLMARK_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')
available_genes = list(set(gene_sets['GOBP_ANGIOGENESIS']) & set(adata_sub.var_names))
gene_sets['GOBP_ANGIOGENESIS']=available_genes
sc.pl.matrixplot(adata_sub, var_names=gene_sets['GOBP_ANGIOGENESIS'],groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_GOBP_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_GOBP_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')


# 计算两组的平均表达
expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['VEGF signaling pathway']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['VEGF signaling pathway']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['VEGF signaling pathway'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_VEGF_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_VEGF_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['HALLMARK_ANGIOGENESIS']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['HALLMARK_ANGIOGENESIS']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['HALLMARK_ANGIOGENESIS'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_HALLMARK_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_HALLMARK_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['GOBP_ANGIOGENESIS']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['GOBP_ANGIOGENESIS']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['GOBP_ANGIOGENESIS'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False).head(30)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_GOBP_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_GOBP_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')



adata_sub = adata_sub[adata_sub.obs['cell_type'].isin(['Trem2 Macrophages']), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/08.Mission/Trem2_Macrophages/"
sc.pl.dotplot(adata_sub, var_names=gene_sets_keys, groupby=['group'], dendrogram=False, show=False,standard_scale="var",swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist_scale.pdf'),bbox_inches='tight')
# gene_sets_keys = ['VEGF signaling pathway','GOBP_ANGIOGENESIS_INVOLVED_IN_WOUND_HEALING']
sc.pl.dotplot(adata_sub, var_names=gene_sets_keys, groupby=['group'], dendrogram=False, show=False,swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist.pdf'),bbox_inches='tight')
sc.pl.matrixplot(adata_sub, var_names=gene_sets['VEGF signaling pathway'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_VEGF_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_VEGF_pathway.pdf'),bbox_inches='tight')
sc.pl.matrixplot(adata_sub, var_names=gene_sets['HALLMARK_ANGIOGENESIS'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_HALLMARK_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_HALLMARK_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')
available_genes = list(set(gene_sets['GOBP_ANGIOGENESIS']) & set(adata_sub.var_names))
gene_sets['GOBP_ANGIOGENESIS']=available_genes
sc.pl.matrixplot(adata_sub, var_names=gene_sets['GOBP_ANGIOGENESIS'],groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_GOBP_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_GOBP_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')


# 计算两组的平均表达
expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['VEGF signaling pathway']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['VEGF signaling pathway']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['VEGF signaling pathway'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_VEGF_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_VEGF_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['HALLMARK_ANGIOGENESIS']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['HALLMARK_ANGIOGENESIS']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['HALLMARK_ANGIOGENESIS'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_HALLMARK_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_HALLMARK_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['GOBP_ANGIOGENESIS']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['GOBP_ANGIOGENESIS']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['GOBP_ANGIOGENESIS'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False).head(30)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_GOBP_ANGIOGENESIS_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_GOBP_ANGIOGENESIS_pathway.pdf'),bbox_inches='tight')



gmt_file = "/mnt/NexData_file/home/sh/Projects/AE2511106372-01/08.Mission/script/trem2.gmt"
gene_sets = gp.read_gmt(gmt_file)
for set_name, genes in gene_sets.items():
    # 检查哪些基因存在
    available = [g for g in genes if g in adata_sub.var_names]
    missing = set(genes) - set(available)
    if missing:
        print(f"Warning: {set_name} missing genes: {missing}")
    if len(available) >= 3:
        sc.tl.score_genes(
            adata_sub, 
            gene_list=available,
            score_name=set_name
        )
        print(f"✓ {set_name}: {len(available)}/{len(genes)} genes used")
    else:
        print(f"✗ {set_name}: insufficient genes ({len(available)}/{len(genes)})")


out_dir="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/08.Mission/Trem2_Macrophages_pathway/"
gene_sets_keys = ['PI3K-Akt signaling pathway', 'VEGF signaling pathway', 'MAPK signaling pathway','GOBP_ERK5_CASCADE']
sc.pl.dotplot(adata_sub, var_names=gene_sets_keys, groupby=['group'], dendrogram=False, show=False,standard_scale="var",swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist_scale.pdf'),bbox_inches='tight')
# gene_sets_keys = ['VEGF signaling pathway','GOBP_ANGIOGENESIS_INVOLVED_IN_WOUND_HEALING']
sc.pl.dotplot(adata_sub, var_names=gene_sets_keys, groupby=['group'], dendrogram=False, show=False,swap_axes=True)
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_group_genelist.pdf'),bbox_inches='tight')
sc.pl.matrixplot(adata_sub, var_names=gene_sets['VEGF signaling pathway'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_VEGF_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_VEGF_pathway.pdf'),bbox_inches='tight')
available_genes = list(set(gene_sets['MAPK signaling pathway']) & set(adata_sub.var_names))
gene_sets['MAPK signaling pathway']=available_genes
sc.pl.matrixplot(adata_sub, var_names=gene_sets['MAPK signaling pathway'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_MAPK_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_MAPK_pathway.pdf'),bbox_inches='tight')
sc.pl.matrixplot(adata_sub, var_names=gene_sets['GOBP_ERK5_CASCADE'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_GOBP_ERK5_CASCADE_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_GOBP_ERK5_CASCADE_pathway.pdf'),bbox_inches='tight')
available_genes = list(set(gene_sets['PI3K-Akt signaling pathway']) & set(adata_sub.var_names))
gene_sets['PI3K-Akt signaling pathway']=available_genes
sc.pl.matrixplot(adata_sub, var_names=gene_sets['PI3K-Akt signaling pathway'], groupby='group',standard_scale="var")
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_PI3K-Akt_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'MatrixPlot_group_PI3K-Akt_pathway.pdf'),bbox_inches='tight')


# 计算两组的平均表达
expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['VEGF signaling pathway']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['VEGF signaling pathway']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['VEGF signaling pathway'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_VEGF_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_VEGF_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['PI3K-Akt signaling pathway']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['PI3K-Akt signaling pathway']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['PI3K-Akt signaling pathway'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False).head(30)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_PI3K-Akt_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_PI3K-Akt_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['MAPK signaling pathway']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['MAPK signaling pathway']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['MAPK signaling pathway'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False).head(30)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_MAPK_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_MAPK_pathway.pdf'),bbox_inches='tight')

expr_A =  np.array(adata_sub[adata_sub.obs['group'] == 'MI_DZ', gene_sets['GOBP_ERK5_CASCADE']].X.mean(axis=0)).flatten()
expr_other = np.array(adata_sub[adata_sub.obs['group'] == 'MI', gene_sets['GOBP_ERK5_CASCADE']].X.mean(axis=0)).flatten()
results = pd.DataFrame({
    'gene': gene_sets['GOBP_ERK5_CASCADE'],
    'mean_A': expr_A,
    'mean_other': expr_other
})
higher_in_MI_DZ = results[results['mean_A'] > results['mean_other']].sort_values('mean_A', ascending=False)
sc.pl.heatmap(adata_sub, var_names=higher_in_MI_DZ['gene'],standard_scale= 'var', groupby='group', cmap='YlOrRd', dendrogram=False, figsize=(8, 12), swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'Heatmap_group_GOBP_ERK5_CASCADE_pathway.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Heatmap_group_GOBP_ERK5_CASCADE_pathway.pdf'),bbox_inches='tight')