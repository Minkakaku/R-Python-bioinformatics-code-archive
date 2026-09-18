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
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/04.Mission/Main_mapping.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages', 'Mast cells', 'Neutrophil', 'Plasma cells', 'Smooth muscle cells', 'T cells'
adata_raw.obs.sample.cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure1/"

with rc_context({'figure.figsize': (7, 7), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':8, 'h':5})


#### 上皮功能热图、气泡图
Angiogenesis=['IHH','NOX1','EDIL3','ACVRL1','SOX17','ITGB3','THSD7A','SHH','ADTRP','ACE','NPR2','TSPAN12','ANPEP','DLL4','ADM2','ANG','SEMA5A','ITGA2B','COL18A1','AMOT','FGF18','SERPINE1','POFUT1','NRCAM','COL27A1','PRKCA','PDGFA','HOXB3','HIF3A','BMP4','CCN3','CCN2','RHOB','EPHA1','FOXC2','RAMP1','RNASE4','VAV2','HOXB13','EFNB2','EPHB4','FZD5','VPS25','HSPG2','SEMA3E','SLIT2','CEMIP2','EPHB2','EFNA1','NAXE']
Extracellular_matrix_organization=['TNFRSF11B','IHH','PRSS1','PRSS2','COL9A1','MMP7','NOX1','SMOC1','VTN','HPN','COMP','ADAMTS6','ADAMTS9','NID2','MIA','COL7A1','COL4A4','PLOD3','WDR72','WASHC1','COL4A3','COL18A1','APLP1','SERPINF2','MATN1','PTX3','SMPD3','PAPLN','COL6A1','SOX9','ADAMTS7','ADAMTS13','COL27A1','COL9A2','CREB3L1','COL28A1','EMID1','MMP1','FOXC2','LAMC1','COL9A3','ADAMTSL5','SCX','TNF','LAMB1','P3H4','COL16A1','KIF9','POMT1','LARGE1','GAS6','MMP11','MMP15','GAS2','ADAMTS1','ELF3','APP','SERPINH1','NOXO1','FLOT1','COL6A2']
adata=sc.read("/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Epithelial/Epithelial.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/"
# adata=adata.raw.to_adata()
sc.tl.score_genes(adata, gene_list=Angiogenesis, score_name='Angiogenesis')
sc.tl.score_genes(adata, gene_list=Extracellular_matrix_organization, score_name='Extracellular_matrix_organization')
sc.pl.umap(adata, color=['Angiogenesis','Extracellular_matrix_organization'], size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_genesets', fig_size={'w': 14, 'h': 6})


import seaborn as sns
from scipy.stats import mannwhitneyu
# 提取数据
df = adata.obs[['group', 'Angiogenesis', 'Extracellular_matrix_organization']].copy()
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
fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'wspace': 0.5})
for idx, feature in enumerate(['Angiogenesis', 'Extracellular_matrix_organization']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.violinplot(data=df, x='group', y=feature, order=groups, ax=ax, 
                   palette='RdBu_r', cut=0, inner=None)
    sns.boxplot(data=df, x='group', y=feature, order=groups, ax=ax, hue='group',
                width=0.2, boxprops={'alpha': 0.5}, showfliers=False)
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
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_genesets_sig', fig_size={'w': 10, 'h': 8})

sc.pl.violin(adata, keys=['Angiogenesis','Extracellular_matrix_organization'], groupby='group', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_genesets', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['Angiogenesis','Extracellular_matrix_organization'], groupby='group', standard_scale='var', cmap='RdBu_r', show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_genesets', fig_size={'w': 6, 'h': 6})

all_genes = {'Angiogenesis':['IHH','NOX1','EDIL3','ACVRL1','SOX17','ITGB3','THSD7A','SHH','ADTRP','ACE','NPR2','TSPAN12','ANPEP','DLL4','ADM2','ANG','SEMA5A','ITGA2B','COL18A1','AMOT','FGF18','SERPINE1','POFUT1','NRCAM','COL27A1','PRKCA','PDGFA','HOXB3','HIF3A','BMP4','CCN3','CCN2','RHOB','EPHA1','FOXC2','RAMP1','RNASE4','VAV2','HOXB13','EFNB2','EPHB4','FZD5','VPS25','HSPG2','SEMA3E','SLIT2','CEMIP2','EPHB2','EFNA1','NAXE'],
'Extracellular_matrix_organization':['TNFRSF11B','IHH','PRSS1','PRSS2','COL9A1','MMP7','NOX1','SMOC1','VTN','HPN','COMP','ADAMTS6','ADAMTS9','NID2','MIA','COL7A1','COL4A4','PLOD3','WDR72','WASHC1','COL4A3','COL18A1','APLP1','SERPINF2','MATN1','PTX3','SMPD3','PAPLN','COL6A1','SOX9','ADAMTS7','ADAMTS13','COL27A1','COL9A2','CREB3L1','COL28A1','EMID1','MMP1','FOXC2','LAMC1','COL9A3','ADAMTSL5','SCX','TNF','LAMB1','P3H4','COL16A1','KIF9','POMT1','LARGE1','GAS6','MMP11','MMP15','GAS2','ADAMTS1','ELF3','APP','SERPINH1','NOXO1','FLOT1','COL6A2']}
sc.pl.heatmap(
    adata,
    var_names=all_genes,
    groupby='group',
    swap_axes=True,
    show_gene_labels=True,
    # standard_scale='var',  # 按基因标准化
    cmap='RdYlBu_r',
    vmin=-0.2,    # 固定范围
    vmax=2
    # figsize=(10, len(available_genes) * 0.3),
    # show=True
)
plt.tight_layout()
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'heatmap_genesets', fig_size={'w': 6, 'h': 20})


#### 组间差异基因提琴图
df_group = adata.obs[['group']].copy()
available_genes = ['KRT20','SST','SMIM24','SPINK1','FXYD2','TRIM54','TNFRSF11B','HHLA2']
expr_matrix = adata.raw[:, available_genes].X.toarray()
df_expr = pd.DataFrame(expr_matrix, 
                           columns=available_genes, 
                           index=adata.obs_names)
# df = adata.obs[['group', 'KRT20','SST','SMIM24','SPINK1','FXYD2','TRIM54','TNFRSF11B','HHLA2']].copy()
df = pd.concat([df_group, df_expr], axis=1)
df = df.dropna()
# 获取分组顺序
groups = df['group'].unique().tolist()

# 绘制小提琴图
fig, axes = plt.subplots(2, 4, figsize=(32, 16), gridspec_kw={'wspace': 0.5})
axes = axes.flatten()
for idx, feature in enumerate(['KRT20','SST','SMIM24','SPINK1','FXYD2','TRIM54','TNFRSF11B','HHLA2']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.violinplot(data=df, x='group', y=feature, order=groups, ax=ax, 
                   palette='RdBu_r', cut=0, inner=None)
    sns.boxplot(data=df, x='group', y=feature, order=groups, ax=ax, hue='group',
                width=0.1, boxprops={'alpha': 0.5}, showfliers=False)
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
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_ADCvsSCC_gene_sig', fig_size={'w': 20, 'h': 16})


#### AAEs的对应结果
adata=sc.read("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/ADEC/AAEs.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/E/"
sc.tl.score_genes(adata, gene_list=Angiogenesis, score_name='Angiogenesis')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='leiden', title='leiden', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_leiden.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_leiden.pdf'))

anno_dict={
    '0':'AAEs',
    '1':'AAEs',
    '2':'AAEs',
    '3':'non-AAEs',
    '4':'non-AAEs',
    '5':'non-AAEs',
    '6':'non-AAEs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':3, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})
sc.tl.score_genes(adata, gene_list=Angiogenesis, score_name='Angiogenesis')
sc.pl.umap(adata, color='Angiogenesis', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_Angiogenesis', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='Angiogenesis', groupby='leiden', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_Angiogenesis', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['Angiogenesis'], groupby='leiden', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_Angiogenesis', fig_size={'w': 6, 'h': 6})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)
adata.write(os.path.join(out_dir, 'AAEs.h5ad'))

out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/F/"
sc.pl.heatmap(
    adata,
    var_names=Angiogenesis,
    groupby='cell_type',
    swap_axes=True,
    show_gene_labels=True,
    # standard_scale='var',  # 按基因标准化
    cmap='RdYlBu_r',
    vmin=-0.2,    # 固定范围
    vmax=2
    # figsize=(10, len(available_genes) * 0.3),
    # show=True
)
plt.tight_layout()

save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'heatmap_Angiogenesis', fig_size={'w': 6, 'h': 20})
sc.pl.dotplot(adata, var_names=Angiogenesis, groupby='cell_type', standard_scale='var', cmap='RdBu_r', show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_Angiogenesis', fig_size={'w': 20, 'h': 6})
sc.pl.matrixplot(adata, var_names=Angiogenesis, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Angiogenesis_heatmap_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Angiogenesis_heatmap_cell_type.pdf'),bbox_inches='tight')
plt.close()


#### AAEs间差异基因提琴图
df_group = adata.obs[['cell_type']].copy()
available_genes = ['GDF15','KRT20','CEACAM6','KRT17','MMP7','SPINK1']
expr_matrix = adata.raw[:, available_genes].X.toarray()
df_expr = pd.DataFrame(expr_matrix, 
                           columns=available_genes, 
                           index=adata.obs_names)
df = pd.concat([df_group, df_expr], axis=1)
df = df.dropna()
# 获取分组顺序
groups = df['cell_type'].unique().tolist()

# 绘制小提琴图
fig, axes = plt.subplots(2, 3, figsize=(32, 12), gridspec_kw={'wspace': 0.5})
axes = axes.flatten()
for idx, feature in enumerate(['GDF15','KRT20','CEACAM6','KRT17','MMP7','SPINK1']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.violinplot(data=df, x='cell_type', y=feature, order=groups, ax=ax, 
                   palette='RdBu_r', cut=0, inner=None)
    sns.boxplot(data=df, x='cell_type', y=feature, order=groups, ax=ax, hue='cell_type',
                width=0.1, boxprops={'alpha': 0.5}, showfliers=False)
    # 计算并添加显著性
    y_max = df[feature].max()
    y_min = df[feature].min()
    y_range = y_max - y_min
    # 组间比较
    comparisons = []
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            data1 = df[df['cell_type'] == groups[i]][feature]
            data2 = df[df['cell_type'] == groups[j]][feature]
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

out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/G/"
plt.tight_layout()
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_AAEs_gene_sig', fig_size={'w': 15, 'h': 16})


#### CAF的对应结果
adata=sc.read("/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Fibroblasts/Fibroblasts.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure3/"
sc.tl.score_genes(adata, gene_list=Angiogenesis, score_name='Angiogenesis')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})
# sc.tl.score_genes(adata, gene_list=Angiogenesis, score_name='Angiogenesis')
sc.pl.umap(adata, color='Angiogenesis', size=25, show=False,cmap='RdBu_r',vmax='p99')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_Angiogenesis', fig_size={'w': 7, 'h': 6})
# sc.pl.violin(adata, keys='Angiogenesis', groupby='group', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
# save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_Angiogenesis', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['Angiogenesis'], groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_Angiogenesis-cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['Angiogenesis'], groupby=['group'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_Angiogenesis', fig_size={'w': 6, 'h': 6})

sc.pl.dotplot(adata, var_names=Angiogenesis, groupby='cell_type', standard_scale='var', cmap='RdBu_r', show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_Angiogenesis', fig_size={'w': 20, 'h': 6})
sc.pl.matrixplot(adata, var_names=Angiogenesis, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'Angiogenesis_heatmap_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Angiogenesis_heatmap_cell_type.pdf'),bbox_inches='tight')
plt.close()
adata.write(os.path.join(out_dir, 'CAF.h5ad'))

sub_genes = {'Antiogenesis−related genes':['THY1','EDNRA','COL8A1','SRPX2','TMEM215','COL27A1','LEF1','VSTM4','GREM1','ECM1','LOXL2','HSPG2','PXDN','UNC5B','TGFBI','PLCD3','COL8A2','ANGPT2','PROK2','NOTCH3','CSPG4','SFRP2','FGF1','VASH1','NRXN1','NPR2'],'Extracellular matrix organization':['MMP11','COL1A1','COL1A2','COL3A1','SPARC','EMILIN1','COL5A1','COL5A2','COL10A1','CREB3L1','MMP7','COL12A1','ADAMTS16','SULF1','SULF2','COMP','BMP1','LOXL1','P3H4','MATN3','KAZALD1','KLK4','OLFML2B','LOXL4','OLFML2A','MMP28','COL11A1','SPINT2','SCARA3','MMP3','MMP10']}
sc.pl.dotplot(adata, var_names=sub_genes, groupby='group', standard_scale='var', cmap='Reds', show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_genesets', fig_size={'w': 20, 'h': 6})

#### CAF组间差异基因提琴图
df_group = adata.obs[['group']].copy()
available_genes = ['FP236383.3','CKB','SFRP4','ALDOA','KCNIP4','GAS6']
expr_matrix = adata.raw[:, available_genes].X.toarray()
df_expr = pd.DataFrame(expr_matrix, 
                           columns=available_genes, 
                           index=adata.obs_names)
df = pd.concat([df_group, df_expr], axis=1)
df = df.dropna()
# 获取分组顺序
groups = df['group'].unique().tolist()

# 绘制小提琴图
fig, axes = plt.subplots(2, 3, figsize=(32, 12), gridspec_kw={'wspace': 0.5})
axes = axes.flatten()
for idx, feature in enumerate(['FP236383.3','CKB','SFRP4','ALDOA','KCNIP4','GAS6']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.violinplot(data=df, x='group', y=feature, order=groups, ax=ax, 
                   palette='RdBu_r', cut=0, inner=None)
    sns.boxplot(data=df, x='group', y=feature, order=groups, ax=ax, hue='group',
                width=0.1, boxprops={'alpha': 0.5}, showfliers=False)
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

out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure3/"
plt.tight_layout()
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_CAF_gene_sig', fig_size={'w': 15, 'h': 16})

#### 做映射
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/04.Mission/Main_mapping.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/E/AAEs.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
adata_raw = adata_raw[adata_raw.obs['cell_type'].isin(['AAEs', 'Fibroblasts', 'Endothelial cells', 'non-AAEs']), :]
adata_raw.write(os.path.join(out_dir, 'Main_sub_cellchat.h5ad'))