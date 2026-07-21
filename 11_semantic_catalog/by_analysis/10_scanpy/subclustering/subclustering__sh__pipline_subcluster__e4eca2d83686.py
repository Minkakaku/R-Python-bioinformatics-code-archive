import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from harmonypy import run_harmony
from matplotlib import font_manager as fm
import seaborn as sns
from scipy.stats import mannwhitneyu

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

# 1. IL16 TRAIL的基因表达在所有细胞类型中的表达水平主要显示IL16 TRAIL在非经典单核细胞的升高，气泡图小提琴图或者热图
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Mast cells', 'Monocytes', 'NK cells', 'Plasmacytoid dendritic cells', 'Platelets', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_1/"
os.makedirs(out_dir, exist_ok=True)
with rc_context({'figure.figsize': (7, 6), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color=['IL16','TNFSF10'], show=False, cmap='coolwarm',vmax='p99')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.pdf'))

sc.pl.dotplot(adata_raw, var_names=['IL16','TNFSF10'], groupby=['cell_type'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'bubble_cell_type_gene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'bubble_cell_type_gene.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.close()
sc.pl.dotplot(adata_raw, var_names=['IL16','TNFSF10'], groupby=['sub_cell_type'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'bubble_sub_cell_type_gene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'bubble_sub_cell_type_gene.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.close()

adata_sub = adata_raw[adata_raw.obs['sub_cell_type'].isin(['Non-classical Monocytes']), :]
sc.pl.dotplot(adata_sub, var_names=['IL16','TNFSF10'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'bubble_NCM_gene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'bubble_NCM_gene.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.close()

#### CAF组间差异基因提琴图
df_group = adata_sub.obs[['group']].copy()
available_genes = ['IL16','TNFSF10']
expr_matrix = adata_sub.raw[:, available_genes].X.toarray()
df_expr = pd.DataFrame(expr_matrix, 
                           columns=available_genes, 
                           index=adata_sub.obs_names)
df = pd.concat([df_group, df_expr], axis=1)
df = df.dropna()
# 获取分组顺序
groups = df['group'].unique().tolist()

# 绘制小提琴图
fig, axes = plt.subplots(1, 2, figsize=(9, 15), gridspec_kw={'wspace': 0.5})
axes = axes.flatten()
for idx, feature in enumerate(['IL16','TNFSF10']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.violinplot(data=df, x='group', y=feature, order=groups, ax=ax, 
               color='black', cut=0, linewidth=1.5,fill=False,inner="quart",common_norm=True,density_norm='width',
               saturation=1, alpha=1)
    sns.stripplot(data=df, x='group', y=feature, order=groups, ax=ax,
              palette='RdBu_r', alpha=0.3, size=5, edgecolor='none', jitter=True)
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
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_gene_sig', fig_size={'w': 10, 'h': 8})


out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_3/"
adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/NK cells/NKcells_anno.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub=adata[adata.obs['cell_type'].isin(['Mature NK c1', 'Mature NK c2']),]
sc.tl.rank_genes_groups(adata_sub,groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata_sub,species='human',output=out_dir)


adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_5/"
method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')


adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/2-result_1/"
adata.obs.cell_type=adata.obs.sub_cell_type
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

method = "sub_cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='sub_cell_type')


keep_cell_type = 'Non-classical Monocytes'
# 创建映射字典
cell_type_map = {ct: 'other' for ct in adata.obs['cell_type'].unique()}
cell_type_map[keep_cell_type] = keep_cell_type
# 应用映射
adata.obs['cell_type'] = adata.obs['cell_type'].map(cell_type_map)
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
# 设置颜色
categories = adata.obs['cell_type'].cat.categories.tolist()
new_colors = ['#D3D3D3'] * len(categories)
if keep_cell_type in categories:
    idx = categories.index(keep_cell_type)
    new_colors[idx] = '#f3e1eb'
adata.uns['cell_type_colors'] = new_colors
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

adata_sub=adata[adata.obs['cell_type'].isin(['Non-classical Monocytes']), :]
group_counts = adata_sub.obs['sample'].value_counts()
group_counts = group_counts.sort_index()
colors = ['#f3e1eb'] * 4 + ['gray'] * 4
plt.figure(figsize=(7, 7))
bars = plt.bar(group_counts.index, group_counts.values, color=colors, edgecolor='black', width=0.6)
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 5,
             f'{int(height)}', ha='center', va='bottom', fontsize=10)

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.ylabel('Total Non-classical Monocytes cell number', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'NCM_number_barplot.png'))
plt.savefig(os.path.join(out_dir, 'NCM_number_barplot.pdf'))