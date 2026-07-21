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


#### CAF的对应结果
adata=sc.read("/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Fibroblasts/Fibroblasts.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/10.Mission/Figure3/"
adata.obs['sample'] = adata.obs['sample'].cat.add_categories('ADC3')
mask = adata.obs['sample'] == 'ADC5'
adata.obs.loc[mask, 'sample'] = 'ADC3'
adata.obs['sample'] = adata.obs['sample'].cat.remove_unused_categories()
adata.obs['sample'] = (
    pd.Categorical(adata.obs['sample'],
                   categories=['ADC1', 'ADC2', 'ADC3', 'ADC4', 'SCC1', 'SCC2', 'SCC3', 'SCC4'],
                   ordered=True)
)
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})

#### CAF组间差异基因提琴图
df_group = adata.obs[['group']].copy()
available_genes = ['ZEB1','SKIL','MEG3','HSPH1','IGF2','SEZ6L2']
expr_matrix = adata.raw[:, available_genes].X.toarray()
df_expr = pd.DataFrame(expr_matrix, 
                           columns=available_genes, 
                           index=adata.obs_names)
df = pd.concat([df_group, df_expr], axis=1)
df = df.dropna()
# 获取分组顺序
groups = df['group'].unique().tolist()

# 绘制小提琴图
fig, axes = plt.subplots(2, 3, figsize=(18, 30), gridspec_kw={'wspace': 0.5})
axes = axes.flatten()
for idx, feature in enumerate(['ZEB1','SKIL','MEG3','HSPH1','IGF2','SEZ6L2']):
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
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_CAF_gene_sig', fig_size={'w': 15, 'h': 16})



fig, axes = plt.subplots(2, 3, figsize=(18, 30), gridspec_kw={'wspace': 0.5})
axes = axes.flatten()
for idx, feature in enumerate(['ZEB1','SKIL','MEG3','HSPH1','IGF2','SEZ6L2']):
    ax = axes[idx]
    # 使用 seaborn 绘制
    sns.barplot(data=df, x='group', y=feature, order=groups, ax=ax, 
               color='black', linewidth=1.5,fill=False,
               saturation=1, alpha=1)
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
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'barplot_CAF_gene_sig', fig_size={'w': 15, 'h': 16})

