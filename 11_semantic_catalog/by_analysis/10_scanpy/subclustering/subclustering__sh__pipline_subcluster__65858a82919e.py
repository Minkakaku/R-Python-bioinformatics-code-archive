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

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/08.Mission/mapping_raw/"

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/Bcells/Bcells.h5ad")
bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/DC/DC.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))


#### 免疫细胞降维
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/08.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/08.Mission/Immune/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['Macrophages','T cells','Neutrophils','Dendritic cells','B cells']), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pl.umap(adata, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_subcelltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_subcelltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.tl.rank_genes_groups(adata,groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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

adata.obs['group'] = (
    pd.Categorical(adata.obs['group'],
                   categories=groups,
                   ordered=True)
)
adata.obs['sample'] = (
    pd.Categorical(adata.obs['sample'],
                   categories=['T01', 'T02', 'T03','H21', 'H22', 'H23', 'H81', 'H82', 'H83', 'H241', 'H243', 'H246','H721', 'H722', 'H723', ],
                   ordered=True)
)
adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Immune.h5ad'))


from scipy.stats import pearsonr
import seaborn as sns

# 1. 计算每个 sample 的 sub_cell_type 比例
sample_proportions = pd.crosstab(adata.obs['sample'], adata.obs['sub_cell_type'], normalize='index') * 100
# 2. 获取每个 sample 对应的 group
sample_group = adata.obs[['sample', 'group']].drop_duplicates().set_index('sample')
# 对齐数据
common_samples = sample_proportions.index.intersection(sample_group.index)
sample_proportions = sample_proportions.loc[common_samples]
sample_group = sample_group.loc[common_samples]
# 3. 计算每个 group 的平均 sub_cell_type 比例
group_proportions = sample_proportions.groupby(sample_group['group']).mean()
# 4. 将 group 名称转换为数值（用于相关性计算）
# 方法：按 group 名称排序后分配数值
groups_sorted = sorted(group_proportions.index)
group_to_numeric = {g: i for i, g in enumerate(groups_sorted)}
groups_numeric = [group_to_numeric[g] for g in group_proportions.index]
# 5. 计算每个 sub_cell_type 比例与 group 数值的相关性
correlations = {}
p_values = {}
for ct in group_proportions.columns:
    proportions_values = group_proportions[ct].values
    corr, p = pearsonr(proportions_values, groups_numeric)
    correlations[ct] = corr
    p_values[ct] = p

# 6. 创建结果 DataFrame
corr_df = pd.DataFrame({
    'sub_cell_type': list(correlations.keys()),
    'correlation': list(correlations.values()),
    'p_value': list(p_values.values())
}).sort_values('correlation', ascending=False)
# 7. 保存结果
corr_df.to_csv(os.path.join(out_dir, 'subcelltype_group_correlation.csv'), index=False)

# 8. 绘制相关性条形图
plt.figure(figsize=(10, max(6, len(corr_df) * 0.3)))
colors = ['red' if c > 0 else 'blue' for c in corr_df['correlation']]
plt.barh(corr_df['sub_cell_type'][::-1], corr_df['correlation'][::-1], color=colors[::-1], alpha=0.7)
plt.axvline(x=0, color='black', linestyle='--', linewidth=1)
plt.xlabel('Correlation with Group')
plt.ylabel('Sub Cell Type')
plt.title('Correlation between Sub Cell Type Proportion and Group')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'subcelltype_group_correlation.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'subcelltype_group_correlation.pdf'), bbox_inches='tight')

# 9. 绘制 group 平均比例热图
plt.figure(figsize=(12, max(6, len(group_proportions) * 0.5)))
sns.heatmap(group_proportions.T, cmap='Blues', annot=True, fmt='.1f', 
            cbar_kws={'label': 'Proportion (%)'})
plt.title('Average Sub Cell Type Proportions by Group')
plt.xlabel('Group')
plt.ylabel('Sub Cell Type')
plt.tight_layout()
plt.savefig('group_avg_proportions_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

# 10. 可选：绘制相关性热图（紧凑版）
corr_matrix = corr_df.set_index('sub_cell_type')[['correlation']].T
plt.figure(figsize=(8, 3))
sns.heatmap(corr_matrix, annot=True, cmap='RdBu_r', center=0, 
            cbar_kws={'label': 'Correlation'})
plt.title('Sub Cell Type Correlation with Group')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'subcelltype_correlation_heatmap.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'subcelltype_correlation_heatmap.pdf'), bbox_inches='tight')