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

# 对Tumor注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/05.Mission/Tumor/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/07.Mission/Tumor/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_raw.obs['raw_cell_type'] = adata_raw.obs['cell_type']
MP = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/07.Mission/Tumor/metaprograms.csv", index_col=0)  # 调整路径
max_cols = MP[['MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6']].idxmax(axis=1)
# 将结果赋值给 adata.obs['cell_type']
adata_raw.obs['cell_type'] = max_cols.values
# 如果需要查看前几行确认
print(adata_raw.obs['cell_type'].head())
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


#### 画评分与比例相关性
adata_raw=sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/mapping_raw/Main_mapping.h5ad")
MP_gene = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/07.Mission/cor/geneNMF_metaprograms_genes_weights.csv", index_col=0)  # 调整路径
mp_genes = {}
for mp, group in MP_gene.groupby('meta_program'):
    mp_genes[mp] = group['gene'].tolist()
sc.tl.score_genes(adata_raw, gene_list=mp_genes['MP1'], score_name='MP1')
sc.tl.score_genes(adata_raw, gene_list=mp_genes['MP2'], score_name='MP2')
sc.tl.score_genes(adata_raw, gene_list=mp_genes['MP3'], score_name='MP3')
sc.tl.score_genes(adata_raw, gene_list=mp_genes['MP4'], score_name='MP4')
sc.tl.score_genes(adata_raw, gene_list=mp_genes['MP5'], score_name='MP5')
sc.tl.score_genes(adata_raw, gene_list=mp_genes['MP6'], score_name='MP6')

# sc.pl.umap(adata, color=['SFTPC','SFTPB','ABCA3','PDPN','HOPX'], size=25, show=False,cmap='RdBu_r',sort_order=True, vmax='p99')
# plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.png'))
# plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.pdf'))
# sc.pl.dotplot(adata, var_names=['SFTPC','SFTPB','ABCA3','PDPN','HOPX'], groupby=['sample'],  cmap='RdBu_r',  show=False)
# save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_marker_sample', fig_size={'w': 6, 'h': 12})
# 快速计算
from scipy.stats import pearsonr
import seaborn as sns
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/07.Mission/cor/"
# 1. 细胞类型比例
prop = pd.crosstab(adata_raw.obs['sample'], adata_raw.obs['cell_type'], normalize='index')
# 2. MP评分（假设已有mp_genes_dict）
mp_scores = adata_raw.obs.groupby('sample')[['MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6']].mean()
# 3. 相关性矩阵
corr_matrix = pd.DataFrame(index=prop.columns, columns=mp_scores.columns)
p_matrix = pd.DataFrame(index=prop.columns, columns=mp_scores.columns)
corr_matrix = corr_matrix.apply(pd.to_numeric, errors='coerce')
p_matrix = p_matrix.apply(pd.to_numeric, errors='coerce')

for ct in prop.columns:
    for mp in mp_scores.columns:
        common = prop.index.intersection(mp_scores.index)
        corr, p = pearsonr(prop.loc[common, ct], mp_scores.loc[common, mp])
        corr_matrix.loc[ct, mp] = corr
        p_matrix.loc[ct, mp] = p

corr_matrix.to_csv(os.path.join(out_dir, 'cell_type_correlation_matrix.csv'))
p_matrix.to_csv(os.path.join(out_dir, 'cell_type_p_value_matrix.csv'))
def plot_correlation_dotplot(corr_matrix, p_matrix, save_path='correlation_dotplot.png',
                             figsize=(10, 8), min_corr=-1, max_corr=1):
    """
    绘制相关性点图，点大小表示显著性，颜色表示相关性
    """
    # 准备数据
    plot_data = corr_matrix.stack().reset_index()
    plot_data.columns = ['Cell_Type', 'MP_Program', 'Correlation']
    
    # 添加p值
    p_values = p_matrix.stack().reset_index()
    plot_data['P_value'] = p_values[0]
    
    # 添加-log10(p-value)用于点大小
    plot_data['neg_log10_p'] = -np.log10(plot_data['P_value'] + 1e-10)
    
    # 添加显著性标记
    plot_data['Significance'] = 'ns'
    plot_data.loc[plot_data['P_value'] < 0.05, 'Significance'] = '*'
    plot_data.loc[plot_data['P_value'] < 0.01, 'Significance'] = '**'
    plot_data.loc[plot_data['P_value'] < 0.001, 'Significance'] = '***'
    
    # 移除缺失值
    plot_data = plot_data.dropna()
    
    # 绘制dotplot
    plt.figure(figsize=figsize)
    
    # 使用散点图模拟dotplot
    scatter = plt.scatter(plot_data['MP_Program'], 
                         plot_data['Cell_Type'],
                         s=plot_data['neg_log10_p'] * 100,  # 点大小与显著性成正比
                         c=plot_data['Correlation'],
                         cmap='RdBu_r',
                         vmin=min_corr,
                         vmax=max_corr,
                         alpha=0.8,
                         edgecolors='black',
                         linewidth=0.5)
    
    # 添加颜色条
    cbar = plt.colorbar(scatter)
    cbar.set_label('Pearson Correlation', rotation=270, labelpad=15)
    
    # 添加显著性文本
    for idx, row in plot_data.iterrows():
        if row['Significance'] != 'ns':
            plt.text(row['MP_Program'], row['Cell_Type'], 
                    row['Significance'], 
                    ha='center', va='center', 
                    fontsize=10, fontweight='bold',
                    color='white' if abs(row['Correlation']) > 0.5 else 'black')
    
    plt.xlabel('MP Programs', fontsize=12)
    plt.ylabel('Cell Types', fontsize=12)
    plt.title('Correlation Dotplot: Cell Type Proportions vs MP Scores', 
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Dotplot已保存至: {save_path}")

# 使用函数
plot_correlation_dotplot(corr_matrix, p_matrix, save_path=os.path.join(out_dir, 'correlation_dotplot.png'))
plot_correlation_dotplot(corr_matrix, p_matrix, save_path=os.path.join(out_dir, 'correlation_dotplot.pdf'))

prop_sub = pd.crosstab(adata_raw.obs['sample'], adata_raw.obs['sub_cell_type'], normalize='index')
# 3. 相关性矩阵
corr_matrix = pd.DataFrame(index=prop.columns, columns=mp_scores.columns)
p_matrix = pd.DataFrame(index=prop.columns, columns=mp_scores.columns)
corr_matrix = corr_matrix.apply(pd.to_numeric, errors='coerce')
p_matrix = p_matrix.apply(pd.to_numeric, errors='coerce')

for ct in prop_sub.columns:
    for mp in mp_scores.columns:
        common = prop_sub.index.intersection(mp_scores.index)
        corr, p = pearsonr(prop_sub.loc[common, ct], mp_scores.loc[common, mp])
        corr_matrix.loc[ct, mp] = corr
        p_matrix.loc[ct, mp] = p

p_matrix_clean = p_matrix.dropna(axis=0, how='any')
corr_matrix_clean = p_matrix.dropna(axis=0, how='any')
p_matrix_clean.to_csv(os.path.join(out_dir, 'sub_cell_type_correlation_matrix.csv'))
corr_matrix_clean.to_csv(os.path.join(out_dir, 'sub_cell_type_p_value_matrix.csv'))
plot_correlation_dotplot(corr_matrix_clean, p_matrix_clean, save_path=os.path.join(out_dir, 'correlation_dotplot_sub.png'),figsize=(10, 16))
plot_correlation_dotplot(corr_matrix_clean, p_matrix_clean, save_path=os.path.join(out_dir, 'correlation_dotplot_sub.pdf'),figsize=(10, 16))





# 对CD8T注释细分 ####后面都没做
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD8T/CD8T.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/04.Mission/CD8T/all_cells_features.xls", index_col=0, sep='\t')  # 调整路径和索引列
umap_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/04.Mission/CD8T/umapreductions.csv", index_col=0)  # 调整路径
adata_raw.obs['cell_type'] = cell_type_df['cell_type']
umap_coords = umap_df.loc[adata_raw.obs.index, ['UMAP_1', 'UMAP_2']].values  # 根据实际列名调整
adata_raw.obsm['X_umap'] = umap_coords

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['LARC', 'N']
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

adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata_raw.write(os.path.join(out_dir, 'CD8T.h5ad'))

sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')



# 映射所有细分 #
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/mapping_raw/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/Monocytes/Macrophage.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD4T/CD4T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/CD8T.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')