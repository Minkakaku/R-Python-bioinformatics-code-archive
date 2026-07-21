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


adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Main_mapping/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/Main_mapping/"
adata_raw.obs['tmp_cell_type']=adata_raw.obs['cell_type']
anno_dict={
    'B cells':'1.B cells', 'Endothelial cells':'2.Endothelial cells', 'Fibroblasts':'3.Fibroblasts', 'Macrophages':'4.Macrophages',
       'Muscle cells':'5.Muscle cells', 'Pericytes':'6.Pericytes', 'Plasma cells':'7.Plasma cells', 'Satelite cells':'8.Satelite cells',
       'T cells':'9.T cells', 'pDC':'10.pDC', 'Myofibroblasts':'11.Myofibroblasts'
}
adata_raw.obs['cell_type'] = [anno_dict[clust] for clust in adata_raw.obs['tmp_cell_type']]
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell type', show=False,legend_loc ='on data',legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell type', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

adata_raw.obs['tmp_cell_type']=adata_raw.obs['sub_cell_type']
anno_dict={
    'B cells':'1.B cells','F13A1 Macrophages':'2.F13A1 Macrophages','Fibroblasts':'3.Fibroblasts','HLA-A Macrophages':'4.HLA-A Macrophages','IFI44L Capillaries':'5.IFI44L Capillaries','MYH1 Myotubes':'6.MYH1 Myotubes','MYH2 Myotubes':'7.MYH2 Myotubes','Myoblasts':'8.Myoblasts','Myocytes':'9.Myocytes','Myofibroblasts':'10.Myofibroblasts','PKD1L1 Capillaries':'11.PKD1L1 Capillaries','Pericytes':'12.Pericytes','Plasma cells':'13.Plasma cells','RGS5 Capillaries':'14RGS5 Capillaries','Satelite cells':'15Satelite cells','T cells':'16T cells','TNNT1 Myotubes':'1.7TNNT1 Myotubes','artery ECs':'18artery ECs','cDC':'19cDC','lymphatic ECs':'20lymphatic ECs','pDC':'21pDC','vein ECs':'22vein ECs'
}
adata_raw.obs['sub_cell_type'] = [anno_dict[clust] for clust in adata_raw.obs['tmp_cell_type']]
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Sub cell type', show=False,legend_loc ='on data',legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'umap_mapping_label.png'))
    plt.savefig(os.path.join(out_dir, 'umap_mapping_label.pdf'))
with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Sub cell type', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'umap_mapping.png'))
    plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'))

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))


adata=sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/Main.h5ad")
bdata=sc.read("/mnt/NexData_file/scRNA/scVision_output/SC2509094832_ERV/result.h5ad")
bdata=bdata.raw.to_adata()
# 1. 指定基因列表
genes = ['LTR/ERVK', 'LTR/ERV1', 'LTR/ERVL', 'LTR/ERVL-MaLR']
# 2. 确保基因在 bdata 中
# genes = [g for g in genes if g in bdata.var_names]
# 3. 提取 bdata 表达矩阵（转为 dense）
b_expr = bdata[:, genes].X.toarray() if hasattr(bdata.X, 'toarray') else bdata[:, genes].X
b_expr = pd.DataFrame(b_expr, index=bdata.obs_names, columns=genes)
# 4. 按 adata 的 barcode 顺序对齐（缺失填 0）
b_expr_aligned = b_expr.reindex(adata.obs_names, fill_value=0)
# 方案 C：逐个基因存入 obs（方便画图）
for g in genes:
    adata.obs[g] = b_expr_aligned[g].values


old_cols = ['LTR/ERVK', 'LTR/ERV1', 'LTR/ERVL', 'LTR/ERVL-MaLR']
new_cols = {old: old.replace('LTR/', '') for old in old_cols if old in adata.obs.columns}
adata.obs.rename(columns=new_cols, inplace=True)
print(adata.obs.columns)
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/ERV_all/"
with rc_context({'figure.figsize': (8, 8), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color=['ERVK', 'ERV1', 'ERVL', 'ERVL-MaLR'], show=False, cmap='coolwarm')
    plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
    plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
with rc_context({'figure.figsize': (8, 8), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color=['ERVK', 'ERV1', 'ERVL', 'ERVL-MaLR'], show=False, cmap='coolwarm',vmax='p99')
    plt.savefig(os.path.join(out_dir, 'umap_genes_p99.png'))
    plt.savefig(os.path.join(out_dir, 'umap_genes_p99.pdf'))

sc.pl.dotplot(adata, var_names=['ERVK', 'ERV1', 'ERVL', 'ERVL-MaLR'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['ERVK', 'ERV1', 'ERVL', 'ERVL-MaLR'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['ERVK', 'ERV1', 'ERVL', 'ERVL-MaLR'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_scale.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['ERVK', 'ERV1', 'ERVL', 'ERVL-MaLR'], groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_cell_type_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_cell_type_scale.pdf'),bbox_inches='tight')


# 'Endothelial cells','Macrophages','Muscle cells','Satelite cells','pDC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/cor_1/"
adata_sub=adata[adata.obs.cell_type=='pDC',]
expr = pd.DataFrame(
    adata_sub.X.toarray() if hasattr(adata_sub.X, 'toarray') else adata_sub.X,
    index=adata_sub.obs_names,
    columns=adata_sub.var_names
)
# erv_genes = ['ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0092_MER57F','ERV_0129_LTR2']
target_genes = ['EIF2AK2','TBK1','IRF7','MX1']#'IKBKE'
all_res = []
for g1, g2 in combinations(target_genes, 2):
    rho, p = spearmanr(expr[g1], expr[g2])
    all_res.append({
        'gene1': g1,
        'gene2': g2,
        'rho': rho,
        'p_val': p,
        'log_p': -np.log10(p + 1e-300)
    })
# 3. 总表导出
pd.DataFrame(all_res).to_csv(os.path.join(out_dir, 'spearman_pDC.csv'), index=False)

res_df = pd.read_csv(os.path.join(out_dir, 'spearman_pDC.csv'))# 获取唯一基因列表
# 获取唯一基因列表
genes = sorted(set(res_df['gene1']) | set(res_df['gene2']))
n = len(genes)
gene_to_idx = {gene: i for i, gene in enumerate(genes)}

# 创建对称矩阵
rho = np.eye(n)
for _, row in res_df.iterrows():
    i, j = gene_to_idx[row['gene1']], gene_to_idx[row['gene2']]
    rho[i, j] = rho[j, i] = row['rho']

# 设置参数
vmin, vmax = -0.5, 0.5
fig, ax = plt.subplots(figsize=(n*0.8+2, n*0.8+1))

# 添加细网格线
ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
ax.grid(True, which='minor', linestyle='-', linewidth=0.5, color='lightgray')
ax.tick_params(which='minor', bottom=False, left=False)

# 绘制下三角
for i in range(n):
    for j in range(i):
        rho_val = np.clip(rho[i, j], vmin, vmax)
        radius = (abs(rho_val) * 800) ** 0.5 / 50
        color = plt.cm.coolwarm((rho_val - vmin) / (vmax - vmin))
        
        ax.add_patch(Circle((j, i), radius, facecolor=color, edgecolor='none', alpha=0.8))
        ax.text(j, i-0.25, f'{rho[i, j]:.2f}', ha='center', va='center', fontsize=8, color='black')

# 坐标轴设置
ax.set(xlim=(-0.5, n-0.5), ylim=(-0.5, n-0.5), aspect='equal')
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(genes, fontsize=10)
ax.invert_yaxis()

# 去掉边框
for spine in ax.spines.values():
    spine.set_visible(False)

# 颜色条
plt.colorbar(plt.cm.ScalarMappable(cmap='coolwarm', norm=mcolors.Normalize(vmin, vmax)), 
             ax=ax, shrink=0.6, pad=0.02, label=f'Spearman Rho')

plt.title('Gene Correlation Matrix', fontsize=12, pad=20)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'correlation_bubble_pDC.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'correlation_bubble_pDC.pdf'), bbox_inches='tight')
plt.show()


# 'Endothelial cells','Macrophages','Muscle cells','Satelite cells','pDC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/04.Mission/cor_2/"
adata_sub=adata[adata.obs.cell_type=='Satelite cells',]
expr = pd.DataFrame(
    adata_sub.X.toarray() if hasattr(adata_sub.X, 'toarray') else adata_sub.X,
    index=adata_sub.obs_names,
    columns=adata_sub.var_names
)
# erv_genes = ['ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0092_MER57F','ERV_0129_LTR2']
target_genes = ['EIF2AK2','IKBKE','IRF7','MX1']#'IKBKE'
all_res = []
for g1, g2 in combinations(target_genes, 2):
    rho, p = spearmanr(expr[g1], expr[g2])
    all_res.append({
        'gene1': g1,
        'gene2': g2,
        'rho': rho,
        'p_val': p,
        'log_p': -np.log10(p + 1e-300)
    })
# 3. 总表导出
pd.DataFrame(all_res).to_csv(os.path.join(out_dir, 'spearman_Satelite.csv'), index=False)

res_df = pd.read_csv(os.path.join(out_dir, 'spearman_Satelite.csv'))# 获取唯一基因列表
# 获取唯一基因列表
genes = sorted(set(res_df['gene1']) | set(res_df['gene2']))
n = len(genes)
gene_to_idx = {gene: i for i, gene in enumerate(genes)}

# 创建对称矩阵
rho = np.eye(n)
for _, row in res_df.iterrows():
    i, j = gene_to_idx[row['gene1']], gene_to_idx[row['gene2']]
    rho[i, j] = rho[j, i] = row['rho']

# 设置参数
vmin, vmax = -0.5, 0.5
fig, ax = plt.subplots(figsize=(n*0.8+2, n*0.8+1))

# 添加细网格线
ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
ax.grid(True, which='minor', linestyle='-', linewidth=0.5, color='lightgray')
ax.tick_params(which='minor', bottom=False, left=False)

# 绘制下三角
for i in range(n):
    for j in range(i):
        rho_val = np.clip(rho[i, j], vmin, vmax)
        radius = (abs(rho_val) * 800) ** 0.5 / 50
        color = plt.cm.coolwarm((rho_val - vmin) / (vmax - vmin))
        
        ax.add_patch(Circle((j, i), radius, facecolor=color, edgecolor='none', alpha=0.8))
        ax.text(j, i-0.25, f'{rho[i, j]:.2f}', ha='center', va='center', fontsize=8, color='black')

# 坐标轴设置
ax.set(xlim=(-0.5, n-0.5), ylim=(-0.5, n-0.5), aspect='equal')
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(genes, fontsize=10)
ax.invert_yaxis()

# 去掉边框
for spine in ax.spines.values():
    spine.set_visible(False)

# 颜色条
plt.colorbar(plt.cm.ScalarMappable(cmap='coolwarm', norm=mcolors.Normalize(vmin, vmax)), 
             ax=ax, shrink=0.6, pad=0.02, label=f'Spearman Rho')

plt.title('Gene Correlation Matrix', fontsize=12, pad=20)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'correlation_bubble_Satelite.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'correlation_bubble_Satelite.pdf'), bbox_inches='tight')
plt.show()