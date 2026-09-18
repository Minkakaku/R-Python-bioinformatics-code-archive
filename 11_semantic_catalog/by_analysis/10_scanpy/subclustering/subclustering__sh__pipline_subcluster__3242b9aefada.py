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
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr

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

# 合并基因表达结果
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Main_mapping/Main_mapping.h5ad")
bdata = sc.read("/mnt/NexData_file/scRNA/scVision_output/SC2509094832_ERV_all/03-Cell_filter/adata_counts.h5ad")
adata = adata.raw.to_adata()

new_genes = ['ERV_0001_MLT1L','ERV_0002_LTR103_Mam','ERV_0003_LTR103_Mam','ERV_0004_MER21B','ERV_0005_MER21B','ERV_0006_MLT1E1A','ERV_0007_MLT2B3','ERV_0008_LTR26','ERV_0009_LTR23','ERV_0010_THE1C-int','ERV_0011_LTR1D1','ERV_0012_MLT1A0','ERV_0013_MLT1G1','ERV_0014_MER39','ERV_0015_MLT1B','ERV_0016_MER72','ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0019_THE1B-int','ERV_0020_LTR8','ERV_0021_LTR8A','ERV_0022_MER4E1','ERV_0023_LTR8A','ERV_0024_MER50','ERV_0025_MLT1C','ERV_0026_MLT1H','ERV_0027_LTR8','ERV_0028_MSTB','ERV_0029_LTR16A','ERV_0030_LTR67B','ERV_0031_ERV3-16A3_I-int','ERV_0032_MER57B1','ERV_0033_MER41D','ERV_0034_MER21C','ERV_0035_MER21C','ERV_0036_MER65D','ERV_0037_MER48','ERV_0038_MER89','ERV_0039_LTR2','ERV_0040_MLT1D','ERV_0041_MamGypLTR3a','ERV_0042_MLT1H','ERV_0043_MLT1H','ERV_0044_LTR2','ERV_0045_MLT1B','ERV_0046_HERVL74-int','ERV_0047_HERVL74-int','ERV_0048_ERVL-B4-int','ERV_0049_MLT1D','ERV_0050_MSTA1','ERV_0051_MER21A','ERV_0052_MLT2A2','ERV_0053_MLT1C','ERV_0054_HUERS-P3-int','ERV_0055_MER101-int','ERV_0056_MamGypLTR2b','ERV_0057_MLT1B','ERV_0058_MLT1B-int','ERV_0059_MLT1H2','ERV_0060_LTR5_Hs','ERV_0061_MLT1H-int','ERV_0062_MLT1H','ERV_0063_LTR49-int','ERV_0064_MLT1M','ERV_0065_MLT1J','ERV_0066_LTR16C','ERV_0067_MLT1A0','ERV_0068_MLT1A0','ERV_0069_MLT1F1','ERV_0070_MER4CL34','ERV_0071_MLT1A0','ERV_0072_MLT1B','ERV_0073_MLT1H2','ERV_0074_MLT1B','ERV_0075_MER21C','ERV_0076_MER4D','ERV_0077_MER51B','ERV_0078_THE1B','ERV_0079_THE1B-int','ERV_0080_THE1B','ERV_0081_MLT1B','ERV_0082_MLT1A0','ERV_0083_LTR48','ERV_0084_MER4-int','ERV_0085_MER4-int','ERV_0086_MER4-int','ERV_0087_MER4B','ERV_0088_MLT2B3','ERV_0089_LTR57','ERV_0090_MER65C','ERV_0091_MER31B','ERV_0092_MER57F','ERV_0093_MLT1A1','ERV_0094_MamGypLTR2c','ERV_0095_LTR5_Hs','ERV_0096_HUERS-P2-int','ERV_0097_MSTB','ERV_0098_MLT1B','ERV_0099_THE1C','ERV_0100_MER51A','ERV_0101_MER51A','ERV_0102_PABL_B','ERV_0103_MLT1I','ERV_0104_LTR12C','ERV_0105_MLT1D','ERV_0106_MLT1A0','ERV_0107_LTR56','ERV_0108_MER21C','ERV_0109_ERV3-16A3_I-int','ERV_0110_MER41D','ERV_0111_LTR78','ERV_0112_LTR16A','ERV_0113_LTR8B','ERV_0114_LTR64','ERV_0115_MSTB','ERV_0116_MLT1B','ERV_0117_ERV3-16A3_I-int','ERV_0118_ERV3-16A3_I-int','ERV_0119_Harlequin-int','ERV_0120_LTR103b_Mam','ERV_0121_MER34B','ERV_0122_MLT1B','ERV_0123_LTR8A','ERV_0124_LTR39','ERV_0125_LTR62','ERV_0126_LTR39','ERV_0127_LTR62','ERV_0128_MLT1F2','ERV_0129_LTR2','ERV_0130_MLT1L','ERV_0131_MER72','ERV_0132_THE1B']
bdata_new = bdata[:, new_genes].copy()

# 2. 取交集 barcode
common_barcodes = adata.obs_names.intersection(bdata_new.obs_names)
# 3. 同时对两个对象切片 → 交集细胞 + 各自基因
adata_common = adata[common_barcodes, :].copy()
bdata_common = bdata_new[common_barcodes, :].copy()
# 4. 沿基因维度拼接（交集细胞，基因并集）
adata_ext = sc.concat([adata_common, bdata_common],
                      axis=1,
                      join='outer',
                      merge='first')
# 5. 验证
print(adata_ext.shape)                       # 行不变，列增加
print(set(new_genes) - set(adata_ext.var_names))  # 应为空
adata_ext.write("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/Main.h5ad")



adata=sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/ERV_all/"
geneset1=['ERV_0001_MLT1L','ERV_0002_LTR103_Mam','ERV_0003_LTR103_Mam','ERV_0004_MER21B','ERV_0005_MER21B','ERV_0006_MLT1E1A','ERV_0007_MLT2B3','ERV_0008_LTR26','ERV_0009_LTR23','ERV_0010_THE1C-int','ERV_0011_LTR1D1','ERV_0012_MLT1A0','ERV_0013_MLT1G1','ERV_0014_MER39','ERV_0015_MLT1B','ERV_0016_MER72','ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0019_THE1B-int','ERV_0020_LTR8','ERV_0021_LTR8A','ERV_0022_MER4E1','ERV_0023_LTR8A','ERV_0024_MER50','ERV_0025_MLT1C','ERV_0026_MLT1H','ERV_0027_LTR8','ERV_0028_MSTB','ERV_0029_LTR16A','ERV_0030_LTR67B','ERV_0031_ERV3-16A3_I-int','ERV_0032_MER57B1','ERV_0033_MER41D']
geneset2=['ERV_0034_MER21C','ERV_0035_MER21C','ERV_0036_MER65D','ERV_0037_MER48','ERV_0038_MER89','ERV_0039_LTR2','ERV_0040_MLT1D','ERV_0041_MamGypLTR3a','ERV_0042_MLT1H','ERV_0043_MLT1H','ERV_0044_LTR2','ERV_0045_MLT1B','ERV_0046_HERVL74-int','ERV_0047_HERVL74-int','ERV_0048_ERVL-B4-int','ERV_0049_MLT1D','ERV_0050_MSTA1','ERV_0051_MER21A','ERV_0052_MLT2A2','ERV_0053_MLT1C','ERV_0054_HUERS-P3-int','ERV_0055_MER101-int','ERV_0056_MamGypLTR2b','ERV_0057_MLT1B','ERV_0058_MLT1B-int','ERV_0059_MLT1H2','ERV_0060_LTR5_Hs','ERV_0061_MLT1H-int','ERV_0062_MLT1H','ERV_0063_LTR49-int','ERV_0064_MLT1M','ERV_0065_MLT1J','ERV_0066_LTR16C']
geneset3=['ERV_0067_MLT1A0','ERV_0068_MLT1A0','ERV_0069_MLT1F1','ERV_0070_MER4CL34','ERV_0071_MLT1A0','ERV_0072_MLT1B','ERV_0073_MLT1H2','ERV_0074_MLT1B','ERV_0075_MER21C','ERV_0076_MER4D','ERV_0077_MER51B','ERV_0078_THE1B','ERV_0079_THE1B-int','ERV_0080_THE1B','ERV_0081_MLT1B','ERV_0082_MLT1A0','ERV_0083_LTR48','ERV_0084_MER4-int','ERV_0085_MER4-int','ERV_0086_MER4-int','ERV_0087_MER4B','ERV_0088_MLT2B3','ERV_0089_LTR57','ERV_0090_MER65C','ERV_0091_MER31B','ERV_0092_MER57F','ERV_0093_MLT1A1','ERV_0094_MamGypLTR2c','ERV_0095_LTR5_Hs','ERV_0096_HUERS-P2-int','ERV_0097_MSTB','ERV_0098_MLT1B','ERV_0099_THE1C']
geneset4=['ERV_0100_MER51A','ERV_0101_MER51A','ERV_0102_PABL_B','ERV_0103_MLT1I','ERV_0104_LTR12C','ERV_0105_MLT1D','ERV_0106_MLT1A0','ERV_0107_LTR56','ERV_0108_MER21C','ERV_0109_ERV3-16A3_I-int','ERV_0110_MER41D','ERV_0111_LTR78','ERV_0112_LTR16A','ERV_0113_LTR8B','ERV_0114_LTR64','ERV_0115_MSTB','ERV_0116_MLT1B','ERV_0117_ERV3-16A3_I-int','ERV_0118_ERV3-16A3_I-int','ERV_0119_Harlequin-int','ERV_0120_LTR103b_Mam','ERV_0121_MER34B','ERV_0122_MLT1B','ERV_0123_LTR8A','ERV_0124_LTR39','ERV_0125_LTR62','ERV_0126_LTR39','ERV_0127_LTR62','ERV_0128_MLT1F2','ERV_0129_LTR2','ERV_0130_MLT1L','ERV_0131_MER72','ERV_0132_THE1B']
sc.pl.umap(adata, color=geneset1, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset1.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset1.pdf'))
sc.pl.dotplot(adata, var_names=geneset1, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset1, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1_cell_type.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=geneset2, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset2.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset2.pdf'))
sc.pl.dotplot(adata, var_names=geneset2, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset2, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2_cell_type.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=geneset3, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset3.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset3.pdf'))
sc.pl.dotplot(adata, var_names=geneset3, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset3, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3_cell_type.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=geneset4, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset4.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset4.pdf'))
sc.pl.dotplot(adata, var_names=geneset4, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset4, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4_cell_type.pdf'),bbox_inches='tight')


genes = geneset1 + geneset2 + geneset3 + geneset4
# 2. 提取并平均
expr = sc.get.obs_df(adata,
                     keys=genes) 
# 2. 把分组信息拼到同一 DataFrame
expr[['cell_type', 'group']] = adata.obs[['cell_type', 'group']]
expr.to_csv(os.path.join(out_dir, 'expression.csv'))
# 3. 按两因素求平均 & 转置
avg_df = expr.groupby(['cell_type', 'group']).mean().T
# 4. 导出
avg_df.to_csv(os.path.join(out_dir, 'average_expression.csv'))


adata=sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/result1/"
geneset=['ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0092_MER57F','ERV_0129_LTR2']
sc.pl.umap(adata, color=geneset, cmap='coolwarm', ncols=2)
plt.savefig(os.path.join(out_dir, 'umap_geneset1.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset1.pdf'))
sc.pl.dotplot(adata, var_names=geneset, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1_cell_type.pdf'),bbox_inches='tight')


out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/result2/"
geneset=['RIGI','IFIH1','MAVS','TBK1','IRF3','EIF2AK2','IKBKE','IRF7','MX1','ISG15']
sc.pl.umap(adata, color=geneset, cmap='coolwarm', ncols=4)
plt.savefig(os.path.join(out_dir, 'umap_geneset.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset.pdf'))
sc.pl.dotplot(adata, var_names=geneset, groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_cell_type.pdf'),bbox_inches='tight')


out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/result3/"
geneset=['ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0092_MER57F','ERV_0129_LTR2','RIGI','IFIH1','MAVS','TBK1','IRF3','EIF2AK2','IKBKE','IRF7','MX1','ISG15']
adata_sub=adata[adata.obs.group=='DM',]
sc.pl.dotplot(adata_sub, var_names=geneset, groupby=['sub_cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset.pdf'),bbox_inches='tight')

adata_sub=adata[adata.obs['cell_type'].isin(['Muscle cells','Macrophages','pDC','Endothelial cells','Fibroblasts']),]
sc.pl.dotplot(adata_sub, var_names=geneset, groupby='sub_cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_sub_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset_sub_cell_type.pdf'),bbox_inches='tight')


out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/result4/"
adata_sub=adata[adata.obs.group=='DM',]
if 'MX1' in adata_sub.var_names:
    gene_idx = list(adata_sub.var_names).index('MX1')
    mx1_expr = adata_sub.X[:, gene_idx]
    mx1_expr = mx1_expr.toarray().flatten() if scipy.sparse.issparse(mx1_expr) else mx1_expr.flatten()
    # 保存表达值
    adata_sub.obs['MX1'] = mx1_expr
    # 按阈值0分类
    adata_sub.obs['MX1_level'] = pd.Series(
        ['High' if x > 0 else 'Low' for x in mx1_expr], 
        index=adata_sub.obs.index
    ).astype('category')
else:
    print("警告: 未找到 MX1 基因")

results = []
for ct in adata_sub.obs['cell_type'].unique():
    sub = adata_sub[adata_sub.obs['cell_type'] == ct, :].copy()
    for g in geneset:
        high_expr = sub[sub.obs['MX1_level'] == 'High', g].X.toarray().flatten()
        low_expr  = sub[sub.obs['MX1_level'] == 'Low',  g].X.toarray().flatten()
        if len(high_expr) < 5 or len(low_expr) < 5:
            continue
        stat, p = mannwhitneyu(high_expr, low_expr, alternative='two-sided')
        results.append({
            'cell_type': ct,
            'gene': g,
            'High_mean': high_expr.mean(),
            'Low_mean': low_expr.mean(),
            'log2FC': np.log2(high_expr.mean() + 1e-9) - np.log2(low_expr.mean() + 1e-9),
            'p_val': p,
            'p_adj': min(p * len(geneset) * len(adata_sub.obs['cell_type'].unique()), 1.0)  # Bonferroni
        })

# 5. 结果表
res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(out_dir, 'wilcoxon_per_subtype_per_gene.csv'), index=False)
del adata_sub.obs['MX1']
# 6. 小提琴图（每组 ≥5 才画）
for ct in adata_sub.obs['cell_type'].unique():
    sc.pl.violin(adata_sub[adata_sub.obs['cell_type'] == ct, :],
                 keys=geneset,
                 groupby='MX1_level',
                 rotation=45,
                 show=False)
    # 用 matplotlib 设置标题
    plt.title(f'{ct}  MX1 High vs Low', fontsize=12)
    plt.savefig(os.path.join(out_dir, f'violin_{ct}.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, f'violin_{ct}.pdf'), bbox_inches='tight')
    plt.close()


out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/result5/"
adata_sub=adata[adata.obs.group=='DM',]
geneset=['ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0092_MER57F','ERV_0129_LTR2']
genes_present = [g for g in geneset if g in adata_sub.var_names]
# 2. 构造 “任意 >0” 指标
#    任意细胞 × 任意基因 >0 → True
any_pos = (
    adata_sub[:, genes_present].X.toarray() > 0      # 布尔矩阵
).any(axis=1)                                       # 按行求任意
# 3. 写回 obs
adata_sub.obs['ERV_gene_level'] = pd.Series(
    np.where(any_pos, 'High', 'Low'),
    index=adata_sub.obs.index,
    dtype='category'
)
print("High 比例:", any_pos.mean())
adata_sub.write(os.path.join('DM_ERV.h5ad'))

geneset=['RIGI','IFIH1','MAVS','TBK1','IRF3','EIF2AK2','IKBKE','IRF7','MX1','ISG15']
results = []
for ct in adata_sub.obs['cell_type'].unique():
    sub = adata_sub[adata_sub.obs['cell_type'] == ct, :].copy()
    for g in geneset:
        high_expr = sub[sub.obs['ERV_gene_level'] == 'High', g].X.toarray().flatten()
        low_expr  = sub[sub.obs['ERV_gene_level'] == 'Low',  g].X.toarray().flatten()
        if len(high_expr) < 5 or len(low_expr) < 5:
            continue
        stat, p = mannwhitneyu(high_expr, low_expr, alternative='two-sided')
        results.append({
            'cell_type': ct,
            'gene': g,
            'High_mean': high_expr.mean(),
            'Low_mean': low_expr.mean(),
            'log2FC': np.log2(high_expr.mean() + 1e-9) - np.log2(low_expr.mean() + 1e-9),
            'p_val': p,
            'p_adj': min(p * len(geneset) * len(adata_sub.obs['cell_type'].unique()), 1.0)  # Bonferroni
        })
# 5. 结果表
res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(out_dir, 'wilcoxon_per_subtype_per_gene.csv'), index=False)
# 6. 小提琴图（每组 ≥5 才画）
for ct in adata_sub.obs['cell_type'].unique():
    sc.pl.violin(adata_sub[adata_sub.obs['cell_type'] == ct, :],
                 keys=geneset,
                 groupby='ERV_gene_level',
                 rotation=45,
                 show=False)
    # 用 matplotlib 设置标题
    plt.title(f'{ct}  ERV_gene High vs Low', fontsize=12)
    plt.savefig(os.path.join(out_dir, f'violin_{ct}.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, f'violin_{ct}.pdf'), bbox_inches='tight')
    plt.close()


out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/03.Mission/result6/"
adata_sub=adata[adata.obs.group=='DM',]
erv_genes = ['ERV_0017_MLT1D','ERV_0018_LTR12','ERV_0092_MER57F','ERV_0129_LTR2']
target_genes = ['RIGI','IFIH1','MAVS','TBK1','IRF3','EIF2AK2','IKBKE','IRF7','MX1','ISG15']
all_res = []
for ct in adata.obs['cell_type'].unique():
    sub = adata[adata.obs['cell_type'] == ct, :].copy()
    if sub.n_obs < 5:  # 太少跳过
        continue
    # 提取表达矩阵
    expr = pd.DataFrame(
        sub.X.toarray() if hasattr(sub.X, 'toarray') else sub.X,
        index=sub.obs_names,
        columns=sub.var_names
    )
    # 逐 ERV 逐目标基因 Spearman
    for erv in erv_genes:
        for tg in geneset:
            rho, p = spearmanr(expr[erv], expr[tg])
            all_res.append({
                'cell_type': ct,
                'ERV': erv,
                'Target': tg,
                'rho': rho,
                'p_val': p
            })
# 3. 总表导出
pd.DataFrame(all_res).to_csv(os.path.join(out_dir, 'spearman_all_celltypes.csv'), index=False)
import seaborn as sns
rho_df = (pd.DataFrame(all_res)
          .pivot_table(index='ERV', columns=['cell_type', 'Target'], values='rho'))
# 2. 清理 NaN / inf
rho_clean = rho_df.replace([np.inf, -np.inf], np.nan)   # 先替换成 NaN
rho_clean = rho_clean.dropna(how='all')                 # 删除全 NaN 行
rho_clean = rho_clean.fillna(0)                         # 余下填 0（或填 p=1）
# 3. 聚类热图
sns.clustermap(rho_clean, cmap='coolwarm', figsize=(20, 8),row_cluster=False,col_cluster=False)
plt.savefig(os.path.join(out_dir, 'rho_clustermap.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'rho_clustermap.pdf'), bbox_inches='tight')