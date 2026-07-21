import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------
# 1. 读取并筛选细胞类型
# ------------------------------
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15.h5ad")

celltypes_to_keep = [
    "CXCL5+ Inf stem tumor Epi", "IGF2+ prolif stem tumor Epi",
    "Hypo Stem tumor Epi", "WNT+ ISC", "ESC-like Epi",
    "Prol Metab Epi", "EMT-like Inv Ent", "Inf Hyp Epi",
    "Early stress Epi", "Metab stress Epi", "Mat Ent",
    "MUC5B+ CGC", "MUC5AC+ MGC", "PCSK1+ AGC",
    "Hypo Metab Epi", "CD8+ T"
]
adata = adata[adata.obs['celltype'].isin(celltypes_to_keep)].copy()

assert 'CD274' in adata.var_names, "CD274 not in var_names"
assert 'PDCD1' in adata.var_names, "PDCD1 not in var_names"

def get_expr(adata, gene):
    x = adata[:, gene].X
    if hasattr(x, 'toarray'):
        return x.toarray().flatten()
    else:
        return np.array(x.flatten()).flatten()

cd274_expr = get_expr(adata, 'CD274')
pdcd1_expr = get_expr(adata, 'PDCD1')

# ------------------------------
# 2. 目标上皮细胞（非CD8+ T 且 CD274>0）
# ------------------------------
is_epi = adata.obs['celltype'] != "CD8+ T"
is_cd274_pos = cd274_expr > 0
epi_idx = np.where(is_epi & is_cd274_pos)[0]

# ------------------------------
# 3. 空间坐标
# ------------------------------
spatial_key = None
for key in ['spatial', 'X_spatial', 'spatial_coords']:
    if key in adata.obsm:
        spatial_key = key
        break
if spatial_key is None:
    raise KeyError("No spatial coordinates found in adata.obsm. Specify the correct key manually.")
coords = adata.obsm[spatial_key]

# ------------------------------
# 4. 所有CD8+ T细胞的空间索引
# ------------------------------
is_cd8t = adata.obs['celltype'] == "CD8+ T"
cd8t_idx = np.where(is_cd8t)[0]
cd8t_coords = coords[cd8t_idx]
cd8t_pdcd1 = pdcd1_expr[cd8t_idx]

tree = cKDTree(cd8t_coords)

# ------------------------------
# 5. 计算每个目标细胞的 PDCD1+ CD8+ T 邻居数量
# ------------------------------
epi_coords = coords[epi_idx]
radius = 200.0
neighbors_list = tree.query_ball_point(epi_coords, r=radius)

# 仅统计 PDCD1 > 0 的数量
pdcd1_pos_count = np.array([np.sum(cd8t_pdcd1[nbrs] > 0) for nbrs in neighbors_list])

# 分类：≥4 为 "High (>=4)"，<4 为 "Low (<4)"
threshold = 4
group_label = np.where(pdcd1_pos_count >= threshold, f'High (>={threshold})', f'Low (<{threshold})')

# 构建结果DataFrame
result_df = pd.DataFrame({
    'cell_index': epi_idx,
    'PDCD1_pos_CD8T_count': pdcd1_pos_count,
    'neighbor_group': group_label,
    'state': adata.obs['state'].iloc[epi_idx].values
})

print("State distribution in target cells:")
print(result_df['state'].value_counts())

# ------------------------------
# 6. 绘制堆积柱状图 + 卡方检验
# ------------------------------
# 准备列联表：行为state，列为邻居组
ct = pd.crosstab(result_df['state'], result_df['neighbor_group'])

# 卡方检验
chi2, p_chi, dof, expected = chi2_contingency(ct)
print(f"\nChi-square test: chi2={chi2:.3f}, p={p_chi:.2e}")

# 绘制堆积柱状图（百分比堆积）
ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100  # 每行百分比

colors = ['#66c2a5', '#fc8d62']  # Low 和 High 的颜色
ax = ct_pct.plot(kind='bar', stacked=True, color=colors, figsize=(6, 5))

# 在柱子上标注百分比数字
for i, (state, row) in enumerate(ct_pct.iterrows()):
    cumsum = 0
    for j, (group, val) in enumerate(row.items()):
        if val > 0:
            ax.text(i, cumsum + val/2, f'{val:.1f}%', ha='center', va='center', fontsize=9, color='black')
            cumsum += val

# 添加计数信息在柱子顶部
total_counts = ct.sum(axis=1)
for i, total in enumerate(total_counts):
    ax.text(i, 101, f'n={total}', ha='center', va='bottom', fontsize=9)

ax.set_ylabel('Percentage of CD274+ epithelial cells')
ax.set_xlabel('State')
ax.set_title('Proportion of cells with ≥4 PDCD1+ CD8+ T neighbors (200 µm)')
ax.legend(title='Neighbor group', loc='upper right')
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig('stacked_bar_state_CD274_PDCD1_CD8T.pdf', dpi=150)
plt.show()

# 输出具体数值
print("\nContingency table (counts):")
print(ct)
print("\nPercentage table:")
print(ct_pct.round(1))
