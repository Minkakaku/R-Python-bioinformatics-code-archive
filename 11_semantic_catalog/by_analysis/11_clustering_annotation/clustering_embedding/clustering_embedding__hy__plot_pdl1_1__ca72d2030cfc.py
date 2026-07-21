import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import kruskal, mannwhitneyu
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
# 5. 计算每个目标细胞的百分比
# ------------------------------
epi_coords = coords[epi_idx]
radius = 200.0
neighbors_list = tree.query_ball_point(epi_coords, r=radius)

total_cd8 = np.array([len(nbrs) for nbrs in neighbors_list])
pdcd1_pos_count = np.array([np.sum(cd8t_pdcd1[nbrs] > 0) for nbrs in neighbors_list])

with np.errstate(divide='ignore', invalid='ignore'):
    percentage = np.where(total_cd8 > 0,
                          pdcd1_pos_count / total_cd8 * 100,
                          np.nan)

result_df = pd.DataFrame({
    'cell_index': epi_idx,
    'total_CD8T': total_cd8,
    'PDCD1_pos_CD8T': pdcd1_pos_count,
    'percentage_PDCD1_CD8T': percentage,
    'state': adata.obs['state'].iloc[epi_idx].values
})

# 剔除周围没有任何CD8+ T细胞的点
valid_df = result_df.dropna(subset=['percentage_PDCD1_CD8T']).copy()
print(f"Total target cells: {len(result_df)}")
print(f"Cells with ≥1 CD8+ T neighbor: {len(valid_df)}")

# ------------------------------
# 6. 小提琴图 + 箱线图 + 不校正显著性
# ------------------------------
states = sorted(valid_df['state'].unique())
if len(states) < 2:
    print("Only one state group, no comparison possible.")
else:
    groups = [valid_df[valid_df['state'] == s]['percentage_PDCD1_CD8T'].values
              for s in states]

    h_stat, p_kruskal = kruskal(*groups)
    print(f"Kruskal-Wallis test: H={h_stat:.3f}, p={p_kruskal:.2e}")

    pairwise_pvals = []
    pairs = []
    for i in range(len(states)):
        for j in range(i+1, len(states)):
            stat, p = mannwhitneyu(groups[i], groups[j], alternative='two-sided')
            pairwise_pvals.append(p)
            pairs.append((states[i], states[j]))

    print("\nPairwise Mann-Whitney U (uncorrected):")
    for (s1, s2), p in zip(pairs, pairwise_pvals):
        sig_mark = '*' if p < 0.05 else 'ns'
        print(f"  {s1} vs {s2}: p={p:.2e} {sig_mark}")

    plt.figure(figsize=(8, 6))
    order = states

    ax = sns.violinplot(x='state', y='percentage_PDCD1_CD8T', data=valid_df,
                        hue='state', order=order, palette='Set2',
                        inner=None, legend=False, cut=0)
    sns.boxplot(x='state', y='percentage_PDCD1_CD8T', data=valid_df,
                hue='state', order=order, palette='Set2',
                width=0.15, showfliers=False, legend=False,
                boxprops={'zorder': 2})

    y_max = valid_df['percentage_PDCD1_CD8T'].max()
    y_min = valid_df['percentage_PDCD1_CD8T'].min()
    y_range = y_max - y_min if y_max != y_min else 1
    y_offset = y_range * 0.05
    annot_height = y_max + y_offset

    for (s1, s2), p in zip(pairs, pairwise_pvals):
        if p < 0.05:
            x1 = order.index(s1)
            x2 = order.index(s2)
            ax.plot([x1, x1, x2, x2],
                    [annot_height, annot_height + y_offset, annot_height + y_offset, annot_height],
                    lw=1.5, color='gray')
            if p < 0.001:
                sig_text = '***'
            elif p < 0.01:
                sig_text = '**'
            else:
                sig_text = '*'
            ax.text((x1 + x2) / 2, annot_height + y_offset, sig_text,
                    ha='center', va='bottom', fontsize=14)
            annot_height += y_offset * 3

    ax.set_ylim(top=annot_height + y_offset * 2)
    ax.set_ylabel('% PDCD1+ among CD8+ T cells within 200 µm')
    ax.set_xlabel('State')
    ax.set_title('PDCD1+ CD8+ T cell percentage near CD274+ epithelial cells')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('CD274_epi_PDCD1_CD8T_percentage_violin.pdf', dpi=150)
    plt.show()
