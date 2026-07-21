# cn10_proportion_boxplot_v6.py
# 修改：箱子按 PRE-R, PRE-NR, POST-R, POST-NR 排列，且箱子紧挨无间隙
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# ==================== 1. 读取数据 ====================
print("Reading h5ad file ...")
adata = sc.read_h5ad("adata_with_CNs_10.h5ad")
obs = adata.obs.copy()

if 'CNs_10' not in obs.columns or 'group_final' not in obs.columns:
    raise ValueError("Columns 'CNs_10' and 'group_final' must be present.")

if 'sample' in obs.columns:
    obs['sample_id'] = obs['sample'].astype(str)
else:
    obs['sample_id'] = obs.index.astype(str)

print(f"Data shape: {obs.shape}")

# ==================== 2. 计算比例 ====================
counts = obs.groupby('sample_id')['CNs_10'].value_counts().unstack(fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)

prop_long = proportions.reset_index().melt(id_vars='sample_id',
                                          var_name='CNs_10',
                                          value_name='proportion')
sample_group = obs[['sample_id', 'group_final']].drop_duplicates()
prop_long = prop_long.merge(sample_group, on='sample_id')

# ==================== 3. 参数（修改顺序） ====================
groups_order = ['PRE-R', 'PRE-NR', 'POST-R', 'POST-NR']
colors = {
    'PRE-R': '#d62728',
    'PRE-NR': '#2ca02c',
    'POST-R': '#ff7f0e',
    'POST-NR': '#1f77b4'
}
cn_categories = sorted(proportions.columns.tolist())
print("CNs_10 categories:", cn_categories)

# ==================== 4. 原生分组箱线图（紧挨着） ====================
fig, ax = plt.subplots(figsize=(len(cn_categories)*4, 7))

cn_positions = np.arange(len(cn_categories))
n_groups = len(groups_order)
total_width = 0.8          # 一组箱子总宽度
box_width = total_width / n_groups   # 每个箱子宽度 = 相邻箱子中心距离

# 计算偏移量：从 -total_width/2 + box_width/2 开始，步长为 box_width
offsets = np.arange(n_groups) * box_width - total_width/2 + box_width/2

box_x_centers = {}

for i, cn in enumerate(cn_categories):
    data_groups = []
    positions_group = []
    for j, grp in enumerate(groups_order):
        vals = prop_long[(prop_long['CNs_10'] == cn) & 
                         (prop_long['group_final'] == grp)]['proportion'].dropna().values
        data_groups.append(vals)
        pos = cn_positions[i] + offsets[j]
        positions_group.append(pos)
    
    bp = ax.boxplot(data_groups, positions=positions_group, widths=box_width,  # 宽度等于间距，无缝衔接
                    patch_artist=True, showfliers=True,
                    flierprops=dict(marker='o', markersize=3),
                    medianprops=dict(color='black', linewidth=1.5),
                    boxprops=dict(linewidth=1))
    
    for j, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[groups_order[j]])
        patch.set_alpha(0.7)
    
    box_x_centers[cn] = {grp: positions_group[j] for j, grp in enumerate(groups_order)}

# x 轴标签
ax.set_xticks(cn_positions)
ax.set_xticklabels(cn_categories, fontsize=12)
ax.set_xlabel('CNs_10 Category', fontsize=14)
ax.set_ylabel('Proportion of cells', fontsize=14)
ax.set_title('Proportion of CNs_10 per sample, by group_final', fontsize=16)

legend_patches = [Patch(facecolor=colors[grp], alpha=0.7, label=grp) for grp in groups_order]
ax.legend(handles=legend_patches, title='Group', bbox_to_anchor=(1.02, 1), loc='upper left')

ax.yaxis.grid(True, linestyle='--', alpha=0.7)
ax.set_axisbelow(True)

# ==================== 5. 显著性检验与标注 ====================
for cn in cn_categories:
    sub = prop_long[prop_long['CNs_10'] == cn]
    xs = box_x_centers[cn]
    
    pairs = list(combinations(groups_order, 2))
    sig_pairs = []
    for g1, g2 in pairs:
        v1 = sub[sub['group_final'] == g1]['proportion'].dropna().values
        v2 = sub[sub['group_final'] == g2]['proportion'].dropna().values
        if len(v1) < 3 or len(v2) < 3:
            continue
        _, p = mannwhitneyu(v1, v2, alternative='two-sided')
        p_adj = p * len(pairs)
        sig = ''
        if p_adj < 0.001:
            sig = '***'
        elif p_adj < 0.01:
            sig = '**'
        elif p_adj < 0.05:
            sig = '*'
        if sig:
            sig_pairs.append((g1, g2, sig, p_adj))
    
    if not sig_pairs:
        continue
    
    sig_pairs.sort(key=lambda x: x[3])
    y_max = sub['proportion'].max()
    y_range = y_max * 0.08
    for i, (g1, g2, sig, _) in enumerate(sig_pairs):
        x1, x2 = xs[g1], xs[g2]
        line_y = y_max + (i+1)*y_range*1.5
        ax.plot([x1, x1, x2, x2],
                [line_y-y_range*0.3, line_y, line_y, line_y-y_range*0.3],
                lw=1.2, color='black')
        ax.text((x1+x2)/2, line_y, sig, ha='center', va='bottom',
                fontsize=11, fontweight='bold')
    top_needed = line_y + y_range*2
    current_ylim = ax.get_ylim()
    if current_ylim[1] < top_needed:
        ax.set_ylim(top=top_needed)

plt.tight_layout()
plt.savefig("CNs_10_proportion_boxplot.pdf", bbox_inches='tight', dpi=150)
plt.savefig("CNs_10_proportion_boxplot.png", bbox_inches='tight', dpi=150)
print("Figure saved successfully.")
