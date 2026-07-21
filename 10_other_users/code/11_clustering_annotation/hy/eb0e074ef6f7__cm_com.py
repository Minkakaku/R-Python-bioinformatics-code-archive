# cn10_proportion_boxplot_v6_points_larger.py
# 箱线图（无离群点）+ 样本点；显著性横线间距更大、图更大
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
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15.h5ad")
obs = adata.obs.copy()

if 'CNs_10' not in obs.columns or 'state' not in obs.columns:
    raise ValueError("Columns 'CNs_10' and 'state' must be present.")

if 'sample_new' in obs.columns:
    obs['sample_new_id'] = obs['sample_new'].astype(str)
else:
    obs['sample_new_id'] = obs.index.astype(str)

print(f"Data shape: {obs.shape}")

# ==================== 2. 计算比例 ====================
counts = obs.groupby('sample_new_id')['CNs_10'].value_counts().unstack(fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)

prop_long = proportions.reset_index().melt(id_vars='sample_new_id',
                                          var_name='CNs_10',
                                          value_name='proportion')
sample_new_group = obs[['sample_new_id', 'state']].drop_duplicates()
prop_long = prop_long.merge(sample_new_group, on='sample_new_id')

# ==================== 3. 参数 ====================
groups_order = ['N', 'A', 'C']
colors = {
    'N': '#d62728',
    'A': '#2ca02c',
    'C': '#ff7f0e'
}
cn_categories = sorted(proportions.columns.tolist())
print("CNs_10 categories:", cn_categories)

# ==================== 4. 箱线图（无离群点） + 叠加样本点 ====================
fig, ax = plt.subplots(figsize=(len(cn_categories)*0.8, 4))   # 更宽、更高

cn_positions = np.arange(len(cn_categories))
n_groups = len(groups_order)
total_width = 0.66          # 一组箱子总宽度
box_width = total_width / n_groups   # 每个箱子宽度

offsets = np.arange(n_groups) * box_width - total_width/2 + box_width/2
box_x_centers = {}

global_top_needed = 0

for i, cn in enumerate(cn_categories):
    data_groups = []
    positions_group = []
    for j, grp in enumerate(groups_order):
        vals = prop_long[(prop_long['CNs_10'] == cn) & 
                         (prop_long['state'] == grp)]['proportion'].dropna().values
        data_groups.append(vals)
        pos = cn_positions[i] + offsets[j]
        positions_group.append(pos)
    
    # 箱线图（无离群点）
    bp = ax.boxplot(data_groups, positions=positions_group, widths=box_width,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.5),
                    boxprops=dict(linewidth=1))
    
    for j, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[groups_order[j]])
        patch.set_alpha(0.7)
    
    # 叠加样本点（抖动）
    for j, grp in enumerate(groups_order):
        vals = prop_long[(prop_long['CNs_10'] == cn) & 
                         (prop_long['state'] == grp)]['proportion'].dropna().values
        x_center = positions_group[j]
        jitter = np.random.normal(0, box_width*0.12, size=len(vals))   # 抖动稍减
        x_jitter = np.full(len(vals), x_center) + jitter
        ax.scatter(x_jitter, vals, s=5, color=colors[grp], edgecolor='black', 
                   linewidth=0.5, alpha=0.8, zorder=3)
    
    box_x_centers[cn] = {grp: positions_group[j] for j, grp in enumerate(groups_order)}
    
    # ==================== 5. 显著性检验与标注（间距增大） ====================
    sub = prop_long[prop_long['CNs_10'] == cn]
    xs = box_x_centers[cn]
    
    pairs = list(combinations(groups_order, 2))
    sig_pairs = []
    for g1, g2 in pairs:
        v1 = sub[sub['state'] == g1]['proportion'].dropna().values
        v2 = sub[sub['state'] == g2]['proportion'].dropna().values
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
    
    if sig_pairs:
        sig_pairs.sort(key=lambda x: x[3])   # 按p值排序
        y_max = sub['proportion'].max()
        y_range = y_max * 0.2               # 基准间距增大（原0.1）
        base_y = y_max + 0.2       # 起始高度（原0.5）
        line_height = 0.15         # 每条线增加的垂直高度（原0.8）
        
        for i_pair, (g1, g2, sig, _) in enumerate(sig_pairs):
            x1, x2 = xs[g1], xs[g2]
            line_y = base_y + i_pair * line_height
            ax.plot([x1, x1, x2, x2],
                    [line_y - y_range*0.2, line_y, line_y, line_y - y_range*0.2],
                    lw=1.2, color='black')
            ax.text((x1+x2)/2, line_y + y_range*0.08, sig, ha='center', va='bottom',
                    fontsize=11, fontweight='bold')
        
        cn_top_needed = base_y + len(sig_pairs) * line_height + y_range * 0.5
        global_top_needed = max(global_top_needed, cn_top_needed)

# 统一y轴上限
current_ylim = ax.get_ylim()
if current_ylim[1] < global_top_needed:
    ax.set_ylim(top=global_top_needed)

ax.set_xticks(cn_positions)
ax.set_xticklabels(cn_categories, fontsize=12)
ax.set_xlabel('CNs_10 Category', fontsize=14)
ax.set_ylabel('Proportion of cells', fontsize=14)
ax.set_title('Proportion of CNs_10 per sample_new, by state', fontsize=16)

legend_patches = [Patch(facecolor=colors[grp], alpha=0.7, label=grp) for grp in groups_order]
ax.legend(handles=legend_patches, title='Group', bbox_to_anchor=(1.02, 1), loc='upper left')

ax.yaxis.grid(True, linestyle='--', alpha=0.7)
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig("CNs_10_proportion_boxplot_larger.pdf", bbox_inches='tight', dpi=150)
plt.savefig("CNs_10_proportion_boxplot_larger.png", bbox_inches='tight', dpi=150)
print("Figure saved successfully (larger spacing and size).")
