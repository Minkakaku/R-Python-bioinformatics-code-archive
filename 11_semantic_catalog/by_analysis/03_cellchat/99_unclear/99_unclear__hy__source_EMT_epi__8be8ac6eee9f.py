#!/usr/bin/env python3
"""
EMT_Epithelium 细胞通讯对比：in_d vs out_a (完整版，同起点，分终点)
- 散点图（通路水平 prob）
- 柱状图（差异通路）
- 网络图：中心 source，周围固定节点，in_d 和 out_a 起点重合，终点在靶标处分开
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyArrowPatch

# ---------- 配置 ----------
in_d_FILE = "in_d.net_lr.txt"    # 改为你的路径
out_a_FILE = "out_a.net_lr.txt"  # 改为你的路径
OUT_DIR = "source_EMT_EPI"
os.makedirs(OUT_DIR, exist_ok=True)

SOURCE_CELL = "EMT_Epithelium"

# ---------- 1. 数据读取与筛选 ----------
mpc_raw = pd.read_csv(in_d_FILE, sep='\t')
tnbc_raw = pd.read_csv(out_a_FILE, sep='\t')

# 筛选 source 为 EMT_Epithelium 的行
mpc_emt = mpc_raw[mpc_raw['source'] == SOURCE_CELL].copy()
tnbc_emt = tnbc_raw[tnbc_raw['source'] == SOURCE_CELL].copy()
print(f"in_d interactions from {SOURCE_CELL}: {len(mpc_emt)}")
print(f"out_a interactions from {SOURCE_CELL}: {len(tnbc_emt)}")

# ---------- 2. 通路标识（用于散点/柱状） ----------
def pathway_id(row):
    if pd.notna(row.get('interaction_name_2')) and row['interaction_name_2'] != '':
        return row['interaction_name_2']
    if pd.notna(row.get('interaction_name')) and row['interaction_name'] != '':
        return row['interaction_name']
    return f"{row['ligand']}_{row['receptor']}"

mpc_emt['pathway'] = mpc_emt.apply(pathway_id, axis=1)
tnbc_emt['pathway'] = tnbc_emt.apply(pathway_id, axis=1)

mpc_path = mpc_emt.groupby('pathway')['prob'].mean().rename('in_d')
tnbc_path = tnbc_emt.groupby('pathway')['prob'].mean().rename('out_a')
pathway_df = pd.concat([mpc_path, tnbc_path], axis=1).fillna(0)
pathway_df['diff'] = pathway_df['in_d'] - pathway_df['out_a']
pathway_df.sort_values('diff', ascending=False, inplace=True)
pathway_df.to_csv(os.path.join(OUT_DIR, "pathway_prob_comparison.csv"))

# ---------- 3. 散点图 ----------
plt.figure(figsize=(8,6))
ax = sns.scatterplot(data=pathway_df, x='out_a', y='in_d')
lims = [0, max(pathway_df[['in_d','out_a']].max().max() * 1.1, 0.1)]
plt.plot(lims, lims, 'k--', alpha=0.5)
plt.xlim(lims); plt.ylim(lims)
plt.xlabel('out_a prob'); plt.ylabel('in_d prob')
plt.title('Pathway activity: in_d vs out_a (EMT_Epithelium source)')
top10 = pathway_df.reindex(pathway_df['diff'].abs().sort_values(ascending=False).index[:10])
for pid, row in top10.iterrows():
    ax.annotate(pid, (row['out_a'], row['in_d']), fontsize=7,
                xytext=(5,5), textcoords='offset points')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "scatter_in_d_vs_out_a.png"), dpi=300)
plt.savefig(os.path.join(OUT_DIR, "scatter_in_d_vs_out_a.pdf"))
plt.close()

# ---------- 4. 柱状图（差异 top 20） ----------
top20 = pathway_df.head(20)
colors_bar = ['#E64B35' if v > 0 else '#4DBBD5' for v in top20['diff']]
plt.figure(figsize=(10,6))
sns.barplot(x=top20['diff'], y=top20.index, palette=colors_bar)
plt.axvline(0, color='black', lw=0.8)
plt.xlabel('in_d prob − out_a prob'); plt.ylabel('Pathway')
plt.title('Top 20 differential pathways')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "barplot_top20_diff.png"), dpi=300)
plt.savefig(os.path.join(OUT_DIR, "barplot_top20_diff.pdf"))
plt.close()

# ---------- 5. 网络图（同起点，分终点） ----------
# 汇总每个 target 的平均 prob
mpc_target = mpc_emt.groupby('target')['prob'].mean()
tnbc_target = tnbc_emt.groupby('target')['prob'].mean()
all_targets = sorted(set(mpc_target.index) | set(tnbc_target.index))
n_targets = len(all_targets)

# 准备数据
mpc_vals = np.array([mpc_target.get(t, 0) for t in all_targets])
tnbc_vals = np.array([tnbc_target.get(t, 0) for t in all_targets])

# --- 设置颜色映射 ---
unique_cell_types = sorted(set(all_targets))
color_palette = sns.color_palette("tab20", n_colors=len(unique_cell_types))
cell_type_color_map = {cell: color_palette[i] for i, cell in enumerate(unique_cell_types)}

# --- 创建画布 ---
fig, ax = plt.subplots(figsize=(12, 10))
ax.set_aspect('equal')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.axis('off')

# --- 中心 Source 圆圈 ---
center_radius = 0.15
source_circle = plt.Circle((0,0), center_radius, color='#FFD700', ec='black', lw=2, zorder=10)
ax.add_patch(source_circle)
ax.text(0, 0, 'EMT\nEpi', ha='center', va='center', fontsize=12, weight='bold', color='black')

# --- Target 节点分布 ---
target_radius = 1.0
node_radius = 0.12
angles = np.linspace(0, 2*np.pi, n_targets, endpoint=False)
target_pos = {}
max_prob = max(mpc_vals.max(), tnbc_vals.max(), 1e-6)

for i, (angle, tname) in enumerate(zip(angles, all_targets)):
    x = target_radius * np.cos(angle)
    y = target_radius * np.sin(angle)
    target_pos[tname] = (x, y)
    cell_color = cell_type_color_map[tname]
    circle = plt.Circle((x, y), node_radius, color=cell_color, ec='black', lw=1, zorder=9)
    ax.add_patch(circle)

# --- 绘制连线（重点修改：同起点，分终点） ---
end_offset = 0.3  # 终点分开的角度，越大分得越开

for i, tname in enumerate(all_targets):
    tx, ty = target_pos[tname]
    angle = angles[i]
    
    # 计算从中心到靶标的方向
    # 起点：统一使用中心圆边缘的同一个点 (角度 = angle)
    sx = center_radius * np.cos(angle)
    sy = center_radius * np.sin(angle)
    
    # in_d 连线 (红色) - 终点沿逆时针方向偏移 end_offset
    if mpc_vals[i] > 0:
        width = 1 + 5 * (mpc_vals[i] / max_prob)
        # in_d 终点角度：angle - end_offset
        end_angle_mpc = angle - end_offset
        ex = tx - node_radius * np.cos(end_angle_mpc)
        ey = ty - node_radius * np.sin(end_angle_mpc)
        
        arrow = FancyArrowPatch((sx, sy), (ex, ey),
                                 arrowstyle='-|>',
                                 mutation_scale=15,
                                 linewidth=width,
                                 color='#E64B35',
                                 alpha=0.7)
        ax.add_patch(arrow)
    
    # out_a 连线 (蓝色) - 终点沿顺时针方向偏移 end_offset
    if tnbc_vals[i] > 0:
        width = 1 + 5 * (tnbc_vals[i] / max_prob)
        # out_a 终点角度：angle + end_offset
        end_angle_tnbc = angle + end_offset
        ex = tx - node_radius * np.cos(end_angle_tnbc)
        ey = ty - node_radius * np.sin(end_angle_tnbc)
        
        arrow = FancyArrowPatch((sx, sy), (ex, ey),
                                 arrowstyle='-|>',
                                 mutation_scale=15,
                                 linewidth=width,
                                 color='#4DBBD5',
                                 alpha=0.7)
        ax.add_patch(arrow)

# --- 图例（右侧） ---
legend_handles = []
cell_legend = [Patch(color=cell_type_color_map[ct], label=ct) for ct in unique_cell_types]
legend_handles.extend(cell_legend)
line_mpc = Line2D([0], [0], color='#E64B35', lw=4, label='in_d')
line_tnbc = Line2D([0], [0], color='#4DBBD5', lw=4, label='out_a')
legend_handles.extend([line_mpc, line_tnbc])

ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, title="Cell Types & Condition", title_fontsize=12)

plt.title('EMT_Epithelium interactions with targets', fontsize=18, pad=20)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "network_fixed.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, "network_fixed.pdf"), bbox_inches='tight')
plt.close()

print(f"✅ 所有结果已保存至文件夹: {OUT_DIR}")
