#!/usr/bin/env python3
"""
Target = EMT_Epithelium 细胞通讯对比：MPC vs TNBC (反向网络图)
- 散点图（通路水平 prob）
- 柱状图（差异通路）
- 网络图：中心 target，外圈 source，箭头由 source 指向 target
- 同 source 的 MPC/TNBC 线起点岔开，终点指向中心同一点
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyArrowPatch

# ---------- 配置 ----------
MPC_FILE = "MPC.net_lr.txt"
TNBC_FILE = "TNBC.net_lr.txt"
OUT_DIR = "target_EMT_Epi_CN7"  # 新文件夹
os.makedirs(OUT_DIR, exist_ok=True)

TARGET_CELL = "EMT_Epithelium"

# ---------- 1. 数据读取与筛选 ----------
mpc_raw = pd.read_csv(MPC_FILE, sep='\t')
tnbc_raw = pd.read_csv(TNBC_FILE, sep='\t')

# 筛选 target 为 EMT_Epithelium 的行
mpc_emt = mpc_raw[mpc_raw['target'] == TARGET_CELL].copy()
tnbc_emt = tnbc_raw[tnbc_raw['target'] == TARGET_CELL].copy()
print(f"MPC interactions targeting {TARGET_CELL}: {len(mpc_emt)}")
print(f"TNBC interactions targeting {TARGET_CELL}: {len(tnbc_emt)}")

# ---------- 2. 通路标识（用于散点/柱状） ----------
def pathway_id(row):
    if pd.notna(row.get('interaction_name_2')) and row['interaction_name_2'] != '':
        return row['interaction_name_2']
    if pd.notna(row.get('interaction_name')) and row['interaction_name'] != '':
        return row['interaction_name']
    return f"{row['ligand']}_{row['receptor']}"

mpc_emt['pathway'] = mpc_emt.apply(pathway_id, axis=1)
tnbc_emt['pathway'] = tnbc_emt.apply(pathway_id, axis=1)

mpc_path = mpc_emt.groupby('pathway')['prob'].mean().rename('MPC')
tnbc_path = tnbc_emt.groupby('pathway')['prob'].mean().rename('TNBC')
pathway_df = pd.concat([mpc_path, tnbc_path], axis=1).fillna(0)
pathway_df['diff'] = pathway_df['MPC'] - pathway_df['TNBC']
pathway_df.sort_values('diff', ascending=False, inplace=True)
pathway_df.to_csv(os.path.join(OUT_DIR, "pathway_prob_comparison.csv"))

# ---------- 3. 散点图 ----------
plt.figure(figsize=(8,6))
ax = sns.scatterplot(data=pathway_df, x='TNBC', y='MPC')
lims = [0, max(pathway_df[['MPC','TNBC']].max().max() * 1.1, 0.1)]
plt.plot(lims, lims, 'k--', alpha=0.5)
plt.xlim(lims); plt.ylim(lims)
plt.xlabel('TNBC prob'); plt.ylabel('MPC prob')
plt.title('Pathway activity: MPC vs TNBC (Target = EMT_Epithelium)')
top10 = pathway_df.reindex(pathway_df['diff'].abs().sort_values(ascending=False).index[:10])
for pid, row in top10.iterrows():
    ax.annotate(pid, (row['TNBC'], row['MPC']), fontsize=7,
                xytext=(5,5), textcoords='offset points')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "scatter_MPC_vs_TNBC.png"), dpi=300)
plt.savefig(os.path.join(OUT_DIR, "scatter_MPC_vs_TNBC.pdf"))
plt.close()

# ---------- 4. 柱状图（差异 top 20） ----------
top20 = pathway_df.head(20)
colors_bar = ['#E64B35' if v > 0 else '#4DBBD5' for v in top20['diff']]
plt.figure(figsize=(10,6))
sns.barplot(x=top20['diff'], y=top20.index, palette=colors_bar)
plt.axvline(0, color='black', lw=0.8)
plt.xlabel('MPC prob − TNBC prob'); plt.ylabel('Pathway')
plt.title('Top 20 differential pathways (Target = EMT_Epithelium)')
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "barplot_top20_diff.png"), dpi=300)
plt.savefig(os.path.join(OUT_DIR, "barplot_top20_diff.pdf"))
plt.close()

# ---------- 5. 网络图（反向：中心 target，外圈 source，箭头指向中心） ----------
# 汇总每个 source 的平均 prob
mpc_source = mpc_emt.groupby('source')['prob'].mean()
tnbc_source = tnbc_emt.groupby('source')['prob'].mean()
all_sources = sorted(set(mpc_source.index) | set(tnbc_source.index))
n_sources = len(all_sources)

# 准备数据
mpc_vals = np.array([mpc_source.get(s, 0) for s in all_sources])
tnbc_vals = np.array([tnbc_source.get(s, 0) for s in all_sources])

# --- 设置颜色映射 ---
unique_cell_types = sorted(set(all_sources))
color_palette = sns.color_palette("tab20", n_colors=len(unique_cell_types))
cell_type_color_map = {cell: color_palette[i] for i, cell in enumerate(unique_cell_types)}

# --- 创建画布 ---
fig, ax = plt.subplots(figsize=(12, 10))
ax.set_aspect('equal')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.axis('off')

# --- 中心 Target 圆圈 ---
center_radius = 0.15
source_circle = plt.Circle((0,0), center_radius, color='#FFD700', ec='black', lw=2, zorder=10)
ax.add_patch(source_circle)
ax.text(0, 0, 'EMT\nEpi', ha='center', va='center', fontsize=12, weight='bold', color='black')

# --- Source 节点分布（外圈） ---
source_radius = 1.0
node_radius = 0.12
angles = np.linspace(0, 2*np.pi, n_sources, endpoint=False)
source_pos = {}
max_prob = max(mpc_vals.max(), tnbc_vals.max(), 1e-6)

for i, (angle, sname) in enumerate(zip(angles, all_sources)):
    x = source_radius * np.cos(angle)
    y = source_radius * np.sin(angle)
    source_pos[sname] = (x, y)
    cell_color = cell_type_color_map[sname]
    circle = plt.Circle((x, y), node_radius, color=cell_color, ec='black', lw=1, zorder=9)
    ax.add_patch(circle)

# --- 绘制连线（重点：同source分起点，箭头指向中心同一点） ---
start_offset = 0.12  # 同一个 source 的 MPC/TNBC 起点岔开的角度

for i, sname in enumerate(all_sources):
    sx, sy = source_pos[sname]
    angle = angles[i]
    
    # 终点：统一指向中心圆圈的边缘（同一点）
    ex = center_radius * np.cos(angle)
    ey = center_radius * np.sin(angle)
    
    # 计算切线方向向量（用于岔开起点）
    tangent_x = -np.sin(angle)
    tangent_y = np.cos(angle)
    
    # MPC 连线 (红色) - 起点沿切线方向偏移 -start_offset
    if mpc_vals[i] > 0:
        width = 1 + 5 * (mpc_vals[i] / max_prob)
        # MPC 起点偏移
        start_x_mpc = sx - start_offset * tangent_x
        start_y_mpc = sy - start_offset * tangent_y
        
        arrow = FancyArrowPatch((start_x_mpc, start_y_mpc), (ex, ey),
                                 arrowstyle='-|>',
                                 mutation_scale=15,
                                 linewidth=width,
                                 color='#E64B35',
                                 alpha=0.7,
                                 shrinkB=0)  # 终点刚好到中心圆边缘
        ax.add_patch(arrow)
    
    # TNBC 连线 (蓝色) - 起点沿切线方向偏移 +start_offset
    if tnbc_vals[i] > 0:
        width = 1 + 5 * (tnbc_vals[i] / max_prob)
        start_x_tnbc = sx + start_offset * tangent_x
        start_y_tnbc = sy + start_offset * tangent_y
        
        arrow = FancyArrowPatch((start_x_tnbc, start_y_tnbc), (ex, ey),
                                 arrowstyle='-|>',
                                 mutation_scale=15,
                                 linewidth=width,
                                 color='#4DBBD5',
                                 alpha=0.7,
                                 shrinkB=0)
        ax.add_patch(arrow)

# --- 图例（右侧） ---
legend_handles = []
# 细胞类型图例（外圈 source）
cell_legend = [Patch(color=cell_type_color_map[ct], label=ct) for ct in unique_cell_types]
legend_handles.extend(cell_legend)
# 线条图例
line_mpc = Line2D([0], [0], color='#E64B35', lw=4, label='MPC')
line_tnbc = Line2D([0], [0], color='#4DBBD5', lw=4, label='TNBC')
legend_handles.extend([line_mpc, line_tnbc])

ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, title="Sources & Condition", title_fontsize=12)

plt.title('Interactions targeting EMT_Epithelium (Source -> Target)', fontsize=18, pad=20)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "network_target_EMT_Epi.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, "network_target_EMT_Epi.pdf"), bbox_inches='tight')
plt.close()

print(f"✅ 所有结果已保存至文件夹: {OUT_DIR}")
