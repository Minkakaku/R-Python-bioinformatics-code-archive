#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CD8+ T cell spatial distance to boundary cells (per sample).
Exclude POST-R group, compare POST-NR, PRE-NR, PRE-R.
Output: sample-level boxplots + cell-level distribution.
Significance bars: horizontal lines only, no vertical connections.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import KDTree
from scipy.stats import kruskal
import scikit_posthocs as sp
from matplotlib import font_manager

# ===================== 0. 字体与输出 =====================
font_path = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'Arial'

output_dir = "CD8T_boundary"
os.makedirs(output_dir, exist_ok=True)

# ===================== 1. 读取数据 =====================
adata = sc.read_h5ad("adata_with_band_by_sample.h5ad")

required_cols = ["sample", "region", "celltype", "is_boundary", "group_final"]
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"Missing column: {col}")
if "spatial" not in adata.obsm:
    raise KeyError("Key 'spatial' not found in adata.obsm")

# ===================== 2. 目标分组（排除 POST-R） =====================
target_groups = ["POST-NR", "PRE-NR", "PRE-R"]

# ===================== 3. 逐样本计算距离 =====================
samples = adata.obs["sample"].unique()
print(f"Total samples: {len(samples)}")

sample_stats_list = []   # 样本级指标
all_cd8_cells = []       # 细胞级数据

for s in samples:
    mask_s = adata.obs["sample"] == s
    ad_s = adata[mask_s].copy()

    # 边界细胞
    boundary_mask = (
        ad_s.obs["is_boundary"]
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )
    boundary_coords = ad_s.obsm["spatial"][boundary_mask]
    if boundary_coords.shape[0] == 0:
        print(f"Sample {s}: no boundary cells, skipped.")
        continue

    # CD8+ T 细胞：region==normal, is_boundary==False, celltype 含 CD8
    cd8_mask = (
        (ad_s.obs["region"] == "normal") &
        (~boundary_mask) &
        (ad_s.obs["celltype"].str.contains("CD8", na=False))
    )
    cd8_coords = ad_s.obsm["spatial"][cd8_mask]
    if cd8_coords.shape[0] == 0:
        print(f"Sample {s}: no CD8+ T cells, skipped.")
        continue

    # 最近边界距离
    tree = KDTree(boundary_coords)
    dists, _ = tree.query(cd8_coords, k=1)

    # 样本的组别
    group = ad_s.obs["group_final"].iloc[0]
    if group not in target_groups:
        continue

    # 样本级指标
    sample_stats_list.append({
        "sample": s,
        "median_distance": np.median(dists),
        "mean_distance": np.mean(dists),
        "group": group
    })

    # 收集每个 CD8+ T 细胞
    for d in dists:
        all_cd8_cells.append({
            "sample": s,
            "distance": d,
            "group": group
        })

    print(f"Sample {s}: median={np.median(dists):.4f}, mean={np.mean(dists):.4f}, group={group}")

# 样本级 DataFrame
df_samples = pd.DataFrame(sample_stats_list)
df_samples["group"] = pd.Categorical(df_samples["group"], categories=target_groups, ordered=True)
df_samples.sort_values("group", inplace=True)
print(f"\nValid samples in target groups: {len(df_samples)}")

# 细胞级 DataFrame
df_cells = pd.DataFrame(all_cd8_cells)
df_cells["group"] = pd.Categorical(df_cells["group"], categories=target_groups, ordered=True)
print(f"Total CD8+ T cells: {len(df_cells)}")

# ===================== 4. 统计检验 =====================
def perform_stats(data, metric_col, label):
    """返回 (p_overall, dunn_df)"""
    if data[metric_col].nunique() <= 1 or data.shape[0] < 3:
        print(f"{label}: insufficient variation for test.")
        return np.nan, None
    groups_data = [data[data["group"] == g][metric_col].values for g in target_groups]
    h_stat, p_overall = kruskal(*groups_data)
    print(f"\n{label} - Kruskal-Wallis: H = {h_stat:.4f}, p = {p_overall:.4e}")
    dunn_df = sp.posthoc_dunn(data, val_col=metric_col, group_col="group", p_adjust="bonferroni")
    print(f"Dunn post-hoc (Bonferroni):\n{dunn_df.round(5)}")
    return p_overall, dunn_df

p_med, dunn_med = perform_stats(df_samples, "median_distance", "Sample-level median")
p_mean, dunn_mean = perform_stats(df_samples, "mean_distance", "Sample-level mean")
p_cell, dunn_cell = perform_stats(df_cells, "distance", "Cell-level distance")

# ===================== 5. 绘图 =====================
palette = {
    "POST-NR": "#E64B35",
    "PRE-NR":  "#00A087",
    "PRE-R":  "#3C5488"
}

def add_significance_bars(ax, dunn_df, data, metric_col, target_groups):
    """仅绘制水平线，无竖线连接"""
    if dunn_df is None:
        return
    y_max = data[metric_col].max()
    y_range = data[metric_col].max() - data[metric_col].min()
    if y_range == 0:
        y_range = y_max * 0.1 if y_max != 0 else 1
    offset = y_range * 0.08

    # 找出显著对 (p < 0.05)
    sig_pairs = []
    for i in range(len(target_groups)):
        for j in range(i+1, len(target_groups)):
            g1, g2 = target_groups[i], target_groups[j]
            p = dunn_df.loc[g1, g2]
            if p < 0.05:
                sig_pairs.append((g1, g2, p))

    for idx, (g1, g2, p) in enumerate(sig_pairs):
        x1 = target_groups.index(g1)
        x2 = target_groups.index(g2)
        y_h = y_max + offset * (idx + 1)
        # 只画水平线
        ax.plot([x1, x2], [y_h, y_h], lw=1.5, color='black')
        # 星号
        if p < 0.001:
            star = '***'
        elif p < 0.01:
            star = '**'
        else:
            star = '*'
        ax.text((x1 + x2) / 2, y_h + offset * 0.15, star, ha='center', va='bottom', fontsize=15, fontweight='bold')

    if sig_pairs:
        ax.set_ylim(top=(y_max + offset * (len(sig_pairs) + 1.5)))

def add_overall_p(ax, p_val):
    if not np.isnan(p_val):
        ax.text(0.5, 0.93, f"Kruskal-Wallis p = {p_val:.2e}",
                transform=ax.transAxes, ha='center', fontsize=12, style='italic')

# --- 5.1 样本级中位数箱线图（高度 6.5 英寸） ---
fig1, ax1 = plt.subplots(figsize=(6, 6.5))
sns.boxplot(data=df_samples, x="group", y="median_distance", order=target_groups,
            palette=palette, width=0.5, linewidth=1.5,
            flierprops=dict(marker='o', markerfacecolor='grey', markersize=5), ax=ax1)
sns.stripplot(data=df_samples, x="group", y="median_distance", order=target_groups,
              color="black", size=7, alpha=0.7, jitter=True, ax=ax1)
ax1.set_ylabel("Median distance to boundary (per sample)", fontweight="bold")
ax1.set_xlabel("Group", fontweight="bold")
ax1.set_title("Sample-level median distance\nCD8+ T cell proximity to boundary", fontweight="bold")
sns.despine(trim=True)
add_significance_bars(ax1, dunn_med, df_samples, "median_distance", target_groups)
add_overall_p(ax1, p_med)
plt.tight_layout()
fig1.savefig(os.path.join(output_dir, "CD8_boundary_median_sample_boxplot.pdf"))
fig1.savefig(os.path.join(output_dir, "CD8_boundary_median_sample_boxplot.png"))
plt.close(fig1)

# --- 5.2 样本级均值箱线图（高度 6.5 英寸） ---
fig2, ax2 = plt.subplots(figsize=(6, 6.5))
sns.boxplot(data=df_samples, x="group", y="mean_distance", order=target_groups,
            palette=palette, width=0.5, linewidth=1.5,
            flierprops=dict(marker='o', markerfacecolor='grey', markersize=5), ax=ax2)
sns.stripplot(data=df_samples, x="group", y="mean_distance", order=target_groups,
              color="black", size=7, alpha=0.7, jitter=True, ax=ax2)
ax2.set_ylabel("Mean distance to boundary (per sample)", fontweight="bold")
ax2.set_xlabel("Group", fontweight="bold")
ax2.set_title("Sample-level mean distance\nCD8+ T cell proximity to boundary", fontweight="bold")
sns.despine(trim=True)
add_significance_bars(ax2, dunn_mean, df_samples, "mean_distance", target_groups)
add_overall_p(ax2, p_mean)
plt.tight_layout()
fig2.savefig(os.path.join(output_dir, "CD8_boundary_mean_sample_boxplot.pdf"))
fig2.savefig(os.path.join(output_dir, "CD8_boundary_mean_sample_boxplot.png"))
plt.close(fig2)

# --- 5.3 细胞级距离分布（小提琴+窄箱线，宽 8 英寸，无散点） ---
fig3, ax3 = plt.subplots(figsize=(8, 6))
sns.violinplot(data=df_cells, x="group", y="distance", order=target_groups,
               palette=palette, inner=None, linewidth=1.2, ax=ax3)
sns.boxplot(data=df_cells, x="group", y="distance", order=target_groups,
            palette=palette, width=0.15, linewidth=1.2,
            boxprops=dict(alpha=0.7), whiskerprops=dict(alpha=0.7),
            capprops=dict(alpha=0.7), medianprops=dict(color='black', linewidth=1.5),
            ax=ax3)
ax3.set_ylabel("Distance to boundary (per cell)", fontweight="bold")
ax3.set_xlabel("Group", fontweight="bold")
ax3.set_title("CD8+ T cell distance distribution\n(normal region, per cell)", fontweight="bold")
sns.despine(trim=True)
add_significance_bars(ax3, dunn_cell, df_cells, "distance", target_groups)
add_overall_p(ax3, p_cell)
plt.tight_layout()
fig3.savefig(os.path.join(output_dir, "CD8_boundary_cell_distribution.pdf"))
fig3.savefig(os.path.join(output_dir, "CD8_boundary_cell_distribution.png"))
plt.close(fig3)

print("All figures saved.")

# ===================== 6. 保存数据 =====================
df_samples.to_csv(os.path.join(output_dir, "sample_distance_metrics.csv"), index=False)
df_cells.to_csv(os.path.join(output_dir, "cell_distance_data.csv"), index=False)
if dunn_med is not None:
    dunn_med.to_csv(os.path.join(output_dir, "dunn_sample_median.csv"), float_format="%.5e")
if dunn_mean is not None:
    dunn_mean.to_csv(os.path.join(output_dir, "dunn_sample_mean.csv"), float_format="%.5e")
if dunn_cell is not None:
    dunn_cell.to_csv(os.path.join(output_dir, "dunn_cell_distance.csv"), float_format="%.5e")

print("\nAll outputs saved in:", output_dir)
