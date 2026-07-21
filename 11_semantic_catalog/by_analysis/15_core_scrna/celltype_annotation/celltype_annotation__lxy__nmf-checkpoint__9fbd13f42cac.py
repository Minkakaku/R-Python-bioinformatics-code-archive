#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NMF 分析细胞亚群频率 — 精简版 + 联合箱线图（FDR 校正） + top10 气泡图
"""

import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from scipy.cluster.hierarchy import linkage, leaves_list
import scipy.stats as stats
from itertools import combinations
from statsmodels.stats.multitest import multipletests  # 新增 FDR
import os
import warnings

# ==================== 字体与全局设置 ====================
import matplotlib
from matplotlib import font_manager
try:
    font_manager.fontManager.addfont(
        '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
    )
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = ['Arial']
except:
    pass

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
warnings.filterwarnings("ignore")

# ==================== 参数 ====================
H5AD_PATH = "merged_annotated_final_group_reamin_new_with_CNs.h5ad"
N_RANK_LIST = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
GROUP_COL = "group"
CELL_TYPE_COL = "singler_best"
SAMPLE_COL = "sample"
OUTPUT_DIR = "nmf_outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==================== 1. 读取与预处理 ====================
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs.copy()
for col in [SAMPLE_COL, CELL_TYPE_COL, GROUP_COL]:
    if col not in obs.columns:
        raise KeyError(f"obs 中缺少必需列: {col}")
obs[SAMPLE_COL] = obs[SAMPLE_COL].astype(str)
obs[CELL_TYPE_COL] = obs[CELL_TYPE_COL].astype(str)
obs[GROUP_COL] = obs[GROUP_COL].astype(str)

freq_raw = (
    obs.groupby([SAMPLE_COL, CELL_TYPE_COL])
    .size()
    .reset_index(name="n")
)
total_per_sample = freq_raw.groupby(SAMPLE_COL)["n"].sum().reset_index(name="n_total")
freq = freq_raw.merge(total_per_sample, on=SAMPLE_COL)
freq["freq"] = freq["n"] / freq["n_total"]
freq_matrix = freq.pivot(index=SAMPLE_COL, columns=CELL_TYPE_COL, values="freq")
freq_matrix = freq_matrix.fillna(0)
non_zero_cols = freq_matrix.columns[freq_matrix.sum(axis=0) > 0]
freq_matrix = freq_matrix[non_zero_cols]

group_info = obs[[SAMPLE_COL, GROUP_COL]].drop_duplicates().set_index(SAMPLE_COL)
group_info = group_info.loc[freq_matrix.index].sort_values(GROUP_COL)
freq_matrix = freq_matrix.loc[group_info.index]

unique_groups = sorted(group_info[GROUP_COL].unique())
print(f"实际 group 列的唯一值: {unique_groups}")   # ⬅️ 关键：请检查输出
palette = sns.color_palette("Set1", len(unique_groups))
group_colors = dict(zip(unique_groups, palette))

# ==================== 2. 单 rank 分析函数 ====================
def run_nmf_rank(n_rank, freq_matrix, group_info, group_colors, unique_groups, rank_dir):
    model = NMF(n_components=n_rank, init='nndsvda', random_state=123, max_iter=2000)
    W = model.fit_transform(freq_matrix.values)
    H = model.components_
    rec_err = model.reconstruction_err_

    cell_types = freq_matrix.columns.tolist()
    samples = freq_matrix.index.tolist()
    mode_names = [f"NMF{i+1}" for i in range(n_rank)]

    # ---------- 2.1 基矩阵热图 ----------
    W_basis = H.T
    W_basis_scaled = (W_basis - W_basis.mean(axis=1, keepdims=True)) / (W_basis.std(axis=1, keepdims=True) + 1e-10)
    dominant_mode_per_cell = np.argmax(H, axis=0)
    order_by_mode = np.argsort(dominant_mode_per_cell)
    cell_types_ordered = [cell_types[i] for i in order_by_mode]
    data_ordered = W_basis_scaled[order_by_mode, :]
    mode_starts = []
    current_mode = -1
    for i, mode_idx in enumerate(dominant_mode_per_cell[order_by_mode]):
        if mode_idx != current_mode:
            current_mode = mode_idx
            mode_starts.append(i)
    mode_starts.append(len(cell_types_ordered))
    fig, ax = plt.subplots(figsize=(6, max(6, len(cell_types_ordered)*0.25)))
    sns.heatmap(data_ordered, xticklabels=mode_names, yticklabels=cell_types_ordered,
                cmap=sns.diverging_palette(250, 10, as_cmap=True), center=0,
                cbar_kws={"label": "Z-score"}, ax=ax)
    ax.set_xlabel("NMF mode"); ax.set_ylabel("Cell type")
    ax.tick_params(axis='x', rotation=0, labelsize=8); ax.tick_params(axis='y', rotation=0, labelsize=6)
    for i in range(len(mode_starts)-1):
        start, end = mode_starts[i], mode_starts[i+1]
        if end-start > 0:
            rect = Rectangle((0, start), data_ordered.shape[1], end-start, linewidth=2, edgecolor='black', facecolor='none')
            ax.add_patch(rect)
    for i in range(len(mode_starts)-1):
        mid = (mode_starts[i] + mode_starts[i+1]) / 2
        ax.text(data_ordered.shape[1]+0.5, mid, mode_names[i], ha='left', va='center', fontsize=9, fontweight='bold')
    plt.title(f"NMF basis with mode blocks (rank={n_rank})")
    plt.tight_layout()
    fig.savefig(os.path.join(rank_dir, f"nmf_heatmap_basis_blocks_rank{n_rank}.pdf"), bbox_inches='tight')
    plt.close()

    # ---------- 2.2 系数矩阵热图 ----------
    H_coef = W.T
    H_coef_scaled = (H_coef - H_coef.mean(axis=1, keepdims=True)) / (H_coef.std(axis=1, keepdims=True) + 1e-10)
    fig_width = max(8, H_coef_scaled.shape[1]*0.12)
    g = sns.clustermap(H_coef_scaled, row_cluster=False, col_cluster=True,
                       cmap=sns.diverging_palette(250, 10, as_cmap=True), center=0,
                       xticklabels=samples, yticklabels=mode_names,
                       figsize=(fig_width, 4), cbar_pos=(0.02, 0.85, 0.03, 0.15))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=5)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
    g.fig.suptitle(f"NMF coefficients with sample clustering (rank={n_rank})", y=1.02, fontsize=12)
    g.fig.subplots_adjust(top=0.95)
    g.savefig(os.path.join(rank_dir, f"nmf_heatmap_coef_clean_rank{n_rank}.pdf"), bbox_inches='tight')
    plt.close()

    # ---------- 2.3 原分组箱线图（保留 Bonferroni） ----------
    weights_df = pd.DataFrame(W, index=samples, columns=mode_names)
    weights_df[GROUP_COL] = group_info.loc[samples, GROUP_COL].values
    fig, axes = plt.subplots(1, n_rank, figsize=(3*n_rank, 4), sharey=False)
    if n_rank == 1: axes = [axes]
    for i, mode in enumerate(mode_names):
        sns.boxplot(x=GROUP_COL, y=mode, data=weights_df, order=unique_groups, palette=group_colors, ax=axes[i])
        sns.stripplot(x=GROUP_COL, y=mode, data=weights_df, order=unique_groups, color='black', size=4, alpha=0.5, ax=axes[i])
        axes[i].set_title(mode); axes[i].set_xlabel(""); axes[i].set_ylabel("Weight" if i==0 else "")
        axes[i].tick_params(axis='x', rotation=45, labelsize=7)
        group_data = {g: weights_df[weights_df[GROUP_COL]==g][mode].values for g in unique_groups}
        pairs = list(combinations(unique_groups, 2))
        pvals = [stats.mannwhitneyu(group_data[g1], group_data[g2], alternative='two-sided')[1] for g1, g2 in pairs]
        pvals_corrected = np.clip(np.array(pvals) * len(pairs), 0, 1)
        y_max = weights_df[mode].max(); y_pos = y_max*1.05; y_step = y_max*0.08
        for idx, (g1, g2) in enumerate(pairs):
            p_corr = pvals_corrected[idx]
            star = '' if p_corr>=0.05 else ('*' if p_corr>=0.01 else ('**' if p_corr>=0.001 else '***'))
            if star:
                x1, x2 = unique_groups.index(g1), unique_groups.index(g2)
                axes[i].plot([x1, x2], [y_pos, y_pos], color='black', lw=1)
                axes[i].text((x1+x2)/2, y_pos+y_step*0.2, star, ha='center', va='bottom', fontsize=10)
                y_pos += y_step
    plt.suptitle(f"NMF mode weights across groups (rank={n_rank})", y=1.05)
    plt.tight_layout()
    fig.savefig(os.path.join(rank_dir, f"nmf_group_boxplot_rank{n_rank}.pdf"), bbox_inches='tight')
    plt.close()

    # ---------- 2.4 原分组小提琴图 ----------
    fig, axes = plt.subplots(1, n_rank, figsize=(3*n_rank, 4), sharey=False)
    if n_rank == 1: axes = [axes]
    for i, mode in enumerate(mode_names):
        sns.violinplot(x=GROUP_COL, y=mode, data=weights_df, order=unique_groups, palette=group_colors,
                       inner='quartile', cut=0, ax=axes[i])
        sns.stripplot(x=GROUP_COL, y=mode, data=weights_df, order=unique_groups,
                      color='black', size=4, alpha=0.5, ax=axes[i])
        axes[i].set_title(mode); axes[i].set_xlabel(""); axes[i].set_ylabel("Weight" if i==0 else "")
        axes[i].tick_params(axis='x', rotation=45, labelsize=7)
    plt.suptitle(f"NMF mode weights (violin) across groups (rank={n_rank})", y=1.02)
    plt.tight_layout()
    fig.savefig(os.path.join(rank_dir, f"nmf_group_violin_rank{n_rank}.pdf"), bbox_inches='tight')
    plt.close()

    # ================== 联合箱线图（FDR 校正） ==================

    # ---- A. NpCR vs pCR ----
    def classify_npcr_pcr(group_val):
        # 根据 group 列的真实值调整，这里保留包含下划线的判断
        # 如果你的 group 列没有下划线，请改为 if group_val == "NpCR": ...
        if "_NpCR_" in group_val or group_val == "NpCR":
            return "NpCR"
        elif "_pCR_" in group_val or group_val == "pCR":
            return "pCR"
        else:
            return "Other"

    group_values = group_info.loc[samples, GROUP_COL].values
    pcr_status = [classify_npcr_pcr(g) for g in group_values]
    valid_mask = [s != "Other" for s in pcr_status]
    if sum(valid_mask) >= 2:
        plot_npcr = weights_df.copy()
        plot_npcr["NpCR_pCR"] = pcr_status
        plot_npcr = plot_npcr[valid_mask]
        plot_long_npcr = pd.melt(plot_npcr, id_vars=[GROUP_COL, "NpCR_pCR"],
                                 value_vars=mode_names, var_name="Mode", value_name="Weight")
        groups_npcr = ["NpCR", "pCR"]
        pvals_raw = []
        for mode in mode_names:
            d1 = plot_npcr[plot_npcr["NpCR_pCR"]==groups_npcr[0]][mode].values
            d2 = plot_npcr[plot_npcr["NpCR_pCR"]==groups_npcr[1]][mode].values
            pvals_raw.append(stats.mannwhitneyu(d1, d2, alternative='two-sided')[1])
        # ===== FDR 校正 =====
        reject, pvals_corrected, _, _ = multipletests(pvals_raw, alpha=0.05, method='fdr_bh')
        fig, ax = plt.subplots(figsize=(max(6, n_rank*1.5), 5))
        sns.boxplot(x="Mode", y="Weight", hue="NpCR_pCR", data=plot_long_npcr,
                    order=mode_names, hue_order=groups_npcr,
                    palette={"NpCR": "#1f77b4", "pCR": "#ff7f0e"}, ax=ax)
        sns.stripplot(x="Mode", y="Weight", hue="NpCR_pCR", data=plot_long_npcr,
                      order=mode_names, hue_order=groups_npcr,
                      color='black', dodge=True, size=4, alpha=0.5, legend=False, ax=ax)
        handles, labels = ax.get_legend_handles_labels()
        new_handles = [h for h in handles if isinstance(h, plt.Rectangle) or isinstance(h, matplotlib.patches.Patch)]
        ax.legend(new_handles, groups_npcr, title="NpCR/pCR")
        y_max_npcr = plot_long_npcr["Weight"].max()
        y_pos = y_max_npcr * 1.05
        y_step = y_max_npcr * 0.06
        for i, mode in enumerate(mode_names):
            p_corr = pvals_corrected[i]
            star = '' if p_corr >= 0.05 else ('*' if p_corr>=0.01 else ('**' if p_corr>=0.001 else '***'))
            if star:
                x1, x2 = i-0.2, i+0.2
                ax.plot([x1, x2], [y_pos, y_pos], color='black', lw=1)
                ax.text((x1+x2)/2, y_pos+y_step*0.2, star, ha='center', va='bottom', fontsize=12)
                y_pos += y_step
        ax.set_title(f"NpCR vs pCR NMF weights (FDR corrected, rank={n_rank})")
        ax.set_ylabel("Weight")
        ax.tick_params(axis='x', rotation=45)
        plt.tight_layout()
        fig.savefig(os.path.join(rank_dir, f"nmf_NpCR_pCR_boxplot_rank{n_rank}.pdf"), bbox_inches='tight')
        plt.close()

    # ---- B. 多组联合箱线图（FDR 校正） ----
    n_groups = len(unique_groups)
    if n_groups >= 2:
        plot_long_group = pd.melt(weights_df, id_vars=[GROUP_COL],
                                  value_vars=mode_names, var_name="Mode", value_name="Weight")
        pairs_group = list(combinations(unique_groups, 2))
        # 收集所有原始 p 值（每个模式 × 每对分组）
        pvals_flat = []
        for i, mode in enumerate(mode_names):
            for j, (g1, g2) in enumerate(pairs_group):
                d1 = weights_df[weights_df[GROUP_COL]==g1][mode].values
                d2 = weights_df[weights_df[GROUP_COL]==g2][mode].values
                pvals_flat.append(stats.mannwhitneyu(d1, d2, alternative='two-sided')[1])
        # FDR 校正
        reject_flat, pvals_corrected_flat, _, _ = multipletests(pvals_flat, alpha=0.05, method='fdr_bh')
        # 重新排列为矩阵
        pval_matrix = np.array(pvals_corrected_flat).reshape(n_rank, len(pairs_group))

        fig, ax = plt.subplots(figsize=(max(8, n_rank*1.5), 5))
        sns.boxplot(x="Mode", y="Weight", hue=GROUP_COL, data=plot_long_group,
                    order=mode_names, hue_order=unique_groups, palette=group_colors, ax=ax)
        sns.stripplot(x="Mode", y="Weight", hue=GROUP_COL, data=plot_long_group,
                      order=mode_names, hue_order=unique_groups,
                      color='black', dodge=True, size=4, alpha=0.5, legend=False, ax=ax)
        handles, labels = ax.get_legend_handles_labels()
        new_handles = [h for h in handles if isinstance(h, plt.Rectangle) or isinstance(h, matplotlib.patches.Patch)]
        ax.legend(new_handles, unique_groups, title=GROUP_COL, bbox_to_anchor=(1.05, 1), loc='upper left')
        y_max_group = plot_long_group["Weight"].max()
        y_offsets = np.linspace(y_max_group*1.05, y_max_group*1.4, n_rank*len(pairs_group))
        offset_idx = 0
        for i, mode in enumerate(mode_names):
            for j, (g1, g2) in enumerate(pairs_group):
                p_corr = pval_matrix[i, j]
                star = '' if p_corr >= 0.05 else ('*' if p_corr>=0.01 else ('**' if p_corr>=0.001 else '***'))
                if star:
                    x1 = i + (unique_groups.index(g1) - (n_groups-1)/2) * 0.25
                    x2 = i + (unique_groups.index(g2) - (n_groups-1)/2) * 0.25
                    y = y_offsets[offset_idx]
                    ax.plot([x1, x2], [y, y], color='black', lw=1)
                    ax.text((x1+x2)/2, y+0.01, star, ha='center', va='bottom', fontsize=8)
                    offset_idx += 1
        ax.set_title(f"Group comparison NMF weights (FDR corrected, rank={n_rank})")
        ax.set_ylabel("Weight")
        ax.tick_params(axis='x', rotation=45)
        plt.tight_layout()
        fig.savefig(os.path.join(rank_dir, f"nmf_multi_group_boxplot_rank{n_rank}.pdf"), bbox_inches='tight')
        plt.close()

    # ---------- top10 气泡图 ----------
    print(f"  Rank {n_rank}: Generating per-mode top10 bubble plot...")
    top_n = min(10, len(cell_types))
    n_cols = min(5, n_rank)
    n_rows = int(np.ceil(n_rank / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rank == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    for idx in range(n_rank):
        ax = axes[idx]
        mode_weights = H[idx, :]
        top_indices = np.argsort(mode_weights)[::-1][:top_n]
        top_cells = [cell_types[i] for i in top_indices]
        top_weights = mode_weights[top_indices]
        df_bubble = pd.DataFrame({'Cell type': top_cells, 'Weight': top_weights})
        sns.scatterplot(data=df_bubble, x='Weight', y='Cell type', size='Weight',
                        sizes=(50, 300), legend=False, ax=ax, color='steelblue')
        ax.set_title(f"{mode_names[idx]} (Top {top_n})", fontsize=10)
        ax.set_xlabel("Weight")
        ax.tick_params(axis='y', labelsize=7)
    for idx in range(n_rank, len(axes)):
        axes[idx].set_visible(False)
    plt.suptitle(f"Top {top_n} cell types per NMF mode (rank={n_rank})", fontsize=14)
    plt.tight_layout()
    fig.savefig(os.path.join(rank_dir, f"nmf_top10_bubble_per_mode_rank{n_rank}.pdf"), bbox_inches='tight')
    plt.close()

    # ---------- 细胞类型相似度热图 ----------
    print(f"  Rank {n_rank}: Generating cell type cosine similarity heatmap...")
    cell_sim = cosine_similarity(H.T)
    order_by_mode = np.argsort(dominant_mode_per_cell)
    cell_types_ordered_sim = [cell_types[i] for i in order_by_mode]
    sim_ordered = cell_sim[order_by_mode, :][:, order_by_mode]
    mode_starts_sim = []; current_mode = -1
    for i, mode_idx in enumerate(dominant_mode_per_cell[order_by_mode]):
        if mode_idx != current_mode:
            current_mode = mode_idx
            mode_starts_sim.append(i)
    mode_starts_sim.append(len(cell_types_ordered_sim))
    mode_names_short = [f"EC{i+1}" for i in range(n_rank)]
    fig, ax = plt.subplots(figsize=(max(8, len(cell_types_ordered_sim)*0.3), max(8, len(cell_types_ordered_sim)*0.3)))
    sns.heatmap(sim_ordered, xticklabels=False, yticklabels=cell_types_ordered_sim,
                cmap='RdBu_r', vmin=0, vmax=1, center=0.5, cbar_kws={'label': 'Cosine similarity'}, ax=ax)
    ax.set_title(f"Cell type × Cell type cosine similarity (rank={n_rank})", fontsize=12, pad=20)
    ax.tick_params(axis='y', rotation=0, labelsize=6)
    for i in range(len(mode_starts_sim)-1):
        start, end = mode_starts_sim[i], mode_starts_sim[i+1]
        if end-start > 0:
            rect = Rectangle((start, start), end-start, end-start, linewidth=2.5, edgecolor='black', facecolor='none')
            ax.add_patch(rect)
    for i in range(len(mode_starts_sim)-1):
        start, end = mode_starts_sim[i], mode_starts_sim[i+1]
        if end-start > 0:
            mid = (start+end)/2
            ax.text(mid, -0.5, mode_names_short[i], ha='center', va='top', fontsize=10, fontweight='bold', color='black')
            ax.text(len(cell_types_ordered_sim)+0.5, mid, mode_names_short[i], ha='left', va='center', fontsize=10, fontweight='bold', color='black')
    plt.tight_layout()
    fig.savefig(os.path.join(rank_dir, f"nmf_celltype_celltype_sim_rank{n_rank}.pdf"), bbox_inches='tight')
    plt.close()

    # 保存CSV
    pd.DataFrame(W, index=samples, columns=mode_names).to_csv(os.path.join(rank_dir, f"nmf_sample_weights_rank{n_rank}.csv"))
    pd.DataFrame(H, index=mode_names, columns=cell_types).to_csv(os.path.join(rank_dir, f"nmf_celltype_weights_rank{n_rank}.csv"))

    return W, H, cell_types, mode_names, rec_err

# ==================== 3. 主循环 ====================
reconstruction_errors = {}
all_results = {}
for n_rank in N_RANK_LIST:
    print(f"\n===== 处理 N_RANK = {n_rank} =====")
    rank_dir = os.path.join(OUTPUT_DIR, f"rank_{n_rank}")
    os.makedirs(rank_dir, exist_ok=True)
    W, H, cell_types, mode_names, rec_err = run_nmf_rank(n_rank, freq_matrix, group_info, group_colors, unique_groups, rank_dir)
    reconstruction_errors[n_rank] = rec_err
    all_results[n_rank] = (W, H, cell_types, mode_names)
    print(f"完成 N_RANK={n_rank}, 重建误差={rec_err:.4f}")

# ==================== 4. 全局汇总图 ====================
ranks = sorted(reconstruction_errors.keys())
errors = [reconstruction_errors[r] for r in ranks]
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(ranks, errors, 'o-', color='steelblue', linewidth=2, markersize=8)
ax.set_xlabel("Number of NMF components (rank)")
ax.set_ylabel("Reconstruction error")
ax.set_title("Reconstruction error vs. NMF rank")
ax.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "nmf_reconstruction_error.pdf"))
plt.close()

all_modes, all_weights = [], []
for n_rank in N_RANK_LIST:
    _, H, cell_types, mode_names = all_results[n_rank]
    for i, mn in enumerate(mode_names):
        all_modes.append(f"R{n_rank}_{mn}")
        all_weights.append(H[i, :])
summary_df = pd.DataFrame(all_weights, index=all_modes, columns=cell_types)
summary_scaled = (summary_df.T - summary_df.T.mean()) / (summary_df.T.std() + 1e-10)
summary_scaled = summary_scaled.loc[summary_scaled.std(axis=1) > 1e-10]
if summary_scaled.shape[0] > 0:
    linkage_rows = linkage(summary_scaled, method='ward')
    order_rows = leaves_list(linkage_rows)
    if len(order_rows) != summary_scaled.shape[0]:
        order_rows = np.arange(summary_scaled.shape[0])
    data_ordered = summary_scaled.iloc[order_rows, :]
    height = max(5, len(cell_types)*0.3)
    fig, ax = plt.subplots(figsize=(max(8, len(all_modes)*0.8), height))
    sns.heatmap(data_ordered, xticklabels=all_modes, yticklabels=data_ordered.index.tolist(),
                cmap=sns.diverging_palette(250, 10, as_cmap=True), center=0, cbar_kws={"label": "Z-score"}, ax=ax)
    ax.set_title("Cross-rank cell type weights (hierarchical clustering)")
    ax.tick_params(axis='both', labelsize=7)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "nmf_cross_rank_heatmap.pdf"), bbox_inches='tight')
    plt.close()

print(f"\n所有结果已保存至: {OUTPUT_DIR}")
