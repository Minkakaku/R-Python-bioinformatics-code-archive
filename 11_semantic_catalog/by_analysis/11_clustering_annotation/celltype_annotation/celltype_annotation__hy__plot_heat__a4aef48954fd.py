#!/usr/bin/env python3
"""
根据 state 分组（N/A/C）绘制 CNs_10 和 CNs_15 的平均样本占比热图，
白色→红色渐变；C vs A、A vs N 秩和检验显著（p<0.05）的 CN 方块上标注大星号。
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.colors import LinearSegmentedColormap
import os
import re

# ---------- 0. 配置 ----------
H5AD_PATH = "All_filter_state_samp_new_with_CNs10-15.h5ad"
OUT_DIR = "./heatmaps"
os.makedirs(OUT_DIR, exist_ok=True)

STATE_ORDER = ["N", "A", "C"]   # 从上到下
CN_COLS = ["CNs_10", "CNs_15"]

# 显著性阈值
P_THRESH = 0.05

# ---------- 1. 读取数据 ----------
print("Reading data...")
adata = sc.read_h5ad(H5AD_PATH)

# 检查必要列
required = ["sample_new", "state"] + CN_COLS
for col in required:
    assert col in adata.obs.columns, f"Missing column: {col}"

# 只保留需要的列，去除缺失值
obs_df = adata.obs[required].dropna().copy()
print(f"Total cells: {len(obs_df)}")

# ---------- 2. 辅助函数 ----------
def sort_cn(cn_list):
    """按末尾数字排序 CN 标签，例如 'CN0','CN1','10' """
    def key(x):
        s = str(x)
        m = re.search(r'(\d+)$', s)
        if m:
            return (0, int(m.group(1)))
        # 如果没有数字，按原始字符串
        return (1, s)
    return sorted(cn_list, key=key)

def compute_sample_cn_fraction(df, cn_col):
    """
    计算每个样本在每个 CN 中的细胞比例。
    返回: pivot 表 (sample x CN)，值为该样本在该 CN 的细胞比例。
    """
    # 每个样本的总细胞数
    sample_counts = df.groupby("sample_new").size()
    # 每个样本在 CN 中的计数
    cn_counts = df.groupby(["sample_new", cn_col]).size().unstack(fill_value=0)
    # 计算比例
    cn_fraction = cn_counts.div(sample_counts, axis=0)
    # 确保所有 CN 列存在（即使某个样本没有某个 CN）
    all_cns = sort_cn(df[cn_col].unique())
    cn_fraction = cn_fraction.reindex(columns=all_cns, fill_value=0.0)
    return cn_fraction

def get_group_mean(cn_fraction_df, state_map):
    """
    根据 state 映射，计算每个 state 分组在每个 CN 上的平均样本比例。
    state_map: series, index=sample, value=state
    返回: DataFrame, index=state (按 STATE_ORDER), columns=CN
    """
    mean_dict = {}
    for state in STATE_ORDER:
        samples = state_map[state_map == state].index
        if len(samples) == 0:
            mean_dict[state] = pd.Series(0, index=cn_fraction_df.columns)
        else:
            mean_dict[state] = cn_fraction_df.loc[samples].mean()
    mean_df = pd.DataFrame(mean_dict).T  # state x CN
    # 只保留存在的 state
    mean_df = mean_df.reindex(STATE_ORDER)
    return mean_df

def annotate_significance(mean_df, cn_fraction_df, state_map, p_thresh=0.05):
    """
    返回与 mean_df 形状相同的字符串矩阵，需要标注 * 的格子为 "*"，否则 ""。
    比较：C vs A 且 C 组显著 -> 在 C 行标 *；
         A vs N 且 A 组显著 -> 在 A 行标 *。
    """
    annot = pd.DataFrame("", index=mean_df.index, columns=mean_df.columns)
    
    # 获取各组样本列表
    samples_C = state_map[state_map == "C"].index
    samples_A = state_map[state_map == "A"].index
    samples_N = state_map[state_map == "N"].index
    
    for cn in mean_df.columns:
        # C vs A
        if len(samples_C) > 0 and len(samples_A) > 0:
            vals_C = cn_fraction_df.loc[samples_C, cn].values
            vals_A = cn_fraction_df.loc[samples_A, cn].values
            try:
                _, p_CA = mannwhitneyu(vals_C, vals_A, alternative='two-sided')
            except ValueError:
                p_CA = 1.0
            if p_CA < p_thresh:
                annot.loc["C", cn] = "*"
        
        # A vs N
        if len(samples_A) > 0 and len(samples_N) > 0:
            vals_A = cn_fraction_df.loc[samples_A, cn].values
            vals_N = cn_fraction_df.loc[samples_N, cn].values
            try:
                _, p_AN = mannwhitneyu(vals_A, vals_N, alternative='two-sided')
            except ValueError:
                p_AN = 1.0
            if p_AN < p_thresh:
                annot.loc["A", cn] = "*"
                
    return annot

def draw_heatmap(mean_df, annot_df, cn_col, out_path):
    """绘制并保存热图：白色→红色，星号加大"""
    plt.figure(figsize=(max(6, len(mean_df.columns)*0.4), 3))
    # 自定义白色到红色的 colormap
    cmap = LinearSegmentedColormap.from_list('white_red', ['white', 'red'])
    
    ax = sns.heatmap(
        mean_df,
        annot=annot_df.values if annot_df is not None else None,
        fmt="",
        linewidths=0.5,
        cmap=cmap,
        cbar_kws={"label": "Mean fraction per sample"},
        vmin=0, vmax=None,
        annot_kws={'size': 16, 'weight': 'bold'}  # 加大星号，加粗
    )
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_title(f"Mean cell fraction per CN ({cn_col})")
    ax.set_xlabel("Cellular Neighborhood (CN)")
    ax.set_ylabel("State")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")

# ---------- 3. 主循环：对两个 CN 列分别处理 ----------
for cn_col in CN_COLS:
    print(f"\nProcessing {cn_col} ...")
    # 每个样本的 CN 比例
    cn_frac = compute_sample_cn_fraction(obs_df, cn_col)
    
    # 获取样本的 state 信息（确保 index 对齐）
    sample_state = obs_df.groupby("sample_new")["state"].first()
    common_samples = cn_frac.index.intersection(sample_state.index)
    cn_frac = cn_frac.loc[common_samples]
    sample_state = sample_state.loc[common_samples]
    
    # 分组平均
    mean_df = get_group_mean(cn_frac, sample_state)
    
    # 显著性标注
    annot_df = annotate_significance(mean_df, cn_frac, sample_state, P_THRESH)
    
    # 输出数值表格（可选）
    print("Mean fraction per sample (state x CN):")
    print(mean_df.round(4))
    print("Significance annotations:")
    print(annot_df)
    
    # 画图
    out_file = os.path.join(OUT_DIR, f"heatmap_{cn_col}.png")
    draw_heatmap(mean_df, annot_df, cn_col, out_file)

print("\nDone.")
