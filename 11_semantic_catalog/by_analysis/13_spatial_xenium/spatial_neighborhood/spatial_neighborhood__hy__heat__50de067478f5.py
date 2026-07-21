#!/usr/bin/env python3
"""
基于 short/long patient 分组绘制 CNs 平均样本占比热图（白色→红色渐变）。
两组间 Mann-Whitney U 检验显著（p<0.05）的 CN 在均值较高组标注星号。
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.colors import LinearSegmentedColormap
import os
import re

# ---------- 0. 配置 ----------
CSV_PATH = "CNs_10_all_cell_features.csv"
OUT_DIR = "./heatmaps"
os.makedirs(OUT_DIR, exist_ok=True)

# 分组列表（patient ID，需与 CSV 中 sample 列匹配）
short_patients = [
    "Z2038176-15","Q2324341-1","Z1530765-11","F1706967-11","F1805049-2",
    "Z1915172-7","Z1954693-9","ZF2207170-8","ZF2202916-8","ZF2208790-1",
    "Z2262137-5","Z2271577-12","ZF2304103-8","Z2323818-8","Z2336500-6",
    "Z2345459-7","Z2332034-5","Z2417954-8","Q2236500-1","Q2406275-2",
    "Z2101519-1","Z2123354-1","Z1642635-14","Z2224197-10","Z1736838-15",
    "Z1944219-14","Z1430967-9","Z2057500-35","Z2171766-3","Z2353385-10",
    "Z1635989-12","Z1814675-2","Z1826868-6","Z1952876-9","Z1802235-10"
]
long_patients = [
    "Z2279427-2","H2510098-1","H2510098-2","Q2243285-1","Z2370717-2",
    "Q2337205-2","Q2233699-2","Q2239242-1","Q2319383-1","Q2322393-1",
    "Q2337350-1"
]

# 自定义分组名称
SHORT_GROUP = "MPC"
LONG_GROUP = "TNBC"
STATE_ORDER = [SHORT_GROUP, LONG_GROUP]

# 显著性阈值
P_THRESH = 0.05

# ---------- 1. 读取数据 ----------
print("Reading data from:", CSV_PATH)
df = pd.read_csv(CSV_PATH,sep="\t")

# 检查必须列
required_base = ['sample']
for col in required_base:
    assert col in df.columns, f"Missing column: {col}"

# 自动检测 CNs_ 列
cn_cols = [c for c in df.columns if c.startswith('CNs_')]
if not cn_cols:
    raise ValueError("No CNs_* columns found in the CSV file.")
print(f"Detected CN columns: {cn_cols}")

# 丢弃 CN 列为空的行
df = df.dropna(subset=cn_cols)
print(f"Total cells after dropping missing CNs: {len(df)}")

# ---------- 2. 辅助函数 ----------
def sort_cn(cn_list):
    """按末尾数字排序 CN 标签，例如 'CN0','CN1','10' """
    def key(x):
        s = str(x)
        m = re.search(r'(\d+)$', s)
        if m:
            return (0, int(m.group(1)))
        return (1, s)
    return sorted(cn_list, key=key)

def compute_sample_cn_fraction(df, cn_col):
    """
    计算每个样本在每个 CN 中的细胞比例。
    返回: pivot 表 (sample x CN)，值为该样本在该 CN 的细胞比例。
    """
    # 每个样本的总细胞数
    sample_counts = df.groupby("sample").size()
    # 每个样本在 CN 中的计数
    cn_counts = df.groupby(["sample", cn_col]).size().unstack(fill_value=0)
    # 计算比例
    cn_fraction = cn_counts.div(sample_counts, axis=0)
    # 确保所有 CN 列存在
    all_cns = sort_cn(df[cn_col].unique())
    cn_fraction = cn_fraction.reindex(columns=all_cns, fill_value=0.0)
    return cn_fraction

def get_group_mean(cn_fraction_df, group_map):
    """
    group_map: Series, index=sample, value=group (MPC/TNBC)
    返回: DataFrame, index=group (按 STATE_ORDER), columns=CN
    """
    mean_dict = {}
    for group in STATE_ORDER:
        samples = group_map[group_map == group].index
        if len(samples) == 0:
            mean_dict[group] = pd.Series(0, index=cn_fraction_df.columns)
        else:
            mean_dict[group] = cn_fraction_df.loc[samples].mean()
    mean_df = pd.DataFrame(mean_dict).T
    mean_df = mean_df.reindex(STATE_ORDER)
    return mean_df

def annotate_significance(mean_df, cn_fraction_df, group_map, p_thresh=0.05):
    """
    比较 MPC vs TNBC，显著则在均值较高的组对应格子标注 "*"。
    """
    annot = pd.DataFrame("", index=mean_df.index, columns=mean_df.columns)
    samples_short = group_map[group_map == SHORT_GROUP].index
    samples_long = group_map[group_map == LONG_GROUP].index

    for cn in mean_df.columns:
        if len(samples_short) > 0 and len(samples_long) > 0:
            vals_short = cn_fraction_df.loc[samples_short, cn].values
            vals_long = cn_fraction_df.loc[samples_long, cn].values
            try:
                _, p = mannwhitneyu(vals_short, vals_long, alternative='two-sided')
            except ValueError:
                p = 1.0
            if p < p_thresh:
                # 在均值较高的组标注 *
                if mean_df.loc[SHORT_GROUP, cn] > mean_df.loc[LONG_GROUP, cn]:
                    annot.loc[SHORT_GROUP, cn] = "*"
                else:
                    annot.loc[LONG_GROUP, cn] = "*"
    return annot

def draw_heatmap(mean_df, annot_df, cn_col, out_path):
    """绘制白色→红色热图，星号加大加粗"""
    plt.figure(figsize=(max(6, len(mean_df.columns)*0.4), 3))
    cmap = LinearSegmentedColormap.from_list('white_red', ['white', 'red'])

    ax = sns.heatmap(
        mean_df,
        annot=annot_df.values if annot_df is not None else None,
        fmt="",
        linewidths=0.5,
        cmap=cmap,
        cbar_kws={"label": "Mean fraction per sample"},
        vmin=0, vmax=None,
        annot_kws={'size': 16, 'weight': 'bold'}
    )
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_title(f"Mean cell fraction per CN ({cn_col})")
    ax.set_xlabel("Cellular Neighborhood (CN)")
    ax.set_ylabel("Group")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")

# ---------- 3. 构建样本 -> 分组映射 ----------
# 根据 sample 列是否在 short/long 列表中分配组别
sample_series = df['sample']
group_map = pd.Series(index=sample_series.unique(), dtype=str)
group_map[group_map.index.isin(short_patients)] = SHORT_GROUP
group_map[group_map.index.isin(long_patients)] = LONG_GROUP
# 移除未匹配的样本
group_map = group_map.dropna()
print(f"Matched samples: {len(group_map)} (MPC: {(group_map==SHORT_GROUP).sum()}, TNBC: {(group_map==LONG_GROUP).sum()})")
if len(group_map) == 0:
    raise ValueError("No samples matched the provided patient lists. Check sample names.")

# ---------- 4. 主循环：对每个 CN 列生成热图 ----------
for cn_col in cn_cols:
    print(f"\nProcessing {cn_col} ...")
    cn_frac = compute_sample_cn_fraction(df, cn_col)

    # 只保留已分组的样本
    common_samples = cn_frac.index.intersection(group_map.index)
    cn_frac = cn_frac.loc[common_samples]
    current_group_map = group_map.loc[common_samples]

    # 分组平均
    mean_df = get_group_mean(cn_frac, current_group_map)

    # 显著性标注
    annot_df = annotate_significance(mean_df, cn_frac, current_group_map, P_THRESH)

    # 输出数值
    print("Mean fraction per sample (group x CN):")
    print(mean_df.round(4))
    print("Significance annotations:")
    print(annot_df)

    # 画图并保存
    out_file = os.path.join(OUT_DIR, f"heatmap_{cn_col}.png")
    draw_heatmap(mean_df, annot_df, cn_col, out_file)

print("\nDone.")
