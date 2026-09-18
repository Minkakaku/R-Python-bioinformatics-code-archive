#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对上皮细胞群体的通路评分进行样本级聚合（中位数/平均值），
比较不同 group 间的差异并绘制柱状图。
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, f_oneway, mannwhitneyu
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# -------------------- 0. 配置 --------------------
DATA_PATH = "all_with_group.csv"

# 目标上皮细胞类型
EPI_TYPES = [
    "CXCL5+ Inf stem tumor Epi",
    "IGF2+ prolif stem tumor Epi",
    "Hypo Stem tumor Epi",
    "WNT+ ISC",
    "ESC-like Epi",
    "Prol Metab Epi",
    "EMT-like Inv Ent",
    "Inf Hyp Epi",
    "Early stress Epi",
    "Metab stress Epi",
    "Mat Ent",
    "MUC5B+ CGC",
    "MUC5AC+ MGC",
    "PCSK1+ AGC",
    "Hypo Metab Epi"
]

# 通路名称与颜色映射
PATHWAY_COLORS = {
    "angiogenesis_score": "#006400",
    "cellular_metabolism_score": "#8B4513",
    "cycling_score": "#FF4500",
    "dedifferentiation_score": "#00008B",
    "DNA_repair_score": "#008B8B",
    "EMT_score": "#DC143C",
    "hypoxia_score": "#800080",
    "inflammation_score": "#006400",      # 与 angiogenesis 相同，保留
    "oncogenic_TFs_score": "#4B0082",
    "stemness_score": "#C71585"
}

IMG_DPI = 300

# -------------------- 1. 数据加载与预处理 --------------------
df = pd.read_csv(DATA_PATH)
print(f"原始数据形状: {df.shape}")

# 获取通路列名列表
pathways = list(PATHWAY_COLORS.keys())

# !!! 关键修复：将通路评分列强制转换为数值类型（避免字符串导致的聚合错误）
for p in pathways:
    df[p] = pd.to_numeric(df[p], errors='coerce')

# 可选：去除那些在所有通路上都为 NaN 的行（保持数据整洁）
df.dropna(subset=pathways, how='all', inplace=True)

# 只保留目标上皮细胞
df_epi = df[df['celltype'].isin(EPI_TYPES)].copy()
print(f"筛选上皮细胞后形状: {df_epi.shape}")
print("Group 类别:", df_epi['group'].unique())

# -------------------- 2. 按 sample + group 聚合：中位数和平均值 --------------------
# 创建每个样本-组别的唯一标识
df_epi['sample_group'] = df_epi['sample'].astype(str) + "|" + df_epi['group'].astype(str)

# 计算每个组合的中位数
agg_median = df_epi.groupby(['sample_group', 'group'])[pathways].median().reset_index()
print(f"中位数聚合后样本-组别数: {agg_median.shape[0]}")

# 计算每个组合的平均值
agg_mean = df_epi.groupby(['sample_group', 'group'])[pathways].mean().reset_index()
print(f"平均值聚合后样本-组别数: {agg_mean.shape[0]}")

# ★★★ 新增：输出所有样本的聚合统计表（CSV）★★★
agg_median.to_csv("pathway_scores_sample_group_median.csv", index=False)
agg_mean.to_csv("pathway_scores_sample_group_mean.csv", index=False)
print("已保存聚合统计表：pathway_scores_sample_group_median.csv 和 pathway_scores_sample_group_mean.csv")

# -------------------- 3. 统计检验函数 --------------------
def compare_groups(data, value_col, group_col='group', test_type='median'):
    """
    对聚合后数据进行组间差异检验。
    test_type: 'median' -> Kruskal-Wallis + 两两 Mann-Whitney U (BH校正)
              'mean'   -> ANOVA + Tukey HSD
    返回: 总体p值, 两两比较结果 [(组1, 组2, p值, 是否显著)]
    """
    groups = data[group_col].unique()
    if len(groups) < 2:
        print("警告: 组别数小于2，无法进行统计检验。")
        return 1.0, []

    group_data = {g: data[data[group_col] == g][value_col].values for g in groups}

    if test_type == 'median':
        # Kruskal-Wallis 检验（非参数）
        stat, overall_p = kruskal(*group_data.values())
        print(f"Kruskal-Wallis H={stat:.3f}, p={overall_p:.4f}")

        pairwise = []
        if overall_p < 0.05:
            combos = list(combinations(groups, 2))
            pvals = []
            for g1, g2 in combos:
                _, p = mannwhitneyu(group_data[g1], group_data[g2], alternative='two-sided')
                pvals.append(p)

            # Benjamini-Hochberg FDR 校正
            from statsmodels.stats.multitest import multipletests
            reject, pvals_corr, _, _ = multipletests(pvals, method='fdr_bh')
            for (g1, g2), p_corr, rej in zip(combos, pvals_corr, reject):
                pairwise.append((g1, g2, p_corr, rej))
        return overall_p, pairwise

    elif test_type == 'mean':
        # 单因素方差分析 (ANOVA)
        stat, overall_p = f_oneway(*group_data.values())
        print(f"ANOVA F={stat:.3f}, p={overall_p:.4f}")

        pairwise = []
        if overall_p < 0.05:
            tukey = pairwise_tukeyhsd(endog=data[value_col],
                                      groups=data[group_col],
                                      alpha=0.05)
            tukey_df = pd.DataFrame(data=tukey.summary().data[1:],
                                    columns=tukey.summary().data[0])
            for _, row in tukey_df.iterrows():
                g1 = row['group1']
                g2 = row['group2']
                p = float(row['p-adj'])
                # 兼容 reject 列可能是 bool 或字符串
                raw_rej = row['reject']
                if isinstance(raw_rej, bool):
                    rej = raw_rej
                else:
                    rej = str(raw_rej).strip().lower() == 'true'
                pairwise.append((g1, g2, p, rej))
        return overall_p, pairwise
    else:
        raise ValueError("test_type 必须为 'median' 或 'mean'")

# -------------------- 4. 绘图函数 --------------------
def plot_pathway_comparison(data, pathway, test_type='median', save_dir='./'):
    """
    绘制单个通路的组间柱状图，添加显著性标注。
    """
    plt.figure(figsize=(6, 5))
    ax = plt.gca()

    group_order = sorted(data['group'].unique())
    color = PATHWAY_COLORS.get(pathway, "#333333")

    # 柱状图 (均值 ± SEM)
    sns.barplot(
        x='group', y=pathway, data=data, order=group_order,
        color=color, alpha=0.7, edgecolor='black',
        capsize=0.1, errwidth=1.5, ax=ax
    )

    # 叠加样本散点
    sns.stripplot(
        x='group', y=pathway, data=data, order=group_order,
        color='black', size=4, alpha=0.7, dodge=False, ax=ax
    )

    # 执行统计检验
    overall_p, pairwise = compare_groups(data, pathway, test_type=test_type)

    # 添加显著性标记
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    sig_pairs = [(g1, g2, p) for (g1, g2, p, rej) in pairwise if rej]
    if sig_pairs:
        step = y_range * 0.08
        for i, (g1, g2, p) in enumerate(sig_pairs):
            x1 = group_order.index(g1)
            x2 = group_order.index(g2)
            y_line = y_max + y_range * 0.05 + i * step
            # 绘制横线
            ax.plot([x1, x1, x2, x2],
                    [y_line, y_line + y_range*0.01, y_line + y_range*0.01, y_line],
                    lw=1.5, c='black')
            # 星号
            if p < 0.001:
                stars = '***'
            elif p < 0.01:
                stars = '**'
            elif p < 0.05:
                stars = '*'
            else:
                stars = 'ns'
            ax.text((x1 + x2)/2, y_line + y_range*0.015, stars,
                    ha='center', va='bottom', fontsize=12, color='black')
        ax.set_ylim(y_min, y_line + step*1.2)

    # 标题与标签
    pathway_name = pathway.replace('_score', '').replace('_', ' ').title()
    ax.set_title(f"{pathway_name} ({test_type.capitalize()} per Sample-Group)", fontsize=12)
    ax.set_xlabel("Group")
    ax.set_ylabel("Score")
    sns.despine()
    plt.tight_layout()

    # 保存图片 —— 同时输出 PNG 和 PDF
    out_name = f"{pathway}_{test_type}_comparison"
    png_path = f"{save_dir}/{out_name}.png"
    pdf_path = f"{save_dir}/{out_name}.pdf"
    plt.savefig(png_path, dpi=IMG_DPI, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')   # ★ 新增 PDF 保存
    plt.close()
    print(f"已保存: {out_name}.png / .pdf")

# -------------------- 5. 批量绘图 --------------------
for test_type, agg_df in zip(['median', 'mean'], [agg_median, agg_mean]):
    print(f"\n{'='*40}")
    print(f"正在处理基于 **{test_type.upper()}** 聚合的数据")
    print(f"{'='*40}")
    for path in pathways:
        print(f"\n通路: {path}")
        plot_pathway_comparison(agg_df, path, test_type=test_type, save_dir='./')

print("\n✅ 所有图形已生成完毕。")
