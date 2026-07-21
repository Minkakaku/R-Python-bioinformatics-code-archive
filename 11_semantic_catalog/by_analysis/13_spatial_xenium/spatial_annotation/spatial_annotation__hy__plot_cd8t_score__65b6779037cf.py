import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np
from statannotations.Annotator import Annotator

# ---------- 1. 读取数据 ----------
df = pd.read_csv("extracted_scores.csv")
print(f"原始数据行数: {len(df)}")

# ---------- 2. 合并 CD8+ T 细胞 ----------
cd8_mask = df['celltype'].isin(['CD8+ T', 'CXCL9+ CD8+ T'])
df_cd8 = df[cd8_mask].copy()
df_cd8['celltype_merged'] = 'CD8+ T'
print(f"CD8+ T 细胞数: {len(df_cd8)}")

# ---------- 3. 定义评分列 ----------
score_cols = [
    'Effector_CD8T_5k_score',
    'Proliferating_CD8T_5k_score',
    'TerminalExhausted_CD8T_5k_score',
    'IFN_Response_5k_score',
    'Hypoxia_Stress_5k_score'
]

# ---------- 4. 按样本聚合 ----------
sample_group = df_cd8.groupby('sample')['group_final'].first().reset_index()

mean_df = df_cd8.groupby('sample')[score_cols].mean().reset_index()
mean_df.columns = ['sample'] + [f"{col}_mean" for col in score_cols]

median_df = df_cd8.groupby('sample')[score_cols].median().reset_index()
median_df.columns = ['sample'] + [f"{col}_median" for col in score_cols]

df_sample = sample_group.merge(mean_df, on='sample').merge(median_df, on='sample')
print(f"样本数量: {df_sample.shape[0]}")
df_sample.to_csv("cd8t_sample_level_scores.csv", index=False)

# ---------- 5. 分组与比较对 ----------
# 自定义顺序，确保箱线图从左到右符合生物学逻辑
groups = ['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R']
pairs = [("POST-NR", "POST-R"),
         ("POST-NR", "PRE-NR"),
         ("POST-R", "PRE-R"),
         ("PRE-NR", "PRE-R")]

# ---------- 6. 输出目录 ----------
out_dir = "score_group"
os.makedirs(out_dir, exist_ok=True)

# ---------- 7. 绘图函数（原始p值，不校正） ----------
def plot_boxplot_with_stats(data, x, y, order, title, filepath, pairs):
    plt.figure(figsize=(7, 5))
    ax = sns.boxplot(x=x, y=y, data=data, order=order, hue=x, palette='Set2', legend=False)
    sns.stripplot(x=x, y=y, data=data, order=order, color='black', size=5, alpha=0.6)

    annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=order)
    # 仅使用 Mann-Whitney U 检验，不校正，星号标注基于原始 p<0.05
    annotator.configure(test='Mann-Whitney', text_format='star')
    annotator.apply_test().annotate()

    plt.title(title)
    plt.xlabel("Group")
    plt.ylabel("Sample Score")
    plt.tight_layout()
    plt.savefig(f"{filepath}.png", dpi=150)
    plt.savefig(f"{filepath}.pdf")
    plt.close()
    print(f"已保存 {filepath}.png / .pdf")

# ---------- 8. 批量绘图 ----------
for col in score_cols:
    mean_col = f"{col}_mean"
    median_col = f"{col}_median"

    plot_boxplot_with_stats(
        data=df_sample,
        x='group_final', y=mean_col,
        order=groups,
        title=f"CD8+ T cells - {col} (mean per sample)",
        filepath=os.path.join(out_dir, f"CD8_{col}_mean"),
        pairs=pairs
    )
    plot_boxplot_with_stats(
        data=df_sample,
        x='group_final', y=median_col,
        order=groups,
        title=f"CD8+ T cells - {col} (median per sample)",
        filepath=os.path.join(out_dir, f"CD8_{col}_median"),
        pairs=pairs
    )

# ---------- 9. 控制台输出原始 p 值（仅指定四对） ----------
def print_comparisons(data, value_col, group_col):
    print(f"\n--- {value_col} ---")
    for g1, g2 in pairs:
        vals1 = data[data[group_col] == g1][value_col].dropna().values
        vals2 = data[data[group_col] == g2][value_col].dropna().values
        if len(vals1) == 0 or len(vals2) == 0:
            print(f"  {g1} vs {g2}: 数据不足，跳过")
            continue
        try:
            _, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
            sig = "*" if p < 0.05 else "ns"
            print(f"  {g1} vs {g2}: p = {p:.4f} ({sig})")
        except ValueError:
            print(f"  {g1} vs {g2}: 检验失败")

print("\n====== 样本平均数原始 p 值 ======")
for col in score_cols:
    print_comparisons(df_sample, f"{col}_mean", 'group_final')

print("\n====== 样本中位数原始 p 值 ======")
for col in score_cols:
    print_comparisons(df_sample, f"{col}_median", 'group_final')
