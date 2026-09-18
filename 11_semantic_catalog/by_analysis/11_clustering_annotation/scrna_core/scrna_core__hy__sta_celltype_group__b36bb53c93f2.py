#!/usr/bin/env python3
"""
按细胞类型分别统计不同 group 间各评分的差异，并绘制指定颜色的箱线图。
需要安装：pip install pandas numpy scipy matplotlib seaborn scikit-posthocs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import os
import sys

try:
    import scikit_posthocs as sp
except ImportError:
    print("请安装 scikit-posthocs: pip install scikit-posthocs")
    sys.exit(1)

# 评分 → 颜色映射
SCORE_COLORS = {
    "angiogenesis_score": "#006400",
    "cellular_metabolism_score": "#8B4513",
    "cycling_score": "#FF4500",
    "dedifferentiation_score": "#00008B",
    "DNA_repair_score": "#008B8B",
    "EMT_score": "#DC143C",
    "hypoxia_score": "#800080",
    "inflammation_score": "#006400",  # 与angiogenesis相同，可保留
    "oncogenic_TFs_score": "#4B0082",
    "stemness_score": "#C71585"
}


def main():
    import argparse
    parser = argparse.ArgumentParser(description='每个细胞类型内 group 差异统计与可视化')
    parser.add_argument('--input', default='all_with_group.csv', help='输入CSV文件路径')
    parser.add_argument('--outdir', default='celltype_group_sta', help='输出根目录')
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    # 前三列: cell_id, sample, celltype；随后10个评分列；最后一列: group
    # 自动识别评分列（以 _score 结尾的列，或通过已知列表）
    known_scores = list(SCORE_COLORS.keys())
    # 确保列名匹配（去除可能空格）
    df.columns = df.columns.str.strip()
    # 筛选存在的评分列
    score_cols = [col for col in known_scores if col in df.columns]
    if not score_cols:
        print("未找到任何评分列，请检查文件格式。")
        sys.exit(1)

    # 将评分列转换为数值
    for col in score_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # 确定分组列（假设为 "group"）
    group_col = 'group'
    if group_col not in df.columns:
        print(f"未找到分组列 '{group_col}'，请确认列名。")
        sys.exit(1)

    celltype_col = 'celltype'
    if celltype_col not in df.columns:
        print(f"未找到细胞类型列 '{celltype_col}'。")
        sys.exit(1)

    # 去掉 celltype 为空的行
    df.dropna(subset=[celltype_col], inplace=True)

    # 按细胞类型分组处理
    celltypes = df[celltype_col].unique()
    os.makedirs(args.outdir, exist_ok=True)

    overall_summary = []  # 保存所有统计结果

    for ct in celltypes:
        ct_df = df[df[celltype_col] == ct].copy()
        if ct_df[group_col].nunique() < 2:
            print(f"[{ct}] 分组数不足（<2），跳过")
            continue

        ct_dir = os.path.join(args.outdir, ct.replace(' ', '_').replace('/', '_'))
        os.makedirs(ct_dir, exist_ok=True)

        for score in score_cols:
            sub = ct_df[[group_col, score]].dropna(subset=[score])
            groups = sorted(sub[group_col].unique())
            if len(groups) < 2:
                continue

            group_data = [sub.loc[sub[group_col] == g, score].values for g in groups]

            # Kruskal-Wallis
            try:
                H, p_kw = kruskal(*group_data)
            except Exception as e:
                print(f"[{ct}] {score}: Kruskal-Wallis 失败 ({e})")
                continue

            overall_summary.append({
                'celltype': ct,
                'score': score,
                'test': 'Kruskal-Wallis',
                'statistic': H,
                'p_value': p_kw,
                'comparison': 'all groups'
            })

            # Dunn 事后检验
            significant_pairs = []
            if p_kw < 0.05 and len(groups) >= 3:
                try:
                    dunn_pval = sp.posthoc_dunn(sub, val_col=score, group_col=group_col,
                                                p_adjust='fdr_bh')
                    for i in range(len(groups)):
                        for j in range(i+1, len(groups)):
                            g1, g2 = groups[i], groups[j]
                            p = dunn_pval.loc[g1, g2]
                            if p < 0.05:
                                significant_pairs.append((g1, g2, p))
                    significant_pairs.sort(key=lambda x: x[2])
                    for g1, g2, p in significant_pairs:
                        overall_summary.append({
                            'celltype': ct,
                            'score': score,
                            'test': 'Dunn',
                            'statistic': np.nan,
                            'p_value': p,
                            'comparison': f'{g1} vs {g2}'
                        })
                except Exception as e:
                    print(f"[{ct}] {score}: Dunn检验失败 ({e})")

            # ---- 绘制箱线图 ----
            color = SCORE_COLORS.get(score, '#888888')
            plt.figure(figsize=(max(6, len(groups) * 1.2), 5.5))
            ax = sns.boxplot(x=group_col, y=score, data=sub, order=groups,
                             color=color, flierprops={'markersize': 3})
            # 可添加散点
            sns.stripplot(x=group_col, y=score, data=sub, order=groups,
                          color='black', size=1.5, alpha=0.4, ax=ax)

            ax.set_title(f"{ct}\n{score}\nKruskal-Wallis p = {p_kw:.4f}", fontsize=11)
            ax.set_xlabel('Group')
            ax.set_ylabel(score)
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()

            # 添加显著星号
            if significant_pairs:
                ymin, ymax = ax.get_ylim()
                y_range = ymax - ymin
                group_to_x = {g: i for i, g in enumerate(groups)}
                drawn = []
                n_show = min(len(significant_pairs), 15)
                y_offset = ymax + y_range * 0.03
                bracket_h = y_range * 0.025

                for idx in range(n_show):
                    g1, g2, p = significant_pairs[idx]
                    if g1 not in group_to_x or g2 not in group_to_x:
                        continue
                    x1, x2 = group_to_x[g1], group_to_x[g2]
                    if x1 > x2:
                        x1, x2 = x2, x1
                    y_top = y_offset
                    # 避免重叠
                    while True:
                        overlap = False
                        for dx1, dx2, dy in drawn:
                            if not (x2 < dx1 or x1 > dx2):
                                if abs(y_top - dy) < bracket_h * 2:
                                    overlap = True
                                    break
                        if overlap:
                            y_top += bracket_h * 3
                        else:
                            break
                    ax.plot([x1, x1, x2, x2],
                            [y_top, y_top + bracket_h, y_top + bracket_h, y_top],
                            lw=1.2, color='black')
                    if p < 0.001:
                        txt = '***'
                    elif p < 0.01:
                        txt = '**'
                    elif p < 0.05:
                        txt = '*'
                    else:
                        txt = f'p={p:.2f}'
                    ax.text((x1 + x2) / 2, y_top + bracket_h * 1.5, txt,
                            ha='center', va='bottom', fontsize=8)
                    drawn.append((x1, x2, y_top))
                    ymax = max(ymax, y_top + bracket_h * 5)
                ax.set_ylim(ymin - y_range * 0.05, ymax + y_range * 0.1)

            # 保存图像
            safe_score = score.replace(' ', '_')
            fig_path = os.path.join(ct_dir, f"{safe_score}.png")
            plt.savefig(fig_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"[{ct}] {score} 图表已保存: {fig_path}")

    # 保存统计总表
    if overall_summary:
        res_df = pd.DataFrame(overall_summary)
        res_csv = os.path.join(args.outdir, 'statistics_summary.csv')
        res_df.to_csv(res_csv, index=False)
        print(f"\n全局统计汇总已保存至: {res_csv}")

if __name__ == '__main__':
    main()
