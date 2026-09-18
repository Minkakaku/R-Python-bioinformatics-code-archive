#!/usr/bin/env python3
"""
统计指定细胞类型（共15种）间各评分的差异，生成带指定颜色的箱线图。
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

TARGET_TYPES = [
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

COLOR_MAP = {
    "CXCL5+ Inf stem tumor Epi": "#FF5733",
    "IGF2+ prolif stem tumor Epi": "#FF7043",
    "Hypo Stem tumor Epi": "#FF8A65",
    "WNT+ ISC": "#FFAB91",
    "ESC-like Epi": "#FFCCBC",
    "Prol Metab Epi": "#1976D2",
    "EMT-like Inv Ent": "#42A5F5",
    "Inf Hyp Epi": "#64B5F6",
    "Early stress Epi": "#90CAF9",
    "Metab stress Epi": "#BBDEFB",
    "Mat Ent": "#388E3C",
    "MUC5B+ CGC": "#4CAF50",
    "MUC5AC+ MGC": "#66BB6A",
    "PCSK1+ AGC": "#81C784",
    "Hypo Metab Epi": "#A5D6A7"
}

def main():
    import argparse
    parser = argparse.ArgumentParser(description='指定细胞类型评分差异统计')
    parser.add_argument('--input', default='all.csv', help='输入CSV文件路径')
    parser.add_argument('--outdir', default='celltype_score_sta', help='输出目录')
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    score_cols = df.columns[3:]

    # 仅保留目标细胞类型
    df_target = df[df['celltype'].isin(TARGET_TYPES)].copy()
    existing_types = [ct for ct in TARGET_TYPES if ct in df_target['celltype'].unique()]
    if not existing_types:
        print("错误：数据中不包含任何指定的细胞类型。")
        sys.exit(1)

    palette = {ct: COLOR_MAP[ct] for ct in existing_types}

    # ★ 关键修复：将所有评分列转换为数值，无法转换的变为 NaN
    for col in score_cols:
        df_target[col] = pd.to_numeric(df_target[col], errors='coerce')

    print(f"分析的细胞类型（{len(existing_types)}种）：{existing_types}")
    summary_rows = []

    for score in score_cols:
        print(f"正在处理: {score}")
        # 删除该评分为 NaN 的行
        sub = df_target[['celltype', score]].dropna(subset=[score])
        # 只保留至少有数据的组
        groups = [ct for ct in existing_types if ct in sub['celltype'].unique()]
        if len(groups) < 2:
            print(f"  跳过（有效组数不足）")
            continue

        group_data = [sub.loc[sub['celltype'] == g, score].values for g in groups]
        try:
            H, p_kw = kruskal(*group_data)
        except Exception as e:
            print(f"  Kruskal-Wallis 失败: {e}")
            continue

        summary_rows.append({
            'score': score,
            'test': 'Kruskal-Wallis',
            'statistic': H,
            'p_value': p_kw,
            'comparison': 'all'
        })

        significant_pairs = []
        if p_kw < 0.05 and len(groups) >= 3:
            try:
                dunn_pval = sp.posthoc_dunn(sub, val_col=score, group_col='celltype',
                                            p_adjust='fdr_bh')
                for i in range(len(groups)):
                    for j in range(i+1, len(groups)):
                        g1, g2 = groups[i], groups[j]
                        p = dunn_pval.loc[g1, g2]
                        if p < 0.05:
                            significant_pairs.append((g1, g2, p))
                significant_pairs.sort(key=lambda x: x[2])
                for g1, g2, p in significant_pairs:
                    summary_rows.append({
                        'score': score,
                        'test': 'Dunn',
                        'statistic': np.nan,
                        'p_value': p,
                        'comparison': f'{g1} vs {g2}'
                    })
            except Exception as e:
                print(f"  Dunn事后检验失败: {e}")

        # 绘制箱线图
        plt.figure(figsize=(max(12, len(groups) * 1.0), 7))
        order = groups
        ax = sns.boxplot(x='celltype', y=score, data=sub, palette=palette, order=order)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        ax.set_xlabel('')
        plt.title(f"{score}\nKruskal-Wallis p = {p_kw:.4f}", fontsize=12)
        plt.tight_layout()

        # 添加显著星号
        if significant_pairs:
            ymin, ymax = ax.get_ylim()
            y_range = ymax - ymin
            group_to_x = {g: i for i, g in enumerate(order)}
            drawn = []
            n_show = min(len(significant_pairs), 20)
            y_offset = ymax + y_range * 0.03
            bracket_height = y_range * 0.025

            for idx in range(n_show):
                g1, g2, p = significant_pairs[idx]
                if g1 not in group_to_x or g2 not in group_to_x:
                    continue
                x1, x2 = group_to_x[g1], group_to_x[g2]
                if x1 > x2:
                    x1, x2 = x2, x1
                y_top = y_offset
                while True:
                    overlap = False
                    for dx1, dx2, dy in drawn:
                        if not (x2 < dx1 or x1 > dx2):
                            if abs(y_top - dy) < bracket_height * 2:
                                overlap = True
                                break
                    if overlap:
                        y_top += bracket_height * 3
                    else:
                        break
                ax.plot([x1, x1, x2, x2],
                        [y_top, y_top + bracket_height, y_top + bracket_height, y_top],
                        lw=1.2, color='black')
                if p < 0.001:
                    text = '***'
                elif p < 0.01:
                    text = '**'
                elif p < 0.05:
                    text = '*'
                else:
                    text = f'p={p:.2f}'
                ax.text((x1 + x2) / 2, y_top + bracket_height * 1.5, text,
                        ha='center', va='bottom', fontsize=8)
                drawn.append((x1, x2, y_top))
                ymax = max(ymax, y_top + bracket_height * 5)
            ax.set_ylim(ymin - y_range * 0.05, ymax + y_range * 0.1)

        fig_path = os.path.join(args.outdir, f"{score}.png")
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  图表已保存: {fig_path}")

    if summary_rows:
        res_df = pd.DataFrame(summary_rows)
        res_csv = os.path.join(args.outdir, 'statistics_summary.csv')
        res_df.to_csv(res_csv, index=False)
        print(f"\n统计结果汇总已保存至: {res_csv}")

if __name__ == '__main__':
    main()
