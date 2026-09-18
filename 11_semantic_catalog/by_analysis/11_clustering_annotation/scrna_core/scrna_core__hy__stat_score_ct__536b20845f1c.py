import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import os
import sys

# 检查并导入scikit_posthocs
try:
    import scikit_posthocs as sp
except ImportError:
    print("需要安装scikit-posthocs: pip install scikit-posthocs")
    sys.exit(1)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='all.csv')
    parser.add_argument('--outdir', default='celltype_score_sta')
    args = parser.parse_args()
    
    df = pd.read_csv(args.input)
    # 前三列是id...
    score_cols = df.columns[3:]  # 假设从第4列开始都是评分
    # 去除可能非数值列，但已知都是数值
    os.makedirs(args.outdir, exist_ok=True)
    
    # 准备记录统计结果
    summary_rows = []
    # 对每个score
    for score in score_cols:
        # 丢弃缺失
        sub = df[['celltype', score]].dropna()
        groups = sub['celltype'].unique()
        if len(groups) < 2:
            continue
        # 分组数据
        group_data = [sub[sub['celltype'] == g][score].values for g in groups]
        try:
            H, p_kw = kruskal(*group_data)
        except Exception as e:
            print(f"Error in Kruskal for {score}: {e}")
            continue
        # 保存整体结果
        summary_rows.append({'score': score, 'test': 'Kruskal-Wallis', 'statistic': H, 'p_value': p_kw, 'comparison': 'all'})
        # 事后检验
        significant_pairs = []
        if p_kw < 0.05:
            # 使用Dunn检验
            try:
                # posthoc_dunn需要数据框或数组，传入一个列表组数据和p_adjust
                dunn_pval = sp.posthoc_dunn(sub, val_col=score, group_col='celltype', p_adjust='fdr_bh')
                # dunn_pval是一个DataFrame，index和columns为组名
                # 提取上三角
                pairs = []
                for i, g1 in enumerate(groups):
                    for j, g2 in enumerate(groups):
                        if i < j:
                            p = dunn_pval.loc[g1, g2]
                            pairs.append((g1, g2, p))
                # 筛选p<0.05
                significant_pairs = [(g1, g2, p) for g1,g2,p in pairs if p < 0.05]
                # 排序
                significant_pairs.sort(key=lambda x: x[2])
                # 保存显著对
                for g1,g2,p in significant_pairs:
                    summary_rows.append({'score': score, 'test': 'Dunn', 'statistic': np.nan, 'p_value': p, 'comparison': f'{g1} vs {g2}'})
            except Exception as e:
                print(f"Posthoc error for {score}: {e}")
                
        # 绘制箱线图
        plt.figure(figsize=(max(6, len(groups)*1.2), 6))
        ax = sns.boxplot(x='celltype', y=score, data=sub, palette='Set2')
        # 旋转x标签
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        plt.title(f"{score}\nKruskal-Wallis p = {p_kw:.4f}")
        # 添加检验星号括号
        if significant_pairs:
            # 获取y轴范围
            ymax = sub[score].max()
            ymin = sub[score].min()
            y_range = ymax - ymin
            # 按组顺序得到索引
            group_to_x = {g: i for i, g in enumerate(groups)}
            # 为每个显著对画线，从低到高偏移
            # 为了避免重叠，计算层级
            # 最多画10个显著对
            n_sig = min(len(significant_pairs), 10)
            # 创建偏移层级：对于x轴上的多个括号，可能需要多层
            # 简单方法：循环，每次增加y位置
            y_offset = ymax + y_range * 0.05
            bracket_height = y_range * 0.03
            drawn = []  # 记录已使用的括号区间，避免重叠 (x1, x2, y_top)
            for idx, (g1, g2, p) in enumerate(significant_pairs[:n_sig]):
                x1 = group_to_x[g1]
                x2 = group_to_x[g2]
                if x1 > x2:
                    x1, x2 = x2, x1
                # 检查与已有括号是否重叠，如果重叠则增加y
                y_top = y_offset
                while True:
                    overlap = False
                    for (dx1, dx2, dy) in drawn:
                        if not (x2 < dx1 or x1 > dx2): # 区间重叠
                            if abs(y_top - dy) < bracket_height: # y接近
                                overlap = True
                                break
                    if overlap:
                        y_top += bracket_height * 2
                    else:
                        break
                # 画线
                ax.plot([x1, x1, x2, x2], [y_top, y_top + bracket_height, y_top + bracket_height, y_top], lw=1.5, color='black')
                # 添加星号
                if p < 0.001:
                    text = '***'
                elif p < 0.01:
                    text = '**'
                elif p < 0.05:
                    text = '*'
                else:
                    text = f'p={p:.2f}'
                ax.text((x1+x2)/2, y_top + bracket_height*1.5, text, ha='center', va='bottom', fontsize=9)
                drawn.append((x1, x2, y_top))
                # 更新整体ymax用于设置ylim
                ymax = max(ymax, y_top + bracket_height*3)
            ax.set_ylim(ymin - y_range*0.1, ymax + y_range*0.1)
        plt.tight_layout()
        # 保存图片
        fig_path = os.path.join(args.outdir, f"{score}.png")
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved {fig_path}")
    
    # 保存统计结果CSV
    if summary_rows:
        res_df = pd.DataFrame(summary_rows)
        res_df.to_csv(os.path.join(args.outdir, 'statistics_summary.csv'), index=False)
        print("统计结果保存至 statistics_summary.csv")
    else:
        print("没有结果。")

if __name__ == '__main__':
    main()
