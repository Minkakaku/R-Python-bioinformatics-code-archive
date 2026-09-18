import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib import rcParams
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")
sns.set_palette("husl")

# 根据图片数据创建DataFrame
data = {
    'patient': [
        'Z2224197-10', 'ZF2208790-1', 'Z2038176-15', 'ZF2207170-8', 
        'F1706967-11', 'Z1915172-7', 'Z1944219-14', 'Z1954693-9',
        'Q2324341-1', 'F1805049-2', 'Z1530765-11', 'Z1802235-10',
        'Z1736838-15', 'Z2101519-1', 'Z1642635-14', 'Z1430967-9',
        'Z2262137-5', 'Q2236500-1', 'Q2406275-2', 'Z2057500-35',
        'ZF2202916-8', 'Z2123354-1', 'Z2171766-3', 'Z1952876-9'
    ],
    'PFS': [
        11.23, 40.23, 15.77, 41.57, 99.5, 81.03, 41.47, 73.93,
        9.4, 91.67, 124.67, 3.73, 7.17, 11.13, 40.97, 1.47,
        38.9, 36.7, 7.67, 36.13, 45.07, 54.2, 47.3, 75.03
    ],
    'percent': [
        0.653240199, 0.650270275, 0.624303959, 0.54034249, 0.504299079,
        0.502780138, 0.499214076, 0.494021075, 0.486289979, 0.437800957,
        0.435040692, 0.431274564, 0.395217345, 0.316687896, 0.315424114,
        0.30690657, 0.290545298, 0.284139584, 0.185051724, 0.164108442,
        0.145352292, 0.108758911, 0.062626548, 0.061972475
    ],
    'group': [
        'T1', 'T1', 'T1', 'T4', 'T2', 'T2', 'T2', 'T1',
        'T4', 'T3', 'T4', 'T1', 'T3', 'T3', 'T1', 'T3',
        'T4', 'T4', 'T4', 'T5', 'T3', 'T4', 'T5', 'T5'
    ],
    'group_new': [
        'T1', 'T1', 'T1', 'T2', 'T1', 'T1', 'T1', 'T1',
        'T2', 'T1', 'T2', 'T1', 'T1', 'T1', 'T1', 'T1',
        'T2', 'T2', 'T2', 'T2', 'T1', 'T2', 'T2', 'T2'
    ]
}

df = pd.DataFrame(data)

print("=== 数据概览 ===")
print(f"总样本数: {len(df)}")
print(f"group_new分组情况: {df['group_new'].unique()}")
print("\n各分组样本数:")
print(df['group_new'].value_counts().sort_index())

print("\n=== 描述性统计 ===")
# 按group_new分组计算PFS的统计量
group_stats = df.groupby('group_new')['PFS'].agg([
    'count', 'mean', 'std', 'median', 'min', 'max', 
    lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75)
]).rename(columns={'<lambda_0>': 'Q1', '<lambda_1>': 'Q3'})

print(group_stats)

print("\n=== 详细分组统计 ===")
for group in sorted(df['group_new'].unique()):
    group_data = df[df['group_new'] == group]['PFS']
    print(f"\n{group}组 (n={len(group_data)}):")
    print(f"  均值: {group_data.mean():.2f} ± {group_data.std():.2f}")
    print(f"  中位数: {group_data.median():.2f}")
    print(f"  范围: {group_data.min():.2f} - {group_data.max():.2f}")
    print(f"  四分位数: Q1={group_data.quantile(0.25):.2f}, Q3={group_data.quantile(0.75):.2f}")
    print(f"  标准差: {group_data.std():.2f}")
    print(f"  变异系数: {(group_data.std()/group_data.mean()*100):.1f}%")

# 1. 正态性检验
print("\n=== 正态性检验 (Shapiro-Wilk检验) ===")
normality_results = {}
for group in df['group_new'].unique():
    group_data = df[df['group_new'] == group]['PFS']
    if len(group_data) >= 3:  # Shapiro-Wilk检验需要至少3个样本
        stat, p_value = stats.shapiro(group_data)
        normality_results[group] = {'statistic': stat, 'p_value': p_value, 'data': group_data}
        print(f"{group}组 (n={len(group_data)}): W = {stat:.3f}, p = {p_value:.3f}")
    else:
        print(f"{group}组 (n={len(group_data)}): 样本量过少，无法进行正态性检验")

# 2. 方差齐性检验
print("\n=== 方差齐性检验 (Levene检验) ===")
levene_groups = [df[df['group_new'] == group]['PFS'].values for group in df['group_new'].unique()]
if len(levene_groups) >= 2:
    levene_stat, levene_p = stats.levene(*levene_groups)
    print(f"Levene检验: 统计量 = {levene_stat:.3f}, p = {levene_p:.3f}")
else:
    print("分组数量不足，无法进行方差齐性检验")

# 3. 组间差异检验
print("\n=== 组间差异检验 ===")

# 检查分组数量
if len(df['group_new'].unique()) < 2:
    print("只有一个分组，无法进行组间比较")
elif len(df['group_new'].unique()) == 2:
    print("两个分组，使用t检验或Mann-Whitney U检验")
    
    # 获取两组数据
    group1_data = df[df['group_new'] == 'T1']['PFS']
    group2_data = df[df['group_new'] == 'T2']['PFS']
    
    # 修复：检查正态性条件
    t1_normal = len(group1_data) >= 3 and ('T1' in normality_results and normality_results['T1']['p_value'] > 0.05)
    t2_normal = len(group2_data) >= 3 and ('T2' in normality_results and normality_results['T2']['p_value'] > 0.05)
    
    print(f"\n正态性检查结果:")
    print(f"  T1组正态性: {'满足' if t1_normal else '不满足'}")
    print(f"  T2组正态性: {'满足' if t2_normal else '不满足'}")
    print(f"  方差齐性: {'满足' if levene_p > 0.05 else '不满足'} (p = {levene_p:.3f})")
    
    if t1_normal and t2_normal and levene_p > 0.05:
        print("\n两组数据均满足正态性且方差齐性，使用独立样本t检验")
        # 独立样本t检验
        t_stat, t_p = stats.ttest_ind(group1_data, group2_data, equal_var=True)
        print(f"独立样本t检验: t = {t_stat:.3f}, p = {t_p:.3f}")
        
        # 计算效应量（Cohen's d）
        pooled_std = np.sqrt(((len(group1_data)-1)*group1_data.var() + (len(group2_data)-1)*group2_data.var()) / (len(group1_data) + len(group2_data) - 2))
        cohens_d = (group1_data.mean() - group2_data.mean()) / pooled_std
        print(f"效应量 (Cohen's d): {cohens_d:.3f}")
        
        # 效应量解读
        if abs(cohens_d) < 0.2:
            effect_size = "极小"
        elif abs(cohens_d) < 0.5:
            effect_size = "小"
        elif abs(cohens_d) < 0.8:
            effect_size = "中等"
        else:
            effect_size = "大"
        print(f"效应量解读: {effect_size}效应")
        
        test_used = "独立样本t检验"
        test_stat = t_stat
        test_p = t_p
        effect_size_value = cohens_d
    else:
        print("\n数据不满足参数检验假设，使用Mann-Whitney U检验")
        # Mann-Whitney U检验（非参数检验）
        u_stat, u_p = stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided')
        print(f"Mann-Whitney U检验: U = {u_stat:.1f}, p = {u_p:.3f}")
        
        # 计算效应量（r = z/√n）
        n1, n2 = len(group1_data), len(group2_data)
        z_stat = stats.norm.ppf(u_p/2) if u_p < 0.5 else -stats.norm.ppf(u_p/2)
        r_effect = abs(z_stat) / np.sqrt(n1 + n2)
        print(f"效应量 (r): {r_effect:.3f}")
        
        # 效应量解读
        if r_effect < 0.1:
            effect_size = "极小"
        elif r_effect < 0.3:
            effect_size = "小"
        elif r_effect < 0.5:
            effect_size = "中等"
        else:
            effect_size = "大"
        print(f"效应量解读: {effect_size}效应")
        
        test_used = "Mann-Whitney U检验"
        test_stat = u_stat
        test_p = u_p
        effect_size_value = r_effect
else:
    print("超过两个分组，使用Kruskal-Wallis H检验（非参数）或ANOVA（参数）")
    
    # 由于样本量小且可能不满足正态性，使用非参数检验
    kruskal_stat, kruskal_p = stats.kruskal(*levene_groups)
    print(f"Kruskal-Wallis H检验: H = {kruskal_stat:.3f}, p = {kruskal_p:.3f}")
    
    test_used = "Kruskal-Wallis H检验"
    test_stat = kruskal_stat
    test_p = kruskal_p
    effect_size_value = None

# 4. 创建可视化图表
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('不同group_new分组的PFS比较分析', fontsize=20, fontweight='bold', y=1.02)

# 4.1 箱线图比较
ax1 = axes[0, 0]
sns.boxplot(x='group_new', y='PFS', data=df, ax=ax1, palette='Set2', order=sorted(df['group_new'].unique()))
sns.stripplot(x='group_new', y='PFS', data=df, ax=ax1, color='black', alpha=0.6, size=6, jitter=True, 
              order=sorted(df['group_new'].unique()))
ax1.set_xlabel('Group_new', fontsize=12, fontweight='bold')
ax1.set_ylabel('PFS (无进展生存期)', fontsize=12, fontweight='bold')
ax1.set_title('不同Group_new的PFS分布比较', fontsize=14, fontweight='bold')

# 在箱线图上添加样本数量
for i, group in enumerate(sorted(df['group_new'].unique())):
    n = len(df[df['group_new'] == group])
    median = df[df['group_new'] == group]['PFS'].median()
    ax1.text(i, df['PFS'].max() * 1.02, f'n={n}', 
             ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax1.text(i, median, f'{median:.1f}', 
             ha='center', va='center', fontsize=10, fontweight='bold', color='white')

# 4.2 小提琴图
ax2 = axes[0, 1]
sns.violinplot(x='group_new', y='PFS', data=df, ax=ax2, palette='Set2', 
               inner='quartile', order=sorted(df['group_new'].unique()))
sns.stripplot(x='group_new', y='PFS', data=df, ax=ax2, color='black', alpha=0.6, 
              size=5, jitter=True, order=sorted(df['group_new'].unique()))
ax2.set_xlabel('Group_new', fontsize=12, fontweight='bold')
ax2.set_ylabel('PFS (无进展生存期)', fontsize=12, fontweight='bold')
ax2.set_title('不同Group_new的PFS小提琴图', fontsize=14, fontweight='bold')

# 4.3 分组条形图（均值±标准差）
ax3 = axes[0, 2]
group_means = df.groupby('group_new')['PFS'].mean().reindex(sorted(df['group_new'].unique()))
group_stds = df.groupby('group_new')['PFS'].std().reindex(sorted(df['group_new'].unique()))
group_counts = df.groupby('group_new')['PFS'].count().reindex(sorted(df['group_new'].unique()))

x_pos = np.arange(len(group_means))
bars = ax3.bar(x_pos, group_means, yerr=group_stds, 
               capsize=8, alpha=0.7, color=sns.color_palette('Set2'),
               edgecolor='black', linewidth=1.5)
ax3.set_xlabel('Group_new', fontsize=12, fontweight='bold')
ax3.set_ylabel('PFS均值±标准差', fontsize=12, fontweight='bold')
ax3.set_title('各Group_new的PFS均值比较', fontsize=14, fontweight='bold')
ax3.set_xticks(x_pos)
ax3.set_xticklabels(group_means.index)

# 在柱子上方添加均值和样本数
for i, (mean, std, count) in enumerate(zip(group_means, group_stds, group_counts)):
    ax3.text(i, mean + std + 2, f'{mean:.1f}±{std:.1f}', 
             ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax3.text(i, mean/2, f'n={count}', 
             ha='center', va='center', color='white', fontsize=10, fontweight='bold')

# 4.4 累积分布函数图
ax4 = axes[1, 0]
for group in sorted(df['group_new'].unique()):
    group_data = df[df['group_new'] == group]['PFS'].sort_values()
    y = np.arange(1, len(group_data) + 1) / len(group_data)
    ax4.plot(group_data, y, marker='o', markersize=4, label=f'{group} (n={len(group_data)})', linewidth=2)
    
ax4.set_xlabel('PFS (无进展生存期)', fontsize=12, fontweight='bold')
ax4.set_ylabel('累积概率', fontsize=12, fontweight='bold')
ax4.set_title('各Group_new的PFS累积分布函数', fontsize=14, fontweight='bold')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 4.5 散点图（按分组）
ax5 = axes[1, 1]
colors = sns.color_palette('Set2', len(df['group_new'].unique()))
markers = ['o', 's']

for i, group in enumerate(sorted(df['group_new'].unique())):
    group_data = df[df['group_new'] == group]
    ax5.scatter(group_data['percent'], group_data['PFS'], 
               label=f'{group} (n={len(group_data)})', 
               color=colors[i], marker=markers[i % len(markers)], 
               s=80, alpha=0.7, edgecolors='black')

ax5.set_xlabel('Percent', fontsize=12, fontweight='bold')
ax5.set_ylabel('PFS', fontsize=12, fontweight='bold')
ax5.set_title('PFS与Percent关系（按Group_new分组）', fontsize=14, fontweight='bold')
ax5.legend()
ax5.grid(True, alpha=0.3)

# 4.6 统计摘要
ax6 = axes[1, 2]
ax6.axis('off')

# 创建统计摘要文本
summary_text = f"""
统计检验结果 (n = {len(df)})

分组信息:
  T1: n = {len(df[df['group_new']=='T1'])} 例
  T2: n = {len(df[df['group_new']=='T2'])} 例

描述性统计:
"""

for group in sorted(df['group_new'].unique()):
    group_data = df[df['group_new'] == group]['PFS']
    summary_text += f"""
  {group}组:
    均值: {group_data.mean():.1f} ± {group_data.std():.1f}
    中位数: {group_data.median():.1f}
    范围: {group_data.min():.1f} - {group_data.max():.1f}
    四分位数: {group_data.quantile(0.25):.1f}, {group_data.quantile(0.75):.1f}
"""

summary_text += f"""
正态性检验 (Shapiro-Wilk):
"""
for group, result in normality_results.items():
    if group in normality_results:  # 确保group存在于normality_results中
        sig = "***" if result['p_value'] < 0.001 else "**" if result['p_value'] < 0.01 else "*" if result['p_value'] < 0.05 else "NS"
        summary_text += f"  {group}: W = {result['statistic']:.3f}, p = {result['p_value']:.3f} {sig}\n"

if 'levene_p' in locals():
    summary_text += f"""
方差齐性检验 (Levene):
  F = {levene_stat:.3f}, p = {levene_p:.3f}
"""

summary_text += f"""
组间差异检验:
  方法: {test_used}
  统计量 = {test_stat:.3f}
  p值 = {test_p:.3f}
"""

if 'effect_size_value' in locals() and effect_size_value is not None:
    summary_text += f"  效应量 = {effect_size_value:.3f}\n"

# 添加统计显著性标记
if test_p < 0.05:
    summary_text += "\n结果: 组间存在统计学显著差异 (p < 0.05)"
    # 添加组间比较
    t1_median = df[df['group_new'] == 'T1']['PFS'].median()
    t2_median = df[df['group_new'] == 'T2']['PFS'].median()
    if t1_median > t2_median:
        summary_text += f"\nT1组PFS中位数({t1_median:.1f}) > T2组PFS中位数({t2_median:.1f})"
    else:
        summary_text += f"\nT2组PFS中位数({t2_median:.1f}) > T1组PFS中位数({t1_median:.1f})"
else:
    summary_text += "\n结果: 组间无统计学显著差异 (p ≥ 0.05)"

ax6.text(0.05, 0.95, summary_text, fontsize=10, 
         verticalalignment='top', horizontalalignment='left',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()

# 保存图表
output_filename = 'Group_new_PFS_comparison_analysis.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"\n综合图表已保存为: {output_filename}")

# 5. 创建单独的高质量箱线图
plt.figure(figsize=(10, 8))
sns.boxplot(x='group_new', y='PFS', data=df, palette='Set2', width=0.6, 
            order=sorted(df['group_new'].unique()))
sns.stripplot(x='group_new', y='PFS', data=df, color='black', alpha=0.7, 
              size=8, jitter=0.2, order=sorted(df['group_new'].unique()))

# 添加统计显著性标记
y_max = df['PFS'].max()
y_height = y_max * 1.1

if test_p < 0.05:
    sig_text = f'统计显著: p = {test_p:.3f}'
    sig_color = 'red'
else:
    sig_text = f'统计不显著: p = {test_p:.3f}'
    sig_color = 'gray'

plt.text(0.5, y_height, sig_text, 
         ha='center', va='bottom', fontsize=12, fontweight='bold', color=sig_color,
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

plt.xlabel('Group_new', fontsize=14, fontweight='bold')
plt.ylabel('PFS (无进展生存期)', fontsize=14, fontweight='bold')
plt.title('不同Group_new分组的PFS比较', fontsize=16, fontweight='bold')
plt.grid(True, alpha=0.3, axis='y')

# 在箱线图上添加样本数量和统计值
for i, group in enumerate(sorted(df['group_new'].unique())):
    group_data = df[df['group_new'] == group]['PFS']
    n = len(group_data)
    median = group_data.median()
    mean = group_data.mean()
    
    # 在箱线图上方添加样本数
    plt.text(i, y_max * 1.02, f'n={n}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # 在箱线图内部添加中位数和均值
    plt.text(i, median, f'{median:.1f}', ha='center', va='center', 
             fontsize=10, fontweight='bold', color='white')
    plt.text(i, mean, f'{mean:.1f}', ha='center', va='bottom', 
             fontsize=9, color='red')

plt.tight_layout()
boxplot_filename = 'Group_new_PFS_boxplot.png'
plt.savefig(boxplot_filename, dpi=300, bbox_inches='tight')
print(f"箱线图已保存为: {boxplot_filename}")

# 6. 详细分析报告
print("\n" + "="*80)
print("不同Group_new分组PFS差异分析报告")
print("="*80)

print("\n一、 样本分布:")
for group in sorted(df['group_new'].unique()):
    group_data = df[df['group_new'] == group]
    print(f"  {group}组: {len(group_data)}例患者")
    print(f"    PFS范围: {group_data['PFS'].min():.2f} - {group_data['PFS'].max():.2f}")
    print(f"    PFS中位数: {group_data['PFS'].median():.2f}")
    print(f"    PFS均值±标准差: {group_data['PFS'].mean():.2f} ± {group_data['PFS'].std():.2f}")

print("\n二、 正态性检验结果:")
normality_passed = all('p_value' in normality_results.get(group, {}) and 
                       normality_results[group]['p_value'] > 0.05 
                       for group in df['group_new'].unique())
if normality_passed:
    print("  两组数据的PFS均符合正态分布 (p > 0.05)")
else:
    print("  部分分组的PFS数据不符合正态分布，建议使用非参数检验")

print("\n三、 方差齐性检验结果:")
if 'levene_p' in locals():
    if levene_p > 0.05:
        print(f"  两组方差齐性 (Levene检验: p = {levene_p:.3f})")
    else:
        print(f"  方差不齐 (Levene检验: p = {levene_p:.3f})，建议使用非参数检验")

print("\n四、 组间差异检验结果:")
print(f"  检验方法: {test_used}")
print(f"  统计量: {test_stat:.3f}")
print(f"  p值: {test_p:.3f}")
if 'effect_size_value' in locals() and effect_size_value is not None:
    print(f"  效应量: {effect_size_value:.3f}")

print("\n五、 临床意义解读:")
if test_p < 0.05:
    print(f"  不同Group_new分组的PFS存在统计学显著差异 (p = {test_p:.3f})")
    
    # 比较两组中位数
    t1_median = df[df['group_new'] == 'T1']['PFS'].median()
    t2_median = df[df['group_new'] == 'T2']['PFS'].median()
    median_diff = abs(t1_median - t2_median)
    
    print(f"\n  临床差异分析:")
    print(f"  - T1组PFS中位数: {t1_median:.1f}")
    print(f"  - T2组PFS中位数: {t2_median:.1f}")
    print(f"  - 中位数差异: {median_diff:.1f}")
    
    if t1_median > t2_median:
        print(f"  - T1组的PFS比T2组长{median_diff:.1f}个月")
        print(f"  - T1组可能具有更好的预后")
    else:
        print(f"  - T2组的PFS比T1组长{median_diff:.1f}个月")
        print(f"  - T2组可能具有更好的预后")
        
    if 'effect_size_value' in locals() and effect_size_value is not None:
        if abs(effect_size_value) < 0.2:
            effect_strength = "极小"
        elif abs(effect_size_value) < 0.5:
            effect_strength = "小"
        elif abs(effect_size_value) < 0.8:
            effect_strength = "中等"
        else:
            effect_strength = "大"
        print(f"  - 效应量大小为: {effect_strength}")
    
else:
    print(f"  不同Group_new分组的PFS没有统计学显著差异 (p = {test_p:.3f})")
    print("  两组间的PFS差异可能由随机变异引起")
    print("  Group_new分组在当前数据中可能不是PFS的重要预测因素")

print("\n六、 研究局限性:")
print("  1. 样本量相对较小，特别是T2组只有10例")
print("  2. PFS数据分布可能存在偏态")
print("  3. 未考虑其他可能的混杂因素（如年龄、性别、治疗方案等）")
print("  4. 分组基于percent阈值，可能存在主观性")

print("\n七、 建议:")
if test_p < 0.05:
    print("  1. 进一步验证Group_new分组在独立队列中的预测价值")
    print("  2. 探索Group_new与其他临床病理特征的关系")
    print("  3. 进行多变量分析，控制可能的混杂因素")
    print("  4. 考虑进行生存分析，计算风险比(HR)")
else:
    print("  1. 增加样本量以提高统计功效")
    print("  2. 考虑其他分组方法或阈值")
    print("  3. 探索其他可能影响PFS的变量")
    print("  4. 结合临床专业知识重新评估分组策略")

print("\n八、 生成文件:")
print(f"  1. {output_filename} - 综合分析图表")
print(f"  2. {boxplot_filename} - 分组箱线图")

# 保存统计结果到CSV
summary_df = df.groupby('group_new')['PFS'].describe()
summary_df['count'] = df.groupby('group_new').size()
summary_df = summary_df[['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']]
summary_df.columns = ['样本数', '均值', '标准差', '最小值', '25分位数', '中位数', '75分位数', '最大值']
summary_df.to_csv('Group_new_PFS_statistics.csv')
print(f"  3. Group_new_PFS_statistics.csv - 详细统计表格")

print("\n" + "="*80)
print("分析完成！")

# 显示图表
plt.show()
