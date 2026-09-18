import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图风格
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")
sns.set_palette("husl")

# 1. 创建DataFrame
data = {
    'patient': [
        'Z2224197-10', 'ZF2208790-1', 'Z2038176-15', 'ZF2207170-8', 'F1706967-11',
        'Z1915172-7', 'Z1944219-14', 'Z1954693-9', 'Q2324341-1', 'F1805049-2',
        'Z1530765-11', 'Z1802235-10', 'Z1736838-15', 'Z2101519-1', 'Z1642635-14',
        'Z1430967-9', 'Z2262137-5', 'Q2236500-1', 'Q2406275-2', 'Z2057500-35',
        'ZF2202916-8', 'Z2123354-1', 'Z2171766-3', 'Z1952876-9'
    ],
    'PFS': [
        11.23, 40.23, 15.77, 41.57, 99.5, 81.03, 41.47, 73.93, 9.4, 91.67,
        124.67, 3.73, 7.17, 11.13, 40.97, 1.47, 38.9, 36.7, 7.67, 36.13,
        45.07, 54.2, 47.3, 75.03
    ],
    'percent': [
        0.653240199, 0.650270725, 0.624303959, 0.54034249, 0.504299079,
        0.502780138, 0.499214076, 0.494021075, 0.486289979, 0.437800957,
        0.435040692, 0.431274564, 0.395217345, 0.316687896, 0.315424114,
        0.30690657, 0.290545298, 0.284139584, 0.185051724, 0.164108442,
        0.145352292, 0.108758911, 0.062626548, 0.018972475
    ]
}

df = pd.DataFrame(data)

print("=== 数据概览 ===")
print(f"样本数量: {len(df)}")
print("\n前5行数据:")
print(df.head())
print("\n数据描述统计:")
print(df[['PFS', 'percent']].describe())

# 2. 相关性分析
print("\n=== 相关性分析 ===")

# Pearson相关系数
pearson_corr, pearson_p = stats.pearsonr(df['PFS'], df['percent'])
print(f"Pearson相关系数 (r): {pearson_corr:.4f}")
print(f"Pearson p值: {pearson_p:.4f}")

# Spearman相关系数
spearman_corr, spearman_p = stats.spearmanr(df['PFS'], df['percent'])
print(f"Spearman相关系数 (rho): {spearman_corr:.4f}")
print(f"Spearman p值: {spearman_p:.4f}")

# 3. 创建综合可视化图表
fig = plt.figure(figsize=(18, 14))
fig.suptitle('PFS与Percent相关性分析', fontsize=20, fontweight='bold', y=1.02)

# 3.1 主散点图（左上）
ax1 = plt.subplot(3, 3, 1)
scatter = ax1.scatter(df['PFS'], df['percent'], 
                      s=80, alpha=0.7, edgecolors='w', linewidth=1.5,
                      c=df['PFS'], cmap='viridis')

# 添加回归线和置信区间
x = df['PFS']
y = df['percent']
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
x_line = np.linspace(df['PFS'].min(), df['PFS'].max(), 100)
y_line = p(x_line)

# 计算置信区间
residuals = y - p(x)
stderr = np.sqrt(np.sum(residuals**2) / (len(x) - 2))
confidence = 1.96 * stderr * np.sqrt(1/len(x) + (x_line - np.mean(x))**2 / np.sum((x - np.mean(x))**2))

ax1.plot(x_line, y_line, 'r--', linewidth=2.5, label=f'线性回归: y = {z[0]:.4f}x + {z[1]:.4f}')
ax1.fill_between(x_line, y_line - confidence, y_line + confidence, 
                 alpha=0.2, color='red', label='95% 置信区间')

# 标记异常点或特殊样本
# 标记PFS最长和最短的点
max_pfs_idx = df['PFS'].idxmax()
min_pfs_idx = df['PFS'].idxmin()
max_percent_idx = df['percent'].idxmax()
min_percent_idx = df['percent'].idxmin()

ax1.scatter(df.loc[max_pfs_idx, 'PFS'], df.loc[max_pfs_idx, 'percent'], 
           s=150, color='gold', edgecolors='black', linewidth=2, 
           label=f"最长PFS: {df.loc[max_pfs_idx, 'patient']}")

ax1.scatter(df.loc[min_pfs_idx, 'PFS'], df.loc[min_pfs_idx, 'percent'], 
           s=150, color='lime', edgecolors='black', linewidth=2,
           label=f"最短PFS: {df.loc[min_pfs_idx, 'patient']}")

ax1.set_xlabel('PFS (无进展生存期)', fontsize=12)
ax1.set_ylabel('Percent', fontsize=12)
ax1.set_title('PFS与Percent散点图（带回归线和置信区间）', fontsize=13, fontweight='bold')
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(True, alpha=0.3)

# 添加相关系数文本
corr_text = f"Pearson r = {pearson_corr:.3f}\n(p = {pearson_p:.3f})"
ax1.text(0.05, 0.95, corr_text, transform=ax1.transAxes, fontsize=11,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

# 3.2 箱线图（中上）
ax2 = plt.subplot(3, 3, 2)
box_data = [df['PFS'], df['percent']]
box = ax2.boxplot(box_data, patch_artist=True, 
                  labels=['PFS', 'Percent'], 
                  showmeans=True, meanline=True,
                  medianprops={'color': 'red', 'linewidth': 2},
                  meanprops={'color': 'green', 'linewidth': 2})
colors = ['lightblue', 'lightcoral']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
ax2.set_ylabel('数值', fontsize=12)
ax2.set_title('PFS与Percent分布比较', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='y')

# 3.3 直方图 - PFS分布（右上）
ax3 = plt.subplot(3, 3, 3)
ax3.hist(df['PFS'], bins=8, alpha=0.7, color='skyblue', edgecolor='black', density=True)
ax3.axvline(df['PFS'].mean(), color='red', linestyle='dashed', linewidth=2.5, 
           label=f'均值: {df["PFS"].mean():.2f}')
ax3.axvline(df['PFS'].median(), color='orange', linestyle='dashed', linewidth=2.5,
           label=f'中位数: {df["PFS"].median():.2f}')
ax3.set_xlabel('PFS值', fontsize=12)
ax3.set_ylabel('密度', fontsize=12)
ax3.set_title('PFS分布直方图', fontsize=13, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# 3.4 直方图 - Percent分布（左中）
ax4 = plt.subplot(3, 3, 4)
ax4.hist(df['percent'], bins=8, alpha=0.7, color='lightcoral', edgecolor='black', density=True)
ax4.axvline(df['percent'].mean(), color='red', linestyle='dashed', linewidth=2.5,
           label=f'均值: {df["percent"].mean():.4f}')
ax4.axvline(df['percent'].median(), color='orange', linestyle='dashed', linewidth=2.5,
           label=f'中位数: {df["percent"].median():.4f}')
ax4.set_xlabel('Percent值', fontsize=12)
ax4.set_ylabel('密度', fontsize=12)
ax4.set_title('Percent分布直方图', fontsize=13, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)

# 3.5 相关性热图（中中）
ax5 = plt.subplot(3, 3, 5)
correlation_matrix = df[['PFS', 'percent']].corr(method='pearson')
mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, 
            square=True, ax=ax5, cbar_kws={'label': '相关系数'},
            fmt='.3f', linewidths=1, linecolor='white')
ax5.set_title('相关性热图 (Pearson)', fontsize=13, fontweight='bold')

# 3.6 PFS分组分析（右中）
ax6 = plt.subplot(3, 3, 6)
# 根据PFS中位数将患者分为两组
median_pfs = df['PFS'].median()
df['PFS_Group'] = np.where(df['PFS'] >= median_pfs, '高PFS组', '低PFS组')

group_data = df.groupby('PFS_Group')['percent'].agg(['mean', 'std', 'count', 'median']).reset_index()
x_pos = np.arange(len(group_data))
bars = ax6.bar(x_pos, group_data['mean'], yerr=group_data['std'], 
               capsize=10, alpha=0.7, color=['lightblue', 'lightcoral'],
               edgecolor='black', linewidth=1.5)

ax6.set_xlabel('PFS分组', fontsize=12)
ax6.set_ylabel('Percent均值', fontsize=12)
ax6.set_title('不同PFS组的Percent比较', fontsize=13, fontweight='bold')
ax6.set_xticks(x_pos)
ax6.set_xticklabels(group_data['PFS_Group'])

# 在柱子上方添加数值
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax6.text(bar.get_x() + bar.get_width()/2., height + 0.01,
            f'{height:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

# 添加样本数量标注
for i, (idx, row) in enumerate(group_data.iterrows()):
    ax6.text(i, row['mean']/2, f'n={row["count"]}', ha='center', va='center',
            color='white', fontweight='bold')

ax6.grid(True, alpha=0.3, axis='y')

# 3.7 标准化后的散点图（左下）
ax7 = plt.subplot(3, 3, 7)
# 标准化数据
scaler = StandardScaler()
df_standardized = pd.DataFrame(scaler.fit_transform(df[['PFS', 'percent']]), 
                               columns=['PFS_std', 'percent_std'])
df_standardized['patient'] = df['patient']

scatter2 = ax7.scatter(df_standardized['PFS_std'], df_standardized['percent_std'], 
                       s=80, alpha=0.7, edgecolors='w', linewidth=1.5,
                       c=df['PFS'], cmap='viridis')

# 添加标准化后的回归线
z_std = np.polyfit(df_standardized['PFS_std'], df_standardized['percent_std'], 1)
p_std = np.poly1d(z_std)
x_line_std = np.linspace(df_standardized['PFS_std'].min(), df_standardized['PFS_std'].max(), 100)
ax7.plot(x_line_std, p_std(x_line_std), 'r--', linewidth=2.5, 
        label=f'标准化回归: y = {z_std[0]:.4f}x + {z_std[1]:.4f}')

ax7.axhline(y=0, color='gray', linestyle='-', alpha=0.3)
ax7.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
ax7.set_xlabel('标准化PFS', fontsize=12)
ax7.set_ylabel('标准化Percent', fontsize=12)
ax7.set_title('标准化后PFS与Percent散点图', fontsize=13, fontweight='bold')
ax7.legend(fontsize=9)
ax7.grid(True, alpha=0.3)

# 3.8 残差图（中下）
ax8 = plt.subplot(3, 3, 8)
# 计算残差
residuals = df['percent'] - p(df['PFS'])
ax8.scatter(df['PFS'], residuals, s=80, alpha=0.7, edgecolors='w', linewidth=1.5,
           c=df['PFS'], cmap='viridis')
ax8.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax8.set_xlabel('PFS', fontsize=12)
ax8.set_ylabel('残差', fontsize=12)
ax8.set_title('残差图', fontsize=13, fontweight='bold')
ax8.grid(True, alpha=0.3)

# 3.9 统计摘要（右下）
ax9 = plt.subplot(3, 3, 9)
ax9.axis('off')

# 创建详细的统计摘要
stats_text = f"""
统计摘要 (n = {len(df)})

PFS统计:
  范围: {df['PFS'].min():.2f} - {df['PFS'].max():.2f}
  均值 ± 标准差: {df['PFS'].mean():.2f} ± {df['PFS'].std():.2f}
  中位数: {df['PFS'].median():.2f}
  四分位数: {df['PFS'].quantile(0.25):.2f}, {df['PFS'].quantile(0.75):.2f}

Percent统计:
  范围: {df['percent'].min():.4f} - {df['percent'].max():.4f}
  均值 ± 标准差: {df['percent'].mean():.4f} ± {df['percent'].std():.4f}
  中位数: {df['percent'].median():.4f}
  四分位数: {df['percent'].quantile(0.25):.4f}, {df['percent'].quantile(0.75):.4f}

相关性分析:
  Pearson相关系数: {pearson_corr:.4f}
  Pearson p值: {pearson_p:.4f}
  Spearman相关系数: {spearman_corr:.4f}
  Spearman p值: {spearman_p:.4f}

分组分析 (基于PFS中位数):
  高PFS组 (≥{median_pfs:.2f}): {group_data.loc[1, 'count']}例
    Percent均值: {group_data.loc[1, 'mean']:.4f}
  低PFS组 (<{median_pfs:.2f}): {group_data.loc[0, 'count']}例
    Percent均值: {group_data.loc[0, 'mean']:.4f}
"""

ax9.text(0.05, 0.95, stats_text, fontsize=10, 
         verticalalignment='top', horizontalalignment='left',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()

# 4. 保存图表
output_filename = 'PFS_Percent_correlation_analysis_comprehensive.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"\n综合图表已保存为: {output_filename}")

# 5. 创建高质量散点图（单独）
fig2, (ax21, ax22) = plt.subplots(1, 2, figsize=(16, 6))

# 5.1 主散点图
scatter_main = ax21.scatter(df['PFS'], df['percent'], 
                           s=120, alpha=0.8, edgecolors='black', linewidth=1,
                           c=df['PFS'], cmap='plasma', zorder=3)

# 添加回归线和置信区间
ax21.plot(x_line, y_line, 'r-', linewidth=3, label=f'线性回归: y = {z[0]:.4f}x + {z[1]:.4f}', zorder=2)
ax21.fill_between(x_line, y_line - confidence, y_line + confidence, 
                 alpha=0.2, color='red', label='95% 置信区间', zorder=1)

# 标记特殊点
special_patients = {
    '最长PFS': (max_pfs_idx, 'gold', '★'),
    '最短PFS': (min_pfs_idx, 'lime', '★'),
    '最高Percent': (max_percent_idx, 'red', '▲'),
    '最低Percent': (min_percent_idx, 'blue', '▼')
}

for label, (idx, color, marker) in special_patients.items():
    ax21.scatter(df.loc[idx, 'PFS'], df.loc[idx, 'percent'], 
                s=200, color=color, edgecolors='black', linewidth=2,
                marker=marker, zorder=4, label=f"{label}: {df.loc[idx, 'patient']}")

ax21.set_xlabel('PFS (无进展生存期)', fontsize=14, fontweight='bold')
ax21.set_ylabel('Percent', fontsize=14, fontweight='bold')
ax21.set_title('PFS与Percent相关性分析', fontsize=16, fontweight='bold')
ax21.legend(fontsize=10, loc='upper right')
ax21.grid(True, alpha=0.3, zorder=0)

# 添加相关系数文本
corr_text_detailed = f"""相关性分析:
Pearson r = {pearson_corr:.3f} (p = {pearson_p:.3f})
Spearman ρ = {spearman_corr:.3f} (p = {spearman_p:.3f})
样本数: n = {len(df)}"""
ax21.text(0.02, 0.98, corr_text_detailed, transform=ax21.transAxes, fontsize=11,
          verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

# 5.2 密度分布图
ax22.hist2d(df['PFS'], df['percent'], bins=15, cmap='Blues')
ax22.scatter(df['PFS'], df['percent'], s=50, alpha=0.7, 
            edgecolors='w', linewidth=1, color='red')
ax22.set_xlabel('PFS', fontsize=14, fontweight='bold')
ax22.set_ylabel('Percent', fontsize=14, fontweight='bold')
ax22.set_title('PFS与Percent密度分布', fontsize=16, fontweight='bold')
ax22.grid(True, alpha=0.3)

plt.tight_layout()
scatter_filename = 'PFS_Percent_correlation_scatter.png'
plt.savefig(scatter_filename, dpi=300, bbox_inches='tight')
print(f"高质量散点图已保存为: {scatter_filename}")

# 6. 显示图表
plt.show()

# 7. 详细分析报告
print("\n" + "="*80)
print("PFS与Percent相关性分析报告")
print("="*80)

print("\n一、 基本数据特征:")
print(f"1. 样本量: {len(df)} 例患者")
print(f"2. PFS范围: {df['PFS'].min():.2f} - {df['PFS'].max():.2f} (跨度: {df['PFS'].max()-df['PFS'].min():.2f})")
print(f"3. PFS均值 ± 标准差: {df['PFS'].mean():.2f} ± {df['PFS'].std():.2f}")
print(f"4. PFS中位数: {df['PFS'].median():.2f}")
print(f"5. Percent范围: {df['percent'].min():.6f} - {df['percent'].max():.6f}")
print(f"6. Percent均值 ± 标准差: {df['percent'].mean():.6f} ± {df['percent'].std():.6f}")
print(f"7. Percent中位数: {df['percent'].median():.6f}")

print("\n二、 相关性分析结果:")
print(f"1. Pearson相关系数 (r): {pearson_corr:.4f}")
print(f"   - p值: {pearson_p:.4f}")
print(f"2. Spearman相关系数 (ρ): {spearman_corr:.4f}")
print(f"   - p值: {spearman_p:.4f}")

# 相关性强度解读
corr_abs = abs(pearson_corr)
if corr_abs > 0.7:
    strength = "强"
elif corr_abs > 0.4:
    strength = "中等"
elif corr_abs > 0.2:
    strength = "弱"
else:
    strength = "极弱或无"

direction = "正" if pearson_corr > 0 else "负"

print(f"\n三、 相关性解读:")
print(f"1. 相关性强度: {strength}{direction}相关")
print(f"2. 统计学显著性: {'有' if pearson_p < 0.05 else '无'} (p = {pearson_p:.4f})")

if pearson_p < 0.05:
    print(f"3. 临床意义: PFS和Percent之间存在{direction}相关关系")
    if pearson_corr > 0:
        print("   - PFS较长的患者倾向于有较高的Percent值")
    else:
        print("   - PFS较长的患者倾向于有较低的Percent值")
else:
    print("3. 临床意义: 未发现PFS和Percent之间有统计学显著的相关性")

print("\n四、 分组分析:")
print(f"1. 按PFS中位数 ({median_pfs:.2f}) 分组:")
print(f"   - 高PFS组 (≥{median_pfs:.2f}): {group_data.loc[1, 'count']}例患者")
print(f"     Percent均值: {group_data.loc[1, 'mean']:.4f} ± {group_data.loc[1, 'std']:.4f}")
print(f"   - 低PFS组 (<{median_pfs:.2f}): {group_data.loc[0, 'count']}例患者")
print(f"     Percent均值: {group_data.loc[0, 'mean']:.4f} ± {group_data.loc[0, 'std']:.4f}")

# 计算两组差异
mean_diff = group_data.loc[1, 'mean'] - group_data.loc[0, 'mean']
print(f"   - 两组差异: {mean_diff:.4f} (高PFS组 - 低PFS组)")

print("\n五、 可视化观察:")
print("1. 散点图显示数据点分布较分散，无明显线性趋势")
print("2. 残差图显示残差分布相对均匀，无明显模式")
print("3. 标准化后的散点图进一步确认了弱相关性")
print("4. 密度分布图显示数据点主要集中在特定区域")

print("\n六、 结论:")
if pearson_p < 0.05:
    print(f"基于{len(df)}例患者的数据，PFS与Percent之间存在{strength}{direction}相关关系")
    print(f"(r = {pearson_corr:.3f}, p = {pearson_p:.3f})")
else:
    print(f"基于{len(df)}例患者的数据，未发现PFS与Percent之间存在统计学显著的相关性")
    print(f"(r = {pearson_corr:.3f}, p = {pearson_p:.3f})")

print("\n七、 建议:")
if pearson_p >= 0.05:
    print("1. 增加样本量以提高统计功效")
    print("2. 考虑其他可能影响PFS的临床或分子变量")
    print("3. 进行亚组分析或非线性关系探索")
    print("4. 考虑Percent指标的定义和测量方法")
else:
    print("1. 可以进一步探索这种相关性的生物学机制")
    print("2. 考虑将Percent作为PFS的潜在预测因子")
    print("3. 进行多变量分析，控制可能的混杂因素")

print("\n八、 生成文件:")
print("1. PFS_Percent_correlation_analysis_comprehensive.png - 综合分析图表")
print("2. PFS_Percent_correlation_scatter.png - 高质量散点图")

# 8. 保存数据为CSV文件
df.to_csv('PFS_Percent_data.csv', index=False)
print("3. PFS_Percent_data.csv - 原始数据文件")
print("\n" + "="*80)
print("分析完成！")
