import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib import rcParams

# 设置中文字体（如果需要显示中文）
rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 1. 创建DataFrame
data = {
    'patient': [
        'Z1642635-14', 'Z1736838-15', 'Z1530765-11', 'Q2324341-1',
        'Z2224197-10', 'Z2038176-15', 'Z1954693-9', 'Z1952876-9',
        'Z1915172-7', 'Z1802235-10', 'Z1944219-14', 'F1805049-2',
        'F1706967-11'
    ],
    'OS': [
        109.47, 11.83, 124.67, 17.83, 30.7, 63.23, 73.93, 75.03,
        81.03, 8.1, 90.77, 91.67, 99.5
    ],
    'Percent': [
        0.315424, 0.395217, 0.435041, 0.48629, 0.65324, 0.624304,
        0.494021, 0.018972, 0.50278, 0.431275, 0.499214, 0.437801,
        0.504299
    ]
}

df = pd.DataFrame(data)

# 2. 基本统计分析
print("=== 基本统计信息 ===")
print(f"样本数量: {len(df)}")
print(f"OS范围: {df['OS'].min():.2f} - {df['OS'].max():.2f}")
print(f"OS均值±标准差: {df['OS'].mean():.2f} ± {df['OS'].std():.2f}")
print(f"Percent范围: {df['Percent'].min():.4f} - {df['Percent'].max():.4f}")
print(f"Percent均值±标准差: {df['Percent'].mean():.4f} ± {df['Percent'].std():.4f}")

# 3. 计算相关性
print("\n=== 相关性分析 ===")

# Pearson相关系数
pearson_corr, pearson_p = stats.pearsonr(df['OS'], df['Percent'])
print(f"Pearson相关系数 (r): {pearson_corr:.4f}")
print(f"Pearson p值: {pearson_p:.4f}")

# Spearman相关系数
spearman_corr, spearman_p = stats.spearmanr(df['OS'], df['Percent'])
print(f"Spearman相关系数 (rho): {spearman_corr:.4f}")
print(f"Spearman p值: {spearman_p:.4f}")

# 4. 创建可视化图表
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('OS与Percent相关性分析', fontsize=16, fontweight='bold', y=1.02)

# 4.1 散点图
ax1 = axes[0, 0]
scatter = ax1.scatter(df['OS'], df['Percent'], alpha=0.7, edgecolors='w', s=100, 
                      c=df['OS'], cmap='viridis')

# 添加回归线
z = np.polyfit(df['OS'], df['Percent'], 1)
p = np.poly1d(z)
x_line = np.linspace(df['OS'].min(), df['OS'].max(), 100)
ax1.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2, label=f'线性拟合: y={z[0]:.6f}x+{z[1]:.4f}')

# 标记每个点
for i, row in df.iterrows():
    ax1.annotate(row['patient'], (row['OS'], row['Percent']), 
                xytext=(5, 5), textcoords='offset points', fontsize=8)

ax1.set_xlabel('OS (总体生存期)', fontsize=12)
ax1.set_ylabel('Percent (百分比)', fontsize=12)
ax1.set_title('散点图与线性回归', fontsize=13)
ax1.legend()
ax1.grid(True, alpha=0.3)

# 添加颜色条
cbar = plt.colorbar(scatter, ax=ax1)
cbar.set_label('OS值', fontsize=10)

# 4.2 箱线图比较
ax2 = axes[0, 1]
data_to_plot = [df['OS'], df['Percent']]
box = ax2.boxplot(data_to_plot, patch_artist=True, 
                  labels=['OS', 'Percent'])
colors = ['lightblue', 'lightcoral']
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
ax2.set_ylabel('值', fontsize=12)
ax2.set_title('变量分布比较', fontsize=13)
ax2.grid(True, alpha=0.3, axis='y')

# 4.3 直方图 - OS分布
ax3 = axes[0, 2]
ax3.hist(df['OS'], bins=6, alpha=0.7, color='skyblue', edgecolor='black')
ax3.axvline(df['OS'].mean(), color='red', linestyle='dashed', linewidth=2, 
           label=f'均值: {df["OS"].mean():.2f}')
ax3.set_xlabel('OS值', fontsize=12)
ax3.set_ylabel('频数', fontsize=12)
ax3.set_title('OS分布直方图', fontsize=13)
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4.4 直方图 - Percent分布
ax4 = axes[1, 0]
ax4.hist(df['Percent'], bins=6, alpha=0.7, color='lightcoral', edgecolor='black')
ax4.axvline(df['Percent'].mean(), color='blue', linestyle='dashed', linewidth=2, 
           label=f'均值: {df["Percent"].mean():.4f}')
ax4.set_xlabel('Percent值', fontsize=12)
ax4.set_ylabel('频数', fontsize=12)
ax4.set_title('Percent分布直方图', fontsize=13)
ax4.legend()
ax4.grid(True, alpha=0.3)

# 4.5 相关性热图
ax5 = axes[1, 1]
correlation_matrix = df[['OS', 'Percent']].corr(method='pearson')
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, 
            square=True, ax=ax5, cbar_kws={'label': '相关系数'})
ax5.set_title('相关性热图 (Pearson)', fontsize=13)

# 4.6 生存期分组分析
ax6 = axes[1, 2]
# 根据OS中位数将患者分为两组
median_os = df['OS'].median()
df['OS_Group'] = np.where(df['OS'] >= median_os, '高OS组', '低OS组')

group_data = df.groupby('OS_Group')['Percent'].agg(['mean', 'std', 'count']).reset_index()
x_pos = np.arange(len(group_data))
bars = ax6.bar(x_pos, group_data['mean'], yerr=group_data['std'], 
               capsize=10, alpha=0.7, color=['lightblue', 'lightcoral'])

ax6.set_xlabel('OS分组', fontsize=12)
ax6.set_ylabel('Percent均值', fontsize=12)
ax6.set_title('不同OS组的Percent比较', fontsize=13)
ax6.set_xticks(x_pos)
ax6.set_xticklabels(group_data['OS_Group'])

# 在柱子上方添加数值
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax6.text(bar.get_x() + bar.get_width()/2., height + 0.01,
            f'{height:.3f}', ha='center', va='bottom', fontsize=10)

ax6.grid(True, alpha=0.3, axis='y')

# 5. 添加统计摘要文本
stats_text = f"""
统计摘要:
样本数: {len(df)}
OS范围: {df['OS'].min():.1f} - {df['OS'].max():.1f}
OS中位数: {df['OS'].median():.1f}
Percent范围: {df['Percent'].min():.3f} - {df['Percent'].max():.3f}
Pearson r = {pearson_corr:.3f} (p = {pearson_p:.3f})
Spearman ρ = {spearman_corr:.3f} (p = {spearman_p:.3f})
"""

fig.text(0.02, 0.02, stats_text, fontsize=9, 
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
         verticalalignment='bottom', horizontalalignment='left')

plt.tight_layout()

# 6. 保存图表
output_filename = 'OS_Percent_correlation_analysis.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"\n图表已保存为: {output_filename}")

# 7. 创建详细的相关性分析报告
print("\n=== 详细分析报告 ===")
print("1. 数据描述:")
print(f"   - OS数据呈现{len(df)}例患者的总体生存期信息")
print(f"   - OS的分布范围较广: {df['OS'].min():.1f}-{df['OS'].max():.1f}")
print(f"   - Percent指标范围: {df['Percent'].min():.1%}-{df['Percent'].max():.1%}")

print("\n2. 相关性分析结果:")
if abs(pearson_corr) > 0.7:
    corr_strength = "强相关"
elif abs(pearson_corr) > 0.4:
    corr_strength = "中等相关"
elif abs(pearson_corr) > 0.2:
    corr_strength = "弱相关"
else:
    corr_strength = "几乎不相关"

if pearson_corr > 0:
    direction = "正相关"
else:
    direction = "负相关"

print(f"   - Pearson相关分析显示{direction}关系 ({corr_strength})")
print(f"   - Pearson相关系数: r = {pearson_corr:.3f}, p = {pearson_p:.3f}")

if pearson_p < 0.05:
    print("   - 结果具有统计学显著性 (p < 0.05)")
else:
    print("   - 结果未达到统计学显著性 (p ≥ 0.05)")

print("\n3. 临床意义分析:")
print("   - 散点图显示数据点分布较为分散，没有明显的线性模式")
print("   - 高OS组和低OS组的Percent均值差异较小")
print("   - 建议进一步分析亚组或考虑其他临床变量的影响")

# 8. 创建单独的高质量散点图
fig2, ax2 = plt.subplots(figsize=(10, 8))

# 使用seaborn的regplot添加回归线和置信区间
sns.regplot(x='OS', y='Percent', data=df, 
            scatter_kws={'s': 100, 'alpha': 0.7, 'edgecolors': 'w'},
            line_kws={'color': 'red', 'linewidth': 2, 'linestyle': '--'},
            ci=95,  # 95%置信区间
            ax=ax2)

# 标记异常点或特殊样本
# 标记Percent最低的点
min_percent_idx = df['Percent'].idxmin()
ax2.scatter(df.loc[min_percent_idx, 'OS'], df.loc[min_percent_idx, 'Percent'], 
           s=150, color='gold', edgecolors='black', linewidth=2, 
           label=f"最低Percent: {df.loc[min_percent_idx, 'patient']}")

# 标记OS最长的点
max_os_idx = df['OS'].idxmax()
ax2.scatter(df.loc[max_os_idx, 'OS'], df.loc[max_os_idx, 'Percent'], 
           s=150, color='lime', edgecolors='black', linewidth=2,
           label=f"最长OS: {df.loc[max_os_idx, 'patient']}")

# 标记OS最短的点
min_os_idx = df['OS'].idxmin()
ax2.scatter(df.loc[min_os_idx, 'OS'], df.loc[min_os_idx, 'Percent'], 
           s=150, color='orange', edgecolors='black', linewidth=2,
           label=f"最短OS: {df.loc[min_os_idx, 'patient']}")

ax2.set_xlabel('总体生存期 (OS)', fontsize=14)
ax2.set_ylabel('Percent指标', fontsize=14)
ax2.set_title('OS与Percent相关性分析 (带95%置信区间)', fontsize=16, fontweight='bold')

# 添加相关系数文本
text_box = f"Pearson r = {pearson_corr:.3f}\np = {pearson_p:.4f}\nn = {len(df)}"
ax2.text(0.05, 0.95, text_box, transform=ax2.transAxes, fontsize=12,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax2.legend(loc='upper right', fontsize=10)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
scatter_filename = 'OS_Percent_scatter_with_CI.png'
plt.savefig(scatter_filename, dpi=300, bbox_inches='tight')
print(f"高质量散点图已保存为: {scatter_filename}")

# 9. 显示图表
plt.show()

# 10. 输出数据摘要
print("\n=== 数据摘要 ===")
print(df.to_string(index=False))

# 11. 保存数据为CSV文件
df.to_csv('OS_Percent_data.csv', index=False)
print(f"\n数据已保存为: OS_Percent_data.csv")

print("\n=== 分析完成 ===")
