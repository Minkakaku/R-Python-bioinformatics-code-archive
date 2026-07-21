import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# 设置图形样式和颜色
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 设置颜色方案
zone_colors = {'inner': '#2E86AB', 'outer': '#A23B72'}
cn_colors = sns.color_palette("husl", 3)

# 创建输出目录
output_dir = "CN_prop_mpc"
os.makedirs(output_dir, exist_ok=True)

# 1. 加载数据
print("Loading data...")
adata = sc.read_h5ad("extracted_samples.h5ad")
print(f"Data shape: {adata.shape}")

# 检查必要的列是否存在
required_cols = ['sample', 'mucus_zone', 'CNs_10']
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"Column '{col}' not found in adata.obs")

# 2. 数据探索
print(f"\nData exploration:")
print(f"Number of samples: {adata.obs['sample'].nunique()}")
print(f"Unique mucus_zones: {sorted(adata.obs['mucus_zone'].dropna().unique())}")
print(f"Unique CNs_10 values: {sorted(adata.obs['CNs_10'].dropna().unique())}")

# 检查是否有inner和outer zones
zones = adata.obs['mucus_zone'].dropna().unique()
if not all(z in zones for z in ['inner', 'outer']):
    print(f"Warning: Expected 'inner' and 'outer' zones, but found: {zones}")
    # 使用实际存在的zone
    zones = [z for z in ['inner', 'outer'] if z in zones]
else:
    zones = ['inner', 'outer']

# 3. 数据准备
# 获取所有CN类型
cn_types = sorted(adata.obs['CNs_10'].dropna().unique())
print(f"CN types found: {cn_types}")

# 创建数据框用于分析
data = adata.obs[['sample', 'mucus_zone', 'CNs_10']].copy()
data = data.dropna(subset=['mucus_zone', 'CNs_10'])

# 只保留inner和outer zones
data = data[data['mucus_zone'].isin(['inner', 'outer'])]

# 4. 计算每个样本每个zone的CN比例
print("\nCalculating CN proportions for each sample and zone...")

# 创建结果数据框
results = []

for sample in data['sample'].unique():
    sample_data = data[data['sample'] == sample]
    
    for zone in zones:
        zone_data = sample_data[sample_data['mucus_zone'] == zone]
        
        if len(zone_data) == 0:
            continue
        
        # 计算每个CN类型的比例
        for cn in cn_types:
            cn_count = (zone_data['CNs_10'] == cn).sum()
            proportion = cn_count / len(zone_data)
            
            results.append({
                'sample': sample,
                'mucus_zone': zone,
                'CN_type': cn,
                'count': cn_count,
                'total_cells': len(zone_data),
                'proportion': proportion
            })

# 转换为DataFrame
results_df = pd.DataFrame(results)

# 5. 为每个CN类型生成箱线图
print(f"\nGenerating boxplots for each CN type...")

for cn in cn_types:
    print(f"  Processing CN type: {cn}")
    
    # 筛选当前CN类型的数据
    cn_data = results_df[results_df['CN_type'] == cn].copy()
    
    if cn_data.empty:
        print(f"    No data for CN type {cn}, skipping...")
        continue
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 创建箱线图
    boxplot_data = []
    boxplot_labels = []
    boxplot_colors = []
    
    for zone in zones:
        zone_cn_data = cn_data[cn_data['mucus_zone'] == zone]
        
        if len(zone_cn_data) > 0:
            boxplot_data.append(zone_cn_data['proportion'].values)
            boxplot_labels.append(zone)
            boxplot_colors.append(zone_colors.get(zone, 'gray'))
    
    # 绘制箱线图
    box = ax.boxplot(boxplot_data, 
                     patch_artist=True,
                     labels=boxplot_labels,
                     widths=0.6)
    
    # 设置箱线图颜色
    for patch, color in zip(box['boxes'], boxplot_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # 设置箱线图其他元素的颜色
    for element in ['whiskers', 'caps', 'medians']:
        for line, color in zip(box[element], boxplot_colors * 2 if element in ['whiskers', 'caps'] else boxplot_colors):
            line.set_color(color)
            if element == 'medians':
                line.set_linewidth(2)
    
    # 添加散点图显示每个样本的数据点
    for i, zone in enumerate(zones, 1):
        zone_cn_data = cn_data[cn_data['mucus_zone'] == zone]
        
        if len(zone_cn_data) > 0:
            # 为每个点添加轻微的水平抖动
            x_jitter = np.random.normal(0, 0.05, len(zone_cn_data))
            x_pos = i + x_jitter
            
            ax.scatter(x_pos, 
                      zone_cn_data['proportion'].values,
                      color=zone_colors.get(zone, 'gray'),
                      alpha=0.6,
                      s=50,
                      edgecolor='white',
                      linewidth=0.5,
                      zorder=3)
    
    # 添加统计信息
    for i, zone in enumerate(zones, 1):
        zone_cn_data = cn_data[cn_data['mucus_zone'] == zone]
        
        if len(zone_cn_data) > 0:
            mean_val = zone_cn_data['proportion'].mean()
            median_val = zone_cn_data['proportion'].median()
            n_samples = len(zone_cn_data)
            
            # 在图上添加统计信息
            ax.text(i, 0.02, f'n={n_samples}', 
                   ha='center', va='bottom', fontsize=9)
            ax.text(i, ax.get_ylim()[1] * 0.98, 
                   f'mean={mean_val:.3f}\nmed={median_val:.3f}', 
                   ha='center', va='top', fontsize=8,
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    
    # 设置图形属性
    ax.set_ylabel(f'Proportion of CN={cn}', fontsize=12)
    ax.set_xlabel('Mucus Zone', fontsize=12)
    ax.set_title(f'Distribution of CN={cn} Proportion in Inner vs Outer Mucus Zones', fontsize=14, pad=15)
    
    # 设置y轴范围
    ax.set_ylim([0, 1.1])
    ax.yaxis.grid(True, alpha=0.3)
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=zone_colors['inner'], alpha=0.7, label='Inner zone'),
                       Patch(facecolor=zone_colors['outer'], alpha=0.7, label='Outer zone')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图形
    png_path = os.path.join(output_dir, f'CN_{cn}_proportion_boxplot.png')
    pdf_path = os.path.join(output_dir, f'CN_{cn}_proportion_boxplot.pdf')
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"    Saved: {png_path}")
    print(f"    Saved: {pdf_path}")
    
    plt.close()

# 6. 创建组合图形（所有CN类型在一个图中）
print("\nCreating combined figure with all CN types...")

fig, axes = plt.subplots(1, len(cn_types), figsize=(5*len(cn_types), 6))

if len(cn_types) == 1:
    axes = [axes]

for idx, (ax, cn) in enumerate(zip(axes, cn_types)):
    # 筛选当前CN类型的数据
    cn_data = results_df[results_df['CN_type'] == cn].copy()
    
    if cn_data.empty:
        ax.text(0.5, 0.5, f"No data for CN={cn}", 
                ha='center', va='center', transform=ax.transAxes)
        continue
    
    # 准备箱线图数据
    boxplot_data = []
    boxplot_labels = []
    
    for zone in zones:
        zone_cn_data = cn_data[cn_data['mucus_zone'] == zone]
        
        if len(zone_cn_data) > 0:
            boxplot_data.append(zone_cn_data['proportion'].values)
            boxplot_labels.append(zone)
    
    # 绘制箱线图
    box = ax.boxplot(boxplot_data, 
                     patch_artist=True,
                     labels=boxplot_labels,
                     widths=0.6)
    
    # 设置颜色
    for patch, color in zip(box['boxes'], [zone_colors.get(z, 'gray') for z in boxplot_labels]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # 添加散点
    for i, zone in enumerate(zones, 1):
        zone_cn_data = cn_data[cn_data['mucus_zone'] == zone]
        
        if len(zone_cn_data) > 0:
            x_jitter = np.random.normal(0, 0.05, len(zone_cn_data))
            x_pos = i + x_jitter
            
            ax.scatter(x_pos, 
                      zone_cn_data['proportion'].values,
                      color=zone_colors.get(zone, 'gray'),
                      alpha=0.6,
                      s=40,
                      edgecolor='white',
                      linewidth=0.5,
                      zorder=3)
    
    # 设置子图属性
    ax.set_ylabel(f'Proportion', fontsize=11)
    ax.set_title(f'CN = {cn}', fontsize=12)
    ax.set_ylim([0, 1.1])
    ax.yaxis.grid(True, alpha=0.3)
    
    # 添加统计信息
    for i, zone in enumerate(zones, 1):
        zone_cn_data = cn_data[cn_data['mucus_zone'] == zone]
        
        if len(zone_cn_data) > 0:
            median_val = zone_cn_data['proportion'].median()
            ax.text(i, 0.05, f'med={median_val:.2f}', 
                   ha='center', va='bottom', fontsize=9)

plt.suptitle('CN Proportion Distribution: Inner vs Outer Mucus Zones', fontsize=14, y=1.02)
plt.tight_layout()

# 保存组合图形
combined_png = os.path.join(output_dir, 'All_CN_proportions_combined.png')
combined_pdf = os.path.join(output_dir, 'All_CN_proportions_combined.pdf')

plt.savefig(combined_png, dpi=300, bbox_inches='tight')
plt.savefig(combined_pdf, bbox_inches='tight')
print(f"Saved combined figure: {combined_png}")
print(f"Saved combined figure: {combined_pdf}")

plt.show()

# 7. 统计分析
print("\n" + "="*50)
print("Statistical Analysis Results")
print("="*50)

# 执行统计检验
from scipy import stats

for cn in cn_types:
    cn_data = results_df[results_df['CN_type'] == cn]
    
    if len(cn_data) < 2:
        print(f"\nCN = {cn}: Insufficient data for statistical test")
        continue
    
    # 分离inner和outer组的数据
    inner_data = cn_data[cn_data['mucus_zone'] == 'inner']['proportion']
    outer_data = cn_data[cn_data['mucus_zone'] == 'outer']['proportion']
    
    if len(inner_data) < 2 or len(outer_data) < 2:
        print(f"\nCN = {cn}: Insufficient data in one or both groups")
        continue
    
    # Shapiro-Wilk正态性检验
    _, inner_p = stats.shapiro(inner_data)
    _, outer_p = stats.shapiro(outer_data)
    
    # 根据正态性选择检验方法
    if inner_p > 0.05 and outer_p > 0.05:
        # 数据正态分布，使用t检验
        t_stat, p_value = stats.ttest_ind(inner_data, outer_data)
        test_name = "Independent t-test"
    else:
        # 数据非正态分布，使用Mann-Whitney U检验
        u_stat, p_value = stats.mannwhitneyu(inner_data, outer_data)
        test_name = "Mann-Whitney U test"
    
    # 计算描述性统计
    inner_mean = inner_data.mean()
    inner_std = inner_data.std()
    outer_mean = outer_data.mean()
    outer_std = outer_data.std()
    
    print(f"\nCN = {cn}:")
    print(f"  Inner zone: n={len(inner_data)}, mean={inner_mean:.3f}, std={inner_std:.3f}")
    print(f"  Outer zone: n={len(outer_data)}, mean={outer_mean:.3f}, std={outer_std:.3f}")
    print(f"  Test: {test_name}, p-value = {p_value:.4f}")
    
    if p_value < 0.05:
        print(f"  Result: SIGNIFICANT difference (p < 0.05)")
    else:
        print(f"  Result: No significant difference")

# 8. 保存结果数据
print("\nSaving result data...")
results_csv = os.path.join(output_dir, 'CN_proportion_results.csv')
results_df.to_csv(results_csv, index=False)
print(f"Saved result data: {results_csv}")

# 9. 生成汇总表格
summary_table = results_df.pivot_table(
    index=['sample', 'mucus_zone'],
    columns='CN_type',
    values='proportion',
    aggfunc='first'
).reset_index()

summary_csv = os.path.join(output_dir, 'CN_proportion_summary.csv')
summary_table.to_csv(summary_csv, index=False)
print(f"Saved summary table: {summary_csv}")

print("\n" + "="*50)
print("Analysis Complete!")
print("="*50)
print(f"Output directory: {output_dir}")
print("\nGenerated files:")
print(f"  1. For each CN type: CN_X_proportion_boxplot.png/pdf")
print(f"  2. Combined figure: All_CN_proportions_combined.png/pdf")
print(f"  3. Raw results: CN_proportion_results.csv")
print(f"  4. Summary table: CN_proportion_summary.csv")
