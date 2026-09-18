import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from matplotlib.patches import Patch
warnings.filterwarnings('ignore')

# 创建输出目录
output_dir = "CD8_analysis"
os.makedirs(output_dir, exist_ok=True)

# 1. 加载数据
adata = sc.read_h5ad("extracted_samples.h5ad")

# 2. 检查数据
print(f"数据形状: {adata.shape}")
print(f"细胞类型列表: {adata.obs['celltype'].unique()}")
print(f"样本列表: {adata.obs['sample'].unique()}")
print(f"组列表: {adata.obs['group'].unique()}")

# 3. 计算每个样本中CD8+ T细胞的比例
cell_type_to_plot = "CD8+ T"

# 检查是否存在CD8+ T细胞类型
if cell_type_to_plot not in adata.obs['celltype'].unique():
    print(f"警告: 数据中没有找到 '{cell_type_to_plot}' 细胞类型")
    print(f"可用的细胞类型: {list(adata.obs['celltype'].unique())}")
    # 尝试查找类似的细胞类型
    similar_types = [ct for ct in adata.obs['celltype'].unique() if 'CD8' in ct or 'T' in ct]
    if similar_types:
        print(f"类似的细胞类型: {similar_types}")
    # 退出脚本或使用第一个类似的细胞类型
    if similar_types:
        cell_type_to_plot = similar_types[0]
        print(f"将使用 '{cell_type_to_plot}' 作为替代")
    else:
        raise ValueError(f"未找到 '{cell_type_to_plot}' 细胞类型")

# 获取所有样本
samples = adata.obs['sample'].unique()
group_for_sample = {}

# 计算每个样本的总细胞数
total_cells_per_sample = {}
for sample in samples:
    sample_mask = adata.obs['sample'] == sample
    total_cells_per_sample[sample] = np.sum(sample_mask)
    # 获取样本的group
    group_for_sample[sample] = adata.obs[sample_mask]['group'].iloc[0]

# 计算每个样本中CD8+ T细胞的数量
cd8_counts = {}
for sample in samples:
    sample_mask = adata.obs['sample'] == sample
    cd8_mask = (adata.obs['celltype'] == cell_type_to_plot) & sample_mask
    cd8_counts[sample] = np.sum(cd8_mask)

# 计算每个样本中CD8+ T细胞的比例
cd8_proportions = {}
for sample in samples:
    if total_cells_per_sample[sample] > 0:
        cd8_proportions[sample] = (cd8_counts[sample] / total_cells_per_sample[sample]) * 100
    else:
        cd8_proportions[sample] = 0

# 4. 创建数据框用于绘图
df_proportions = pd.DataFrame({
    'Sample': list(cd8_proportions.keys()),
    'Proportion': list(cd8_proportions.values()),
    'Group': [group_for_sample[sample] for sample in cd8_proportions.keys()]
})

# 按照组和比例排序
df_proportions = df_proportions.sort_values(['Group', 'Proportion'], ascending=[True, False])

# 5. 计算每个组中others的比例
# 首先计算每个样本中非CD8+ T细胞的比例
others_proportions = {}
for sample in samples:
    if total_cells_per_sample[sample] > 0:
        others_proportions[sample] = 100 - cd8_proportions[sample]
    else:
        others_proportions[sample] = 0

df_others = pd.DataFrame({
    'Sample': list(others_proportions.keys()),
    'Proportion': list(others_proportions.values()),
    'Celltype': 'Others',
    'Group': [group_for_sample[sample] for sample in others_proportions.keys()]
})

# 6. 创建合并的数据框用于堆积柱状图
df_cd8 = pd.DataFrame({
    'Sample': list(cd8_proportions.keys()),
    'Proportion': list(cd8_proportions.values()),
    'Celltype': cell_type_to_plot,
    'Group': [group_for_sample[sample] for sample in cd8_proportions.keys()]
})

df_stacked = pd.concat([df_cd8, df_others])
df_stacked = df_stacked.merge(df_proportions[['Sample', 'Group']], on=['Sample', 'Group'])

# 按照样本在df_proportions中的顺序排序
sample_order = df_proportions['Sample'].tolist()
df_stacked['Sample'] = pd.Categorical(df_stacked['Sample'], categories=sample_order, ordered=True)
df_stacked = df_stacked.sort_values('Sample')

# 7. 创建df_cd8_only数据框（提前定义，避免后面引用错误）
df_cd8_only = pd.DataFrame({
    'Sample': sample_order,
    'Group': df_proportions['Group'].tolist(),
    'CD8_T_proportion': [cd8_proportions[sample] for sample in sample_order],
    'Others_proportion': [others_proportions[sample] for sample in sample_order],
    'Total_cells': [total_cells_per_sample[sample] for sample in sample_order],
    'CD8_T_count': [cd8_counts[sample] for sample in sample_order]
})

# 8. 设置绘图样式
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_style("whitegrid")

# 9. 为不同的group定义颜色
unique_groups = df_proportions['Group'].unique()
group_colors = {}
cmap = plt.cm.tab10
for i, group in enumerate(unique_groups):
    group_colors[group] = cmap(i % 10)

# 获取样本和组信息
samples_ordered = df_proportions['Sample'].tolist()
groups_ordered = df_proportions['Group'].tolist()

# 创建x轴位置
x_positions = np.arange(len(samples_ordered))

# 获取CD8+ T细胞比例
cd8_props = [cd8_proportions[sample] for sample in samples_ordered]
others_props = [others_proportions[sample] for sample in samples_ordered]

# 10. 绘制柱状图
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), height_ratios=[3, 1])

# 10.1 绘制堆积柱状图显示CD8+ T细胞和others的比例
# 为CD8+ T细胞和others定义颜色
colors = {cell_type_to_plot: '#E41A1C', 'Others': '#377EB8'}

# 绘制堆积柱状图
# 先绘制Others
ax1.bar(x_positions, others_props, color=colors['Others'], label='Others', 
        edgecolor='black', linewidth=0.5)

# 再绘制CD8+ T细胞
ax1.bar(x_positions, cd8_props, bottom=others_props, color=colors[cell_type_to_plot], 
        label=cell_type_to_plot, edgecolor='black', linewidth=0.5)

# 设置图表属性
ax1.set_ylabel('Proportion (%)', fontsize=12, fontweight='bold')
ax1.set_title(f'Proportion of {cell_type_to_plot} Cells in Each Sample', fontsize=14, fontweight='bold')
ax1.set_ylim(0, 100)

# 添加网格
ax1.grid(True, linestyle='--', alpha=0.5, axis='y')

# 添加图例 - 移到右侧外部
ax1.legend(title='Cell Type', fontsize=10, title_fontsize=11,
           frameon=True, fancybox=True, shadow=True,
           bbox_to_anchor=(1.05, 0.5), loc='center left')

# 10.2 在第二个子图中显示样本名称，按照group颜色区分
# 设置x轴
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(0, 1)
ax2.axis('off')

# 在底部添加样本名称，按照group颜色
for i, (sample, group) in enumerate(zip(samples_ordered, groups_ordered)):
    color = group_colors[group]
    # 旋转45度，并调整位置靠近X轴
    ax2.text(i, 0.2, sample, rotation=45, fontsize=9, 
             ha='right', va='top', color=color, fontweight='bold')

# 添加group颜色图例
legend_elements = [Patch(facecolor=group_colors[group], edgecolor='black', 
                         label=group) for group in unique_groups]
ax2.legend(handles=legend_elements, title='Group', fontsize=9, 
           title_fontsize=10, frameon=True, fancybox=True, 
           shadow=True, bbox_to_anchor=(1.05, 0.5), loc='center left')

# 调整布局
plt.tight_layout()

# 保存图表
fig.savefig(f'{output_dir}/CD8_T_cell_proportion.png', 
            dpi=300, bbox_inches='tight')
fig.savefig(f'{output_dir}/CD8_T_cell_proportion.pdf', 
            bbox_inches='tight')

plt.close()

# 11. 绘制简单的柱状图（只显示CD8+ T细胞比例）
fig, ax = plt.subplots(figsize=(12, 6))

# 为不同group的样本使用不同颜色
bar_colors = [group_colors[group] for group in groups_ordered]

# 绘制柱状图
bars = ax.bar(x_positions, cd8_props, color=bar_colors, edgecolor='black', linewidth=0.5)

# 设置图表属性
ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
ax.set_ylabel(f'{cell_type_to_plot} Cell Proportion (%)', fontsize=12, fontweight='bold')
ax.set_title(f'Proportion of {cell_type_to_plot} Cells in Each Sample', fontsize=14, fontweight='bold')
ax.set_xticks(x_positions)

# 修改X轴标签：旋转45度，并调整位置
ax.set_xticklabels(samples_ordered, rotation=45, fontsize=9, ha='right', va='top')

# 调整X轴标签与X轴的距离
ax.tick_params(axis='x', which='major', pad=2)

# 设置Y轴范围
ax.set_ylim(0, max(cd8_props) * 1.1 if max(cd8_props) > 0 else 10)

# 添加网格
ax.grid(True, linestyle='--', alpha=0.5, axis='y')

# 添加group颜色图例 - 移到右侧外部
legend_elements = [Patch(facecolor=group_colors[group], edgecolor='black', 
                         label=group) for group in unique_groups]
ax.legend(handles=legend_elements, title='Group', fontsize=10, 
          title_fontsize=11, frameon=True, fancybox=True, shadow=True,
          bbox_to_anchor=(1.05, 0.5), loc='center left')

# 调整布局，为右侧图例留出空间
plt.subplots_adjust(right=0.85)

# 保存图表
fig.savefig(f'{output_dir}/CD8_T_cell_proportion_simple.png', 
            dpi=300, bbox_inches='tight')
fig.savefig(f'{output_dir}/CD8_T_cell_proportion_simple.pdf', 
            bbox_inches='tight')

plt.close()

# 12. 绘制简单柱状图（垂直版本，标签不旋转）
fig, ax = plt.subplots(figsize=(12, 6))

# 为不同group的样本使用不同颜色
bar_colors = [group_colors[group] for group in groups_ordered]

# 绘制柱状图
bars = ax.bar(x_positions, cd8_props, color=bar_colors, edgecolor='black', linewidth=0.5)

# 设置图表属性
ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
ax.set_ylabel(f'{cell_type_to_plot} Cell Proportion (%)', fontsize=12, fontweight='bold')
ax.set_title(f'Proportion of {cell_type_to_plot} Cells in Each Sample', fontsize=14, fontweight='bold')
ax.set_xticks(x_positions)

# 垂直标签，靠近X轴
ax.set_xticklabels(samples_ordered, rotation=90, fontsize=9, ha='center', va='top')

# 调整X轴标签与X轴的距离
ax.tick_params(axis='x', which='major', pad=0)

# 设置Y轴范围
ax.set_ylim(0, max(cd8_props) * 1.1 if max(cd8_props) > 0 else 10)

# 添加网格
ax.grid(True, linestyle='--', alpha=0.5, axis='y')

# 添加group颜色图例 - 移到右侧外部
legend_elements = [Patch(facecolor=group_colors[group], edgecolor='black', 
                         label=group) for group in unique_groups]
ax.legend(handles=legend_elements, title='Group', fontsize=10, 
          title_fontsize=11, frameon=True, fancybox=True, shadow=True,
          bbox_to_anchor=(1.05, 0.5), loc='center left')

# 调整布局，为右侧图例留出空间
plt.subplots_adjust(right=0.85, bottom=0.2)  # 增加底部边距以容纳垂直标签

# 保存图表
fig.savefig(f'{output_dir}/CD8_T_cell_proportion_simple_vertical.png', 
            dpi=300, bbox_inches='tight')
fig.savefig(f'{output_dir}/CD8_T_cell_proportion_simple_vertical.pdf', 
            bbox_inches='tight')

plt.close()

# 13. 绘制分组统计图
fig, ax = plt.subplots(figsize=(8, 6))

# 为每个组绘制点图
for i, group in enumerate(unique_groups):
    group_data = df_cd8_only[df_cd8_only['Group'] == group]
    x_pos = i + 1
    y_vals = group_data['CD8_T_proportion'].values
    
    # 绘制散点
    ax.scatter([x_pos] * len(y_vals), y_vals, color=group_colors[group], 
               s=50, alpha=0.6, label=group)
    
    # 绘制均值和误差线
    mean_val = np.mean(y_vals)
    std_val = np.std(y_vals)
    ax.errorbar(x_pos, mean_val, yerr=std_val, color=group_colors[group], 
                capsize=5, capthick=2, elinewidth=2, fmt='o', markersize=8)

# 设置图表属性
ax.set_xlabel('Group', fontsize=12, fontweight='bold')
ax.set_ylabel(f'{cell_type_to_plot} Cell Proportion (%)', fontsize=12, fontweight='bold')
ax.set_title(f'Proportion of {cell_type_to_plot} Cells by Group', fontsize=14, fontweight='bold')
ax.set_xticks(range(1, len(unique_groups) + 1))
ax.set_xticklabels(unique_groups, fontsize=11)
ax.set_ylim(0, max(df_cd8_only['CD8_T_proportion']) * 1.1 if max(df_cd8_only['CD8_T_proportion']) > 0 else 10)

# 添加网格
ax.grid(True, linestyle='--', alpha=0.5, axis='y')

# 调整布局
plt.tight_layout()

# 保存图表
fig.savefig(f'{output_dir}/CD8_T_cell_proportion_by_group.png', 
            dpi=300, bbox_inches='tight')
fig.savefig(f'{output_dir}/CD8_T_cell_proportion_by_group.pdf', 
            bbox_inches='tight')

plt.close()

# 14. 保存分组统计
group_stats = []
for group in unique_groups:
    group_samples = df_cd8_only[df_cd8_only['Group'] == group]
    if len(group_samples) > 0:
        mean_prop = group_samples['CD8_T_proportion'].mean()
        std_prop = group_samples['CD8_T_proportion'].std()
        group_stats.append({
            'Group': group,
            'N_samples': len(group_samples),
            'Mean_CD8_T_proportion': mean_prop,
            'Std_CD8_T_proportion': std_prop,
            'Min_CD8_T_proportion': group_samples['CD8_T_proportion'].min(),
            'Max_CD8_T_proportion': group_samples['CD8_T_proportion'].max()
        })

df_group_stats = pd.DataFrame(group_stats)
df_group_stats.to_csv(f'{output_dir}/CD8_T_cell_proportions_group_stats.csv', index=False)

# 15. 保存CD8+ T细胞比例数据
df_cd8_only.to_csv(f'{output_dir}/CD8_T_cell_proportions.csv', index=False)

print(f"分析完成！结果已保存到 {output_dir} 目录")
print(f"已生成 4 张图表")
print(f"数据汇总已保存为 CSV 文件")
print(f"总样本数: {len(samples_ordered)}")
print(f"CD8+ T细胞平均比例: {np.mean(cd8_props):.2f}%")
print(f"CD8+ T细胞比例范围: {min(cd8_props):.2f}% - {max(cd8_props):.2f}%")
