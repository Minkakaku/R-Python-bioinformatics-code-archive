import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 1. 读取h5ad文件
adata = sc.read_h5ad("merged_celltype_stage.h5ad")  # 替换为您的实际文件路径

# 2. 从adata.obs中提取所需信息
df = adata.obs[['sample', 'stage', 'celltype']].copy()

# 3. 关键修改：去除stage中的"others"组（确保大小写正确）
print("原始数据中的stage类型:", df['stage'].unique().tolist())

# 使用更健壮的方式去除"others"（不区分大小写）
df = df[~df['stage'].str.lower().str.contains('other', na=False)]

print("去除'others'后的stage类型:", df['stage'].unique().tolist())
print(f"过滤后剩余细胞数: {len(df)}")

# 4. 定义目标免疫细胞类型
target_immune_cells = ['T cells', 'B cells', 'Macrophages', 'Myeloids', 'Plasma']
available_immune_cells = [ct for ct in target_immune_cells if ct in df['celltype'].unique()]

print(f"\n目标分析的免疫细胞类型: {available_immune_cells}")

if not available_immune_cells:
    raise ValueError("未找到任何目标免疫细胞类型，请检查细胞类型名称")

# 5. 计算每个(sample, stage, celltype)组合的细胞数量
count_df = df.groupby(['sample', 'stage', 'celltype']).size().reset_index(name='count')

# 6. 计算每个(sample, stage)分组内各细胞类型的比例
count_df['total_per_group'] = count_df.groupby(['sample', 'stage'])['count'].transform('sum')
count_df['proportion'] = count_df['count'] / count_df['total_per_group']

print("\n数据概览（已去除'others'）：")
print(f"样本数量: {df['sample'].nunique()}")
print(f"Stage数量: {df['stage'].nunique()}")
print(f"细胞类型数量: {df['celltype'].nunique()}")

# 关键修改：手动设置stage顺序为 Normal, Adenoma, Tumor
stage_order = ['Normal', 'Adenoma', 'Tumor']
print(f"分析用的stage顺序: {stage_order}")

# 7. 计算每个stage中免疫细胞的比例
def calculate_immune_proportions(df, immune_cells, stage_order):
    """计算每个stage中免疫细胞的比例"""
    
    results = []
    
    for stage in stage_order:
        # 筛选当前stage的数据
        stage_data = df[df['stage'] == stage]
        total_cells_in_stage = len(stage_data)
        
        if total_cells_in_stage == 0:
            # 如果该stage没有数据，添加0值记录
            for immune_cell in immune_cells:
                results.append({
                    'stage': stage,
                    'celltype': immune_cell,
                    'count': 0,
                    'proportion': 0,
                    'total_cells': 0
                })
            continue
            
        # 计算每种免疫细胞的数量和比例
        for immune_cell in immune_cells:
            immune_cell_count = len(stage_data[stage_data['celltype'] == immune_cell])
            proportion = immune_cell_count / total_cells_in_stage
            
            results.append({
                'stage': stage,
                'celltype': immune_cell,
                'count': immune_cell_count,
                'proportion': proportion * 100,  # 转换为百分比
                'total_cells': total_cells_in_stage
            })
    
    return pd.DataFrame(results)

# 计算比例
immune_prop_df = calculate_immune_proportions(df, available_immune_cells, stage_order)

print("\n免疫细胞比例统计（已去除'others'）:")
print(immune_prop_df.round(3))

# 8. 可视化（基于去除others后的数据）
print("\n正在生成可视化图表（已去除'others'）...")

# 设置颜色
immune_colors = {
    'T cells': '#E31A1C',      # 红色
    'B cells': '#1F78B4',     # 蓝色  
    'Macrophages': '#33A02C',  # 绿色
    'Myeloids': '#FF7F00',    # 橙色
    'Plasma': '#6A3D9A'       # 紫色
}

# 只保留存在的细胞类型的颜色
plot_colors = {k: v for k, v in immune_colors.items() if k in available_immune_cells}

# 创建综合可视化图表
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Immune_cell_proportions', fontsize=16, fontweight='bold')

# 图1: 堆积条形图
ax1 = axes[0, 0]
if not immune_prop_df.empty:
    # 创建数据透视表
    pivot_data = immune_prop_df.pivot(index='stage', columns='celltype', values='proportion').fillna(0)
    
    # 关键修改：确保数据按照stage_order顺序排列
    pivot_data = pivot_data.reindex(stage_order)
    
    # 绘制堆积条形图
    bottom = np.zeros(len(pivot_data))
    for celltype in available_immune_cells:
        if celltype in pivot_data.columns:
            proportions = pivot_data[celltype].values
            color = plot_colors.get(celltype, 'gray')
            ax1.bar(range(len(pivot_data)), proportions, bottom=bottom, 
                   label=celltype, color=color, alpha=0.8)
            bottom += proportions
    
    ax1.set_xticks(range(len(pivot_data)))
    ax1.set_xticklabels(pivot_data.index, rotation=45)
    ax1.set_ylabel('Percentages(%)')
    ax1.set_title('Stages')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

# 图2: 分组条形图
ax2 = axes[0, 1]
if not immune_prop_df.empty:
    # 关键修改：确保stage顺序正确
    immune_prop_df['stage'] = pd.Categorical(immune_prop_df['stage'], categories=stage_order, ordered=True)
    immune_prop_df = immune_prop_df.sort_values('stage')
    
    sns.barplot(data=immune_prop_df, x='stage', y='proportion', hue='celltype',
                palette=plot_colors, ax=ax2, errorbar=None, order=stage_order)  # 添加order参数
    ax2.set_ylabel('Percentages(%)')
    ax2.set_title('Immune_cell')
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend(title='Celltypes')

# 图3: 热图
ax3 = axes[1, 0]
if not immune_prop_df.empty:
    heatmap_data = immune_prop_df.pivot(index='celltype', columns='stage', values='proportion').fillna(0)
    heatmap_data = heatmap_data[stage_order]  # 确保列顺序
    
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd', 
                cbar_kws={'label': 'Proportion (%)'}, ax=ax3)
    ax3.set_title('Immune_cell heatmap')
    ax3.set_ylabel('Celltypes')

# 图4: 趋势线图
ax4 = axes[1, 1]
if not immune_prop_df.empty:
    for celltype in available_immune_cells:
        celltype_data = immune_prop_df[immune_prop_df['celltype'] == celltype]
        if not celltype_data.empty:
            # 关键修改：确保数据按stage_order顺序排列
            celltype_data = celltype_data.set_index('stage').reindex(stage_order).reset_index()
            proportions = celltype_data['proportion'].fillna(0).values
            
            color = plot_colors.get(celltype, 'gray')
            ax4.plot(range(len(stage_order)), proportions, marker='o', label=celltype,
                    color=color, linewidth=2, markersize=6)
    
    ax4.set_xticks(range(len(stage_order)))
    ax4.set_xticklabels(stage_order, rotation=45)
    ax4.set_ylabel('Percentages(%)')
    ax4.set_xlabel('Stage')
    ax4.set_title('Immune_cell')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('immune_cell_proportions_no_others.png', dpi=300, bbox_inches='tight')
plt.savefig('immune_cell_proportions_no_others.pdf', dpi=300, bbox_inches='tight')
plt.show()

print("可视化图表已保存为: immune_cell_proportions_no_others.png")

# 9. 统计检验（基于去除others后的数据）
print("\n正在进行统计检验（已去除'others'）...")

# 使用原始数据而非汇总数据进行统计检验
stats_results = []

for immune_cell in available_immune_cells:
    # 准备ANOVA数据：创建二进制向量（1=目标细胞，0=其他细胞）
    anova_data = []
    groups = []
    
    for stage in stage_order:  # 按照指定顺序处理stage
        stage_data = df[df['stage'] == stage]
        if len(stage_data) > 0:
            target_cells = len(stage_data[stage_data['celltype'] == immune_cell])
            other_cells = len(stage_data) - target_cells
            
            # 创建二进制数据
            anova_data.extend([1] * target_cells + [0] * other_cells)
            groups.extend([stage] * len(stage_data))
    
    if len(set(groups)) >= 2:  # 至少需要两个组
        try:
            # 分组进行ANOVA
            group_data = []
            for stage in set(groups):
                group_values = [anova_data[i] for i, g in enumerate(groups) if g == stage]
                group_data.append(group_values)
            
            if all(len(g) >= 2 for g in group_data):  # 每组至少2个样本
                f_stat, p_value = stats.f_oneway(*group_data)
                
                stats_results.append({
                    'celltype': immune_cell,
                    'f_statistic': f_stat,
                    'p_value': p_value,
                    'n_stages': len(group_data)
                })
        except Exception as e:
            print(f"警告: 对 {immune_cell} 进行ANOVA检验时出错: {str(e)}")

# 显示统计结果
if stats_results:
    stats_df = pd.DataFrame(stats_results)
    stats_df['significant'] = stats_df['p_value'] < 0.05
    
    print("\n统计检验结果 (ANOVA):")
    print(stats_df.round(4))
    
    # 保存统计结果
    stats_df.to_csv('immune_cell_stats_no_others.csv', index=False)
    print("统计检验结果已保存为: immune_cell_stats_no_others.csv")
    
    significant_cells = stats_df[stats_df['significant']]['celltype'].tolist()
    if significant_cells:
        print(f"在p<0.05水平上显著的细胞类型: {significant_cells}")
    else:
        print("没有发现显著差异的细胞类型")
else:
    print("样本量不足，无法进行统计检验")

# 10. 保存处理后的数据
immune_prop_df.to_csv('immune_cell_proportions_no_others.csv', index=False)
print("\n免疫细胞比例数据已保存为: immune_cell_proportions_no_others.csv")

# 11. 生成总结报告
print("\n" + "="*60)
print("免疫细胞比例分析总结（已去除'others'阶段）")
print("="*60)

# 计算总体统计
total_cells_analyzed = len(df)
total_samples = df['sample'].nunique()

print(f"分析的总细胞数: {total_cells_analyzed}")
print(f"涉及的样本数: {total_samples}")
print(f"分析的免疫细胞类型数: {len(available_immune_cells)}")
print(f"分析的stage数量: {len(stage_order)}")
print(f"Stage顺序: {stage_order}")

# 显示每个免疫细胞类型的统计
print("\n各免疫细胞类型在不同stage中的平均比例 (%):")
summary_table = immune_prop_df.groupby(['celltype', 'stage'])['proportion'].mean().unstack().fillna(0)
summary_table = summary_table[stage_order]  # 确保顺序
print(summary_table.round(2))

print("\n分析完成！所有结果基于去除'others'阶段后的数据。")
