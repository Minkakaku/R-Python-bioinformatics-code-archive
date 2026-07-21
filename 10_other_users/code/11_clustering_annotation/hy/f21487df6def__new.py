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
adata = sc.read_h5ad("adata_celltype_stage.h5ad")  # 替换为实际路径

# 2. 从adata.obs中提取所需信息
df = adata.obs[['sample', 'stage', 'celltype']].copy()

# 3. 去除stage中的"Others"组
df = df[df['stage'] != 'others']

# 4. 计算每个(sample, stage, celltype)组合的细胞数量
count_df = df.groupby(['sample', 'stage', 'celltype']).size().reset_index(name='count')

# 5. 计算每个(sample, stage)分组内各细胞类型的比例
count_df['total_per_group'] = count_df.groupby(['sample', 'stage'])['count'].transform('sum')
count_df['proportion'] = count_df['count'] / count_df['total_per_group']

print("数据概览：")
print(f"样本数量: {df['sample'].nunique()}")
print(f"Stage数量: {df['stage'].nunique()}")
print(f"细胞类型数量: {df['celltype'].nunique()}")

# 第一部分：每个样本的细胞类型比例可视化（图例在外）
print("\n正在生成每个样本的细胞类型比例图...")

# 设置stage顺序
stage_order = ['Normal', 'Adenoma', 'Tumor']

# 为每个样本创建单独的图表
samples = count_df['sample'].unique()
celltypes = count_df['celltype'].unique()

# 定义颜色映射
palette = sns.color_palette('Set3', len(celltypes))
color_dict = dict(zip(celltypes, palette))

for sample in samples:
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 6))
    
    sample_data = count_df[count_df['sample'] == sample]
    
    # 为每个stage创建堆积柱状图
    x_positions = np.arange(len(stage_order))
    bar_width = 0.6
    bottom = np.zeros(len(stage_order))
    
    # 收集图例信息
    legend_handles = []
    legend_labels = []
    
    for celltype in celltypes:
        proportions = []
        
        for stage in stage_order:
            stage_data = sample_data[sample_data['stage'] == stage]
            celltype_data = stage_data[stage_data['celltype'] == celltype]
            
            if not celltype_data.empty:
                proportion = celltype_data['proportion'].iloc[0]
            else:
                proportion = 0
            proportions.append(proportion)
        
        # 绘制柱状图
        bars = ax.bar(x_positions, proportions, bottom=bottom, 
                     color=color_dict[celltype], label=celltype, 
                     width=bar_width, edgecolor='white', linewidth=0.5)
        
        # 更新底部位置
        bottom += np.array(proportions)
        
        # 添加到图例（只添加一次）
        if celltype not in legend_labels:
            legend_handles.append(bars[0])
            legend_labels.append(celltype)
    
    # 设置图表属性
    ax.set_xlabel('Stage', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title(f'Cell Type Proportions in {sample}', fontsize=14, fontweight='bold')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(stage_order)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis='y')
    
    # 将图例放在图表外部右侧
    ax.legend(legend_handles, legend_labels, 
              title='Cell Type', 
              bbox_to_anchor=(1.05, 1), 
              loc='upper left',
              frameon=True,
              fancybox=True,
              shadow=True)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表
    filename = f'{sample}_celltype_proportions.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"已生成图表: {filename}")

# 第二部分：不同stage分组间的细胞类型比例变化统计
print("\n正在分析不同stage间的细胞类型比例变化...")

# 6. 统计不同stage间的比例变化
# 计算每个stage下各细胞类型的平均比例
stage_celltype_stats = count_df.groupby(['stage', 'celltype'])['proportion'].agg(['mean', 'std', 'count']).reset_index()
stage_celltype_stats.rename(columns={'mean': 'mean_proportion', 'std': 'std_proportion', 'count': 'n_samples'}, inplace=True)

# 确保stage顺序
stage_celltype_stats['stage'] = pd.Categorical(stage_celltype_stats['stage'], categories=stage_order, ordered=True)
stage_celltype_stats = stage_celltype_stats.sort_values('stage')

print("\n各Stage细胞类型平均比例:")
print(stage_celltype_stats.round(3))

# 7. 可视化stage间的比例变化
# 创建热图
plt.figure(figsize=(12, 8))

# 创建数据透视表用于热图
heatmap_data = stage_celltype_stats.pivot(index='celltype', columns='stage', values='mean_proportion').fillna(0)

# 确保列的顺序为Normal, Adenoma, Tumor
heatmap_data = heatmap_data[stage_order]

# 绘制热图
plt.subplot(1, 2, 1)
sns.heatmap(heatmap_data, annot=True, cmap='YlOrRd', fmt='.3f', 
            cbar_kws={'label': 'Proportion'}, linewidths=0.5)
plt.title('Cell Type Proportions Heatmap by Stage', fontsize=12, fontweight='bold')
plt.ylabel('Cell Type')
plt.xlabel('Stage')

# 创建折线图
plt.subplot(1, 2, 2)
# 为每个细胞类型绘制一条线
for celltype in celltypes:
    celltype_data = stage_celltype_stats[stage_celltype_stats['celltype'] == celltype]
    # 确保数据按stage_order顺序排列
    celltype_data = celltype_data.set_index('stage').reindex(stage_order).reset_index()
    proportions = celltype_data['mean_proportion'].values
    plt.plot(stage_order, proportions, marker='o', label=celltype, 
             linewidth=2, markersize=8, alpha=0.8)

plt.xlabel('Stage', fontsize=12)
plt.ylabel('Mean Proportion', fontsize=12)
plt.title('Cell Type Proportions Across Stages', fontsize=12, fontweight='bold')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig('stage_comparison_combined.png', dpi=300, bbox_inches='tight')
plt.show()

print("Stage比较图已保存为: stage_comparison_combined.png")

# 8. 统计检验 - 比较不同stage间的显著性差异
print("\n正在进行统计检验...")

# 选择样本量足够的细胞类型进行检验
celltype_counts = count_df['celltype'].value_counts()
valid_celltypes = celltype_counts[celltype_counts >= 3].index  # 至少3个样本

stats_results = []

for celltype in valid_celltypes:
    celltype_data = count_df[count_df['celltype'] == celltype]
    
    if len(celltype_data['stage'].unique()) >= 2:  # 至少两个stage才能比较
        # 准备ANOVA数据
        stage_groups = [group['proportion'].values for name, group in celltype_data.groupby('stage')]
        
        if len(stage_groups) >= 2 and all(len(g) >= 2 for g in stage_groups):  # 每组至少2个样本
            # 单因素ANOVA检验
            f_stat, p_value = stats.f_oneway(*stage_groups)
            
            stats_results.append({
                'celltype': celltype,
                'f_statistic': f_stat,
                'p_value': p_value,
                'n_stages': len(stage_groups),
                'total_samples': len(celltype_data)
            })

# 转换为DataFrame并显示结果
if stats_results:
    stats_df = pd.DataFrame(stats_results)
    stats_df['significant'] = stats_df['p_value'] < 0.05
    print("\n统计检验结果 (ANOVA):")
    print(stats_df.round(4))
    
    # 保存统计结果
    stats_df.to_csv('statistical_test_results.csv', index=False)
    print("\n统计检验结果已保存为: statistical_test_results.csv")
    
    # 标记显著的细胞类型
    significant_celltypes = stats_df[stats_df['significant']]['celltype'].tolist()
    if significant_celltypes:
        print(f"\n在p<0.05水平上显著的细胞类型: {significant_celltypes}")
    else:
        print("\n没有发现显著差异的细胞类型")
else:
    print("\n样本量不足，无法进行统计检验")

# 9. 保存比例表格
proportion_table = count_df.pivot_table(
    index=['sample', 'stage'], 
    columns='celltype', 
    values='proportion', 
    fill_value=0
)

# 保存比例数据
proportion_table.to_csv('celltype_proportions.csv')
print(f"\n比例表格已保存，形状: {proportion_table.shape}")

print("\n分析完成！所有图表和结果已保存。")
