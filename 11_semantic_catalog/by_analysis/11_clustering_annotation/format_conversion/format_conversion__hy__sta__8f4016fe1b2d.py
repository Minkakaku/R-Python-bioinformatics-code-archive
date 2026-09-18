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

# 检查数据基本信息
print("数据概览：")
print(f"总细胞数: {len(df)}")
print(f"样本数量: {df['sample'].nunique()}")
print(f"Stage类型: {df['stage'].unique().tolist()}")
print(f"所有细胞类型: {df['celltype'].unique().tolist()}")

# 3. 定义目标免疫细胞类型
target_immune_cells = ['T cells', 'B cells', 'Macrophages', 'Myeloids', 'Plasma']
available_immune_cells = [ct for ct in target_immune_cells if ct in df['celltype'].unique()]

print(f"\n目标分析的免疫细胞类型: {available_immune_cells}")
print(f"未找到的细胞类型: {set(target_immune_cells) - set(available_immune_cells)}")

if not available_immune_cells:
    raise ValueError("未找到任何目标免疫细胞类型，请检查细胞类型名称")

# 4. 计算每个stage中免疫细胞的比例
def calculate_immune_proportions(df, immune_cells):
    """计算每个stage中免疫细胞的比例"""
    
    # 获取所有stage
    stages = df['stage'].unique()
    
    results = []
    
    for stage in stages:
        # 筛选当前stage的数据
        stage_data = df[df['stage'] == stage]
        total_cells_in_stage = len(stage_data)
        
        # 计算每种免疫细胞的数量和比例
        for immune_cell in immune_cells:
            immune_cell_count = len(stage_data[stage_data['celltype'] == immune_cell])
            proportion = immune_cell_count / total_cells_in_stage if total_cells_in_stage > 0 else 0
            
            results.append({
                'stage': stage,
                'celltype': immune_cell,
                'count': immune_cell_count,
                'proportion': proportion,
                'total_cells': total_cells_in_stage
            })
    
    return pd.DataFrame(results)

# 计算比例
immune_prop_df = calculate_immune_proportions(df, available_immune_cells)

print("\n免疫细胞比例统计:")
print(immune_prop_df)

# 5. 可视化
print("\n正在生成可视化图表...")

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

# 图1: 堆积条形图显示各stage免疫细胞组成
plt.figure(figsize=(12, 8))

# 创建数据透视表用于堆积条形图
pivot_data = immune_prop_df.pivot(index='stage', columns='celltype', values='proportion').fillna(0)

# 计算其他细胞类型的比例
for stage in pivot_data.index:
    immune_total = pivot_data.loc[stage].sum()
    pivot_data.loc[stage, 'Other'] = 1 - immune_total

# 绘制堆积条形图
ax1 = plt.subplot(2, 2, 1)
bottom = np.zeros(len(pivot_data))

# 先绘制免疫细胞
for celltype in available_immune_cells + ['Other']:
    if celltype in pivot_data.columns:
        proportions = pivot_data[celltype].values
        color = plot_colors.get(celltype, 'lightgray')
        ax1.bar(pivot_data.index, proportions, bottom=bottom, 
                label=celltype, color=color, alpha=0.8, edgecolor='white', linewidth=0.5)
        bottom += proportions

ax1.set_ylabel('比例')
ax1.set_title('各Stage免疫细胞组成（堆积图）')
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax1.tick_params(axis='x', rotation=45)

# 图2: 分组条形图比较各免疫细胞在不同stage的比例
ax2 = plt.subplot(2, 2, 2)
sns.barplot(data=immune_prop_df, x='stage', y='proportion', hue='celltype',
            palette=plot_colors, ax=ax2, errwidth=1.5, capsize=0.05)
ax2.set_ylabel('比例')
ax2.set_title('免疫细胞比例比较（分组条形图）')
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.tick_params(axis='x', rotation=45)

# 图3: 热图显示比例
ax3 = plt.subplot(2, 2, 3)
heatmap_data = immune_prop_df.pivot(index='celltype', columns='stage', values='proportion').fillna(0)
sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='YlOrRd', 
            cbar_kws={'label': '比例'}, ax=ax3)
ax3.set_title('免疫细胞比例热图')
ax3.set_ylabel('细胞类型')

# 图4: 点图显示趋势
ax4 = plt.subplot(2, 2, 4)
for celltype in available_immune_cells:
    cell_data = immune_prop_df[immune_prop_df['celltype'] == celltype]
    ax4.plot(cell_data['stage'], cell_data['proportion'], 
             marker='o', label=celltype, linewidth=2, markersize=8,
             color=plot_colors.get(celltype, 'gray'))

ax4.set_ylabel('比例')
ax4.set_xlabel('Stage')
ax4.set_title('免疫细胞比例变化趋势')
ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax4.tick_params(axis='x', rotation=45)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('immune_cell_proportions_by_stage.png', dpi=300, bbox_inches='tight')
plt.show()

# 6. 统计检验 - 比较不同stage间免疫细胞比例的差异
print("\n正在进行统计检验...")

stats_results = []

for immune_cell in available_immune_cells:
    cell_data = immune_prop_df[immune_prop_df['celltype'] == immune_cell]
    
    if len(cell_data) >= 2:  # 至少两个stage才能比较
        # 准备ANOVA数据（需要原始数据而不是汇总数据）
        immune_cell_raw_data = []
        stages_for_test = []
        
        # 使用原始数据而不是汇总比例
        for stage in df['stage'].unique():
            stage_data = df[df['stage'] == stage]
            total_cells = len(stage_data)
            immune_cells = len(stage_data[stage_data['celltype'] == immune_cell])
            
            # 创建二进制向量（1表示目标细胞，0表示其他）
            if total_cells > 0:
                immune_vector = [1] * immune_cells + [0] * (total_cells - immune_cells)
                immune_cell_raw_data.extend(immune_vector)
                stages_for_test.extend([stage] * total_cells)
        
        if len(set(stages_for_test)) >= 2:  # 确保至少两个stage
            # 分组数据
            groups = []
            for stage in set(stages_for_test):
                group_data = [immune_cell_raw_data[i] for i, s in enumerate(stages_for_test) if s == stage]
                groups.append(group_data)
            
            # 单因素ANOVA检验
            try:
                f_stat, p_value = stats.f_oneway(*groups)
                
                stats_results.append({
                    'celltype': immune_cell,
                    'f_statistic': f_stat,
                    'p_value': p_value,
                    'n_stages': len(groups)
                })
            except:
                print(f"警告: 无法对 {immune_cell} 进行ANOVA检验")

# 显示统计结果
if stats_results:
    stats_df = pd.DataFrame(stats_results)
    stats_df['significant'] = stats_df['p_value'] < 0.05
    stats_df['p_value_adj'] = stats_df['p_value']  # 可以添加多重检验校正
    
    print("\n统计检验结果 (ANOVA):")
    print(stats_df.round(4))
    
    # 保存统计结果
    stats_df.to_csv('immune_cell_statistical_test.csv', index=False)
    print("统计检验结果已保存为: immune_cell_statistical_test.csv")
    
    # 标记显著的细胞类型
    significant_cells = stats_df[stats_df['significant']]['celltype'].tolist()
    if significant_cells:
        print(f"在p<0.05水平上显著的免疫细胞类型: {significant_cells}")
        
        # 对显著细胞进行事后检验
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
        
        print("\n显著细胞类型的事后检验 (Tukey HSD):")
        for celltype in significant_cells:
            print(f"\n{celltype}:")
            # 准备数据
            tukey_data = []
            groups = []
            
            for stage in df['stage'].unique():
                stage_data = df[df['stage'] == stage]
                total_cells = len(stage_data)
                immune_cells = len(stage_data[stage_data['celltype'] == celltype])
                
                if total_cells > 0:
                    # 创建二进制向量
                    immune_vector = [1] * immune_cells + [0] * (total_cells - immune_cells)
                    tukey_data.extend(immune_vector)
                    groups.extend([stage] * total_cells)
            
            if len(set(groups)) >= 2:
                try:
                    tukey_result = pairwise_tukeyhsd(tukey_data, groups, alpha=0.05)
                    print(tukey_result)
                except:
                    print("  无法进行事后检验")
    else:
        print("没有发现显著差异的免疫细胞类型")
else:
    print("样本量不足，无法进行统计检验")

# 7. 生成详细报告
print("\n" + "="*60)
print("免疫细胞比例分析报告")
print("="*60)

# 计算每个stage的免疫细胞总比例
stage_immune_summary = immune_prop_df.groupby('stage').agg({
    'proportion': 'sum',
    'total_cells': 'first'
}).reset_index()

stage_immune_summary.rename(columns={'proportion': 'total_immune_proportion'}, inplace=True)
stage_immune_summary['total_immune_proportion'] = stage_immune_summary['total_immune_proportion'].round(4)

print("\n各Stage免疫细胞总比例:")
for _, row in stage_immune_summary.iterrows():
    print(f"{row['stage']}: {row['total_immune_proportion']*100:.2f}% (总细胞数: {row['total_cells']})")

# 显示每个免疫细胞类型的最大比例阶段
print("\n各免疫细胞类型比例最高的Stage:")
for celltype in available_immune_cells:
    cell_data = immune_prop_df[immune_prop_df['celltype'] == celltype]
    if not cell_data.empty:
        max_row = cell_data.loc[cell_data['proportion'].idxmax()]
        print(f"{celltype}: {max_row['stage']} ({max_row['proportion']*100:.2f}%)")

# 8. 保存比例数据
immune_prop_df.to_csv('immune_cell_proportions.csv', index=False)
print(f"\n免疫细胞比例数据已保存为: immune_cell_proportions.csv")

print("\n分析完成！所有图表和结果已生成。")
