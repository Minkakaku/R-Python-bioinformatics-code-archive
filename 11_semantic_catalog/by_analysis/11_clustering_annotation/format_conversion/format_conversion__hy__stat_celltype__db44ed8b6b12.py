import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体（如果标签包含中文）
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置绘图样式
sns.set_style("whitegrid")
sns.set_palette("Set2")

# 1. 读取数据
print("正在读取数据...")
adata = sc.read_h5ad('adata_celltype.h5ad')
print(f"数据形状: {adata.shape}")
print(f"样本数量: {adata.obs['sample'].nunique()}")
print(f"细胞类型数量: {adata.obs['celltype'].nunique()}")
print(f"组别: {sorted(adata.obs['group'].unique())}")

# 检查必要列是否存在
required_cols = ['celltype', 'group', 'sample']
missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    raise ValueError(f"缺失必要的列: {missing_cols}")

# 2. 统计每个group中每个celltype的数量
print("\n1. 每个group中每个celltype的数量统计:")
group_celltype_counts = pd.crosstab(adata.obs['group'], adata.obs['celltype'])
print(group_celltype_counts)

# 保存统计结果到CSV
group_celltype_counts.to_csv('group_celltype_counts.csv')
print("已保存到: group_celltype_counts.csv")

# 3. 统计每个sample中每种celltype的数量
print("\n2. 每个sample中每种celltype的数量统计:")
# 先创建计数表
sample_celltype_counts = pd.crosstab(adata.obs['sample'], adata.obs['celltype'])
print(f"计数表形状: {sample_celltype_counts.shape}")
print(sample_celltype_counts.head())

# 获取每个sample的group信息
sample_info = adata.obs[['sample', 'group']].drop_duplicates().set_index('sample')

# 计算每个样本的总细胞数
sample_celltype_counts['total_cells'] = sample_celltype_counts.sum(axis=1)

# 添加group信息
sample_celltype_counts['group'] = sample_info.loc[sample_celltype_counts.index, 'group'].values

# 保存统计结果到CSV
sample_celltype_counts.to_csv('sample_celltype_counts.csv')
print("已保存到: sample_celltype_counts.csv")

# 4. 计算每个sample中每种celltype的比例
print("\n3. 计算每个sample中每种celltype的比例...")
# 获取细胞类型列（排除'total_cells'和'group'）
celltype_columns = [col for col in sample_celltype_counts.columns 
                    if col not in ['total_cells', 'group']]

# 创建比例数据框
sample_celltype_proportions = sample_celltype_counts[['group']].copy()

for col in celltype_columns:
    sample_celltype_proportions[col] = sample_celltype_counts[col] / sample_celltype_counts['total_cells'] * 100
    sample_celltype_proportions[f'{col}_count'] = sample_celltype_counts[col]  # 保存原始计数

# 保存比例数据
sample_celltype_proportions.to_csv('sample_celltype_proportions.csv')
print("已保存到: sample_celltype_proportions.csv")

# 5. 准备箱线图数据
print("\n4. 准备箱线图数据...")
# 创建长格式数据用于绘图
plot_data_list = []

for celltype in celltype_columns:
    for idx, row in sample_celltype_counts.iterrows():
        sample = idx
        group = row['group']
        count = row[celltype]
        total = row['total_cells']
        prop = (count / total * 100) if total > 0 else 0
        
        plot_data_list.append({
            'CellType': celltype,
            'Sample': sample,
            'Group': group,
            'Count': count,
            'Proportion': prop
        })

plot_data = pd.DataFrame(plot_data_list)

# 6. 绘制箱线图（按数量）
print("\n5. 绘制每个细胞类型的箱线图（按数量）...")
celltypes = plot_data['CellType'].unique()
n_celltypes = len(celltypes)

# 计算子图布局
n_cols = 4
n_rows = int(np.ceil(n_celltypes / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 4), squeeze=False)
axes = axes.flatten()

for idx, celltype in enumerate(celltypes):
    if idx < len(axes):
        ax = axes[idx]
        celltype_data = plot_data[plot_data['CellType'] == celltype]
        
        if not celltype_data.empty:
            # 绘制箱线图
            sns.boxplot(data=celltype_data, x='Group', y='Count', ax=ax, 
                       hue='Group', dodge=False, palette='Set2')
            
            # 添加散点显示个体样本
            sns.stripplot(data=celltype_data, x='Group', y='Count', ax=ax,
                         color='black', size=6, jitter=0.2, alpha=0.7)
            
            # 计算并添加组间显著性
            groups = celltype_data['Group'].dropna().unique()
            if len(groups) >= 2:
                anova_data = []
                valid_groups = []
                
                for g in groups:
                    group_counts = celltype_data[celltype_data['Group'] == g]['Count'].values
                    if len(group_counts) > 0:
                        anova_data.append(group_counts)
                        valid_groups.append(g)
                
                if len(anova_data) >= 2:
                    try:
                        h_stat, p_value = stats.kruskal(*anova_data)
                        
                        # 添加显著性标注
                        if p_value < 0.05:
                            y_max = celltype_data['Count'].max()
                            y_pos = y_max + 0.1 * y_max
                            
                            # 添加显著性星号
                            stars = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*'
                            ax.text(0.5, y_pos, f'p={p_value:.3e} {stars}', 
                                   ha='center', va='bottom', fontsize=9, fontweight='bold')
                    except Exception as e:
                        print(f"在{celltype}的显著性检验中出错: {e}")
            
            ax.set_title(f'{celltype}', fontsize=12, fontweight='bold')
            ax.set_xlabel('Group', fontsize=10)
            ax.set_ylabel('Cell Count', fontsize=10)
            ax.tick_params(axis='x', rotation=45)
            if ax.get_legend() is not None:
                ax.get_legend().remove()

# 隐藏多余的子图
for idx in range(len(celltypes), len(axes)):
    axes[idx].set_visible(False)

plt.suptitle('Cell Type Counts by Group (Boxplots with Individual Samples)', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('celltype_counts_boxplots.png', dpi=300, bbox_inches='tight')
plt.show()

# 7. 绘制箱线图（按比例）
print("\n6. 绘制每个细胞类型的箱线图（按比例）...")
fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 4), squeeze=False)
axes2 = axes2.flatten()

for idx, celltype in enumerate(celltypes):
    if idx < len(axes2):
        ax = axes2[idx]
        celltype_data = plot_data[plot_data['CellType'] == celltype]
        
        if not celltype_data.empty:
            # 绘制箱线图
            sns.boxplot(data=celltype_data, x='Group', y='Proportion', ax=ax, 
                       hue='Group', dodge=False, palette='Set2')
            
            # 添加散点显示个体样本
            sns.stripplot(data=celltype_data, x='Group', y='Proportion', ax=ax,
                         color='black', size=6, jitter=0.2, alpha=0.7)
            
            # 计算并添加组间显著性
            groups = celltype_data['Group'].dropna().unique()
            if len(groups) >= 2:
                anova_data = []
                valid_groups = []
                
                for g in groups:
                    group_props = celltype_data[celltype_data['Group'] == g]['Proportion'].values
                    if len(group_props) > 0:
                        anova_data.append(group_props)
                        valid_groups.append(g)
                
                if len(anova_data) >= 2:
                    try:
                        h_stat, p_value = stats.kruskal(*anova_data)
                        
                        if p_value < 0.05:
                            y_max = celltype_data['Proportion'].max()
                            y_pos = y_max + 0.1 * y_max
                            
                            stars = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*'
                            ax.text(0.5, y_pos, f'p={p_value:.3e} {stars}', 
                                   ha='center', va='bottom', fontsize=9, fontweight='bold')
                    except Exception as e:
                        print(f"在{celltype}的显著性检验中出错: {e}")
            
            ax.set_title(f'{celltype}', fontsize=12, fontweight='bold')
            ax.set_xlabel('Group', fontsize=10)
            ax.set_ylabel('Proportion (%)', fontsize=10)
            ax.tick_params(axis='x', rotation=45)
            if ax.get_legend() is not None:
                ax.get_legend().remove()

# 隐藏多余的子图
for idx in range(len(celltypes), len(axes2)):
    axes2[idx].set_visible(False)

plt.suptitle('Cell Type Proportions by Group (Boxplots with Individual Samples)', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('celltype_proportions_boxplots.png', dpi=300, bbox_inches='tight')
plt.show()

# 8. 创建汇总统计表
print("\n7. 生成汇总统计表...")
summary_stats = []

for celltype in celltypes:
    for group in adata.obs['group'].unique():
        celltype_data = plot_data[(plot_data['CellType'] == celltype) & 
                                  (plot_data['Group'] == group)]
        
        if not celltype_data.empty and len(celltype_data) > 0:
            counts = celltype_data['Count'].values
            props = celltype_data['Proportion'].values
            
            summary_stats.append({
                'CellType': celltype,
                'Group': group,
                'N_samples': len(celltype_data),
                'Mean_Count': np.mean(counts),
                'Std_Count': np.std(counts),
                'Median_Count': np.median(counts),
                'Min_Count': np.min(counts),
                'Max_Count': np.max(counts),
                'Mean_Prop': np.mean(props),
                'Std_Prop': np.std(props),
                'Median_Prop': np.median(props),
                'Min_Prop': np.min(props),
                'Max_Prop': np.max(props)
            })

if summary_stats:
    summary_df = pd.DataFrame(summary_stats)
    summary_df = summary_df.sort_values(['CellType', 'Group'])
    summary_df.to_csv('celltype_summary_statistics.csv', index=False)
    print("汇总统计已保存到: celltype_summary_statistics.csv")
else:
    print("警告: 汇总统计表为空!")

# 9. 可选：绘制热图展示细胞类型组成
print("\n8. 绘制细胞类型组成热图...")
if 'summary_df' in locals() and not summary_df.empty:
    # 计算每个group的平均比例
    heatmap_data = summary_df.pivot(index='CellType', columns='Group', values='Mean_Prop')
    heatmap_data = heatmap_data.fillna(0)
    
    plt.figure(figsize=(10, max(8, len(celltypes) * 0.5)))
    sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='YlOrRd', 
               cbar_kws={'label': 'Mean Proportion (%)'}, linewidths=0.5)
    plt.title('Mean Cell Type Proportions by Group', fontsize=14, fontweight='bold')
    plt.xlabel('Group', fontsize=12)
    plt.ylabel('Cell Type', fontsize=12)
    plt.tight_layout()
    plt.savefig('celltype_composition_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()
else:
    print("无法绘制热图: 汇总数据缺失")

# 10. 可选：绘制堆叠柱状图
print("\n9. 绘制细胞类型组成堆叠柱状图...")
if 'summary_df' in locals() and not summary_df.empty:
    group_means = summary_df.pivot(index='Group', columns='CellType', values='Mean_Prop')
    
    plt.figure(figsize=(12, 6))
    group_means.plot(kind='bar', stacked=True, colormap='tab20c', edgecolor='black')
    plt.title('Cell Type Composition by Group (Mean Proportions)', fontsize=14, fontweight='bold')
    plt.xlabel('Group', fontsize=12)
    plt.ylabel('Proportion (%)', fontsize=12)
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('celltype_composition_stacked_bar.png', dpi=300, bbox_inches='tight')
    plt.show()
else:
    print("无法绘制堆叠柱状图: 汇总数据缺失")

print("\n分析完成！")
print("生成的文件:")
print("1. group_celltype_counts.csv - 每组细胞类型数量")
print("2. sample_celltype_counts.csv - 每个样本细胞类型数量")
print("3. sample_celltype_proportions.csv - 每个样本细胞类型比例")
print("4. celltype_summary_statistics.csv - 汇总统计")
print("5. celltype_counts_boxplots.png - 数量箱线图")
print("6. celltype_proportions_boxplots.png - 比例箱线图")
if 'summary_df' in locals() and not summary_df.empty:
    print("7. celltype_composition_heatmap.png - 组成热图")
    print("8. celltype_composition_stacked_bar.png - 堆叠柱状图")
