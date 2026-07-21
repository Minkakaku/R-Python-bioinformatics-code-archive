import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 1. 读取h5ad数据
print("读取h5ad数据...")
adata = sc.read_h5ad("adata_celltype.h5ad")  # 请替换为您的实际文件名
print(f"数据形状: {adata.shape}")
print(f"细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
print(f"obs列: {adata.obs.columns.tolist()}")

# 检查必需的列
if 'sample' not in adata.obs.columns:
    print("错误: 数据中缺少 'sample' 列")
    sys.exit(1)
    
if 'celltype' not in adata.obs.columns:
    print("错误: 数据中缺少 'celltype' 列")
    sys.exit(1)

# 2. 添加group_final列
print("\n添加group_final列...")

# 定义映射字典
sample_to_group_final = {
    "PRE_5-1": "PRE-R", "PRE_5-2": "PRE-R", "PRE_6-1": "PRE-R", 
    "PRE_6-2": "PRE-R", "PRE_6-3": "PRE-R", "PRE_8": "PRE-R", 
    "PRE_15-1": "PRE-R", "PRE_15-2": "PRE-R", "PRE_15-3": "PRE-R", 
    "PRE_15-4": "PRE-R", "PRE_16": "PRE-R", "PRE_17-1": "PRE-R", 
    "PRE_17-2": "PRE-R", "PRE_17-3": "PRE-R", "PRE_18-1": "PRE-R", 
    "PRE_18-2": "PRE-R", "PRE_19": "PRE-R", "PRE_3-1": "PRE-NR", 
    "PRE_3-2": "PRE-NR", "PRE_10-1": "PRE-NR", "PRE_10-2": "PRE-NR", 
    "PRE_20-1": "PRE-NR", "PRE_20-2": "PRE-NR", "PRE_20-3": "PRE-NR", 
    "PRE_20-4": "PRE-NR", "PRE_21": "PRE-NR", "PRE_22-1": "PRE-NR", 
    "PRE_22-2": "PRE-NR", "PRE_23-1": "PRE-NR", "PRE_23-2": "PRE-NR", 
    "PRE_25-1": "PRE-NR", "PRE_25-2": "PRE-NR", "PRE_25-3": "PRE-NR", 
    "PRE_25-4": "PRE-NR", "PRE_26-1": "PRE-NR", "PRE_26-3": "PRE-NR", 
    "PRE_27-1": "PRE-NR", "PRE_27-2": "PRE-NR", "POST_5": "POST-R", 
    "POST_6": "POST-R", "POST_8": "POST-R", "POST_15": "POST-R", 
    "POST_16": "POST-R", "POST_18": "POST-R", "POST_19": "POST-R", 
    "POST_3": "POST-NR", "POST_20": "POST-NR", "POST_21": "POST-NR", 
    "POST_22-1": "POST-NR", "POST_25": "POST-NR", "POST_26": "POST-NR", 
    "POST_27": "POST-NR"
}

# 应用映射
adata.obs['group_final'] = adata.obs['sample'].map(sample_to_group_final)

# 检查未映射的样本
unmapped_samples = adata.obs[adata.obs['group_final'].isna()]['sample'].unique()
if len(unmapped_samples) > 0:
    print(f"警告: 以下样本未找到group_final映射关系: {unmapped_samples.tolist()}")
    print(f"未映射的细胞数: {adata.obs['group_final'].isna().sum()}")
    # 移除未映射的样本
    adata = adata[~adata.obs['group_final'].isna(), :].copy()

# 统计group_final分布
print(f"\ngroup_final分布:")
group_counts = adata.obs['group_final'].value_counts()
for group, count in group_counts.items():
    print(f"  {group}: {count} 个细胞 ({count/len(adata)*100:.2f}%)")

# 3. 以样本为维度计算细胞类型比例
print("\n计算每个样本的细胞类型比例...")

# 创建数据副本
df = adata.obs[['sample', 'celltype', 'group_final']].copy()

# 计算每个样本中每个细胞类型的细胞数
sample_celltype_counts = df.groupby(['sample', 'group_final', 'celltype']).size().reset_index(name='count')

# 计算每个样本的总细胞数
sample_totals = df.groupby(['sample', 'group_final']).size().reset_index(name='total')

# 合并数据并计算比例
sample_data = pd.merge(sample_celltype_counts, sample_totals, on=['sample', 'group_final'])
sample_data['proportion'] = (sample_data['count'] / sample_data['total'] * 100).round(4)

# 显示样本数
unique_samples = sample_data['sample'].nunique()
print(f"样本数: {unique_samples}")

# 显示细胞类型数
unique_celltypes = sample_data['celltype'].nunique()
print(f"细胞类型数: {unique_celltypes}")

# 保存样本维度的数据
sample_data.to_csv("sample_celltype_proportions.csv", index=False, encoding='utf-8-sig')
print(f"样本细胞类型比例数据已保存到: sample_celltype_proportions.csv")

# 4. 创建输出目录
output_dir = "celltype_proportion_boxplots"
os.makedirs(output_dir, exist_ok=True)

# 5. 统计检验函数
def perform_statistical_test(data, celltype, group_col='group_final'):
    """执行统计检验，返回p值和检验方法"""
    # 按组提取数据
    groups = data[group_col].unique()
    group_data = []
    
    for group in groups:
        subset = data[(data[group_col] == group) & (data['celltype'] == celltype)]
        group_data.append(subset['proportion'].values)
    
    # 检查每组是否有足够的数据
    valid_groups = [g for g in group_data if len(g) >= 2]
    if len(valid_groups) < 2:
        return None, None, False
    
    # 如果有两组，使用Mann-Whitney U检验
    if len(valid_groups) == 2:
        stat, p_value = stats.mannwhitneyu(valid_groups[0], valid_groups[1], alternative='two-sided')
        test_method = "Mann-Whitney U"
        return test_method, p_value, p_value < 0.05
    # 如果有多组，使用Kruskal-Wallis检验
    else:
        stat, p_value = stats.kruskal(*valid_groups)
        test_method = "Kruskal-Wallis"
        return test_method, p_value, p_value < 0.05

# 6. 绘制每个细胞类型的箱线图
print(f"\n绘制细胞类型比例的箱线图...")
celltypes = df['celltype'].unique()

# 创建汇总表存储统计结果
stats_results = []

for i, celltype in enumerate(celltypes, 1):
    print(f"处理细胞类型 {i}/{len(celltypes)}: {celltype}")
    
    # 过滤当前细胞类型的数据
    celltype_data = sample_data[sample_data['celltype'] == celltype].copy()
    
    if len(celltype_data) == 0:
        continue
    
    # 检查有多少个组
    groups = celltype_data['group_final'].unique()
    if len(groups) < 2:
        print(f"  跳过: 只有{len(groups)}个组的数据")
        continue
    
    # 执行统计检验
    test_method, p_value, is_significant = perform_statistical_test(sample_data, celltype)
    
    # 保存统计结果
    if p_value is not None:
        stats_results.append({
            'celltype': celltype,
            'test_method': test_method,
            'p_value': p_value,
            'significant': is_significant
        })
    
    # 创建图形
    plt.figure(figsize=(max(8, len(groups)*1.5), 6))
    
    # 绘制箱线图
    ax = sns.boxplot(x='group_final', y='proportion', data=celltype_data, 
                     palette='Set2', width=0.6, fliersize=3)
    
    # 添加散点图显示每个样本
    sns.stripplot(x='group_final', y='proportion', data=celltype_data, 
                  color='black', size=5, alpha=0.7, jitter=True)
    
    # 添加标题和标签
    title = f"Cell Type: {celltype}"
    if p_value is not None and is_significant:
        title += f"\n{test_method} test: p = {p_value:.3e}"
    
    plt.title(title, fontsize=12, fontweight='bold')
    plt.xlabel('Group', fontsize=11)
    plt.ylabel('Proportion (%)', fontsize=11)
    
    # 添加样本数量标签
    for j, group in enumerate(groups):
        group_subset = celltype_data[celltype_data['group_final'] == group]
        n_samples = len(group_subset)
        median_prop = group_subset['proportion'].median()
        
        # 在箱线图上方添加标签
        y_max = celltype_data['proportion'].max()
        plt.text(j, y_max * 1.1, f'n={n_samples}\nmedian={median_prop:.2f}%', 
                ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    # 保存图形
    safe_celltype = celltype.replace('/', '_').replace('\\', '_').replace(':', '_')
    output_path = os.path.join(output_dir, f"{safe_celltype}_proportion.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

# 7. 保存统计检验结果
if stats_results:
    stats_df = pd.DataFrame(stats_results)
    # 按p值排序
    stats_df = stats_df.sort_values('p_value')
    stats_df.to_csv("celltype_proportion_stats.csv", index=False, encoding='utf-8-sig')
    print(f"\n统计检验结果已保存到: celltype_proportion_stats.csv")
    
    # 显示显著差异的细胞类型
    significant_celltypes = stats_df[stats_df['significant']]
    if len(significant_celltypes) > 0:
        print(f"\n显著差异的细胞类型 ({len(significant_celltypes)}个):")
        for _, row in significant_celltypes.iterrows():
            print(f"  {row['celltype']}: p = {row['p_value']:.3e} ({row['test_method']})")

# 8. 生成汇总报告
print(f"\n汇总报告:")
print(f"  总细胞类型数: {len(celltypes)}")
print(f"  分析的细胞类型数: {len([ct for ct in celltypes if len(sample_data[sample_data['celltype'] == ct]) > 0])}")
print(f"  样本数: {unique_samples}")
print(f"  分组: {sorted(df['group_final'].unique().tolist())}")
print(f"  图形保存到目录: {output_dir}/")
print(f"  每个细胞类型生成一个PNG文件")

# 9. 保存更新后的h5ad数据
output_h5ad = "data_with_group_final.h5ad"
adata.write_h5ad(output_h5ad, compression='gzip')
print(f"\n更新后的h5ad数据已保存到: {output_h5ad}")

print("\n分析完成!")
