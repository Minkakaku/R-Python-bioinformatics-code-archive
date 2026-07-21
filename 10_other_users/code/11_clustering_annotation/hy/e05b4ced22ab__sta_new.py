import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置字体（使用英文字体，避免中文显示问题）
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 创建输出目录
sample_dir = "sample_com"
patient_dir = "patient_com"
os.makedirs(sample_dir, exist_ok=True)
os.makedirs(patient_dir, exist_ok=True)

# 1. 读取h5ad文件
print("读取h5ad文件...")
adata = sc.read_h5ad("All_filter_state_samp_new.h5ad")
print(f"数据形状: {adata.shape}")
print(f"obs列: {adata.obs.columns.tolist()}")

# 检查必需列
required_cols = ['sample_new', 'state', 'celltype', 'sample']
for col in required_cols:
    if col not in adata.obs.columns:
        print(f"错误: 缺少必需列 '{col}'")
        print(f"可用列: {adata.obs.columns.tolist()}")
        exit(1)

# 创建数据副本
df = adata.obs[required_cols].copy()

# 检查数据
print(f"\n数据基本信息:")
print(f"总细胞数: {len(df)}")
print(f"sample_new 唯一值: {df['sample_new'].nunique()}")
print(f"state 唯一值: {df['state'].nunique()}")
print(f"sample 唯一值: {df['sample'].nunique()}")
print(f"celltype 唯一值: {df['celltype'].nunique()}")

# 2. 统计所有细胞类型的数量和占比
print("\n" + "="*50)
print("细胞类型总体统计")
print("="*50)

total_cells = len(df)
celltype_stats = df['celltype'].value_counts().reset_index()
celltype_stats.columns = ['celltype', 'count']
celltype_stats['percentage'] = (celltype_stats['count'] / total_cells * 100).round(2)

print(f"总细胞数: {total_cells}")
print(f"细胞类型数: {len(celltype_stats)}")
print(f"state分组: {sorted(df['state'].dropna().unique().tolist())}")
print(f"sample数量: {len(df['sample'].unique())}")
print("\n各细胞类型统计:")
for i, row in celltype_stats.iterrows():
    print(f"  {i+1:2d}. {row['celltype']:20s}: {row['count']:6d} 个细胞 ({row['percentage']:6.2f}%)")

# 保存总体统计
celltype_stats.to_csv("celltype_overall_statistics.csv", index=False, encoding='utf-8-sig')
print(f"\n总体统计已保存到: celltype_overall_statistics.csv")

# 3. 以样本为维度的分析（使用sample_new和state）
print("\n" + "="*50)
print("以样本为维度的分析 (使用sample_new和state)")
print("="*50)

# 3.1 获取sample_new到state的映射
sample_state_mapping = df.drop_duplicates(['sample_new', 'state'])[['sample_new', 'state']]

# 3.2 计算每个样本中每个细胞类型的比例
sample_cell_counts = df.groupby(['sample_new', 'celltype']).size().unstack(fill_value=0)
sample_totals = sample_cell_counts.sum(axis=1)
sample_proportions = sample_cell_counts.div(sample_totals, axis=0) * 100  # 转换为百分比

# 3.3 将state信息添加到sample_proportions
sample_proportions = sample_proportions.reset_index()
# 合并state信息
sample_proportions = sample_proportions.merge(sample_state_mapping, on='sample_new', how='left')

# 3.4 重新组织数据用于绘图
# 获取细胞类型列
celltype_cols = [col for col in sample_proportions.columns if col not in ['sample_new', 'state']]

sample_long = sample_proportions.melt(
    id_vars=['sample_new', 'state'], 
    value_vars=celltype_cols,
    var_name='celltype', 
    value_name='proportion'
).dropna(subset=['state', 'proportion'])  # 移除state或proportion为NaN的行

# 移除比例为0的行
sample_long = sample_long[sample_long['proportion'] > 0]

# 3.5 保存样本维度的数据
sample_summary = sample_proportions.copy()
sample_summary['total_cells'] = sample_summary[celltype_cols].sum(axis=1)
sample_summary.to_csv(os.path.join(sample_dir, "sample_celltype_proportions.csv"), 
                     index=False, encoding='utf-8-sig')
print(f"样本维度数据已保存到: {sample_dir}/sample_celltype_proportions.csv")

# 4. 以患者为维度的分析（使用sample和state）
print("\n" + "="*50)
print("以患者为维度的分析 (使用sample和state)")
print("="*50)

# 4.1 获取sample到state的映射
patient_state_mapping = df.drop_duplicates(['sample', 'state'])[['sample', 'state']]

# 4.2 计算每个患者中每个细胞类型的比例
patient_cell_counts = df.groupby(['sample', 'celltype']).size().unstack(fill_value=0)
patient_totals = patient_cell_counts.sum(axis=1)
patient_proportions = patient_cell_counts.div(patient_totals, axis=0) * 100  # 转换为百分比

# 4.3 将state信息添加到patient_proportions
patient_proportions = patient_proportions.reset_index()
# 合并state信息
patient_proportions = patient_proportions.merge(patient_state_mapping, on='sample', how='left')

# 4.4 重新组织数据用于绘图
# 获取细胞类型列
patient_celltype_cols = [col for col in patient_proportions.columns if col not in ['sample', 'state']]

patient_long = patient_proportions.melt(
    id_vars=['sample', 'state'], 
    value_vars=patient_celltype_cols,
    var_name='celltype', 
    value_name='proportion'
).dropna(subset=['state', 'proportion'])  # 移除state或proportion为NaN的行

# 移除比例为0的行
patient_long = patient_long[patient_long['proportion'] > 0]

# 4.5 保存患者维度的数据
patient_summary = patient_proportions.copy()
patient_summary['total_cells'] = patient_summary[patient_celltype_cols].sum(axis=1)
patient_summary.to_csv(os.path.join(patient_dir, "patient_celltype_proportions.csv"), 
                      index=False, encoding='utf-8-sig')
print(f"患者维度数据已保存到: {patient_dir}/patient_celltype_proportions.csv")

# 5. 统计检验函数
def perform_statistical_tests(data_long, celltype, group_col='state', value_col='proportion'):
    """
    执行统计检验
    返回: 检验方法, p值, 是否显著, 组间比较结果
    """
    # 按组提取数据
    groups = data_long[group_col].unique()
    group_data = {}
    
    for group in groups:
        subset = data_long[(data_long[group_col] == group) & (data_long['celltype'] == celltype)]
        group_data[group] = subset[value_col].values
    
    # 检查是否有足够的数据
    valid_groups = {g: d for g, d in group_data.items() if len(d) >= 2}
    if len(valid_groups) < 2:
        return None, None, False, "数据不足"
    
    # 提取有效组的数据
    data_list = [group_data[g] for g in sorted(valid_groups.keys())]
    group_names = sorted(valid_groups.keys())
    
    # 检查方差齐性
    if len(data_list) == 2:
        # 两组比较：Mann-Whitney U检验（非参数）
        stat, p_value = stats.mannwhitneyu(data_list[0], data_list[1], alternative='two-sided')
        test_method = "Mann-Whitney U"
        
        # 计算效应量
        n1, n2 = len(data_list[0]), len(data_list[1])
        effect_size = 2 * stat / (n1 * n2) - 1
        
        return test_method, p_value, p_value < 0.05, {
            'groups': group_names,
            'statistic': stat,
            'effect_size': effect_size
        }
    else:
        # 多组比较：Kruskal-Wallis检验
        stat, p_value = stats.kruskal(*data_list)
        test_method = "Kruskal-Wallis"
        
        if p_value < 0.05:
            # 事后两两比较
            pairwise_results = []
            for i in range(len(group_names)):
                for j in range(i+1, len(group_names)):
                    stat_pair, p_pair = stats.mannwhitneyu(data_list[i], data_list[j], alternative='two-sided')
                    pairwise_results.append({
                        'group1': group_names[i],
                        'group2': group_names[j],
                        'p_value': p_pair
                    })
            return test_method, p_value, p_value < 0.05, {
                'groups': group_names,
                'pairwise': pairwise_results
            }
        else:
            return test_method, p_value, False, {'groups': group_names}

# 6. 绘图函数
def plot_celltype_proportions(data_long, celltype, output_dir, level="sample", p_value_info=None, group_col='state'):
    """
    绘制单个细胞类型的箱线图
    """
    # 过滤特定细胞类型的数据
    celltype_data = data_long[data_long['celltype'] == celltype].copy()
    
    if len(celltype_data) == 0:
        print(f"  警告: 细胞类型 '{celltype}' 在{level}维度无数据")
        return None
    
    # 获取组名并排序
    groups = sorted(celltype_data[group_col].unique())
    
    if len(groups) < 2:
        print(f"  警告: 细胞类型 '{celltype}' 在{level}维度只有{len(groups)}个组，跳过绘图")
        return None
    
    # 创建图形
    plt.figure(figsize=(max(6, len(groups)*1.5), 6))
    
    # 绘制箱线图和小提琴图
    ax = sns.boxplot(x=group_col, y='proportion', data=celltype_data, 
                     palette='Set2', width=0.5, fliersize=3)
    sns.stripplot(x=group_col, y='proportion', data=celltype_data, 
                  color='black', size=4, alpha=0.6, jitter=True)
    
    # 添加统计信息
    title_parts = [f"Cell Type: {celltype}"]  # 使用英文标题
    if p_value_info and p_value_info.get('significant', False):
        test_method = p_value_info.get('test_method', '')
        p_val = p_value_info.get('p_value', 1)
        title_parts.append(f"{test_method} test: p={p_val:.3e}")
    
    plt.title('\n'.join(title_parts), fontsize=12, fontweight='bold')
    plt.xlabel('State', fontsize=11)  # 使用英文标签
    plt.ylabel('Cell Proportion (%)', fontsize=11)  # 使用英文标签
    
    # 调整y轴范围
    y_max = celltype_data['proportion'].max()
    plt.ylim(0, y_max * 1.2)
    
    # 添加样本/患者数量标签
    for i, group in enumerate(groups):
        group_data_subset = celltype_data[celltype_data[group_col] == group]
        n = len(group_data_subset)
        median_val = group_data_subset['proportion'].median()
        
        # 在x轴上方添加数量标签
        plt.text(i, y_max * 1.1, f'n={n}\nmedian={median_val:.2f}%', 
                ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    # 保存图形
    safe_celltype = celltype.replace('/', '_').replace('\\', '_').replace(':', '_')
    output_file = os.path.join(output_dir, f"{level}_{safe_celltype}_proportion.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

# 7. 执行样本维度的分析和绘图
print("\n" + "="*50)
print("生成样本维度箱线图 (基于state分组)")
print("="*50)

sample_stats_results = []
celltypes = df['celltype'].unique()

for i, celltype in enumerate(celltypes, 1):
    print(f"处理细胞类型 {i}/{len(celltypes)}: {celltype}")
    
    # 执行统计检验
    test_method, p_value, is_sig, details = perform_statistical_tests(
        sample_long, celltype, 'state', 'proportion'
    )
    
    # 保存统计结果
    if p_value is not None:
        sample_stats_results.append({
            'celltype': celltype,
            'test_method': test_method,
            'p_value': p_value,
            'significant': is_sig,
            'details': str(details)
        })
    
    # 绘制图形
    p_value_info = {
        'test_method': test_method,
        'p_value': p_value,
        'significant': is_sig
    } if p_value is not None else None
    
    plot_celltype_proportions(sample_long, celltype, sample_dir, "sample", p_value_info, 'state')

# 保存样本维度统计结果
if sample_stats_results:
    sample_stats_df = pd.DataFrame(sample_stats_results)
    sample_stats_df.to_csv(os.path.join(sample_dir, "sample_statistical_tests.csv"), 
                          index=False, encoding='utf-8-sig')
    print(f"样本维度统计检验结果已保存到: {sample_dir}/sample_statistical_tests.csv")

# 8. 执行患者维度的分析和绘图
print("\n" + "="*50)
print("生成患者维度箱线图 (基于state分组)")
print("="*50)

patient_stats_results = []

for i, celltype in enumerate(celltypes, 1):
    print(f"处理细胞类型 {i}/{len(celltypes)}: {celltype}")
    
    # 执行统计检验
    test_method, p_value, is_sig, details = perform_statistical_tests(
        patient_long, celltype, 'state', 'proportion'
    )
    
    # 保存统计结果
    if p_value is not None:
        patient_stats_results.append({
            'celltype': celltype,
            'test_method': test_method,
            'p_value': p_value,
            'significant': is_sig,
            'details': str(details)
        })
    
    # 绘制图形
    p_value_info = {
        'test_method': test_method,
        'p_value': p_value,
        'significant': is_sig
    } if p_value is not None else None
    
    plot_celltype_proportions(patient_long, celltype, patient_dir, "patient", p_value_info, 'state')

# 保存患者维度统计结果
if patient_stats_results:
    patient_stats_df = pd.DataFrame(patient_stats_results)
    patient_stats_df.to_csv(os.path.join(patient_dir, "patient_statistical_tests.csv"), 
                           index=False, encoding='utf-8-sig')
    print(f"患者维度统计检验结果已保存到: {patient_dir}/patient_statistical_tests.csv")

# 9. 生成汇总报告
print("\n" + "="*50)
print("生成汇总报告")
print("="*50)

with open("celltype_analysis_summary.txt", "w", encoding='utf-8') as f:
    f.write("细胞类型比例分析汇总报告\n")
    f.write("="*60 + "\n\n")
    
    f.write("1. 数据概览\n")
    f.write("-"*40 + "\n")
    f.write(f"总细胞数: {total_cells}\n")
    f.write(f"细胞类型数: {len(celltype_stats)}\n")
    f.write(f"样本数 (sample_new): {sample_proportions.shape[0]}\n")
    f.write(f"患者数 (sample): {patient_proportions.shape[0]}\n")
    f.write(f"state分组: {sorted(df['state'].dropna().unique().tolist())}\n")
    f.write(f"样本数量: {len(df['sample_new'].unique())}\n")
    f.write(f"患者数量: {len(df['sample'].unique())}\n\n")
    
    f.write("2. 细胞类型总体分布\n")
    f.write("-"*40 + "\n")
    for _, row in celltype_stats.iterrows():
        f.write(f"{row['celltype']:25s}: {row['count']:6d} 个细胞 ({row['percentage']:6.2f}%)\n")
    
    f.write("\n3. 样本维度统计 (基于state)\n")
    f.write("-"*40 + "\n")
    if sample_stats_results:
        sig_count = sum(1 for r in sample_stats_results if r['significant'])
        f.write(f"显著差异的细胞类型数: {sig_count}/{len(sample_stats_results)}\n")
        for result in sample_stats_results:
            if result['significant']:
                f.write(f"  {result['celltype']}: p={result['p_value']:.3e} ({result['test_method']})\n")
    
    f.write("\n4. 患者维度统计 (基于state)\n")
    f.write("-"*40 + "\n")
    if patient_stats_results:
        sig_count = sum(1 for r in patient_stats_results if r['significant'])
        f.write(f"显著差异的细胞类型数: {sig_count}/{len(patient_stats_results)}\n")
        for result in patient_stats_results:
            if result['significant']:
                f.write(f"  {result['celltype']}: p={result['p_value']:.3e} ({result['test_method']})\n")
    
    f.write("\n5. 输出文件列表\n")
    f.write("-"*40 + "\n")
    f.write("总体统计:\n")
    f.write("  celltype_overall_statistics.csv\n\n")
    f.write("样本维度输出:\n")
    f.write(f"  {sample_dir}/sample_celltype_proportions.csv\n")
    f.write(f"  {sample_dir}/sample_statistical_tests.csv\n")
    f.write(f"  {sample_dir}/sample_*.png (箱线图)\n\n")
    f.write("患者维度输出:\n")
    f.write(f"  {patient_dir}/patient_celltype_proportions.csv\n")
    f.write(f"  {patient_dir}/patient_statistical_tests.csv\n")
    f.write(f"  {patient_dir}/patient_*.png (箱线图)\n")

print(f"汇总报告已保存到: celltype_analysis_summary.txt")
print("\n分析完成!")
