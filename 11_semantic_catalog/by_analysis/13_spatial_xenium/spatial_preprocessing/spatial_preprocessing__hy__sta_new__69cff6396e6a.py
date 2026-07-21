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
adata = sc.read_h5ad("all_celltype.h5ad")
print(f"数据形状: {adata.shape}")
print(f"obs列: {adata.obs.columns.tolist()}")

# 检查必需列
required_cols = ['sample', 'group_new', 'celltype', 'patient']
for col in required_cols:
    if col not in adata.obs.columns:
        print(f"错误: 缺少必需列 '{col}'")
        print(f"可用列: {adata.obs.columns.tolist()}")
        exit(1)

# 创建数据副本
df = adata.obs[required_cols].copy()

# 2. 根据映射关系重新赋值group_final
print("\n" + "="*50)
print("重新映射group_final和patient_new")
print("="*50)

# 定义样本到group_final的映射关系
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

# 定义样本到patient_new的映射关系
sample_to_patient_new = {
    "POST_3": "P3", "POST_8": "P8", "PRE_8": "P8", "POST_6": "P6", 
    "PRE_3-1": "P3", "PRE_3-2": "P3", "PRE_6-1": "P6", "PRE_6-2": "P6", 
    "PRE_5-1": "P5", "PRE_5-2": "P5", "PRE_10-1": "P10", "PRE_10-2": "P10", 
    "POST_5": "P5", "POST_16": "P16", "POST_19": "P19", "POST_22-1": "P22", 
    "POST_22-2": "P22", "PRE_15-1": "P15", "PRE_15-2": "P15", "PRE_15-3": "P15", 
    "PRE_15-4": "P15", "PRE_16": "P16", "PRE_17-1": "P17", "PRE_17-2": "P17", 
    "PRE_17-3": "P17", "PRE_20-1": "P20", "PRE_20-2": "P20", "PRE_21": "P21", 
    "PRE_22-1": "P22", "PRE_22-2": "P22", "PRE_23-1": "P23", "PRE_23-2": "P23", 
    "PRE_24-1": "P24", "PRE_24-2": "P24", "PRE_26-1": "P26", "PRE_26-2": "P26", 
    "PRE_26-3": "P26", "POST_15": "P15", "POST_18": "P18", "POST_20": "P20", 
    "POST_21": "P21", "POST_25": "P25", "POST_26": "P26", "POST_27": "P27", 
    "PRE_18-1": "P18", "PRE_18-2": "P18", "PRE_19": "P19", "PRE_25-1": "P25", 
    "PRE_25-2": "P25", "PRE_25-3": "P25", "PRE_25-4": "P25", "PRE_27-1": "P27", 
    "PRE_27-2": "P27"
}

# 记录原始group_new分布
print("原始group_new分布:")
print(df['group_new'].value_counts())

# 应用映射
original_sample_count = len(df)
df['group_final'] = df['sample'].map(sample_to_group_final)
df['patient_new'] = df['sample'].map(sample_to_patient_new)

# 检查未映射的样本
unmapped_group = df[df['group_final'].isna()]['sample'].unique()
unmapped_patient = df[df['patient_new'].isna()]['sample'].unique()

if len(unmapped_group) > 0:
    print(f"\n警告: 以下样本未找到group_final映射关系: {unmapped_group.tolist()}")
    print("将保留原始group_new值作为group_final")
    df.loc[df['group_final'].isna(), 'group_final'] = df.loc[df['group_final'].isna(), 'group_new']

if len(unmapped_patient) > 0:
    print(f"\n警告: 以下样本未找到patient_new映射关系: {unmapped_patient.tolist()}")
    print("将保留原始patient值作为patient_new")
    df.loc[df['patient_new'].isna(), 'patient_new'] = df.loc[df['patient_new'].isna(), 'patient']

# 检查映射后的结果
print(f"\n映射后group_final分布:")
print(df['group_final'].value_counts())
print(f"\n映射后patient_new分布:")
print(df['patient_new'].value_counts())

# 检查是否有空的group_final或patient_new
if df['group_final'].isna().any():
    print(f"警告: 仍有 {df['group_final'].isna().sum()} 个细胞的group_final为空")
    df = df.dropna(subset=['group_final'])

if df['patient_new'].isna().any():
    print(f"警告: 仍有 {df['patient_new'].isna().sum()} 个细胞的patient_new为空")
    df = df.dropna(subset=['patient_new'])

print(f"处理后细胞数: {len(df)} (原始: {original_sample_count})")

# 3. 统计所有细胞类型的数量和占比
print("\n" + "="*50)
print("细胞类型总体统计")
print("="*50)

total_cells = len(df)
celltype_stats = df['celltype'].value_counts().reset_index()
celltype_stats.columns = ['celltype', 'count']
celltype_stats['percentage'] = (celltype_stats['count'] / total_cells * 100).round(2)

print(f"总细胞数: {total_cells}")
print(f"细胞类型数: {len(celltype_stats)}")
print(f"group_final分组: {sorted(df['group_final'].unique().tolist())}")
print(f"patient_new数量: {len(df['patient_new'].unique())}")
print("\n各细胞类型统计:")
for i, row in celltype_stats.iterrows():
    print(f"  {i+1:2d}. {row['celltype']:20s}: {row['count']:6d} 个细胞 ({row['percentage']:6.2f}%)")

# 保存总体统计
celltype_stats.to_csv("celltype_overall_statistics.csv", index=False, encoding='utf-8-sig')
print(f"\n总体统计已保存到: celltype_overall_statistics.csv")

# 4. 以样本为维度的分析（使用group_final作为分组）
print("\n" + "="*50)
print("以样本为维度的分析 (使用group_final分组)")
print("="*50)

# 4.1 计算每个样本中每个细胞类型的比例
sample_cell_counts = df.groupby(['sample', 'celltype']).size().unstack(fill_value=0)
sample_totals = sample_cell_counts.sum(axis=1)
sample_proportions = sample_cell_counts.div(sample_totals, axis=0) * 100  # 转换为百分比

# 4.2 添加group_final信息
# 确保每个sample只有一个group_final
sample_groups = df[['sample', 'group_final']].drop_duplicates()
# 检查是否有重复的sample
if sample_groups['sample'].duplicated().any():
    print("警告: sample存在重复的group_final，使用第一个出现的group_final")
    sample_groups = sample_groups.drop_duplicates(subset='sample', keep='first')

# 创建映射字典
group_mapping = dict(zip(sample_groups['sample'], sample_groups['group_final']))
sample_proportions['group_final'] = sample_proportions.index.map(group_mapping)

print(f"样本数: {sample_proportions.shape[0]}")
print(f"group_final分组: {sorted(sample_proportions['group_final'].dropna().unique().tolist())}")

# 4.3 重新组织数据用于绘图
sample_long = sample_proportions.reset_index().melt(
    id_vars=['sample', 'group_final'], 
    value_vars=[col for col in sample_proportions.columns if col not in ['sample', 'group_final']],
    var_name='celltype', 
    value_name='proportion'
).dropna(subset=['group_final'])  # 移除group_final为NaN的行

# 移除比例为0的行
sample_long = sample_long[sample_long['proportion'] > 0]

# 4.4 保存样本维度的数据
sample_summary = sample_proportions.copy()
sample_summary['total_cells'] = sample_totals
sample_summary.to_csv(os.path.join(sample_dir, "sample_celltype_proportions.csv"), encoding='utf-8-sig')
print(f"样本维度数据已保存到: {sample_dir}/sample_celltype_proportions.csv")

# 5. 以患者为维度的分析（使用patient_new和group_final）
print("\n" + "="*50)
print("以患者为维度的分析 (使用patient_new和group_final)")
print("="*50)

# 5.1 计算每个患者中每个细胞类型的比例
patient_cell_counts = df.groupby(['patient_new', 'celltype']).size().unstack(fill_value=0)
patient_totals = patient_cell_counts.sum(axis=1)
patient_proportions = patient_cell_counts.div(patient_totals, axis=0) * 100  # 转换为百分比

# 5.2 添加group_final信息
# 确保每个patient_new只有一个group_final
patient_groups = df[['patient_new', 'group_final']].drop_duplicates()
# 检查是否有重复的patient_new
if patient_groups['patient_new'].duplicated().any():
    print("警告: patient_new存在重复的group_final，使用第一个出现的group_final")
    patient_groups = patient_groups.drop_duplicates(subset='patient_new', keep='first')

# 创建映射字典
patient_group_mapping = dict(zip(patient_groups['patient_new'], patient_groups['group_final']))
patient_proportions['group_final'] = patient_proportions.index.map(patient_group_mapping)

print(f"患者数: {patient_proportions.shape[0]}")
print(f"group_final分组: {sorted(patient_proportions['group_final'].dropna().unique().tolist())}")

# 5.3 重新组织数据用于绘图
patient_long = patient_proportions.reset_index().melt(
    id_vars=['patient_new', 'group_final'], 
    value_vars=[col for col in patient_proportions.columns if col not in ['patient_new', 'group_final']],
    var_name='celltype', 
    value_name='proportion'
).dropna(subset=['group_final'])  # 移除group_final为NaN的行

# 移除比例为0的行
patient_long = patient_long[patient_long['proportion'] > 0]

# 5.4 保存患者维度的数据
patient_summary = patient_proportions.copy()
patient_summary['total_cells'] = patient_totals
patient_summary.to_csv(os.path.join(patient_dir, "patient_celltype_proportions.csv"), encoding='utf-8-sig')
print(f"患者维度数据已保存到: {patient_dir}/patient_celltype_proportions.csv")

# 6. 统计检验函数
def perform_statistical_tests(data_long, celltype, group_col='group_final', value_col='proportion'):
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

# 7. 绘图函数
def plot_celltype_proportions(data_long, celltype, output_dir, level="sample", p_value_info=None, group_col='group_final'):
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
    plt.xlabel('Group (group_final)', fontsize=11)  # 使用英文标签
    plt.ylabel('Cell Proportion (%)', fontsize=11)  # 使用英文标签
    
    # 调整y轴范围
    y_max = celltype_data['proportion'].max()
    plt.ylim(0, y_max * 1.2)
    
    # 添加样本/患者数量标签
    for i, group in enumerate(groups):
        group_data = celltype_data[celltype_data[group_col] == group]
        n = len(group_data)
        median_val = group_data['proportion'].median()
        
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

# 8. 执行样本维度的分析和绘图
print("\n" + "="*50)
print("生成样本维度箱线图 (基于group_final分组)")
print("="*50)

sample_stats_results = []
celltypes = df['celltype'].unique()

for i, celltype in enumerate(celltypes, 1):
    print(f"处理细胞类型 {i}/{len(celltypes)}: {celltype}")
    
    # 执行统计检验
    test_method, p_value, is_sig, details = perform_statistical_tests(
        sample_long, celltype, 'group_final', 'proportion'
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
    
    plot_celltype_proportions(sample_long, celltype, sample_dir, "sample", p_value_info, 'group_final')

# 保存样本维度统计结果
if sample_stats_results:
    sample_stats_df = pd.DataFrame(sample_stats_results)
    sample_stats_df.to_csv(os.path.join(sample_dir, "sample_statistical_tests.csv"), 
                          index=False, encoding='utf-8-sig')
    print(f"样本维度统计检验结果已保存到: {sample_dir}/sample_statistical_tests.csv")

# 9. 执行患者维度的分析和绘图
print("\n" + "="*50)
print("生成患者维度箱线图 (基于group_final分组)")
print("="*50)

patient_stats_results = []

for i, celltype in enumerate(celltypes, 1):
    print(f"处理细胞类型 {i}/{len(celltypes)}: {celltype}")
    
    # 执行统计检验
    test_method, p_value, is_sig, details = perform_statistical_tests(
        patient_long, celltype, 'group_final', 'proportion'
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
    
    plot_celltype_proportions(patient_long, celltype, patient_dir, "patient", p_value_info, 'group_final')

# 保存患者维度统计结果
if patient_stats_results:
    patient_stats_df = pd.DataFrame(patient_stats_results)
    patient_stats_df.to_csv(os.path.join(patient_dir, "patient_statistical_tests.csv"), 
                           index=False, encoding='utf-8-sig')
    print(f"患者维度统计检验结果已保存到: {patient_dir}/patient_statistical_tests.csv")

# 10. 生成汇总报告
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
    f.write(f"样本数: {sample_proportions.shape[0]}\n")
    f.write(f"患者数: {patient_proportions.shape[0]}\n")
    f.write(f"group_final分组: {sorted(df['group_final'].unique().tolist())}\n")
    f.write(f"patient_new数量: {len(df['patient_new'].unique())}\n\n")
    
    f.write("2. 细胞类型总体分布\n")
    f.write("-"*40 + "\n")
    for _, row in celltype_stats.iterrows():
        f.write(f"{row['celltype']:25s}: {row['count']:6d} 个细胞 ({row['percentage']:6.2f}%)\n")
    
    f.write("\n3. 样本维度统计 (基于group_final)\n")
    f.write("-"*40 + "\n")
    if sample_stats_results:
        sig_count = sum(1 for r in sample_stats_results if r['significant'])
        f.write(f"显著差异的细胞类型数: {sig_count}/{len(sample_stats_results)}\n")
        for result in sample_stats_results:
            if result['significant']:
                f.write(f"  {result['celltype']}: p={result['p_value']:.3e} ({result['test_method']})\n")
    
    f.write("\n4. 患者维度统计 (基于group_final和patient_new)\n")
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
    
    f.write("\n6. 样本映射信息\n")
    f.write("-"*40 + "\n")
    f.write("group_final映射:\n")
    for sample, group in sample_to_group_final.items():
        f.write(f"  {sample}: {group}\n")
    
    f.write("\npatient_new映射:\n")
    for sample, patient in sample_to_patient_new.items():
        f.write(f"  {sample}: {patient}\n")

print(f"汇总报告已保存到: celltype_analysis_summary.txt")
print("\n分析完成!")
