import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置字体和字号 - 增加基础字号
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 12  # 增加全局字体大小
plt.rcParams['axes.titlesize'] = 16  # 标题字体大小
plt.rcParams['axes.labelsize'] = 14  # 坐标轴标签字体大小
plt.rcParams['xtick.labelsize'] = 12  # X轴刻度标签字体大小
plt.rcParams['ytick.labelsize'] = 12  # Y轴刻度标签字体大小
plt.rcParams['legend.fontsize'] = 12  # 图例字体大小

# 创建输出目录
output_dir = "cn15_com_T"
os.makedirs(output_dir, exist_ok=True)

# 1. 读取h5ad文件
print("读取h5ad文件...")
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15_T.h5ad")
print(f"数据形状: {adata.shape}")
print(f"obs列: {adata.obs.columns.tolist()}")

# 检查必需列
required_cols = ['sample_new', 'T', 'CNs_15']
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
print(f"T 唯一值: {df['T'].nunique()}")
print(f"CNs_15 唯一值: {df['CNs_15'].nunique()}")

# 检查T的取值
T_values = sorted(df['T'].dropna().unique().tolist())
print(f"T分组: {T_values}")

# 2. 统计所有细胞类型的数量和占比
print("\n" + "="*50)
print("细胞类型总体统计")
print("="*50)

total_cells = len(df)
CNs_15_stats = df['CNs_15'].value_counts().reset_index()
CNs_15_stats.columns = ['CNs_15', 'count']
CNs_15_stats['percentage'] = (CNs_15_stats['count'] / total_cells * 100).round(2)

print(f"总细胞数: {total_cells}")
print(f"细胞类型数: {len(CNs_15_stats)}")
print(f"样本数量: {len(df['sample_new'].unique())}")
print("\n各细胞类型统计:")
for i, row in CNs_15_stats.iterrows():
    print(f"  {i+1:2d}. {row['CNs_15']:20s}: {row['count']:6d} 个细胞 ({row['percentage']:6.2f}%)")

# 保存总体统计
CNs_15_stats.to_csv(os.path.join(output_dir, "CNs_15_overall_statistics.csv"), index=False, encoding='utf-8-sig')
print(f"\n总体统计已保存到: {output_dir}/CNs_15_overall_statistics.csv")

# 3. 以样本为维度的分析（使用sample_new和T）
print("\n" + "="*50)
print("以样本为维度的分析 (使用sample_new和T)")
print("="*50)

# 3.1 获取sample_new到T的映射
sample_T_mapping = df.drop_duplicates(['sample_new', 'T'])[['sample_new', 'T']]

# 3.2 计算每个样本中每个细胞类型的比例
sample_cell_counts = df.groupby(['sample_new', 'CNs_15']).size().unstack(fill_value=0)
sample_totals = sample_cell_counts.sum(axis=1)
sample_proportions = sample_cell_counts.div(sample_totals, axis=0) * 100  # 转换为百分比

# 3.3 将T信息添加到sample_proportions
sample_proportions = sample_proportions.reset_index()
# 合并T信息
sample_proportions = sample_proportions.merge(sample_T_mapping, on='sample_new', how='left')

# 3.4 重新组织数据用于绘图
# 获取细胞类型列
CNs_15_cols = [col for col in sample_proportions.columns if col not in ['sample_new', 'T']]

sample_long = sample_proportions.melt(
    id_vars=['sample_new', 'T'], 
    value_vars=CNs_15_cols,
    var_name='CNs_15', 
    value_name='proportion'
).dropna(subset=['T', 'proportion'])  # 移除T或proportion为NaN的行

# 移除比例为0的行
sample_long = sample_long[sample_long['proportion'] > 0]

# 3.5 保存样本维度的数据
sample_summary = sample_proportions.copy()
sample_summary['total_cells'] = sample_summary[CNs_15_cols].sum(axis=1)
sample_summary.to_csv(os.path.join(output_dir, "sample_CNs_15_proportions.csv"), 
                     index=False, encoding='utf-8-sig')
print(f"样本维度数据已保存到: {output_dir}/sample_CNs_15_proportions.csv")

# 4. 统计检验函数 - 简化版本，只处理两组
def perform_statistical_tests(data_long, CNs_15, group_col='T', value_col='proportion'):
    """
    执行统计检验 - 专门处理两组比较
    返回: 检验方法, p值, 是否显著, 两两比较结果
    """
    # 只考虑T2和T3组
    target_groups = ['T2', 'T3']
    group_data = {}
    
    for group in target_groups:
        subset = data_long[(data_long[group_col] == group) & (data_long['CNs_15'] == CNs_15)]
        group_data[group] = subset[value_col].values
    
    # 检查是否有足够的数据
    valid_groups = {g: d for g, d in group_data.items() if len(d) >= 2}
    if len(valid_groups) < 2:
        return None, None, False, None, "数据不足"
    
    # 提取有效组的数据
    data_list = [group_data['T2'], group_data['T3']]
    
    # 两组比较：Mann-Whitney U检验（非参数）
    stat, p_value = stats.mannwhitneyu(data_list[0], data_list[1], alternative='two-sided')
    test_method = "Mann-Whitney U"
    
    # 计算效应量
    n1, n2 = len(data_list[0]), len(data_list[1])
    effect_size = 2 * stat / (n1 * n2) - 1
    
    # 两两比较结果
    pairwise_results = [{
        'group1': 'T2',
        'group2': 'T3',
        'p_value': p_value
    }]
    
    return test_method, p_value, p_value < 0.05, pairwise_results, {
        'groups': ['T2', 'T3'],
        'statistic': stat,
        'effect_size': effect_size
    }

# 5. 绘图函数 - 简化版本，只处理T2和T3
def plot_CNs_15_proportions_with_pvalues(data_long, CNs_15, output_dir, level="sample", 
                                          p_value_info=None, pairwise_results=None, 
                                          group_col='T', group_order=['T2', 'T3']):
    """
    绘制单个细胞类型的箱线图，并添加两组间比较的P值
    """
    # 过滤特定细胞类型的数据
    CNs_15_data = data_long[data_long['CNs_15'] == CNs_15].copy()
    
    if len(CNs_15_data) == 0:
        print(f"  警告: 细胞类型 '{CNs_15}' 在{level}维度无数据")
        return None
    
    # 只保留存在的组
    existing_groups = [g for g in group_order if g in CNs_15_data[group_col].unique()]
    if len(existing_groups) < 2:
        print(f"  警告: 细胞类型 '{CNs_15}' 在{level}维度只有{len(existing_groups)}个组，跳过绘图")
        return None
    
    # 计算两两比较的数量
    pairwise_count = len(pairwise_results) if pairwise_results else 0
    
    # 根据两两比较的数量调整图形高度
    base_height = 9  # 增加基础高度
    fig_height = base_height
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(max(10, len(existing_groups)*3), fig_height))  # 增加宽度
    
    # 绘制箱线图和小提琴图
    sns.boxplot(x=group_col, y='proportion', data=CNs_15_data, 
                order=existing_groups, palette='Set2', width=0.6, 
                fliersize=4, ax=ax)  # 增加异常点大小
    sns.stripplot(x=group_col, y='proportion', data=CNs_15_data, 
                  order=existing_groups, color='black', size=5, alpha=0.7, 
                  jitter=True, ax=ax)  # 增加点的大小
    
    # 添加统计信息标题
    title_parts = [f"Cell Type: {CNs_15}"]  # 使用英文标题
    if p_value_info and p_value_info.get('significant', False):
        test_method = p_value_info.get('test_method', '')
        p_val = p_value_info.get('p_value', 1)
        title_parts.append(f"Mann-Whitney U test: p={p_val:.3e}")
    elif p_value_info:
        p_val = p_value_info.get('p_value', 1)
        title_parts.append(f"Mann-Whitney U test: p={p_val:.3e}")
    
    plt.title('\n'.join(title_parts), fontsize=18, fontweight='bold', pad=20)  # 增大标题字号
    plt.xlabel('T Group', fontsize=16)  # 增大X轴标签字号
    plt.ylabel('Cell Proportion (%)', fontsize=16)  # 增大Y轴标签字号
    
    # 设置刻度标签大小
    ax.tick_params(axis='x', labelsize=14)  # 增大X轴刻度标签大小
    ax.tick_params(axis='y', labelsize=14)  # 增大Y轴刻度标签大小
    
    # 获取数据的范围
    y_max = CNs_15_data['proportion'].max()
    y_min = CNs_15_data['proportion'].min()
    y_range = y_max - y_min
    
    # 动态计算Y轴范围
    # 基础比例：为连接线留出40%的空间（增加0.4）
    y_top = y_max + y_range * 0.4
    
    # 确保至少有一点空间
    if y_top <= y_max:
        y_top = y_max + y_range * 0.4
    
    plt.ylim(0, y_top)  # 从0开始，确保所有数据和标注都在图内
    
    # 添加样本数量标签
    for i, group in enumerate(existing_groups):
        group_data_subset = CNs_15_data[CNs_15_data[group_col] == group]
        n = len(group_data_subset)
        median_val = group_data_subset['proportion'].median()
        
        # 在x轴上方添加数量标签
        label_y = y_max + y_range * 0.1
        plt.text(i, label_y, f'n={n}\nmedian={median_val:.2f}%', 
                ha='center', va='bottom', fontsize=12, fontweight='bold')  # 增大标签字号
    
    # 添加两两组间比较的P值
    if pairwise_results and len(pairwise_results) > 0:
        # 计算图形中的位置
        y_start = y_max + y_range * 0.2
        
        # 为每对比较添加P值
        for i, comp in enumerate(pairwise_results):
            group1 = comp['group1']
            group2 = comp['group2']
            p_val = comp['p_value']
            
            # 确保组在现有组中
            if group1 in existing_groups and group2 in existing_groups:
                # 获取组在x轴上的位置
                x1 = existing_groups.index(group1)
                x2 = existing_groups.index(group2)
                
                # 计算P值星号表示
                if p_val < 0.0001:
                    p_text = '****'
                elif p_val < 0.001:
                    p_text = '***'
                elif p_val < 0.01:
                    p_text = '**'
                elif p_val < 0.05:
                    p_text = '*'
                else:
                    p_text = f'p={p_val:.3f}'
                
                # 绘制连接线
                y_pos = y_start
                
                # 绘制水平连接线
                plt.plot([x1, x2], [y_pos, y_pos], lw=2, color='black')  # 增加线宽
                
                # 绘制垂直线端
                vertical_line_length = 0.008 * y_range
                plt.plot([x1, x1], [y_pos-vertical_line_length, y_pos], lw=2, color='black')
                plt.plot([x2, x2], [y_pos-vertical_line_length, y_pos], lw=2, color='black')
                
                # 添加P值文本
                text_y_offset = 0.015 * y_range
                plt.text((x1 + x2) * 0.5, y_pos + text_y_offset, p_text, 
                        ha='center', va='bottom', fontsize=12, fontweight='bold')  # 增大P值文本字号
    
    # 添加边距，确保所有元素都在图内
    plt.tight_layout(pad=3.0)  # 增加padding
    
    # 保存图形
    safe_CNs_15 = CNs_15.replace('/', '_').replace('\\', '_').replace(':', '_')
    output_file = os.path.join(output_dir, f"{level}_{safe_CNs_15}_proportion.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

# 6. 执行样本维度的分析和绘图
print("\n" + "="*50)
print("生成样本维度箱线图 (基于T分组)")
print("="*50)

# 设置组顺序
group_order = ['T2', 'T3']
print(f"组顺序: {group_order}")

sample_stats_results = []
CNs_15s = df['CNs_15'].unique()

for i, CNs_15 in enumerate(CNs_15s, 1):
    print(f"处理细胞类型 {i}/{len(CNs_15s)}: {CNs_15}")
    
    # 执行统计检验
    test_method, p_value, is_sig, pairwise_results, details = perform_statistical_tests(
        sample_long, CNs_15, 'T', 'proportion'
    )
    
    # 保存统计结果
    if p_value is not None:
        sample_stats_results.append({
            'CNs_15': CNs_15,
            'test_method': test_method,
            'p_value': p_value,
            'significant': is_sig,
            'pairwise_results': pairwise_results,
            'details': str(details)
        })
    
    # 绘制图形
    p_value_info = {
        'test_method': test_method,
        'p_value': p_value,
        'significant': is_sig
    } if p_value is not None else None
    
    # 计算两两比较的数量
    pairwise_count = len(pairwise_results) if pairwise_results else 0
    print(f"  两两比较数量: {pairwise_count}")
    
    plot_CNs_15_proportions_with_pvalues(
        sample_long, CNs_15, output_dir, "sample", 
        p_value_info, pairwise_results, 'T', group_order
    )

# 保存样本维度统计结果
if sample_stats_results:
    # 创建详细的统计结果DataFrame
    detailed_stats = []
    
    for result in sample_stats_results:
        CNs_15 = result['CNs_15']
        test_method = result['test_method']
        p_value = result['p_value']
        is_sig = result['significant']
        pairwise_results = result.get('pairwise_results', [])
        
        # 添加总体检验结果
        detailed_stats.append({
            'CNs_15': CNs_15,
            'comparison': 'T2_vs_T3',
            'test_method': test_method,
            'p_value': p_value,
            'significant': is_sig
        })
    
    detailed_stats_df = pd.DataFrame(detailed_stats)
    detailed_stats_df.to_csv(os.path.join(output_dir, "sample_statistical_tests_detailed.csv"), 
                          index=False, encoding='utf-8-sig')
    
    # 也保存简化版的结果
    simple_stats_df = pd.DataFrame(sample_stats_results)
    simple_stats_df.to_csv(os.path.join(output_dir, "sample_statistical_tests.csv"), 
                          index=False, encoding='utf-8-sig')
    
    print(f"样本维度统计检验结果已保存到: {output_dir}/sample_statistical_tests.csv")
    print(f"详细两两比较结果已保存到: {output_dir}/sample_statistical_tests_detailed.csv")

# 7. 生成汇总报告
print("\n" + "="*50)
print("生成汇总报告")
print("="*50)

with open(os.path.join(output_dir, "CNs_15_analysis_summary.txt"), "w", encoding='utf-8') as f:
    f.write("细胞类型比例分析汇总报告 (基于T分组)\n")
    f.write("="*60 + "\n\n")
    
    f.write("1. 数据概览\n")
    f.write("-"*40 + "\n")
    f.write(f"总细胞数: {total_cells}\n")
    f.write(f"细胞类型数: {len(CNs_15_stats)}\n")
    f.write(f"样本数 (sample_new): {sample_proportions.shape[0]}\n")
    f.write(f"T分组: {T_values}\n")
    f.write(f"样本数量: {len(df['sample_new'].unique())}\n")
    f.write(f"组顺序: {group_order}\n\n")
    
    f.write("2. 细胞类型总体分布\n")
    f.write("-"*40 + "\n")
    for _, row in CNs_15_stats.iterrows():
        f.write(f"{row['CNs_15']:25s}: {row['count']:6d} 个细胞 ({row['percentage']:6.2f}%)\n")
    
    f.write("\n3. 样本维度统计 (基于T分组)\n")
    f.write("-"*40 + "\n")
    if sample_stats_results:
        sig_count = sum(1 for r in sample_stats_results if r['significant'])
        f.write(f"显著差异的细胞类型数 (Mann-Whitney U检验): {sig_count}/{len(sample_stats_results)}\n\n")
        
        f.write("显著差异的细胞类型详情:\n")
        for result in sample_stats_results:
            if result['significant']:
                f.write(f"  {result['CNs_15']}: p={result['p_value']:.3e} ({result['test_method']})\n")
                f.write("\n")
    
    f.write("\n4. 输出文件列表\n")
    f.write("-"*40 + "\n")
    f.write("总体统计:\n")
    f.write(f"  {output_dir}/CNs_15_overall_statistics.csv\n\n")
    f.write("样本维度数据:\n")
    f.write(f"  {output_dir}/sample_CNs_15_proportions.csv\n")
    f.write(f"  {output_dir}/sample_statistical_tests.csv\n")
    f.write(f"  {output_dir}/sample_statistical_tests_detailed.csv\n")
    f.write(f"  {output_dir}/sample_*.png (箱线图)\n")

print(f"汇总报告已保存到: {output_dir}/CNs_15_analysis_summary.txt")

# 8. 创建T分布统计
print("\n" + "="*50)
print("T分布统计")
print("="*50)

T_counts = df['T'].value_counts().reset_index()
T_counts.columns = ['T', 'count']
T_counts['percentage'] = (T_counts['count'] / total_cells * 100).round(2)

print(f"\nT分布:")
for _, row in T_counts.iterrows():
    print(f"  {row['T']}: {row['count']} 个细胞 ({row['percentage']:.2f}%)")

T_counts.to_csv(os.path.join(output_dir, "T_distribution.csv"), 
                   index=False, encoding='utf-8-sig')
print(f"T分布统计已保存到: {output_dir}/T_distribution.csv")

print("\n分析完成!")
