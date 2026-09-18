import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

def analyze_cns_proportions(file_path, short_patients, long_patients, 
                           short_group_name="short", long_group_name="long"):
    """
    分析CNs_10列在两组间的比例差异
    
    Parameters:
    -----------
    file_path : str
        数据文件路径
    short_patients : list
        short组patient列表
    long_patients : list
        long组patient列表
    short_group_name : str
        short组名称
    long_group_name : str
        long组名称
    """
    
    # 1. 读取数据
    print("正在读取数据...")
    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"数据形状: {df.shape}")
        print(f"数据列: {list(df.columns)}")
    except Exception as e:
        print(f"读取文件错误: {e}")
        return None
    
    # 检查必要的列是否存在
    if 'patient' not in df.columns:
        print("错误: 数据中未找到'patient'列")
        return None
    if 'CNs_10' not in df.columns:
        print("错误: 数据中未找到'CNs_10'列")
        return None
    
    # 2. 添加分组信息
    print("\n正在添加分组信息...")
    df['group'] = None
    
    # 为short组添加标签
    short_mask = df['patient'].isin(short_patients)
    df.loc[short_mask, 'group'] = short_group_name
    
    # 为long组添加标签
    long_mask = df['patient'].isin(long_patients)
    df.loc[long_mask, 'group'] = long_group_name
    
    # 过滤掉不在两组中的patient
    df_filtered = df[df['group'].notna()].copy()
    print(f"过滤后数据形状: {df_filtered.shape}")
    print(f"分组分布:\n{df_filtered['group'].value_counts()}")
    
    # 检查是否有patient未被匹配
    matched_short = df_filtered[df_filtered['group'] == short_group_name]['patient'].nunique()
    matched_long = df_filtered[df_filtered['group'] == long_group_name]['patient'].nunique()
    
    print(f"short组匹配patient数: {matched_short}/{len(short_patients)}")
    print(f"long组匹配patient数: {matched_long}/{len(long_patients)}")
    
    # 3. 计算每个patient的CNs_10比例
    print("\n正在计算CNs_10比例...")
    
    # 获取所有唯一的CNs_10值
    cns_values = df_filtered['CNs_10'].unique()
    print(f"CNs_10的唯一值: {sorted(cns_values)}")
    
    # 按patient和group分组，计算每个CNs_10值的比例
    patient_proportions = []
    
    for patient in df_filtered['patient'].unique():
        patient_data = df_filtered[df_filtered['patient'] == patient]
        group = patient_data['group'].iloc[0]
        
        # 计算该patient中每个CNs_10值的计数
        cns_counts = patient_data['CNs_10'].value_counts()
        total_cells = len(patient_data)
        
        # 计算比例
        proportions = {}
        for cns_value in cns_values:
            count = cns_counts.get(cns_value, 0)
            proportions[cns_value] = count / total_cells if total_cells > 0 else 0
        
        proportions['patient'] = patient
        proportions['group'] = group
        proportions['total_cells'] = total_cells
        
        patient_proportions.append(proportions)
    
    # 转换为DataFrame
    proportion_df = pd.DataFrame(patient_proportions)
    
    # 4. 统计学检验
    print("\n正在进行统计学检验...")
    results = {}
    
    for cns_value in cns_values:
        short_props = proportion_df[proportion_df['group'] == short_group_name][cns_value]
        long_props = proportion_df[proportion_df['group'] == long_group_name][cns_value]
        
        # Mann-Whitney U检验（比较比例分布）
        try:
            stat, p_value = mannwhitneyu(short_props, long_props, alternative='two-sided')
            results[cns_value] = {
                'short_mean': short_props.mean(),
                'long_mean': long_props.mean(),
                'statistic': stat,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
        except Exception as e:
            print(f"CNs_10值 {cns_value} 检验错误: {e}")
            results[cns_value] = {
                'short_mean': short_props.mean() if len(short_props) > 0 else 0,
                'long_mean': long_props.mean() if len(long_props) > 0 else 0,
                'statistic': np.nan,
                'p_value': np.nan,
                'significant': False
            }
    
    # 5. 可视化
    print("\n正在生成可视化图表...")
    create_visualizations(proportion_df, cns_values, results, 
                         short_group_name, long_group_name)
    
    # 6. 保存结果
    output_filename = f"CNs_10_analysis_{short_group_name}_vs_{long_group_name}.xlsx"
    with pd.ExcelWriter(output_filename) as writer:
        proportion_df.to_excel(writer, sheet_name='Patient_Proportions', index=False)
        
        # 创建结果摘要
        summary_data = []
        for cns_value, result in results.items():
            summary_data.append({
                'CNs_10': cns_value,
                f'{short_group_name}_mean': result['short_mean'],
                f'{long_group_name}_mean': result['long_mean'],
                'statistic': result['statistic'],
                'p_value': result['p_value'],
                'significant': result['significant']
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Statistical_Results', index=False)
    
    print(f"\n结果已保存到: {output_filename}")
    
    return proportion_df, results

def create_visualizations(proportion_df, cns_values, results, short_name, long_name):
    """
    创建可视化图表
    """
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # 准备数据用于绘图
    plot_data = []
    for cns_value in cns_values:
        short_mean = proportion_df[proportion_df['group'] == short_name][cns_value].mean()
        long_mean = proportion_df[proportion_df['group'] == long_name][cns_value].mean()
        
        plot_data.append({'CNs_10': cns_value, 'Group': short_name, 'Proportion': short_mean})
        plot_data.append({'CNs_10': cns_value, 'Group': long_name, 'Proportion': long_mean})
    
    plot_df = pd.DataFrame(plot_data)
    
    # 1. 柱状图 - 两组CNs_10比例比较
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 柱状图
    sns.barplot(data=plot_df, x='CNs_10', y='Proportion', hue='Group', ax=ax1)
    ax1.set_title(f'CNs_10 : {short_name} vs {long_name}')
    ax1.set_xlabel('CNs_10')
    ax1.set_ylabel('Average Proportion')
    ax1.tick_params(axis='x', rotation=45)
    
    # 添加显著性标记
    for i, cns_value in enumerate(cns_values):
        if results[cns_value]['significant']:
            p_val = results[cns_value]['p_value']
            star = '*' if p_val < 0.05 else 'ns'
            ax1.text(i, max(plot_df[plot_df['CNs_10'] == cns_value]['Proportion']) + 0.02, 
                    star, ha='center', fontsize=12, fontweight='bold')
    
    # 2. 箱线图 - 显示分布情况
    melt_df = proportion_df.melt(id_vars=['patient', 'group'], 
                                value_vars=cns_values,
                                var_name='CNs_10', value_name='Proportion')
    
    sns.boxplot(data=melt_df, x='CNs_10', y='Proportion', hue='group', ax=ax2)
    ax2.set_title(f'CNs_10 Proportion: {short_name} vs {long_name}')
    ax2.set_xlabel('CNs_10')
    ax2.set_ylabel('Proportion')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f'CNs_10_comparison_{short_name}_vs_{long_name}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. 热图 - 显示比例差异
    plt.figure(figsize=(10, 6))
    
    # 创建热图数据
    heatmap_data = []
    for cns_value in cns_values:
        short_prop = proportion_df[proportion_df['group'] == short_name][cns_value].mean()
        long_prop = proportion_df[proportion_df['group'] == long_name][cns_value].mean()
        diff = long_prop - short_prop
        heatmap_data.append({'CNs_10': cns_value, short_name: short_prop, 
                           long_name: long_prop, 'Difference': diff})
    
    heatmap_df = pd.DataFrame(heatmap_data)
    heatmap_df = heatmap_df.set_index('CNs_10')
    
    # 绘制热图
    sns.heatmap(heatmap_df[['Difference']], annot=True, cmap='RdBu_r', center=0,
                cbar_kws={'label': f'{long_name} - {short_name} 比例差异'})
    plt.title(f'CNs_10 Proportion Heatmap ({long_name} - {short_name})')
    plt.tight_layout()
    plt.savefig(f'CNs_10_heatmap_{short_name}_vs_{long_name}.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    # 配置参数
    file_path = "CNs_10_all_cell_features.xls"
    
    # 定义分组
    short_patients = ["Z1954693-9","Z2101519-1","Z1642635-14","Z1430967-9","Z1814675-2","Z1952876-9","Z2123354-1"]
    long_patients = ["Z2345459-7","Q2236500-1","Q2406275-2","Z1736838-15","Z2171766-3","Z1635989-12","Z2038176-15","Q2324341-1","Z1530765-11","F1706967-11","F1805049-2","Z1915172-7","ZF2207170-8","ZF2202916-8","ZF2208790-1","Z2262137-5","Z2271577-12","ZF2304103-8","Z2323818-8","Z2336500-6","Z2332034-5","Z2417954-8","Z2224197-10","Z1944219-14","Z2057500-35","Z2353385-10","Z1826868-6","Z1802235-10"]
    # 自定义分组名称（可以修改）
    short_group_name = "young"
    long_group_name = "old"
    
    # 执行分析
    result_df, stats_results = analyze_cns_proportions(
        file_path=file_path,
        short_patients=short_patients,
        long_patients=long_patients,
        short_group_name=short_group_name,
        long_group_name=long_group_name
    )
    
    if result_df is not None:
        # 打印摘要结果
        print("\n" + "="*60)
        print("分析结果摘要")
        print("="*60)
        
        significant_count = sum(1 for r in stats_results.values() if r['significant'])
        print(f"分析的CNs_10值数量: {len(stats_results)}")
        print(f"具有显著差异的CNs_10值数量: {significant_count}")
        
        if significant_count > 0:
            print("\n具有显著差异的CNs_10值:")
            for cns_value, result in stats_results.items():
                if result['significant']:
                    p_val = result['p_value']
                    print(f"  {cns_value}: p = {p_val:.4f}")
        
        print(f"\n可视化图表已保存为:")
        print(f"  - CNs_10_comparison_{short_group_name}_vs_{long_group_name}.png")
        print(f"  - CNs_10_heatmap_{short_group_name}_vs_{long_group_name}.png")

if __name__ == "__main__":
    main()
