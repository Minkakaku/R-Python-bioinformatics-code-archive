import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu
from scipy.stats import f_oneway
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和图形样式
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def read_and_merge_data(cns_file, releate_file):
    """
    读取并合并数据文件
    """
    print("正在读取数据文件...")
    
    # 读取CNs数据文件
    try:
        cns_df = pd.read_csv(cns_file, sep='\t')
        print(f"CNs数据形状: {cns_df.shape}")
        print(f"CNs数据列名: {list(cns_df.columns)}")
    except Exception as e:
        print(f"读取CNs文件错误: {e}")
        return None
    
    # 检查必要的列
    if 'patient' not in cns_df.columns or 'CNs_10' not in cns_df.columns:
        print("错误: CNs数据中缺少'patient'或'CNs_10'列")
        return None
    
    # 读取releate文件（尝试多种分隔符）
    try:
        # 尝试制表符分隔
        try:
            releate_df = pd.read_csv(releate_file, sep='\t', header=None)
        except:
            # 尝试空格分隔
            try:
                releate_df = pd.read_csv(releate_file, sep='\s+', header=None)
            except:
                # 尝试逗号分隔
                releate_df = pd.read_csv(releate_file, sep=',', header=None)
        
        # 设置列名
        if releate_df.shape[1] >= 2:
            releate_df = releate_df.iloc[:, :2]  # 只取前两列
            releate_df.columns = ['patient', 'T_group']
        else:
            print("错误: releate文件列数不足")
            return None
            
        print(f"releate数据形状: {releate_df.shape}")
        print(f"releate数据预览:\n{releate_df.head()}")
    except Exception as e:
        print(f"读取releate文件错误: {e}")
        return None
    
    # 确保patient列类型一致
    cns_df['patient'] = cns_df['patient'].astype(str)
    releate_df['patient'] = releate_df['patient'].astype(str)
    
    # 合并数据
    print("\n正在合并数据...")
    merged_df = cns_df.merge(releate_df, on='patient', how='inner')
    
    # 重命名分组列
    merged_df.rename(columns={'T_group': 'MPCpercent'}, inplace=True)
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"MPCpercent分组分布:\n{merged_df['MPCpercent'].value_counts()}")
    
    if merged_df.shape[0] == 0:
        print("错误: 合并后无数据，请检查patient列是否匹配")
        return None
    
    return merged_df

def calculate_proportions(merged_df):
    """
    计算不同分组下各CNs_10的比例
    """
    print("\n正在计算CNs_10比例...")
    
    # 方法1: 直接计算每个分组中各CNs_10的比例
    contingency_table = pd.crosstab(merged_df['MPCpercent'], merged_df['CNs_10'])
    proportions = contingency_table.div(contingency_table.sum(axis=1), axis=0)
    
    # 方法2: 按patient计算比例，然后按分组求平均（更稳健）
    patient_proportions = []
    patient_groups = {}
    
    for patient in merged_df['patient'].unique():
        patient_data = merged_df[merged_df['patient'] == patient]
        group = patient_data['MPCpercent'].iloc[0]
        
        # 计算该patient中各CNs_10的比例
        cns_counts = patient_data['CNs_10'].value_counts()
        total_cells = len(patient_data)
        
        proportions_dict = {'patient': patient, 'MPCpercent': group, 'total_cells': total_cells}
        for cns_value in merged_df['CNs_10'].unique():
            count = cns_counts.get(cns_value, 0)
            proportions_dict[f'CNs_{cns_value}'] = count / total_cells if total_cells > 0 else 0
        
        patient_proportions.append(proportions_dict)
        patient_groups[patient] = group
    
    proportion_df = pd.DataFrame(patient_proportions)
    
    print("各分组CNs_10平均比例:")
    cns_columns = [col for col in proportion_df.columns if col.startswith('CNs_')]
    mean_proportions = proportion_df.groupby('MPCpercent')[cns_columns].mean()
    print(mean_proportions.round(4))
    
    return proportion_df, contingency_table

def statistical_analysis(proportion_df, contingency_table):
    """
    进行统计学差异分析
    """
    print("\n正在进行差异分析...")
    
    results = {}
    cns_columns = [col for col in proportion_df.columns if col.startswith('CNs_')]
    groups = proportion_df['MPCpercent'].unique()
    
    if len(groups) < 2:
        print("警告: 分组数量少于2，无法进行差异分析")
        return None
    
    print(f"发现的分组: {list(groups)}")
    
    # 1. 卡方检验 - 检验整体分布差异
    try:
        chi2, p_value_chi2, dof, expected = chi2_contingency(contingency_table)
        results['overall_chi2'] = {
            'test': 'Chi-square',
            'statistic': chi2,
            'p_value': p_value_chi2,
            'dof': dof,
            'significant': p_value_chi2 < 0.05
        }
        print(f"卡方检验 - 卡方值: {chi2:.4f}, P值: {p_value_chi2:.4f}, 自由度: {dof}")
    except Exception as e:
        print(f"卡方检验错误: {e}")
    
    # 2. 对每个CNs_10值进行组间比较
    for cns_col in cns_columns:
        group_data = []
        valid_groups = []
        
        for group in groups:
            group_values = proportion_df[proportion_df['MPCpercent'] == group][cns_col].dropna()
            if len(group_values) > 0:
                group_data.append(group_values)
                valid_groups.append(group)
        
        if len(group_data) < 2:
            results[cns_col] = {'test': '无法检验', 'p_value': np.nan, 'significant': False}
            continue
        
        # 选择检验方法
        if len(group_data) == 2:
            # Mann-Whitney U检验
            try:
                stat, p_value = mannwhitneyu(group_data[0], group_data[1], alternative='two-sided')
                test_used = 'Mann-Whitney U'
            except:
                p_value = np.nan
                test_used = '检验失败'
        else:
            # Kruskal-Wallis检验
            try:
                stat, p_value = kruskal(*group_data)
                test_used = 'Kruskal-Wallis'
            except:
                p_value = np.nan
                test_used = '检验失败'
        
        significant = p_value < 0.05 if not np.isnan(p_value) else False
        results[cns_col] = {
            'test': test_used,
            'p_value': p_value,
            'significant': significant
        }
    
    return results

def create_visualizations(proportion_df, contingency_table, results, output_prefix):
    """
    创建可视化图表
    """
    print("\n正在生成可视化图表...")
    
    cns_columns = [col for col in proportion_df.columns if col.startswith('CNs_')]
    cns_values = [col.replace('CNs_', '') for col in cns_columns]
    
    # 1. 热图 - 显示各分组CNs_10比例
    plt.figure(figsize=(12, 8))
    
    # 计算比例热图数据
    heatmap_data = contingency_table.div(contingency_table.sum(axis=1), axis=0)
    
    plt.subplot(2, 2, 1)
    sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='YlOrRd', 
                cbar_kws={'label': '比例'})
    plt.title('各分组CNs_10比例热图')
    plt.xlabel('CNs_10值')
    plt.ylabel('MPCpercent分组')
    
    # 2. 堆叠柱状图
    plt.subplot(2, 2, 2)
    heatmap_data.T.plot(kind='bar', stacked=True, ax=plt.gca())
    plt.title('CNs_10比例堆叠图')
    plt.xlabel('CNs_10值')
    plt.ylabel('比例')
    plt.legend(title='MPCpercent分组', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    
    # 3. 箱线图 - 显示显著差异的CNs_10
    significant_cns = [col for col in cns_columns 
                      if col in results and results[col].get('significant', False)]
    
    if significant_cns:
        # 选择前4个显著差异的进行展示
        display_cns = significant_cns[:4] if len(significant_cns) > 4 else significant_cns
        
        for i, cns_col in enumerate(display_cns, 3):
            if i > 4:  # 最多显示2个箱线图
                break
                
            plt.subplot(2, 2, i)
            plot_data = []
            for group in proportion_df['MPCpercent'].unique():
                group_data = proportion_df[proportion_df['MPCpercent'] == group]
                for value in group_data[cns_col]:
                    plot_data.append({'Group': group, 'Proportion': value, 'CNs': cns_col})
            
            plot_df = pd.DataFrame(plot_data)
            sns.boxplot(data=plot_df, x='Group', y='Proportion')
            p_val = results[cns_col].get('p_value', np.nan)
            plt.title(f'{cns_col}\n(p={p_val:.4f})')
            plt.xlabel('MPCpercent分组')
            plt.ylabel('比例')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_analysis_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. 显著性星号图
    if 'overall_chi2' in results:
        plt.figure(figsize=(10, 6))
        overall_p = results['overall_chi2']['p_value']
        
        # 创建差异对比图
        group_means = proportion_df.groupby('MPCpercent')[cns_columns].mean()
        
        # 选择差异最大的几个CNs值
        if len(cns_columns) > 5:
            # 计算组间差异
            variances = []
            for col in cns_columns:
                if col in results and not np.isnan(results[col]['p_value']):
                    variances.append(results[col]['p_value'])
                else:
                    variances.append(1.0)
            
            # 选择p值最小的5个
            top_indices = np.argsort(variances)[:5]
            display_cols = [cns_columns[i] for i in top_indices]
        else:
            display_cols = cns_columns
        
        # 绘制柱状图
        group_means[display_cols].plot(kind='bar', figsize=(12, 6))
        plt.title(f'CNs_10比例比较 (整体卡方检验p={overall_p:.4f})')
        plt.ylabel('平均比例')
        plt.xlabel('MPCpercent分组')
        plt.legend(title='CNs_10值', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_significance_chart.png', dpi=300, bbox_inches='tight')
        plt.show()

def main():
    """
    主函数
    """
    # 文件路径配置
    cns_file = "CNs_10_all_cell_features.xls"  # 替换为您的文件路径
    releate_file = "releate"                   # 替换为您的文件路径
    output_prefix = "CNs_10_MPCpercent_analysis"
    
    try:
        # 1. 读取和合并数据
        merged_df = read_and_merge_data(cns_file, releate_file)
        if merged_df is None:
            return
        
        # 2. 计算比例
        proportion_df, contingency_table = calculate_proportions(merged_df)
        
        # 3. 统计学分析
        results = statistical_analysis(proportion_df, contingency_table)
        
        if results is None:
            return
        
        # 4. 创建可视化
        create_visualizations(proportion_df, contingency_table, results, output_prefix)
        
        # 5. 保存结果（避免保存过大数据）
        print("\n正在保存结果...")
        
        # 只保存汇总数据，避免Excel大小限制
        with pd.ExcelWriter(f'{output_prefix}_summary.xlsx') as writer:
            # 保存比例汇总
            summary_proportions = proportion_df.groupby('MPCpercent').mean(numeric_only=True)
            summary_proportions.to_excel(writer, sheet_name='比例汇总')
            
            # 保存列联表
            contingency_table.to_excel(writer, sheet_name='细胞计数列联表')
            
            # 保存统计学结果
            stats_data = []
            for key, result in results.items():
                if key == 'overall_chi2':
                    stats_data.append({
                        'CNs_10': '整体分布',
                        '检验方法': result['test'],
                        '统计量': result['statistic'],
                        'P值': result['p_value'],
                        '自由度': result.get('dof', ''),
                        '显著性': result['significant']
                    })
                else:
                    stats_data.append({
                        'CNs_10': key.replace('CNs_', ''),
                        '检验方法': result['test'],
                        '统计量': '',
                        'P值': result['p_value'],
                        '自由度': '',
                        '显著性': result['significant']
                    })
            
            stats_df = pd.DataFrame(stats_data)
            stats_df.to_excel(writer, sheet_name='统计学结果', index=False)
        
        # 6. 打印摘要结果
        print("\n" + "="*60)
        print("分析结果摘要")
        print("="*60)
        
        if 'overall_chi2' in results:
            overall_p = results['overall_chi2']['p_value']
            print(f"整体卡方检验P值: {overall_p:.4f}")
            if overall_p < 0.05:
                print("✅ 不同MPCpercent分组间的CNs_10分布存在显著差异")
            else:
                print("❌ 不同MPCpercent分组间的CNs_10分布无显著差异")
        
        significant_count = sum(1 for key, result in results.items() 
                             if key != 'overall_chi2' and result.get('significant', False))
        
        print(f"\n分析的CNs_10值数量: {len([k for k in results.keys() if k != 'overall_chi2'])}")
        print(f"具有显著差异的CNs_10值数量: {significant_count}")
        
        if significant_count > 0:
            print("\n具有显著差异的CNs_10值:")
            for key, result in results.items():
                if key != 'overall_chi2' and result.get('significant', False):
                    cns_value = key.replace('CNs_', '')
                    p_val = result['p_value']
                    print(f"  {cns_value}: {result['test']}, p = {p_val:.4f}")
        
        print(f"\n结果文件已保存:")
        print(f"  - {output_prefix}_summary.xlsx (汇总数据)")
        print(f"  - {output_prefix}_analysis_results.png (分析图表)")
        print(f"  - {output_prefix}_significance_chart.png (显著性图表)")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
