import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# 设置样式和字体
plt.style.use('default')
sns.set_palette("Set2")

def load_and_process_data(main_file, relate_file):
    """
    加载并处理主数据和分级数据
    """
    print("Loading main data file...")
    df = pd.read_csv(main_file,sep='\t')
    print(f"Main data shape: {df.shape}")
    
    print("Loading relate file...")
    # 假设relate文件是CSV格式，包含patient和分级两列
    relate_df = pd.read_csv(relate_file,sep='\t')
    print(f"Relate data shape: {relate_df.shape}")
    print(f"Relate file columns: {list(relate_df.columns)}")
    
    # 检查列名
    if len(relate_df.columns) < 2:
        raise ValueError("Relate file should have at least two columns: patient and grading")
    
    # 使用第一列作为patient，第二列作为分级
    patient_col = relate_df.columns[0]
    grading_col = relate_df.columns[1]
    
    print(f"Using column '{patient_col}' as patient ID")
    print(f"Using column '{grading_col}' as grading")
    
    # 合并数据，只保留主数据中在分级文件中存在的patient
    merged_df = df.merge(relate_df[[patient_col, grading_col]], 
                        left_on='patient', right_on=patient_col, how='inner')
    
    # 重命名分级列为tie
    merged_df = merged_df.rename(columns={grading_col: 'tie'})
    
    print(f"Merged data shape: {merged_df.shape}")
    print(f"Tie distribution:\n{merged_df['tie'].value_counts()}")
    
    return merged_df

def calculate_cn_proportions(df):
    """
    计算每个patient中不同CNs_10的比例
    """
    print("Calculating CN proportions per patient...")
    
    # 检查必要的列
    if 'patient' not in df.columns:
        raise ValueError("'patient' column not found in data")
    if 'CNs_10' not in df.columns:
        raise ValueError("'CNs_10' column not found in data")
    if 'tie' not in df.columns:
        raise ValueError("'tie' column not found in data")
    
    # 计算每个patient的总细胞数
    patient_totals = df.groupby('patient').size().reset_index(name='total_cells')
    
    # 计算每个patient中每个CNs_10的细胞数量
    cn_counts = df.groupby(['patient', 'CNs_10']).size().reset_index(name='cn_count')
    
    # 合并总细胞数和CN计数
    cn_data = pd.merge(cn_counts, patient_totals, on='patient')
    
    # 计算每个CNs_10在patient中的比例
    cn_data['proportion'] = cn_data['cn_count'] / cn_data['total_cells'] * 100
    
    # 添加tie信息
    patient_tie = df[['patient', 'tie']].drop_duplicates()
    cn_data = pd.merge(cn_data, patient_tie, on='patient')
    
    return cn_data

def statistical_analysis(cn_data):
    """
    进行统计学差异分析
    """
    print("Performing statistical analysis...")
    
    results = {}
    unique_cns = cn_data['CNs_10'].unique()
    ties = cn_data['tie'].unique()
    
    print(f"Found {len(unique_cns)} unique CNs_10 categories")
    print(f"Found {len(ties)} tie groups: {list(ties)}")
    
    if len(ties) < 2:
        print("Warning: Need at least 2 tie groups for comparison")
        return results
    
    for cn in unique_cns:
        cn_subset = cn_data[cn_data['CNs_10'] == cn]
        
        if len(cn_subset) == 0:
            results[cn] = {'test': 'No data', 'p_value': np.nan, 'significant': False}
            continue
        
        # 准备分组数据
        group_data = []
        for tie in ties:
            tie_data = cn_subset[cn_subset['tie'] == tie]['proportion'].values
            if len(tie_data) > 0:
                group_data.append(tie_data)
        
        # 检查是否有足够的数据进行检验
        if len(group_data) < 2:
            results[cn] = {'test': 'Insufficient data', 'p_value': np.nan, 'significant': False}
            continue
        
        # 根据分组数量选择检验方法
        if len(ties) == 2:
            # Mann-Whitney U检验（两组比较）
            try:
                group1_data = cn_subset[cn_subset['tie'] == ties[0]]['proportion'].dropna()
                group2_data = cn_subset[cn_subset['tie'] == ties[1]]['proportion'].dropna()
                
                if len(group1_data) > 0 and len(group2_data) > 0:
                    stat, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
                    test_used = 'Mann-Whitney U'
                else:
                    p_value = np.nan
                    test_used = 'Insufficient data'
            except Exception as e:
                print(f"Error in Mann-Whitney U test for {cn}: {e}")
                p_value = np.nan
                test_used = 'Test failed'
        else:
            # Kruskal-Wallis检验（多组比较）
            try:
                group_lists = []
                for tie in ties:
                    tie_data = cn_subset[cn_subset['tie'] == tie]['proportion'].dropna()
                    if len(tie_data) > 0:
                        group_lists.append(tie_data)
                
                if len(group_lists) >= 2:
                    stat, p_value = kruskal(*group_lists)
                    test_used = 'Kruskal-Wallis'
                else:
                    p_value = np.nan
                    test_used = 'Insufficient data'
            except Exception as e:
                print(f"Error in Kruskal-Wallis test for {cn}: {e}")
                p_value = np.nan
                test_used = 'Test failed'
        
        significant = p_value < 0.05 if not np.isnan(p_value) else False
        results[cn] = {
            'test': test_used,
            'p_value': p_value,
            'significant': significant
        }
    
    return results

def create_visualizations(cn_data, results, output_prefix):
    """
    创建可视化图表
    """
    print("Creating visualizations...")
    
    unique_cns = cn_data['CNs_10'].unique()
    ties = cn_data['tie'].unique()
    
    # 1. 为每个CNs_10创建箱线图
    n_cns = len(unique_cns)
    n_cols = 4
    n_rows = (n_cns + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5*n_rows))
    if n_rows == 1:
        if n_cols == 1:
            axes = [axes]
        else:
            axes = axes
    else:
        axes = axes.ravel()
    
    for i, cn in enumerate(unique_cns):
        if i >= len(axes):
            break
            
        ax = axes[i] if n_rows == 1 and n_cols > 1 else axes[i]
        
        # 准备数据
        cn_subset = cn_data[cn_data['CNs_10'] == cn]
        
        # 绘制箱线图
        sns.boxplot(data=cn_subset, x='tie', y='proportion', ax=ax)
        
        # 添加统计信息
        p_value = results.get(cn, {}).get('p_value', np.nan)
        test_used = results.get(cn, {}).get('test', '')
        
        title = f'{cn}\n'
        if not np.isnan(p_value):
            title += f'p-value: {p_value:.4f}'
            if results[cn]['significant']:
                title += ' *'
                ax.patch.set_facecolor('#FFF0F0')  # 浅红色背景标记显著差异
        
        ax.set_title(title)
        ax.set_xlabel('Tie')
        ax.set_ylabel('Proportion (%)')
        ax.tick_params(axis='x', rotation=45)
    
    # 隐藏多余的子图
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_cn_proportions_by_tie.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. 创建热图显示平均比例
    mean_proportions = cn_data.groupby(['tie', 'CNs_10'])['proportion'].mean().unstack()
    
    plt.figure(figsize=(max(10, len(unique_cns)), max(6, len(ties))))
    sns.heatmap(mean_proportions, annot=True, cmap='YlOrRd', 
                cbar_kws={'label': 'Mean Proportion (%)'}, fmt='.2f')
    plt.title('Mean CNs_10 Proportions by Tie')
    plt.xlabel('CNs_10')
    plt.ylabel('Tie')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 创建显著差异的CNs_10图表
    significant_cns = [cn for cn, result in results.items() 
                      if result.get('significant', False)]
    
    if significant_cns:
        plt.figure(figsize=(12, 8))
        
        # 准备数据
        sig_data = cn_data[cn_data['CNs_10'].isin(significant_cns)]
        
        # 绘制点图
        sns.pointplot(data=sig_data, x='CNs_10', y='proportion', hue='tie',
                     capsize=0.1, errwidth=1.5, markers='o', linestyles='-')
        plt.title('Significant CNs_10 Proportions by Tie')
        plt.xlabel('CNs_10')
        plt.ylabel('Proportion (%)')
        plt.xticks(rotation=45)
        plt.legend(title='Tie')
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_significant_cns.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print("Visualizations created successfully")

def main():
    # 配置参数
    main_file = "CNs_10_all_cell_features.xls"
    relate_file = "releate.csv"  # 假设为CSV格式
    output_prefix = "tie_cn_analysis"
    
    try:
        # 1. 加载和处理数据
        df = load_and_process_data(main_file, relate_file)
        
        # 2. 计算CNs_10比例
        cn_data = calculate_cn_proportions(df)
        
        print(f"CN proportion data shape: {cn_data.shape}")
        print(f"Unique CNs_10: {cn_data['CNs_10'].unique()}")
        
        # 3. 统计学分析
        results = statistical_analysis(cn_data)
        
        # 4. 创建可视化
        create_visualizations(cn_data, results, output_prefix)
        
        # 5. 保存结果
        # 描述性统计
        desc_stats = cn_data.groupby(['tie', 'CNs_10'])['proportion'].describe()
        
        # 结果汇总
        results_df = pd.DataFrame.from_dict(results, orient='index')
        results_df['CNs_10'] = results_df.index
        
        # 保存到Excel
        with pd.ExcelWriter(f'{output_prefix}_results.xlsx') as writer:
            cn_data.to_excel(writer, sheet_name='CN Proportion Data', index=False)
            desc_stats.to_excel(writer, sheet_name='Descriptive Statistics')
            results_df.to_excel(writer, sheet_name='Statistical Results', index=False)
        
        # 6. 打印摘要结果
        print("\n" + "="*60)
        print("Analysis Summary")
        print("="*60)
        
        significant_count = sum(1 for r in results.values() if r.get('significant', False))
        total_cns = len(results)
        
        print(f"Total CNs_10 analyzed: {total_cns}")
        print(f"CNs_10 with significant differences: {significant_count}")
        
        if significant_count > 0:
            print("\nCNs_10 with significant differences:")
            for cn, result in results.items():
                if result.get('significant', False):
                    p_val = result.get('p_value', np.nan)
                    test_used = result.get('test', 'N/A')
                    print(f"  {cn}: {test_used}, p = {p_val:.4f}")
        
        print(f"\nResults saved to:")
        print(f"  - {output_prefix}_results.xlsx")
        print(f"  - {output_prefix}_cn_proportions_by_tie.png")
        print(f"  - {output_prefix}_heatmap.png")
        if significant_count > 0:
            print(f"  - {output_prefix}_significant_cns.png")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
