import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal
from scipy.stats import f_oneway
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

def load_and_merge_data(main_file, relate_file):
    """
    加载主数据文件和分级文件，并合并数据
    """
    print("正在读取主数据文件...")
    # 读取主Excel文件[1,2](@ref)
    df = pd.read_csv(main_file,sep="\t")
    print(f"主数据形状: {df.shape}")
    print(f"主数据列名: {list(df.columns)}")
    
    print("\n正在读取分级文件...")
    # 读取分级文件（假设为CSV格式）[1](@ref)
    relate_df = pd.read_csv(relate_file,sep="\t")
    print(f"分级数据形状: {relate_df.shape}")
    print(f"分级文件列名: {list(relate_df.columns)}")
    
    # 检查分级文件列名
    if len(relate_df.columns) < 2:
        raise ValueError("分级文件应包含至少两列：patient和分级")
    
    # 假设第一列是patient，第二列是分级
    patient_col = relate_df.columns[0]
    grading_col = relate_df.columns[1]
    
    print(f"使用列: {patient_col} 作为患者ID, {grading_col} 作为分级")
    
    # 合并数据，只保留主数据中在分级文件中存在的patient[2](@ref)
    merged_df = df.merge(relate_df[[patient_col, grading_col]], 
                        left_on='patient', right_on=patient_col, how='inner')
    
    # 重命名分级列为tie
    merged_df = merged_df.rename(columns={grading_col: 'tie'})
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"分级分布:\n{merged_df['tie'].value_counts()}")
    
    return merged_df

def check_variable_type(df, column_name):
    """
    检查变量类型：数值型还是分类型
    """
    if df[column_name].dtype in ['int64', 'float64']:
        # 检查唯一值的数量来判断是否是分类变量
        unique_count = df[column_name].nunique()
        total_count = len(df[column_name])
        
        # 如果唯一值数量少于总样本的20%或少于20个，认为是分类变量
        if unique_count < 20 or unique_count / total_count < 0.2:
            return 'categorical'
        else:
            return 'numerical'
    else:
        return 'categorical'

def statistical_analysis(df, group_var, target_var):
    """
    进行统计学差异分析
    根据变量类型自动选择合适的检验方法[9,10,11](@ref)
    """
    print(f"\n正在进行统计学分析: {group_var} 对 {target_var}")
    
    # 检查变量类型
    target_type = check_variable_type(df, target_var)
    group_categories = df[group_var].nunique()
    
    print(f"{target_var} 变量类型: {target_type}")
    print(f"{group_var} 分组数量: {group_categories}")
    
    results = {}
    
    if target_type == 'numerical':
        # 数值型变量：使用方差分析或Kruskal-Wallis检验[9](@ref)
        groups = df[group_var].unique()
        group_data = [df[df[group_var] == group][target_var].dropna() for group in groups]
        
        # 检查正态性和方差齐性（简化处理）
        if len(groups) == 2:
            # 两组比较：t检验或Mann-Whitney U检验
            from scipy.stats import ttest_ind, mannwhitneyu
            
            # 简单检查方差齐性
            var_ratio = np.var(group_data[0]) / np.var(group_data[1])
            if 0.5 < var_ratio < 2:  # 方差相对齐性
                stat, p_value = ttest_ind(group_data[0], group_data[1])
                test_used = '独立样本t检验'
            else:
                stat, p_value = mannwhitneyu(group_data[0], group_data[1])
                test_used = 'Mann-Whitney U检验'
        else:
            # 多组比较：ANOVA或Kruskal-Wallis
            try:
                stat, p_value = f_oneway(*group_data)
                test_used = '单因素方差分析(ANOVA)'
            except:
                stat, p_value = kruskal(*group_data)
                test_used = 'Kruskal-Wallis H检验'
        
        results = {
            'test': test_used,
            'statistic': stat,
            'p_value': p_value,
            'significant': p_value < 0.05 if not np.isnan(p_value) else False
        }
        
    else:
        # 分类型变量：使用卡方检验[9,11](@ref)
        # 创建列联表
        contingency_table = pd.crosstab(df[group_var], df[target_var])
        
        # 确保列联表有效
        if contingency_table.size > 0 and contingency_table.sum().sum() > 0:
            try:
                chi2, p_value, dof, expected = chi2_contingency(contingency_table)
                results = {
                    'test': '卡方检验',
                    'chi2': chi2,
                    'p_value': p_value,
                    'dof': dof,
                    'significant': p_value < 0.05 if not np.isnan(p_value) else False
                }
            except Exception as e:
                print(f"卡方检验出错: {e}")
                results = {'test': '检验失败', 'p_value': np.nan, 'significant': False}
        else:
            results = {'test': '数据不足', 'p_value': np.nan, 'significant': False}
    
    return results

def create_visualizations(df, group_var, target_var, results, output_prefix):
    """
    创建可视化图表[1,12](@ref)
    """
    print("\n正在生成可视化图表...")
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # 1. 根据变量类型创建不同的可视化
    target_type = check_variable_type(df, target_var)
    
    if target_type == 'numerical':
        # 数值型变量：箱线图和小提琴图
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # 箱线图
        sns.boxplot(data=df, x=group_var, y=target_var, ax=ax1)
        ax1.set_title(f'{target_var} {group_var} proportion')
        ax1.tick_params(axis='x', rotation=45)
        
        # 小提琴图
        sns.violinplot(data=df, x=group_var, y=target_var, ax=ax2)
        ax2.set_title(f'{target_var} {group_var} proportion')
        ax2.tick_params(axis='x', rotation=45)
        
        # 添加统计结果
        p_value = results.get('p_value', np.nan)
        test_used = results.get('test', '')
        fig.suptitle(f'{test_used} p value = {p_value:.4f} {"(sig)" if p_value < 0.05 else "(non-sig)"}', 
                    fontsize=14, fontweight='bold')
        
    else:
        # 分类型变量：堆叠条形图和百分比堆叠条形图
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # 创建列联表
        contingency_table = pd.crosstab(df[group_var], df[target_var])
        
        # 堆叠条形图（计数）
        contingency_table.plot(kind='bar', stacked=True, ax=ax1)
        ax1.set_title(f'{target_var} {group_var} proportion')
        ax1.set_ylabel('cell number')
        ax1.legend(title=target_var, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.tick_params(axis='x', rotation=45)
        
        # 百分比堆叠条形图
        contingency_table_perc = contingency_table.div(contingency_table.sum(axis=1), axis=0) * 100
        contingency_table_perc.plot(kind='bar', stacked=True, ax=ax2)
        ax2.set_title(f'{target_var} {group_var} proportion')
        ax2.set_ylabel('(%)')
        ax2.set_ylim(0, 100)
        ax2.legend(title=target_var, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.tick_params(axis='x', rotation=45)
        
        # 添加统计结果
        p_value = results.get('p_value', np.nan)
        test_used = results.get('test', '')
        fig.suptitle(f'{test_used}: p value = {p_value:.4f} {"(sig)" if p_value < 0.05 else "(non-sig)"}', 
                    fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_{group_var}_vs_{target_var}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. 热图显示分布情况（特别是对于分类变量）
    if target_type == 'categorical':
        plt.figure(figsize=(12, 8))
        
        # 创建百分比列联表
        contingency_table = pd.crosstab(df[group_var], df[target_var])
        contingency_table_perc = contingency_table.div(contingency_table.sum(axis=1), axis=0) * 100
        
        sns.heatmap(contingency_table_perc.T, annot=True, cmap='YlOrRd', 
                   cbar_kws={'label': '百分比 (%)'}, fmt='.1f')
        plt.title(f'{target_var} {group_var} proportion heatmap')
        plt.xlabel(group_var)
        plt.ylabel(target_var)
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_{group_var}_vs_{target_var}_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"可视化图表已保存: {output_prefix}_*.png")

def main():
    """
    主函数
    """
    # 配置参数
    main_file = "CNs_10_all_cell_features.xls"
    relate_file = "releate.csv"  # 假设为CSV文件
    output_prefix = "tie_CNs10_analysis"
    
    try:
        # 1. 读取和合并数据
        df = load_and_merge_data(main_file, relate_file)
        
        # 检查必要的列是否存在
        if 'tie' not in df.columns:
            raise ValueError("合并后数据中未找到'tie'列")
        if 'CNs_10' not in df.columns:
            raise ValueError("数据中未找到'CNs_10'列")
        
        print(f"\n数据基本信息:")
        print(f"总样本数: {len(df)}")
        print(f"分级(tie)的唯一值: {df['tie'].unique()}")
        print(f"CNs_10的唯一值数量: {df['CNs_10'].nunique()}")
        
        # 2. 统计学分析
        results = statistical_analysis(df, 'tie', 'CNs_10')
        
        # 3. 生成可视化
        create_visualizations(df, 'tie', 'CNs_10', results, output_prefix)
        
        # 4. 保存详细结果
        # 描述性统计
        if check_variable_type(df, 'CNs_10') == 'numerical':
            desc_stats = df.groupby('tie')['CNs_10'].describe()
        else:
            # 对于分类变量，显示频数统计
            desc_stats = pd.crosstab(df['tie'], df['CNs_10'], margins=True)
        
        # 保存结果到Excel
        with pd.ExcelWriter(f'{output_prefix}_results.xlsx') as writer:
            #df.to_excel(writer, sheet_name='合并数据', index=False)
            #desc_stats.to_excel(writer, sheet_name='描述性统计')
            
            # 保存检验结果
            results_df = pd.DataFrame([results])
            results_df.to_excel(writer, sheet_name='统计检验结果', index=False)
        
        # 5. 打印摘要结果
        print("\n" + "="*60)
        print("分析结果摘要")
        print("="*60)
        
        print(f"分组变量: tie")
        print(f"目标变量: CNs_10")
        print(f"使用的统计方法: {results.get('test', 'N/A')}")
        
        p_value = results.get('p_value', np.nan)
        significant = results.get('significant', False)
        
        if not np.isnan(p_value):
            print(f"p值: {p_value:.6f}")
            print(f"统计显著性: {'显著 (p < 0.05)' if significant else '不显著 (p ≥ 0.05)'}")
            
            if significant:
                print("结论: 不同分级(tie)之间的CNs_10分布存在显著差异")
            else:
                print("结论: 未发现不同分级(tie)之间的CNs_10分布存在显著差异")
        else:
            print("警告: 无法计算有效的p值")
        
        # 显示详细的分布情况
        print(f"\n分布情况:")
        if check_variable_type(df, 'CNs_10') == 'categorical':
            contingency_table = pd.crosstab(df['tie'], df['CNs_10'])
            print("列联表（计数）:")
            print(contingency_table)
            
            print("\n列联表（百分比）:")
            contingency_table_perc = contingency_table.div(contingency_table.sum(axis=1), axis=0) * 100
            print(contingency_table_perc.round(2))
        else:
            group_stats = df.groupby('tie')['CNs_10'].agg(['count', 'mean', 'std', 'min', 'max'])
            print(group_stats)
        
        print(f"\n分析完成！结果文件:")
        print(f"  - {output_prefix}_results.xlsx (详细数据)")
        print(f"  - {output_prefix}_tie_vs_CNs_10.png (主要可视化)")
        print(f"  - {output_prefix}_tie_vs_CNs_10_heatmap.png (热图)")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
