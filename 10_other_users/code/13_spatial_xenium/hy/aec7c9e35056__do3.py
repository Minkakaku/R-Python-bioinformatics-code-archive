import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def load_and_merge_data(main_file, relate_file):
    """
    读取主数据文件和分级文件，并合并数据
    """
    print("正在读取主数据文件...")
    # 读取主Excel文件
    df = pd.read_csv(main_file,sep="\t")
    print(f"主数据形状: {df.shape}")
    
    print("正在读取分级文件...")
    # 读取分级文件（假设为CSV格式）
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
    
    # 合并数据
    merged_df = df.merge(relate_df[[patient_col, grading_col]], 
                        left_on='patient', right_on=patient_col, how='inner')
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"分级分布:\n{merged_df[grading_col].value_counts()}")
    
    return merged_df, grading_col

def calculate_celltype_proportions(df, grading_col, celltypes_of_interest):
    """
    计算每个patient中指定细胞类型的比例
    """
    print("正在计算细胞类型比例...")
    
    # 检查必要的列是否存在
    if 'patient' not in df.columns:
        raise ValueError("数据中未找到'patient'列，请检查数据")
    
    if grading_col not in df.columns:
        raise ValueError(f"数据中未找到'{grading_col}'列，请检查数据")
    
    # 使用celltype作为细胞类型列（根据用户说明）
    celltype_col = 'celltype'
    if celltype_col not in df.columns:
        # 如果celltype不存在，尝试使用annotation
        celltype_col = 'annotation'
        if celltype_col not in df.columns:
            raise ValueError("数据中未找到'celltype'或'annotation'列，请检查数据")
    
    print(f"使用 '{celltype_col}' 作为细胞类型列")
    
    # 获取数据中实际存在的细胞类型
    actual_celltypes = df[celltype_col].unique()
    print(f"数据中存在的细胞类型: {list(actual_celltypes)}")
    
    # 筛选出我们感兴趣的且实际存在的细胞类型
    valid_celltypes = [ct for ct in celltypes_of_interest if ct in actual_celltypes]
    missing_celltypes = [ct for ct in celltypes_of_interest if ct not in actual_celltypes]
    
    if missing_celltypes:
        print(f"警告: 以下细胞类型在数据中不存在: {missing_celltypes}")
    
    if not valid_celltypes:
        raise ValueError("没有找到任何有效的细胞类型进行分析")
    
    print(f"将分析的细胞类型: {valid_celltypes}")
    
    # 计算每个patient的总细胞数
    patient_totals = df['patient'].value_counts().sort_index()
    
    # 为每个patient计算每种细胞类型的数量
    celltype_counts = pd.crosstab(df['patient'], df[celltype_col])
    
    # 只保留我们感兴趣的细胞类型
    celltype_counts = celltype_counts.reindex(columns=valid_celltypes, fill_value=0)
    
    # 计算比例
    proportions = celltype_counts.div(patient_totals, axis=0)
    
    # 添加分级信息
    patient_grades = df[['patient', grading_col]].drop_duplicates().set_index('patient')
    proportions[grading_col] = patient_grades[grading_col]
    
    return proportions, valid_celltypes, grading_col

def statistical_analysis(proportions, celltypes, grading_col):
    """
    进行统计学差异分析
    """
    print("\n正在进行统计学差异分析...")
    
    results = {}
    groups = proportions[grading_col].unique()
    
    if len(groups) < 2:
        print("警告: 分组数量少于2，无法进行差异分析")
        return None
    
    print(f"发现的分组: {list(groups)}")
    
    for celltype in celltypes:
        group_data = []
        for group in groups:
            group_values = proportions[proportions[grading_col] == group][celltype].values
            group_data.append(group_values)
        
        # 检查是否有足够的数据进行检验
        valid_groups = [data for data in group_data if len(data) > 0 and not np.all(np.isnan(data))]
        
        if len(valid_groups) < 2:
            results[celltype] = {'test': '无法检验', 'p_value': np.nan, 'significant': False}
            continue
        
        # 根据分组数量选择检验方法
        if len(groups) == 2:
            # Mann-Whitney U检验（非参数检验）
            try:
                group1_data = proportions[proportions[grading_col] == groups[0]][celltype].dropna()
                group2_data = proportions[proportions[grading_col] == groups[1]][celltype].dropna()
                
                if len(group1_data) > 0 and len(group2_data) > 0:
                    stat, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
                    test_used = 'Mann-Whitney U'
                else:
                    p_value = np.nan
                    test_used = '数据不足'
            except Exception as e:
                print(f"进行Mann-Whitney U检验时出错: {e}")
                p_value = np.nan
                test_used = '检验失败'
        else:
            # Kruskal-Wallis检验（多组比较）
            try:
                group_lists = []
                for group in groups:
                    group_data = proportions[proportions[grading_col] == group][celltype].dropna()
                    if len(group_data) > 0:
                        group_lists.append(group_data)
                
                if len(group_lists) >= 2:
                    stat, p_value = kruskal(*group_lists)
                    test_used = 'Kruskal-Wallis'
                else:
                    p_value = np.nan
                    test_used = '数据不足'
            except Exception as e:
                print(f"进行Kruskal-Wallis检验时出错: {e}")
                p_value = np.nan
                test_used = '检验失败'
        
        significant = p_value < 0.05 if not np.isnan(p_value) else False
        results[celltype] = {
            'test': test_used,
            'p_value': p_value,
            'significant': significant
        }
    
    return results

def multiple_testing_correction(results):
    """
    进行多重检验校正
    """
    p_values = [result['p_value'] for result in results.values() if not np.isnan(result['p_value'])]
    celltypes_with_pvals = [ct for ct, result in results.items() if not np.isnan(result['p_value'])]
    
    if len(p_values) == 0:
        return results
    
    # Benjamini-Hochberg FDR校正
    rejected, corrected_pvals, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    
    # 更新结果
    for i, celltype in enumerate(celltypes_with_pvals):
        results[celltype]['corrected_p_value'] = corrected_pvals[i]
        results[celltype]['fdr_significant'] = rejected[i]
    
    return results

def create_visualizations(proportions, celltypes, results, grading_col, output_prefix):
    """
    创建可视化图表
    """
    print("\n正在生成可视化图表...")
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # 1. 细胞类型比例箱线图
    n_celltypes = len(celltypes)
    n_cols = 4
    n_rows = (n_celltypes + n_cols - 1) // n_cols  # 向上取整
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5*n_rows))
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.ravel()
    
    for i, celltype in enumerate(celltypes):
        if i >= len(axes):
            break
            
        # 准备数据
        plot_data = []
        for group in proportions[grading_col].unique():
            group_data = proportions[proportions[grading_col] == group][celltype]
            for value in group_data:
                plot_data.append({'Group': group, 'Proportion': value, 'CellType': celltype})
        
        plot_df = pd.DataFrame(plot_data)
        
        # 绘制箱线图
        sns.boxplot(data=plot_df, x='Group', y='Proportion', ax=axes[i])
        axes[i].set_title(f'{celltype}\np-value: {results[celltype].get("p_value", np.nan):.4f}')
        axes[i].tick_params(axis='x', rotation=45)
        
        # 标记显著性
        if results[celltype].get('significant', False):
            axes[i].patch.set_facecolor('#FFF0F0')  # 浅红色背景标记显著差异
    
    # 隐藏多余的子图
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_celltype_proportions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. 热图显示各细胞类型在不同分组中的平均比例
    mean_proportions = proportions.groupby(grading_col)[celltypes].mean()
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(mean_proportions.T, annot=True, cmap='YlOrRd', 
                cbar_kws={'label': 'Proportion'}, fmt='.3f')
    plt.title(f'Average Cell Type Proportions by {grading_col}')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 显著性星号图
    significant_celltypes = [ct for ct in celltypes 
                           if results[ct].get('significant', False)]
    
    if significant_celltypes:
        plt.figure(figsize=(12, 6))
        significant_means = proportions.groupby(grading_col)[significant_celltypes].mean()
        
        significant_means.T.plot(kind='bar', figsize=(12, 6))
        plt.title('Significant Cell Type Proportions by Grading')
        plt.ylabel('Proportion')
        plt.xlabel('Cell Type')
        plt.xticks(rotation=45)
        plt.legend(title=grading_col)
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_significant_celltypes.png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    # 配置参数
    main_file = "CNs_10_all_cell_features.xls"
    relate_file = "releate.csv"  # 假设为CSV文件
    output_prefix = "grading_analysis"
    
    # 感兴趣的细胞类型列表
    celltypes_of_interest = ["Stress-responsive_glandular_Epithelium",
        "Proliferative_Epithelium", 
        "Luminal_Epithelium",
        "FGFR2+_Basal_Epithelium",
        "Basal_Epithelium", 
        "EMT_Epithelium",
        "Renal_Epithelium"
    ]
    
    try:
        # 1. 读取和合并数据
        df, grading_col = load_and_merge_data(main_file, relate_file)
        
        # 2. 计算细胞类型比例
        proportions, valid_celltypes, grading_col = calculate_celltype_proportions(
            df, grading_col, celltypes_of_interest)
        
        print(f"比例数据形状: {proportions.shape}")
        print(f"有效细胞类型: {valid_celltypes}")
        
        # 3. 统计学分析
        results = statistical_analysis(proportions, valid_celltypes, grading_col)
        
        if results:
            # 4. 多重检验校正
            results = multiple_testing_correction(results)
            
            # 5. 保存结果
            # 描述性统计
            desc_stats = proportions.groupby(grading_col)[valid_celltypes].describe()
            
            # 结果汇总
            results_df = pd.DataFrame.from_dict(results, orient='index')
            results_df['celltype'] = results_df.index
            
            # 保存到Excel
            with pd.ExcelWriter(f'{output_prefix}_results.xlsx') as writer:
                proportions.to_excel(writer, sheet_name='原始比例数据')
                desc_stats.to_excel(writer, sheet_name='描述性统计')
                results_df.to_excel(writer, sheet_name='差异分析结果')
            
            # 6. 生成可视化
            create_visualizations(proportions, valid_celltypes, results, grading_col, output_prefix)
            
            # 7. 打印摘要结果
            print("\n" + "="*60)
            print("分析结果摘要")
            print("="*60)
            
            significant_count = sum(1 for r in results.values() if r.get('significant', False))
            fdr_significant_count = sum(1 for r in results.values() if r.get('fdr_significant', False))
            
            print(f"分析的细胞类型数量: {len(valid_celltypes)}")
            print(f"具有显著差异的细胞类型数量: {significant_count}")
            print(f"FDR校正后显著的细胞类型数量: {fdr_significant_count}")
            
            if significant_count > 0:
                print("\n具有显著差异的细胞类型:")
                for celltype, result in results.items():
                    if result.get('significant', False):
                        p_val = result.get('p_value', np.nan)
                        corr_p = result.get('corrected_p_value', 'N/A')
                        test_used = result.get('test', 'N/A')
                        print(f"  {celltype}: {test_used}, p = {p_val:.4f}, 校正后p = {corr_p}")
            
            print(f"\n结果已保存到:")
            print(f"  - {output_prefix}_results.xlsx")
            print(f"  - {output_prefix}_celltype_proportions.png")
            print(f"  - {output_prefix}_heatmap.png")
            if significant_count > 0:
                print(f"  - {output_prefix}_significant_celltypes.png")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
