import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, chi2_contingency
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def combine_m2_celltypes(adata):
    """
    合并M2和EMT_M2为M2
    """
    print("正在合并M2和EMT_M2细胞类型...")
    
    # 检查annotation列是否存在
    if 'annotation' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'annotation'列，请检查数据")
    
    # 获取原始的细胞类型分布
    original_counts = adata.obs['annotation'].value_counts()
    print("原始细胞类型分布:")
    print(original_counts)
    
    # 创建新的细胞类型列
    adata.obs['annotation_combined'] = adata.obs['annotation'].copy()
    
    # 将EMT_M2改为M2
    m2_mask = (adata.obs['annotation'] == 'M2') | (adata.obs['annotation'] == 'EMT_M2')
    adata.obs.loc[m2_mask, 'annotation_combined'] = 'M2'
    
    # 获取合并后的细胞类型分布
    combined_counts = adata.obs['annotation_combined'].value_counts()
    print("\n合并后的细胞类型分布:")
    print(combined_counts)
    
    # 计算M2合并前后数量的变化
    if 'M2' in original_counts.index and 'EMT_M2' in original_counts.index:
        original_m2 = original_counts.get('M2', 0)
        original_emt_m2 = original_counts.get('EMT_M2', 0)
        combined_m2 = combined_counts.get('M2', 0)
        print(f"\nM2合并统计:")
        print(f"  原始M2: {original_m2}")
        print(f"  原始EMT_M2: {original_emt_m2}")
        print(f"  合并后M2: {combined_m2}")
        print(f"  合计: {original_m2 + original_emt_m2}")
    
    return adata

def calculate_celltype_proportions(adata, celltype_of_interest):
    """
    计算每个patient中指定细胞类型的比例
    """
    print("正在计算细胞类型比例...")
    
    # 检查细胞类型列是否存在
    if 'annotation_combined' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'annotation_combined'列，请先运行combine_m2_celltypes函数")
    
    # 检查MPCpercent列是否存在
    if 'MPCpercent' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'MPCpercent'列，请检查数据")
    
    # 检查patient列是否存在
    if 'patient' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'patient'列，请检查数据")
    
    # 检查指定的细胞类型是否存在
    if celltype_of_interest not in adata.obs['annotation_combined'].unique():
        print(f"警告: 细胞类型 '{celltype_of_interest}' 在数据中不存在")
        print(f"数据中存在的细胞类型: {list(adata.obs['annotation_combined'].unique())}")
        return None, None
    
    print(f"将分析的细胞类型: {celltype_of_interest}")
    
    # 计算每个patient的总细胞数
    patient_totals = adata.obs['patient'].value_counts().sort_index()
    
    # 为每个patient计算M2细胞类型的数量
    m2_counts = adata.obs[adata.obs['annotation_combined'] == celltype_of_interest].groupby('patient').size()
    
    # 确保所有patient都在m2_counts中（即使计数为0）
    m2_counts = m2_counts.reindex(patient_totals.index, fill_value=0)
    
    # 计算M2比例
    proportions = pd.DataFrame({
        'patient': patient_totals.index,
        'total_cells': patient_totals.values,
        f'{celltype_of_interest}_count': m2_counts.values,
        f'{celltype_of_interest}_proportion': m2_counts.values / patient_totals.values
    })
    proportions.set_index('patient', inplace=True)
    
    # 添加分组信息
    patient_groups = adata.obs[['patient', 'MPCpercent']].drop_duplicates().set_index('patient')
    proportions['MPCpercent'] = patient_groups['MPCpercent']
    
    # 打印一些统计信息
    print(f"\n细胞类型比例计算完成:")
    print(f"  分析的病人数量: {len(proportions)}")
    print(f"  M2细胞总数: {m2_counts.sum()}")
    print(f"  M2平均比例: {proportions[f'{celltype_of_interest}_proportion'].mean():.4f}")
    print(f"  M2比例范围: {proportions[f'{celltype_of_interest}_proportion'].min():.4f} - {proportions[f'{celltype_of_interest}_proportion'].max():.4f}")
    
    return proportions, celltype_of_interest

def statistical_analysis(proportions, celltype):
    """
    进行统计学差异分析
    """
    print("\n正在进行统计学差异分析...")
    
    # 创建结果字典
    results = {}
    
    # 获取分组
    groups = proportions['MPCpercent'].unique()
    
    if len(groups) < 2:
        print("警告: 分组数量少于2，无法进行差异分析")
        return None
    
    print(f"发现的分组: {list(groups)}")
    
    # 准备数据
    group_data = []
    for group in groups:
        group_values = proportions[proportions['MPCpercent'] == group][f'{celltype}_proportion'].values
        group_data.append(group_values)
        print(f"  分组 {group}: {len(group_values)} 个样本, 平均M2比例: {np.mean(group_values):.4f}")
    
    # 检查是否有足够的数据进行检验
    valid_groups = [data for data in group_data if len(data) > 0 and not np.all(np.isnan(data))]
    
    if len(valid_groups) < 2:
        results[celltype] = {'test': '无法检验', 'p_value': np.nan, 'significant': False}
        return results
    
    # 根据分组数量选择检验方法
    if len(groups) == 2:
        # Mann-Whitney U检验（非参数检验）
        try:
            stat, p_value = mannwhitneyu(
                proportions[proportions['MPCpercent'] == groups[0]][f'{celltype}_proportion'].dropna(),
                proportions[proportions['MPCpercent'] == groups[1]][f'{celltype}_proportion'].dropna(),
                alternative='two-sided'
            )
            test_used = 'Mann-Whitney U'
            print(f"  使用{test_used}检验")
        except Exception as e:
            print(f"  Mann-Whitney U检验失败: {e}")
            p_value = np.nan
            test_used = '检验失败'
    else:
        # Kruskal-Wallis检验（多组比较）
        try:
            group_lists = []
            for group in groups:
                group_data = proportions[proportions['MPCpercent'] == group][f'{celltype}_proportion'].dropna()
                if len(group_data) > 0:
                    group_lists.append(group_data)
            
            if len(group_lists) >= 2:
                stat, p_value = kruskal(*group_lists)
                test_used = 'Kruskal-Wallis'
                print(f"  使用{test_used}检验")
            else:
                p_value = np.nan
                test_used = '无法检验'
        except Exception as e:
            print(f"  Kruskal-Wallis检验失败: {e}")
            p_value = np.nan
            test_used = '检验失败'
    
    significant = p_value < 0.05 if not np.isnan(p_value) else False
    results[celltype] = {
        'test': test_used,
        'p_value': p_value,
        'significant': significant
    }
    
    print(f"  p值: {p_value:.4f}")
    print(f"  是否显著(p<0.05): {significant}")
    
    return results

def create_visualizations(proportions, celltype, results, output_prefix):
    """
    创建可视化图表
    """
    print("\n正在生成可视化图表...")
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_style("whitegrid")
    
    # 1. M2比例箱线图
    plt.figure(figsize=(10, 8))
    
    # 准备数据
    plot_data = []
    for group in proportions['MPCpercent'].unique():
        group_data = proportions[proportions['MPCpercent'] == group][f'{celltype}_proportion']
        for value in group_data:
            plot_data.append({'MPCpercent': group, 'Proportion': value})
    
    plot_df = pd.DataFrame(plot_data)
    
    # 绘制箱线图
    ax = sns.boxplot(data=plot_df, x='MPCpercent', y='Proportion', palette='Set2')
    
    # 添加散点显示个体数据点
    sns.stripplot(data=plot_df, x='MPCpercent', y='Proportion', 
                  color='black', alpha=0.5, size=6, jitter=True)
    
    # 设置标题和标签
    plt.title(f'{celltype}细胞比例在不同MPCpercent分组中的分布\np-value: {results[celltype].get("p_value", np.nan):.4f}', 
              fontsize=16, pad=20)
    plt.xlabel('MPCpercent分组', fontsize=14)
    plt.ylabel(f'{celltype}细胞比例', fontsize=14)
    
    # 标记显著性
    if results[celltype].get('significant', False):
        # 在图表顶部添加显著星号
        y_max = plot_df['Proportion'].max()
        plt.text(0.5, y_max + 0.01, '*', fontsize=24, ha='center', color='red')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(f'{output_prefix}_{celltype}_proportion_boxplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_{celltype}_proportion_boxplot.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. M2比例小提琴图
    plt.figure(figsize=(10, 8))
    
    # 绘制小提琴图
    ax = sns.violinplot(data=plot_df, x='MPCpercent', y='Proportion', palette='Set2', inner='quartile')
    
    # 添加散点
    sns.stripplot(data=plot_df, x='MPCpercent', y='Proportion', 
                  color='black', alpha=0.5, size=6, jitter=True)
    
    # 设置标题和标签
    plt.title(f'{celltype}细胞比例分布（小提琴图）', fontsize=16, pad=20)
    plt.xlabel('MPCpercent分组', fontsize=14)
    plt.ylabel(f'{celltype}细胞比例', fontsize=14)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(f'{output_prefix}_{celltype}_proportion_violin.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_{celltype}_proportion_violin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 分组平均M2比例条形图
    plt.figure(figsize=(8, 6))
    
    # 计算每组平均比例
    mean_proportions = proportions.groupby('MPCpercent')[f'{celltype}_proportion'].mean()
    std_proportions = proportions.groupby('MPCpercent')[f'{celltype}_proportion'].std()
    
    # 绘制条形图
    ax = mean_proportions.plot(kind='bar', yerr=std_proportions, 
                                capsize=5, color=sns.color_palette('Set2'), 
                                edgecolor='black', linewidth=1.5)
    
    # 设置标题和标签
    plt.title(f'不同MPCpercent分组的平均{celltype}细胞比例', fontsize=16, pad=20)
    plt.xlabel('MPCpercent分组', fontsize=14)
    plt.ylabel(f'平均{celltype}细胞比例', fontsize=14)
    plt.xticks(rotation=0)
    
    # 在每个柱子上方添加数值
    for i, v in enumerate(mean_proportions):
        plt.text(i, v + 0.001, f'{v:.3f}', ha='center', fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(f'{output_prefix}_{celltype}_mean_proportion_bar.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_{celltype}_mean_proportion_bar.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  可视化图表已保存为PNG和PDF格式")

def main():
    # 配置参数
    h5ad_file = "adata_filtered_with_MPCpercent.h5ad"
    output_prefix = "MPCpercent_M2_analysis"
    
    # 要分析的细胞类型（现在是合并后的M2）
    celltype_of_interest = "M2"
    
    try:
        # 1. 读取数据
        print("正在读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        print(f"样本数量: {adata.obs['patient'].nunique()}")
        print(f"MPCpercent分组分布:\n{adata.obs['MPCpercent'].value_counts()}")
        
        # 2. 合并M2和EMT_M2为M2
        adata = combine_m2_celltypes(adata)
        
        # 3. 计算细胞类型比例
        proportions, celltype = calculate_celltype_proportions(adata, celltype_of_interest)
        
        if proportions is None:
            print("无法计算细胞类型比例，程序退出")
            return
        
        # 4. 统计学分析
        results = statistical_analysis(proportions, celltype)
        
        if results:
            # 5. 保存结果
            # 创建结果数据框
            results_df = pd.DataFrame.from_dict(results, orient='index')
            results_df['celltype'] = results_df.index
            
            # 描述性统计
            desc_stats = proportions.groupby('MPCpercent')[f'{celltype}_proportion'].describe()
            
            # 保存到Excel
            with pd.ExcelWriter(f'{output_prefix}_results.xlsx') as writer:
                proportions.to_excel(writer, sheet_name='原始数据')
                desc_stats.to_excel(writer, sheet_name='描述性统计')
                results_df.to_excel(writer, sheet_name='差异分析结果')
            
            # 6. 生成可视化
            create_visualizations(proportions, celltype, results, output_prefix)
            
            # 7. 打印摘要结果
            print("\n" + "="*60)
            print("分析结果摘要")
            print("="*60)
            
            p_value = results[celltype].get('p_value', np.nan)
            significant = results[celltype].get('significant', False)
            
            print(f"分析的细胞类型: {celltype}")
            print(f"统计检验方法: {results[celltype].get('test', 'N/A')}")
            print(f"p值: {p_value:.6f}")
            print(f"是否显著(p<0.05): {significant}")
            
            if significant:
                print(f"结论: 不同MPCpercent分组间的{celltype}细胞比例存在显著差异")
            else:
                print(f"结论: 不同MPCpercent分组间的{celltype}细胞比例无显著差异")
            
            # 打印每个分组的平均M2比例
            print(f"\n各分组{celltype}细胞平均比例:")
            for group in proportions['MPCpercent'].unique():
                mean_prop = proportions[proportions['MPCpercent'] == group][f'{celltype}_proportion'].mean()
                n_samples = len(proportions[proportions['MPCpercent'] == group])
                print(f"  {group}: {mean_prop:.4f} (n={n_samples})")
            
            print(f"\n结果已保存到:")
            print(f"  - {output_prefix}_results.xlsx")
            print(f"  - {output_prefix}_{celltype}_proportion_boxplot.png/pdf")
            print(f"  - {output_prefix}_{celltype}_proportion_violin.png/pdf")
            print(f"  - {output_prefix}_{celltype}_mean_proportion_bar.png/pdf")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
