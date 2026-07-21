import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, chi2_contingency
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def calculate_celltype_proportions(adata, celltypes_of_interest):
    """
    计算每个patient中指定细胞类型的比例
    """
    print("正在计算细胞类型比例...")
    
    # 检查细胞类型列是否存在
    if 'annotation' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'cell_type'列，请检查数据")
    
    # 检查MPCpercent列是否存在
    if 'MPCpercent' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'MPCpercent'列，请检查数据")
    
    # 检查patient列是否存在
    if 'patient' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'patient'列，请检查数据")
    
    # 获取数据中实际存在的细胞类型
    actual_celltypes = adata.obs['annotation'].unique()
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
    patient_totals = adata.obs['patient'].value_counts().sort_index()
    
    # 为每个patient计算每种细胞类型的数量
    celltype_counts = pd.crosstab(adata.obs['patient'], adata.obs['annotation'])
    
    # 只保留我们感兴趣的细胞类型
    celltype_counts = celltype_counts[valid_celltypes]
    
    # 计算比例
    proportions = celltype_counts.div(patient_totals, axis=0)
    
    # 添加分组信息
    patient_groups = adata.obs[['patient', 'MPCpercent']].drop_duplicates().set_index('patient')
    proportions['MPCpercent'] = patient_groups['MPCpercent']
    
    return proportions, valid_celltypes

def statistical_analysis(proportions, celltypes):
    """
    进行统计学差异分析
    """
    print("\n正在进行统计学差异分析...")
    
    results = {}
    groups = proportions['MPCpercent'].unique()
    
    if len(groups) < 2:
        print("警告: 分组数量少于2，无法进行差异分析")
        return None
    
    print(f"发现的分组: {list(groups)}")
    
    for celltype in celltypes:
        group_data = []
        for group in groups:
            group_values = proportions[proportions['MPCpercent'] == group][celltype].values
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
                stat, p_value = mannwhitneyu(
                    proportions[proportions['MPCpercent'] == groups[0]][celltype].dropna(),
                    proportions[proportions['MPCpercent'] == groups[1]][celltype].dropna(),
                    alternative='two-sided'
                )
                test_used = 'Mann-Whitney U'
            except:
                p_value = np.nan
                test_used = '检验失败'
        else:
            # Kruskal-Wallis检验（多组比较）
            try:
                group_lists = []
                for group in groups:
                    group_data = proportions[proportions['MPCpercent'] == group][celltype].dropna()
                    if len(group_data) > 0:
                        group_lists.append(group_data)
                
                if len(group_lists) >= 2:
                    stat, p_value = kruskal(*group_lists)
                    test_used = 'Kruskal-Wallis'
                else:
                    p_value = np.nan
                    test_used = '无法检验'
            except:
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

def create_visualizations(proportions, celltypes, results, output_prefix):
    """
    创建可视化图表
    """
    print("\n正在生成可视化图表...")
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # 1. 细胞类型比例箱线图
    fig, axes = plt.subplots(4, 4, figsize=(20, 16))
    axes = axes.ravel()
    
    for i, celltype in enumerate(celltypes):
        if i >= len(axes):
            break
            
        # 准备数据
        plot_data = []
        for group in proportions['MPCpercent'].unique():
            group_data = proportions[proportions['MPCpercent'] == group][celltype]
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
    mean_proportions = proportions.groupby('MPCpercent')[celltypes].mean()
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(mean_proportions.T, annot=True, cmap='YlOrRd', 
                cbar_kws={'label': 'Proportion'}, fmt='.3f')
    plt.title('MPCpercent group')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. 显著性星号图
    significant_celltypes = [ct for ct in celltypes 
                           if results[ct].get('significant', False)]
    
    if significant_celltypes:
        plt.figure(figsize=(10, 6))
        significant_means = proportions.groupby('MPCpercent')[significant_celltypes].mean()
        
        significant_means.T.plot(kind='bar', figsize=(12, 6))
        plt.title('Sig')
        plt.ylabel('proportion')
        plt.xticks(rotation=45)
        plt.legend(title='MPCpercent分组')
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_significant_celltypes.png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    # 配置参数
    h5ad_file = "adata_filtered_with_MPCpercent.h5ad"
    output_prefix = "MPCpercent_immu"
    
    # 感兴趣的细胞类型列表
    celltypes_of_interest = [
        "CD8+_Tem", "M2", "MKI67+_T", "B_follicular", "CD4+_Tex",
        "LYVE1+_TAM", "CD8+_Tcm", "NK", "IFN-T", "B_germinal_center",
        "EMT_M2", "Treg", "CD8+_Tex", "CD4+_Tem", "B_memory", "TAM","Plasma","CD4+_Tcm","MKI67+_Plasma"
    ]
    
    try:
        # 1. 读取数据
        print("正在读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        print(f"MPCpercent分组分布:\n{adata.obs['MPCpercent'].value_counts()}")
        
        # 2. 计算细胞类型比例
        proportions, valid_celltypes = calculate_celltype_proportions(adata, celltypes_of_interest)
        
        # 3. 统计学分析
        results = statistical_analysis(proportions, valid_celltypes)
        
        if results:
            # 4. 多重检验校正
            results = multiple_testing_correction(results)
            
            # 5. 保存结果
            # 描述性统计
            desc_stats = proportions.groupby('MPCpercent')[valid_celltypes].describe()
            
            # 结果汇总
            results_df = pd.DataFrame.from_dict(results, orient='index')
            results_df['celltype'] = results_df.index
            
            # 保存到Excel
            with pd.ExcelWriter(f'{output_prefix}_results.xlsx') as writer:
                proportions.to_excel(writer, sheet_name='原始比例数据')
                desc_stats.to_excel(writer, sheet_name='描述性统计')
                results_df.to_excel(writer, sheet_name='差异分析结果')
            
            # 6. 生成可视化
            create_visualizations(proportions, valid_celltypes, results, output_prefix)
            
            # 7. 打印摘要结果
            print("\n" + "="*60)
            print("分析结果摘要")
            print("="*60)
            
            significant_count = sum(1 for r in results.values() if r.get('significant', False))
            print(f"分析的细胞类型数量: {len(valid_celltypes)}")
            print(f"具有显著差异的细胞类型数量: {significant_count}")
            
            if significant_count > 0:
                print("\n具有显著差异的细胞类型:")
                for celltype, result in results.items():
                    if result.get('significant', False):
                        p_val = result.get('p_value', np.nan)
                        corr_p = result.get('corrected_p_value', 'N/A')
                        print(f"  {celltype}: p = {p_val:.4f}, 校正后p = {corr_p}")
            
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
