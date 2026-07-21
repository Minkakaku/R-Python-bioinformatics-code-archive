import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal, chi2_contingency
from statsmodels.stats.multitest import multipletests
import warnings
import os
warnings.filterwarnings('ignore')

def calculate_celltype_proportions(adata):
    """
    计算每个sample中所有细胞类型的比例
    """
    print("正在计算细胞类型比例...")
    
    # 检查细胞类型列是否存在
    if 'annotation' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'annotation'列，请检查数据")
    
    # 检查group列是否存在
    if 'mucus_zone' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'mucus_zone'列，请检查数据")
    
    # 检查sample列是否存在
    if 'sample' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'sample'列，请检查数据")
    
    # 获取数据中所有的细胞类型
    all_celltypes = adata.obs['annotation'].unique()
    print(f"数据中存在的细胞类型数量: {len(all_celltypes)}")
    print(f"细胞类型列表: {list(all_celltypes)}")
    
    # 计算每个sample的总细胞数
    sample_totals = adata.obs['sample'].value_counts().sort_index()
    
    # 为每个sample计算每种细胞类型的数量
    celltype_counts = pd.crosstab(adata.obs['sample'], adata.obs['annotation'])
    
    # 确保列名是字符串
    celltype_counts.columns = [str(col) for col in celltype_counts.columns]
    
    # 计算比例
    proportions = celltype_counts.div(sample_totals, axis=0)
    
    # 添加分组信息
    sample_groups = adata.obs[['sample', 'mucus_zone']].drop_duplicates().set_index('sample')
    proportions['mucus_zone'] = sample_groups['mucus_zone']
    
    return proportions, all_celltypes

def statistical_analysis(proportions, celltypes):
    """
    进行统计学差异分析
    """
    print("\n正在进行统计学差异分析...")
    
    results = {}
    groups = proportions['mucus_zone'].unique()
    
    if len(groups) < 2:
        print("警告: 分组数量少于2，无法进行差异分析")
        return None
    
    print(f"发现的分组: {list(groups)}")
    
    for celltype in celltypes:
        celltype_str = str(celltype)  # 确保细胞类型名称是字符串
        
        group_data = []
        for group in groups:
            group_values = proportions[proportions['mucus_zone'] == group][celltype_str].values
            group_data.append(group_values)
        
        # 检查是否有足够的数据进行检验
        valid_groups = [data for data in group_data if len(data) > 0 and not np.all(np.isnan(data))]
        
        if len(valid_groups) < 2:
            results[celltype_str] = {'test': '无法检验', 'p_value': np.nan, 'significant': False}
            continue
        
        # 根据分组数量选择检验方法
        if len(groups) == 2:
            # Mann-Whitney U检验（非参数检验）
            try:
                stat, p_value = mannwhitneyu(
                    proportions[proportions['mucus_zone'] == groups[0]][celltype_str].dropna(),
                    proportions[proportions['mucus_zone'] == groups[1]][celltype_str].dropna(),
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
                    group_data = proportions[proportions['mucus_zone'] == group][celltype_str].dropna()
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
        results[celltype_str] = {
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

def create_individual_boxplots(proportions, celltypes, results, output_dir):
    """
    为每个细胞类型创建单独的箱线图
    """
    print("正在为每个细胞类型生成单独的箱线图...")
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    for celltype in celltypes:
        celltype_str = str(celltype)  # 确保细胞类型名称是字符串
        
        # 检查这个细胞类型是否在proportions的列中
        if celltype_str not in proportions.columns:
            print(f"警告: 细胞类型 '{celltype_str}' 不在proportions的列中，跳过")
            continue
        
        # 准备数据
        plot_data = []
        for group in proportions['mucus_zone'].unique():
            group_data = proportions[proportions['mucus_zone'] == group][celltype_str]
            for value in group_data:
                plot_data.append({'Group': group, 'Proportion': value, 'CellType': celltype_str})
        
        plot_df = pd.DataFrame(plot_data)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # 绘制箱线图
        sns.boxplot(data=plot_df, x='Group', y='Proportion', ax=ax, palette="Set2")
        
        # 添加散点显示个体数据点
        sns.stripplot(data=plot_df, x='Group', y='Proportion', ax=ax, 
                     color='black', alpha=0.6, size=4, jitter=True)
        
        # 设置标题和标签
        p_value = results[celltype_str].get('p_value', np.nan) if celltype_str in results else np.nan
        significance = "***" if (celltype_str in results and results[celltype_str].get('significant', False)) else "ns"
        
        ax.set_title(f'{celltype_str}\n(p-value: {p_value:.4f} {significance})', fontsize=14, fontweight='bold')
        ax.set_xlabel('group Group', fontsize=12)
        ax.set_ylabel('Cell Proportion', fontsize=12)
        ax.tick_params(axis='x', rotation=45)
        
        # 设置y轴从0开始
        ax.set_ylim(bottom=0)
        
        # 如果显著，添加背景色
        if celltype_str in results and results[celltype_str].get('significant', False):
            ax.patch.set_facecolor('#FFF0F0')
            ax.patch.set_alpha(0.3)
        
        # 添加统计信息注释
        stats_text = ""
        if celltype_str in results:
            stats_text = f"Test: {results[celltype_str].get('test', 'N/A')}\n"
            if 'corrected_p_value' in results[celltype_str]:
                stats_text += f"FDR p-value: {results[celltype_str]['corrected_p_value']:.4f}"
        
        if stats_text:
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                   fontsize=10)
        
        plt.tight_layout()
        
        # 保存为PNG和PDF
        safe_celltype_name = str(celltype_str).replace('/', '_').replace(' ', '_')
        plt.savefig(f'{output_dir}/{safe_celltype_name}_boxplot.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_dir}/{safe_celltype_name}_boxplot.pdf', bbox_inches='tight')
        plt.close()
        
    print(f"已为 {len(celltypes)} 个细胞类型生成单独的箱线图")

def create_summary_visualizations(proportions, celltypes, results, output_prefix):
    """
    创建汇总可视化图表
    """
    print("\n正在生成汇总可视化图表...")
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # 确保celltypes是字符串列表
    celltypes_str = [str(ct) for ct in celltypes]
    
    # 1. 所有细胞类型比例热图
    # 只选择在proportions中存在的细胞类型
    available_celltypes = [ct for ct in celltypes_str if ct in proportions.columns]
    
    if available_celltypes:
        mean_proportions = proportions.groupby('mucus_zone')[available_celltypes].mean()
        
        plt.figure(figsize=(max(12, len(available_celltypes)*0.3), max(8, len(proportions['mucus_zone'].unique())*0.5)))
        sns.heatmap(mean_proportions.T, annot=True, cmap='YlOrRd', 
                    cbar_kws={'label': 'Proportion'}, fmt='.3f')
        plt.title('Average Cell Type Proportions by group Group')
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_heatmap.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_prefix}_heatmap.pdf', bbox_inches='tight')
        plt.close()
    else:
        print("警告: 没有可用的细胞类型用于生成热图")
    
    # 2. 显著性细胞类型条形图
    significant_celltypes = [ct for ct in celltypes_str 
                           if ct in results and results[ct].get('significant', False)]
    
    if significant_celltypes:
        # 确保这些细胞类型在proportions中存在
        valid_sig_celltypes = [ct for ct in significant_celltypes if ct in proportions.columns]
        
        if valid_sig_celltypes:
            significant_means = proportions.groupby('mucus_zone')[valid_sig_celltypes].mean()
            
            plt.figure(figsize=(14, 8))
            significant_means.T.plot(kind='bar', figsize=(14, 8))
            plt.title('Significantly Different Cell Types by group Group')
            plt.ylabel('Average Proportion')
            plt.xlabel('Cell Type')
            plt.xticks(rotation=45)
            plt.legend(title='group Group')
            plt.tight_layout()
            plt.savefig(f'{output_prefix}_significant_celltypes.png', dpi=300, bbox_inches='tight')
            plt.savefig(f'{output_prefix}_significant_celltypes.pdf', bbox_inches='tight')
            plt.close()
    else:
        print("没有显著差异的细胞类型，跳过生成显著性条形图")
    
    # 3. p值分布图
    p_values = [results[ct].get('p_value', np.nan) for ct in celltypes_str if ct in results]
    valid_p_values = [p for p in p_values if not np.isnan(p)]
    
    if valid_p_values:
        plt.figure(figsize=(10, 6))
        plt.hist(valid_p_values, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        plt.axvline(x=0.05, color='red', linestyle='--', label='p=0.05')
        plt.xlabel('p-value')
        plt.ylabel('Frequency')
        plt.title('Distribution of p-values across all cell types')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_pvalue_distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_prefix}_pvalue_distribution.pdf', bbox_inches='tight')
        plt.close()

def main():
    # 配置参数
    h5ad_file = "extracted_samples.h5ad"
    output_prefix = "mz_sample"
    individual_plots_dir = "plot_all_mz_sample"
    
    try:
        # 1. 读取数据
        print("正在读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        print(f"观察到的变量: {adata.obs.columns.tolist()}")
        print(f"group分组分布:\n{adata.obs['mucus_zone'].value_counts()}")
        
        # 2. 计算所有细胞类型的比例
        proportions, all_celltypes = calculate_celltype_proportions(adata)
        print(f"成功计算了 {len(all_celltypes)} 个细胞类型的比例")
        
        # 3. 统计学分析
        results = statistical_analysis(proportions, all_celltypes)
        
        if results:
            # 4. 多重检验校正
            results = multiple_testing_correction(results)
            
            # 5. 为每个细胞类型生成单独的箱线图
            create_individual_boxplots(proportions, all_celltypes, results, individual_plots_dir)
            
            # 6. 生成汇总可视化
            create_summary_visualizations(proportions, all_celltypes, results, output_prefix)
            
            # 7. 保存结果到Excel
            # 描述性统计
            # 确保细胞类型列是字符串
            str_celltypes = [str(ct) for ct in all_celltypes if str(ct) in proportions.columns]
            desc_stats = proportions.groupby('mucus_zone')[str_celltypes].describe()
            
            # 结果汇总
            results_df = pd.DataFrame.from_dict(results, orient='index')
            results_df['annotation'] = results_df.index
            
            # 保存到Excel
            with pd.ExcelWriter(f'{output_prefix}_results.xlsx') as writer:
                proportions.to_excel(writer, sheet_name='原始比例数据')
                desc_stats.to_excel(writer, sheet_name='描述性统计')
                results_df.to_excel(writer, sheet_name='差异分析结果')
            
            # 8. 打印摘要结果
            print("\n" + "="*60)
            print("分析结果摘要")
            print("="*60)
            
            significant_count = sum(1 for r in results.values() if r.get('significant', False))
            fdr_significant_count = sum(1 for r in results.values() if r.get('fdr_significant', False))
            
            print(f"分析的细胞类型数量: {len(all_celltypes)}")
            print(f"具有显著差异的细胞类型数量 (p < 0.05): {significant_count}")
            print(f"FDR校正后显著的细胞类型数量 (FDR < 0.05): {fdr_significant_count}")
            
            if significant_count > 0:
                print("\n具有显著差异的细胞类型 (p < 0.05):")
                for celltype, result in results.items():
                    if result.get('significant', False):
                        p_val = result.get('p_value', np.nan)
                        corr_p = result.get('corrected_p_value', 'N/A')
                        fdr_sig = "是" if result.get('fdr_significant', False) else "否"
                        print(f"  {celltype}: p = {p_val:.4f}, 校正后p = {corr_p}, FDR显著 = {fdr_sig}")
            
            print(f"\n结果已保存到:")
            print(f"  - {output_prefix}_results.xlsx")
            print(f"  - {individual_plots_dir}/ (包含{len(all_celltypes)}个单独的箱线图)")
            print(f"  - {output_prefix}_heatmap.png/pdf")
            print(f"  - {output_prefix}_pvalue_distribution.png/pdf")
            if significant_count > 0:
                print(f"  - {output_prefix}_significant_celltypes.png/pdf")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
