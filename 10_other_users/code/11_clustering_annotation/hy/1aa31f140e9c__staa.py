import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats import f_oneway, kruskal
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图风格
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def load_and_prepare_data(file_path):
    """加载和准备数据"""
    print("=== 加载数据 ===")
    adata = sc.read_h5ad(file_path)
    print(f"数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列
    required_cols = ['celltype', 'group_new']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"错误: 缺少必要的列: {missing_cols}")
        print(f"可用的列: {list(adata.obs.columns)}")
        return None
    
    print(f"细胞类型数量: {adata.obs['celltype'].nunique()}")
    print(f"分组数量: {adata.obs['group_new'].nunique()}")
    print(f"分组详情: {list(adata.obs['group_new'].unique())}")
    
    return adata

def create_composition_table(adata):
    """创建细胞组成表"""
    print("\n=== 创建细胞组成表 ===")
    
    # 查找样本标识符列
    sample_col_candidates = ['sample', 'Sample', 'sample_id', 'SampleID', 'patient', 'Patient', 'orig.ident']
    sample_col = None
    
    for col in sample_col_candidates:
        if col in adata.obs.columns:
            sample_col = col
            break
    
    if sample_col is None:
        print("警告: 未找到样本标识符列，使用group_new作为样本标识符")
        print("注意: 这可能导致每个分组只有一个样本，统计检验不可靠")
        sample_col = 'group_new'
    
    print(f"使用样本标识符列: {sample_col}")
    
    # 创建样本-细胞类型计数表
    cross_table = pd.crosstab(adata.obs[sample_col], adata.obs['celltype'])
    
    # 添加分组信息
    sample_groups = adata.obs.groupby(sample_col)['group_new'].first()
    cross_table['group_new'] = sample_groups
    
    print(f"样本数量: {len(cross_table)}")
    print(f"分组分布: {cross_table['group_new'].value_counts().to_dict()}")
    
    return cross_table, sample_col

def manual_multiple_test_correction(p_values, method='fdr_bh'):
    """手动实现多重比较校正"""
    print(f"\n=== 执行多重比较校正 ({method}) ===")
    
    n = len(p_values)
    if n == 0:
        return np.array([]), np.array([])
    
    # 排序p值
    sorted_indices = np.argsort(p_values)
    sorted_p = np.array(p_values)[sorted_indices]
    
    if method == 'bonferroni':
        # Bonferroni校正
        corrected_p = np.minimum(sorted_p * n, 1.0)
    
    elif method == 'fdr_bh':
        # Benjamini-Hochberg FDR校正
        ranks = np.arange(1, n + 1)
        corrected_p = sorted_p * n / ranks
        # 确保校正后的p值单调不减
        for i in range(n-2, -1, -1):
            corrected_p[i] = min(corrected_p[i], corrected_p[i+1])
        corrected_p = np.minimum(corrected_p, 1.0)
    
    elif method == 'fdr_by':
        # Benjamini-Yekutieli FDR校正
        c = np.sum(1.0 / np.arange(1, n + 1))
        corrected_p = sorted_p * n * c / np.arange(1, n + 1)
        for i in range(n-2, -1, -1):
            corrected_p[i] = min(corrected_p[i], corrected_p[i+1])
        corrected_p = np.minimum(corrected_p, 1.0)
    
    else:
        raise ValueError(f"不支持的校正方法: {method}")
    
    # 恢复原始顺序
    final_corrected = np.zeros(n)
    final_corrected[sorted_indices] = corrected_p
    
    # 判断显著性
    significant = final_corrected < 0.05
    
    return final_corrected, significant

def perform_chi_square_test(cross_table):
    """执行卡方检验"""
    print("\n=== 卡方检验 - 整体分布差异 ===")
    
    # 提取细胞类型计数（排除group_new列）
    cell_types = [col for col in cross_table.columns if col != 'group_new']
    count_data = cross_table[cell_types]
    
    # 执行卡方检验
    chi2, p_value, dof, expected = chi2_contingency(count_data)
    
    print(f"卡方值: {chi2:.4f}")
    print(f"P值: {p_value:.4e}")
    print(f"自由度: {dof}")
    
    if p_value < 0.05:
        print("结论: 细胞类型分布在组间存在显著差异 (p < 0.05)")
    else:
        print("结论: 细胞类型分布在组间无显著差异 (p ≥ 0.05)")
    
    return chi2, p_value, dof

def perform_pairwise_comparisons(cross_table, sample_col):
    """执行成对比较分析"""
    print("\n=== 成对比较分析 ===")
    
    groups = cross_table['group_new'].unique()
    cell_types = [col for col in cross_table.columns if col != 'group_new' and col != sample_col]
    
    all_results = []
    
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            group1, group2 = groups[i], groups[j]
            
            print(f"\n比较: {group1} vs {group2}")
            
            # 提取两个组的计数数据
            group1_data = cross_table[cross_table['group_new'] == group1][cell_types].sum()
            group2_data = cross_table[cross_table['group_new'] == group2][cell_types].sum()
            
            group_results = []
            
            for cell_type in cell_types:
                # 创建2x2列联表
                table = np.array([
                    [group1_data[cell_type], group1_data.sum() - group1_data[cell_type]],
                    [group2_data[cell_type], group2_data.sum() - group2_data[cell_type]]
                ])
                
                # 执行Fisher精确检验
                try:
                    odds_ratio, p_value = fisher_exact(table)
                except:
                    odds_ratio, p_value = 1.0, 1.0
                
                # 计算比例
                prop1 = group1_data[cell_type] / group1_data.sum() if group1_data.sum() > 0 else 0
                prop2 = group2_data[cell_type] / group2_data.sum() if group2_data.sum() > 0 else 0
                
                group_results.append({
                    'celltype': cell_type,
                    'odds_ratio': odds_ratio,
                    'p_value': p_value,
                    'proportion1': prop1,
                    'proportion2': prop2,
                    'proportion_diff': prop1 - prop2
                })
                
                # 显示显著结果
                if p_value < 0.05:
                    significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*"
                    print(f"  {cell_type}: OR = {odds_ratio:.3f}, p = {p_value:.4f} {significance}")
            
            # 提取p值进行多重比较校正
            p_values = [result['p_value'] for result in group_results]
            corrected_p, significant = manual_multiple_test_correction(p_values, method='fdr_bh')
            
            # 添加校正后的p值和显著性标记
            for k, result in enumerate(group_results):
                result['p_value_corrected'] = corrected_p[k]
                result['significant'] = significant[k]
                result['group1'] = group1
                result['group2'] = group2
                all_results.append(result)
    
    results_df = pd.DataFrame(all_results)
    return results_df

def perform_anova_analysis(cross_table):
    """对每种细胞类型执行ANOVA分析"""
    print("\n=== 细胞类型特异性ANOVA分析 ===")
    
    # 计算每个样本中每种细胞类型的比例
    cell_types = [col for col in cross_table.columns if col != 'group_new']
    proportion_data = []
    
    for idx, row in cross_table.iterrows():
        total = row[cell_types].sum()
        for cell_type in cell_types:
            proportion = row[cell_type] / total if total > 0 else 0
            proportion_data.append({
                'sample': idx,
                'group': row['group_new'],
                'celltype': cell_type,
                'proportion': proportion
            })
    
    prop_df = pd.DataFrame(proportion_data)
    
    # 对每种细胞类型执行ANOVA
    anova_results = []
    for cell_type in cell_types:
        group_data = []
        groups = cross_table['group_new'].unique()
        
        for group in groups:
            proportions = prop_df[(prop_df['celltype'] == cell_type) & 
                                (prop_df['group'] == group)]['proportion'].values
            if len(proportions) > 0:
                group_data.append(proportions)
        
        if len(group_data) >= 2:  # 至少需要两个组
            # 尝试ANOVA，如果数据不满足假设则使用Kruskal-Wallis检验
            try:
                # 检查方差齐性
                variances = [np.var(g) for g in group_data if len(g) > 1]
                if len(variances) >= 2 and max(variances) / min(variances) < 10:  # 方差比阈值
                    f_stat, p_value = f_oneway(*group_data)
                    test_type = 'ANOVA'
                else:
                    # 使用非参数检验
                    f_stat, p_value = kruskal(*group_data)
                    test_type = 'Kruskal-Wallis'
                
                anova_results.append({
                    'celltype': cell_type,
                    'f_statistic': f_stat,
                    'p_value': p_value,
                    'test_type': test_type,
                    'significant': p_value < 0.05
                })
            except:
                anova_results.append({
                    'celltype': cell_type,
                    'f_statistic': np.nan,
                    'p_value': 1.0,
                    'test_type': 'Failed',
                    'significant': False
                })
    
    anova_df = pd.DataFrame(anova_results)
    
    # 多重比较校正
    p_values = anova_df['p_value'].fillna(1.0).values
    corrected_p, significant = manual_multiple_test_correction(p_values, method='fdr_bh')
    anova_df['p_value_corrected'] = corrected_p
    anova_df['significant_corrected'] = significant
    
    print("ANOVA分析结果:")
    print(anova_df.to_string(index=False))
    
    return anova_df

def create_visualizations(cross_table, pairwise_results, anova_results, sample_col):
    """创建可视化图表"""
    print("\n=== 生成可视化图表 ===")
    
    # 设置图形布局
    fig = plt.figure(figsize=(20, 15))
    
    # 提取细胞类型
    cell_types = [col for col in cross_table.columns if col != 'group_new' and col != sample_col]
    
    # 1. 堆叠柱状图 - 细胞类型比例
    ax1 = plt.subplot(2, 3, 1)
    
    # 计算各分组比例
    groups = cross_table['group_new'].unique()
    proportion_data = []
    
    for group in groups:
        group_data = cross_table[cross_table['group_new'] == group][cell_types].sum()
        total = group_data.sum()
        proportions = (group_data / total * 100).round(2)
        
        for cell_type, prop in proportions.items():
            proportion_data.append({
                'group': group,
                'celltype': cell_type,
                'proportion': prop
            })
    
    proportion_df = pd.DataFrame(proportion_data)
    pivot_df = proportion_df.pivot(index='group', columns='celltype', values='proportion')
    
    pivot_df.plot(kind='bar', stacked=True, ax=ax1, colormap='tab20', alpha=0.8)
    ax1.set_title('各分组细胞类型比例分布', fontsize=14, fontweight='bold')
    ax1.set_ylabel('比例 (%)')
    ax1.set_xlabel('分组')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. 热图 - 细胞类型比例
    ax2 = plt.subplot(2, 3, 2)
    sns.heatmap(pivot_df.T, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': '比例 (%)'}, ax=ax2)
    ax2.set_title('细胞类型比例热图', fontsize=14, fontweight='bold')
    ax2.set_xlabel('分组')
    ax2.set_ylabel('细胞类型')
    
    # 3. 显著性热图
    ax3 = plt.subplot(2, 3, 3)
    
    if len(pairwise_results) > 0:
        # 创建显著性矩阵
        sig_matrix = []
        unique_comparisons = pairwise_results[['group1', 'group2']].drop_duplicates()
        
        # 选择前10个细胞类型
        top_celltypes = cross_table[cell_types].sum().nlargest(10).index
        
        for celltype in top_celltypes:
            row = []
            for _, comp in unique_comparisons.iterrows():
                mask = ((pairwise_results['group1'] == comp['group1']) & 
                       (pairwise_results['group2'] == comp['group2']) & 
                       (pairwise_results['celltype'] == celltype))
                if mask.any():
                    p_val = pairwise_results[mask]['p_value_corrected'].iloc[0]
                    # 转换为显著性分数
                    sig_score = -np.log10(p_val + 1e-10)
                    row.append(sig_score)
                else:
                    row.append(0)
            
            if row:
                sig_matrix.append(row)
        
        if sig_matrix:
            im = ax3.imshow(sig_matrix, cmap='Reds', aspect='auto')
            ax3.set_title('差异显著性热图', fontsize=14, fontweight='bold')
            ax3.set_xlabel('组间比较')
            ax3.set_ylabel('细胞类型')
            
            # 设置刻度标签
            if len(unique_comparisons) <= 10:
                comp_labels = [f"{row['group1']}vs{row['group2']}" for _, row in unique_comparisons.iterrows()]
                ax3.set_xticks(range(len(comp_labels)))
                ax3.set_xticklabels(comp_labels, rotation=45, ha='right')
            
            ax3.set_yticks(range(len(top_celltypes)))
            ax3.set_yticklabels(top_celltypes)
            plt.colorbar(im, ax=ax3, label='-log10(校正p值)')
    
    # 4. 优势比分布
    ax4 = plt.subplot(2, 3, 4)
    if len(pairwise_results) > 0:
        sig_results = pairwise_results[pairwise_results['significant']]
        if len(sig_results) > 0:
            # 对优势比取log2
            sig_results['log2_OR'] = np.log2(sig_results['odds_ratio'].replace([0, np.inf, -np.inf], np.nan).fillna(1))
            
            scatter = ax4.scatter(sig_results['proportion_diff'], sig_results['log2_OR'],
                                c=-np.log10(sig_results['p_value_corrected'] + 1e-10), 
                                cmap='viridis', s=50, alpha=0.7)
            ax4.axhline(y=0, color='red', linestyle='--', alpha=0.5)
            ax4.axvline(x=0, color='red', linestyle='--', alpha=0.5)
            ax4.set_xlabel('比例差异 (组1 - 组2)')
            ax4.set_ylabel('log2(优势比)')
            ax4.set_title('显著差异结果分布', fontsize=14, fontweight='bold')
            plt.colorbar(scatter, ax=ax4, label='-log10(校正p值)')
    
    # 5. 细胞类型特异性比例箱线图
    ax5 = plt.subplot(2, 3, 5)
    
    # 选择前6个细胞类型
    top_celltypes = cross_table[cell_types].sum().nlargest(6).index
    
    plot_data = []
    for celltype in top_celltypes:
        for group in groups:
            group_data = cross_table[cross_table['group_new'] == group]
            if len(group_data) > 0:
                proportions = group_data[cell_types].apply(lambda x: x[celltype] / x.sum() * 100, axis=1)
                for prop in proportions:
                    plot_data.append({
                        'celltype': celltype,
                        'group': group,
                        'proportion': prop
                    })
    
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        sns.boxplot(data=plot_df, x='celltype', y='proportion', hue='group', ax=ax5)
        ax5.set_title('主要细胞类型比例分布', fontsize=14, fontweight='bold')
        ax5.set_ylabel('比例 (%)')
        ax5.set_xlabel('细胞类型')
        ax5.tick_params(axis='x', rotation=45)
        ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 6. 统计摘要
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    # 执行卡方检验获取整体p值
    chi2, p_value, dof = perform_chi_square_test(cross_table)
    
    summary_text = "统计分析摘要\n\n"
    summary_text += f"总样本数: {len(cross_table)}\n"
    summary_text += f"分组数量: {len(groups)}\n"
    summary_text += f"细胞类型数量: {len(cell_types)}\n"
    summary_text += f"卡方检验p值: {p_value:.4e}\n"
    
    if len(pairwise_results) > 0:
        sig_comparisons = pairwise_results['significant'].sum()
        summary_text += f"显著成对比较: {sig_comparisons}\n"
    
    if len(anova_results) > 0:
        sig_anova = anova_results['significant_corrected'].sum()
        summary_text += f"ANOVA显著细胞类型: {sig_anova}\n"
    
    if p_value < 0.05:
        summary_text += "\n✓ 细胞类型分布在组间存在显著差异\n"
    else:
        summary_text += "\n○ 细胞类型分布在组间无显著差异\n"
    
    ax6.text(0.1, 0.95, summary_text, transform=ax6.transAxes, fontsize=12,
            verticalalignment='top', linespacing=1.5,
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    plt.tight_layout()
    plt.show()

def export_results(cross_table, pairwise_results, anova_results, sample_col, output_prefix='celltype_composition_analysis'):
    """导出分析结果"""
    import os
    
    output_dir = 'composition_analysis_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存计数表
    count_table = cross_table.copy()
    count_table.to_csv(os.path.join(output_dir, f'{output_prefix}_counts.csv'))
    
    # 保存比例表
    cell_types = [col for col in cross_table.columns if col != 'group_new' and col != sample_col]
    groups = cross_table['group_new'].unique()
    proportion_data = []
    
    for group in groups:
        group_data = cross_table[cross_table['group_new'] == group][cell_types].sum()
        total = group_data.sum()
        for cell_type in cell_types:
            proportion = group_data[cell_type] / total * 100 if total > 0 else 0
            proportion_data.append({
                'group': group,
                'celltype': cell_type,
                'proportion': proportion,
                'count': group_data[cell_type]
            })
    
    proportion_df = pd.DataFrame(proportion_data)
    proportion_df.to_csv(os.path.join(output_dir, f'{output_prefix}_proportions.csv'), index=False)
    
    # 保存成对检验结果
    if len(pairwise_results) > 0:
        pairwise_results.to_csv(os.path.join(output_dir, f'{output_prefix}_pairwise_tests.csv'), index=False)
    
    # 保存ANOVA结果
    if len(anova_results) > 0:
        anova_results.to_csv(os.path.join(output_dir, f'{output_prefix}_anova_results.csv'), index=False)
    
    # 保存汇总统计 - 修复这里：使用concat替代append
    chi2, p_value, dof = perform_chi_square_test(cross_table)
    
    # 创建基础汇总统计
    summary_data = {
        'metric': ['总样本数', '分组数', '细胞类型数', '卡方检验p值', '显著差异(卡方)'],
        'value': [
            len(cross_table),
            len(groups),
            len(cell_types),
            p_value,
            p_value < 0.05
        ]
    }
    
    # 使用列表收集所有行数据，然后一次性创建DataFrame
    summary_rows = []
    
    # 添加基础数据
    for i in range(len(summary_data['metric'])):
        summary_rows.append({
            'metric': summary_data['metric'][i],
            'value': summary_data['value'][i]
        })
    
    # 添加成对比较结果
    if len(pairwise_results) > 0:
        summary_rows.append({
            'metric': '显著成对比较数',
            'value': pairwise_results['significant'].sum()
        })
    
    # 添加ANOVA结果
    if len(anova_results) > 0:
        summary_rows.append({
            'metric': 'ANOVA显著细胞类型数',
            'value': anova_results['significant_corrected'].sum()
        })
    
    # 一次性创建最终的DataFrame
    summary_stats = pd.DataFrame(summary_rows)
    summary_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_summary.csv'), index=False)
    
    print(f"\n分析结果已导出到目录: {output_dir}/")

def main():
    """主执行函数"""
    # 设置文件路径
    input_file = "merged_single_cell_data_with_group_new.h5ad"
    
    try:
        # 加载数据
        adata = load_and_prepare_data(input_file)
        if adata is None:
            return
        
        # 创建细胞组成表
        cross_table, sample_col = create_composition_table(adata)
        
        if len(cross_table) < 2:
            print("错误: 样本数量不足，无法进行统计检验")
            return
        
        # 执行卡方检验
        chi2, p_value, dof = perform_chi_square_test(cross_table)
        
        # 执行成对比较分析
        pairwise_results = perform_pairwise_comparisons(cross_table, sample_col)
        
        # 执行ANOVA分析
        anova_results = perform_anova_analysis(cross_table)
        
        # 创建可视化
        create_visualizations(cross_table, pairwise_results, anova_results, sample_col)
        
        # 导出结果
        export_results(cross_table, pairwise_results, anova_results, sample_col)
        
        # 总结报告
        print("\n" + "="*60)
        print("分析完成!")
        print("="*60)
        
        # 显示显著结果
        if len(pairwise_results) > 0:
            sig_results = pairwise_results[pairwise_results['significant']]
            if len(sig_results) > 0:
                print("显著差异结果 (前5个):")
                top5 = sig_results.nsmallest(5, 'p_value_corrected')
                for _, row in top5.iterrows():
                    print(f"{row['celltype']}: {row['group1']} vs {row['group2']}")
                    print(f"  优势比: {row['odds_ratio']:.3f}, 校正p值: {row['p_value_corrected']:.4f}")
            else:
                print("未发现显著差异")
        
    except Exception as e:
        print(f"分析过程中出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
