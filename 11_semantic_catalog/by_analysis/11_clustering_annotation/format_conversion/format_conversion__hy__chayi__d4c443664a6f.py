import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from scipy.stats import chi2_contingency, fisher_exact
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图风格
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'Arial']
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
    
    return adata

def create_celltype_count_table(adata):
    """创建细胞类型计数表"""
    print("\n=== 创建细胞类型计数表 ===")
    
    # 创建样本-细胞类型交叉表
    cross_table = pd.crosstab(adata.obs['group_new'], adata.obs['celltype'])
    print("细胞类型计数表:")
    print(cross_table)
    
    # 计算比例表
    proportion_table = cross_table.div(cross_table.sum(axis=1), axis=0)
    print("\n细胞类型比例表:")
    print(proportion_table.round(3))
    
    return cross_table, proportion_table

def perform_chi_square_test(cross_table):
    """执行卡方检验"""
    print("\n=== 卡方检验 - 整体差异 ===")
    
    # 执行卡方检验
    chi2, p_value, dof, expected = chi2_contingency(cross_table)
    
    print(f"卡方值: {chi2:.4f}")
    print(f"P值: {p_value:.4e}")
    print(f"自由度: {dof}")
    
    if p_value < 0.05:
        print("结论: 细胞类型分布在组间存在显著差异 (p < 0.05)")
    else:
        print("结论: 细胞类型分布在组间无显著差异 (p ≥ 0.05)")
    
    return chi2, p_value, dof

def perform_pairwise_fisher_tests(cross_table):
    """执行成对Fisher精确检验"""
    print("\n=== 成对Fisher精确检验 ===")
    
    groups = cross_table.index.tolist()
    results = []
    
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            group1, group2 = groups[i], groups[j]
            
            # 提取两个组的计数数据
            data_subset = cross_table.loc[[group1, group2]]
            
            # 执行Fisher精确检验
            oddsratio, p_value = fisher_exact(data_subset)
            
            results.append({
                'group1': group1,
                'group2': group2,
                'odds_ratio': oddsratio,
                'p_value': p_value
            })
            
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            print(f"{group1} vs {group2}: OR = {oddsratio:.3f}, p = {p_value:.4f} {significance}")
    
    return pd.DataFrame(results)

def calculate_odds_ratios(cross_table):
    """计算优势比(Odds Ratios)"""
    print("\n=== 优势比分析 ===")
    
    # 计算每个细胞类型在每个组中的比例
    proportion_table = cross_table.div(cross_table.sum(axis=1), axis=0)
    
    # 计算相对于所有组平均值的优势比
    overall_means = proportion_table.mean(axis=0)
    odds_ratios = proportion_table.div(overall_means, axis=1)
    
    print("优势比矩阵 (值>1表示在该组中富集):")
    print(odds_ratios.round(3))
    
    return odds_ratios

def create_visualizations(adata, cross_table, proportion_table, odds_ratios):
    """创建可视化图表"""
    print("\n=== 生成可视化图表 ===")
    
    # 设置图形布局
    fig = plt.figure(figsize=(20, 16))
    
    # 1. 堆叠柱状图 - 细胞类型比例
    ax1 = plt.subplot(2, 3, 1)
    proportion_table.plot(kind='bar', stacked=True, ax=ax1, 
                         colormap='tab20', alpha=0.8)
    ax1.set_title('各分组细胞类型比例分布', fontsize=14, fontweight='bold')
    ax1.set_ylabel('比例')
    ax1.set_xlabel('分组')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. 热图 - 细胞类型计数
    ax2 = plt.subplot(2, 3, 2)
    # 对计数进行log转换以便更好可视化
    log_counts = np.log10(cross_table + 1)
    sns.heatmap(log_counts, annot=cross_table, fmt='d', cmap='YlOrRd', 
                cbar_kws={'label': 'Log10(细胞数量+1)'}, ax=ax2)
    ax2.set_title('各分组细胞类型计数热图', fontsize=14, fontweight='bold')
    
    # 3. 优势比热图
    ax3 = plt.subplot(2, 3, 3)
    # 创建对称的颜色映射，以1为中心
    vmax = max(2, odds_ratios.max().max())
    vmin = min(0.5, odds_ratios.min().min())
    sns.heatmap(odds_ratios, annot=True, fmt='.2f', cmap='RdBu_r', 
                center=1, vmin=vmin, vmax=vmax, ax=ax3)
    ax3.set_title('优势比热图 (值>1表示富集)', fontsize=14, fontweight='bold')
    
    # 4. 箱线图 - 细胞类型分布
    ax4 = plt.subplot(2, 3, 4)
    # 准备数据
    plot_data = []
    for group in cross_table.index:
        for celltype in cross_table.columns:
            count = cross_table.loc[group, celltype]
            total = cross_table.loc[group].sum()
            proportion = count / total if total > 0 else 0
            plot_data.append({
                'group': group,
                'celltype': celltype,
                'proportion': proportion * 100  # 转换为百分比
            })
    
    plot_df = pd.DataFrame(plot_data)
    # 选择前9个最多的细胞类型进行可视化
    top_celltypes = cross_table.sum().nlargest(9).index
    plot_df_top = plot_df[plot_df['celltype'].isin(top_celltypes)]
    
    sns.boxplot(data=plot_df_top, x='celltype', y='proportion', hue='group', ax=ax4)
    ax4.set_title('主要细胞类型在各分组中的比例分布', fontsize=14, fontweight='bold')
    ax4.set_ylabel('比例 (%)')
    ax4.set_xlabel('细胞类型')
    ax4.tick_params(axis='x', rotation=45)
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 5. 分组差异雷达图
    ax5 = plt.subplot(2, 3, 5, polar=True)
    
    # 准备雷达图数据
    angles = np.linspace(0, 2 * np.pi, len(cross_table.columns), endpoint=False).tolist()
    angles += angles[:1]  # 闭合图形
    
    # 绘制每个组的数据
    for i, group in enumerate(cross_table.index):
        values = proportion_table.loc[group].values.tolist()
        values += values[:1]  # 闭合图形
        ax5.plot(angles, values, 'o-', linewidth=2, label=group)
        ax5.fill(angles, values, alpha=0.1)
    
    ax5.set_xticks(angles[:-1])
    ax5.set_xticklabels(cross_table.columns)
    ax5.set_title('细胞类型比例雷达图', size=14, fontweight='bold')
    ax5.legend(bbox_to_anchor=(1.3, 1), loc='upper left')
    
    # 6. 统计摘要
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    # 执行卡方检验获取整体p值
    chi2, p_value, dof = perform_chi_square_test(cross_table)
    
    summary_text = "统计分析摘要\n\n"
    summary_text += f"总细胞数: {cross_table.sum().sum():,}\n"
    summary_text += f"分组数量: {len(cross_table)}\n"
    summary_text += f"细胞类型数量: {len(cross_table.columns)}\n"
    summary_text += f"卡方检验p值: {p_value:.4e}\n\n"
    
    if p_value < 0.05:
        summary_text += "✓ 细胞类型分布在组间存在显著差异\n"
    else:
        summary_text += "○ 细胞类型分布在组间无显著差异\n"
    
    # 添加每个组的细胞数量信息
    summary_text += "\n各分组细胞数量:\n"
    for group in cross_table.index:
        count = cross_table.loc[group].sum()
        summary_text += f"  {group}: {count:,} 细胞\n"
    
    ax6.text(0.1, 0.95, summary_text, transform=ax6.transAxes, fontsize=12,
            verticalalignment='top', linespacing=1.5,
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    plt.tight_layout()
    plt.show()
    
    return fig

def perform_detailed_celltype_analysis(adata, cross_table):
    """对每种细胞类型进行详细的组间比较"""
    print("\n=== 细胞类型特异性分析 ===")
    
    results = []
    celltypes = cross_table.columns
    
    for celltype in celltypes:
        print(f"\n分析细胞类型: {celltype}")
        
        # 提取该细胞类型在各组中的计数
        celltype_counts = cross_table[celltype]
        group_totals = cross_table.sum(axis=1)
        
        # 计算比例
        proportions = (celltype_counts / group_totals * 100).round(2)
        
        print(f"在各组中的比例: {proportions.to_dict()}")
        
        # 如果只有两组，执行t检验
        if len(celltype_counts) == 2:
            groups = celltype_counts.index.tolist()
            prop1 = proportions[groups[0]]
            prop2 = proportions[groups[1]]
            
            # 执行比例差异的z检验
            from statsmodels.stats.proportion import proportions_ztest
            
            count = [celltype_counts[groups[0]], celltype_counts[groups[1]]]
            nobs = [group_totals[groups[0]], group_totals[groups[1]]]
            
            z_stat, p_value = proportions_ztest(count, nobs)
            
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            print(f"比例差异检验: z = {z_stat:.3f}, p = {p_value:.4f} {significance}")
            
            results.append({
                'celltype': celltype,
                'group1': groups[0],
                'group2': groups[1],
                'prop1': prop1,
                'prop2': prop2,
                'z_stat': z_stat,
                'p_value': p_value,
                'significance': significance
            })
    
    return pd.DataFrame(results)

def export_results(cross_table, proportion_table, odds_ratios, pairwise_results, output_prefix='celltype_group_analysis'):
    """导出分析结果"""
    import os
    
    output_dir = 'celltype_composition_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存计数表
    cross_table.to_csv(os.path.join(output_dir, f'{output_prefix}_counts.csv'))
    
    # 保存比例表
    proportion_table.to_csv(os.path.join(output_dir, f'{output_prefix}_proportions.csv'))
    
    # 保存优势比
    odds_ratios.to_csv(os.path.join(output_dir, f'{output_prefix}_odds_ratios.csv'))
    
    # 保存成对检验结果
    if pairwise_results is not None and len(pairwise_results) > 0:
        pairwise_results.to_csv(os.path.join(output_dir, f'{output_prefix}_pairwise_tests.csv'), index=False)
    
    # 保存汇总统计
    summary_stats = pd.DataFrame({
        '统计量': ['总细胞数', '分组数', '细胞类型数', '卡方检验p值'],
        '值': [cross_table.sum().sum(), len(cross_table), len(cross_table.columns), 
              chi2_contingency(cross_table)[1]]
    })
    summary_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_summary.csv'), index=False)
    
    print(f"\n分析结果已导出到目录: {output_dir}/")

def main():
    """主执行函数"""
    # 设置文件路径
    input_file = "merged_single_cell_data_with_group_new.h5ad"
    
    # 加载数据
    adata = load_and_prepare_data(input_file)
    if adata is None:
        return
    
    # 创建细胞类型计数表
    cross_table, proportion_table = create_celltype_count_table(adata)
    
    # 执行卡方检验
    chi2, p_value, dof = perform_chi_square_test(cross_table)
    
    # 执行成对Fisher检验
    pairwise_results = perform_pairwise_fisher_tests(cross_table)
    
    # 计算优势比
    odds_ratios = calculate_odds_ratios(cross_table)
    
    # 创建可视化
    create_visualizations(adata, cross_table, proportion_table, odds_ratios)
    
    # 执行详细的细胞类型分析
    detailed_results = perform_detailed_celltype_analysis(adata, cross_table)
    
    # 导出结果
    export_results(cross_table, proportion_table, odds_ratios, pairwise_results)
    
    print("\n" + "="*60)
    print("分析完成!")
    print("="*60)
    print("关键发现:")
    print(f"- 分析了 {cross_table.sum().sum():,} 个细胞")
    print(f"- 比较 {len(cross_table)} 个分组")
    print(f"- 涵盖 {len(cross_table.columns)} 种细胞类型")
    
    if p_value < 0.05:
        print("- ✓ 细胞类型分布在组间存在显著差异")
    else:
        print("- ○ 细胞类型分布在组间无显著差异")

if __name__ == "__main__":
    main()
