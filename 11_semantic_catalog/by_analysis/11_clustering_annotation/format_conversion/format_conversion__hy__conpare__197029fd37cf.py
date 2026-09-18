import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def create_compare_directory():
    """创建保存结果的目录"""
    if not os.path.exists('compare'):
        os.makedirs('compare')
    print("创建 compare 目录成功")

def load_and_validate_data(h5ad_file):
    """加载h5ad文件并验证必要的列存在"""
    try:
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据加载成功: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
        
        # 验证必要的列存在
        if 'group' not in adata.obs.columns:
            raise ValueError("obs中未找到'group'列")
        if 'celltype' not in adata.obs.columns:
            raise ValueError("obs中未找到'celltype'列")
            
        print(f"分组信息: {adata.obs['group'].unique()}")
        print(f"细胞类型: {adata.obs['celltype'].unique()}")
        
        return adata
    except Exception as e:
        print(f"数据加载失败: {e}")
        return None

def calculate_celltype_proportions(adata):
    """计算各分组中细胞类型的比例"""
    # 创建交叉表计算细胞数量
    count_table = pd.crosstab(adata.obs['celltype'], adata.obs['group'])
    
    # 计算比例（按列归一化）
    proportion_table = count_table.div(count_table.sum(axis=0), axis=1) * 100
    
    # 转置以便于绘图
    proportion_table_T = proportion_table.T
    
    # 保存比例数据
    count_table.to_csv('compare/celltype_counts.csv')
    proportion_table_T.to_csv('compare/celltype_proportions.csv')
    
    print("细胞类型数量和比例计算完成")
    return count_table, proportion_table_T

def statistical_analysis(adata, proportion_table):
    """进行统计学分析"""
    results = []
    celltypes = proportion_table.columns
    groups = proportion_table.index
    
    # 对每个细胞类型进行组间比较
    for celltype in celltypes:
        # 准备数据：每个细胞的细胞类型标签（0/1）和分组信息
        celltype_data = []
        group_labels = []
        
        for group in groups:
            # 筛选该组该细胞类型的细胞
            mask = (adata.obs['group'] == group) & (adata.obs['celltype'] == celltype)
            celltype_cells = mask.sum()
            total_cells = (adata.obs['group'] == group).sum()
            
            celltype_data.extend([1] * celltype_cells + [0] * (total_cells - celltype_cells))
            group_labels.extend([group] * total_cells)
        
        # 转换为DataFrame
        df_test = pd.DataFrame({
            'celltype_presence': celltype_data,
            'group': group_labels
        })
        
        # 单因素ANOVA或Kruskal-Wallis检验
        if len(groups) > 2:
            # 准备每组的数据
            group_data = []
            for group in groups:
                group_mask = df_test['group'] == group
                group_data.append(df_test['celltype_presence'][group_mask].values)
            
            # 方差齐性检验
            _, levene_p = stats.levene(*group_data)
            
            if levene_p > 0.05:  # 方差齐性，使用ANOVA
                f_stat, anova_p = stats.f_oneway(*group_data)
                test_used = "ANOVA"
            else:  # 方差不齐，使用Kruskal-Wallis
                h_stat, kruskal_p = stats.kruskal(*group_data)
                anova_p = kruskal_p
                test_used = "Kruskal-Wallis"
                f_stat = h_stat
        else:
            # 两组比较使用Mann-Whitney U检验
            group1_data = df_test[df_test['group'] == groups[0]]['celltype_presence']
            group2_data = df_test[df_test['group'] == groups[1]]['celltype_presence']
            u_stat, anova_p = stats.mannwhitneyu(group1_data, group2_data)
            test_used = "Mann-Whitney U"
            f_stat = u_stat
        
        results.append({
            'celltype': celltype,
            'test_used': test_used,
            'statistic': f_stat,
            'p_value': anova_p
        })
    
    # 转换为DataFrame并校正p值
    stats_df = pd.DataFrame(results)
    if not stats_df.empty:
        _, stats_df['p_adjusted'], _, _ = multipletests(stats_df['p_value'], method='fdr_bh')
    
    stats_df.to_csv('compare/statistical_analysis.csv', index=False)
    print("统计学分析完成")
    return stats_df

def plot_stacked_bar(proportion_table, adata):
    """绘制堆叠条形图"""
    plt.figure(figsize=(12, 8))
    
    # 创建颜色映射
    n_celltypes = len(proportion_table.columns)
    colors = plt.cm.Set3(np.linspace(0, 1, n_celltypes))
    
    # 绘制堆叠条形图
    ax = proportion_table.plot(kind='bar', stacked=True, color=colors, 
                              edgecolor='black', linewidth=0.5, figsize=(12, 8))
    
    plt.title('各分组细胞类型比例分布', fontsize=16, fontweight='bold')
    plt.xlabel('分组', fontsize=12)
    plt.ylabel('细胞比例 (%)', fontsize=12)
    plt.xticks(rotation=45)
    plt.legend(title='细胞类型', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # 在柱子上添加百分比标签
    for container in ax.containers:
        ax.bar_label(container, labels=[f'{v.get_height():.1f}%' if v.get_height() > 5 else '' 
                                        for v in container], fontsize=8)
    
    plt.savefig('compare/stacked_barplot.png', dpi=300, bbox_inches='tight')
    plt.savefig('compare/stacked_barplot.pdf', bbox_inches='tight')
    plt.close()

def plot_heatmap(proportion_table):
    """绘制热图展示比例差异"""
    plt.figure(figsize=(10, 8))
    
    # 创建热图
    sns.heatmap(proportion_table, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': '细胞比例 (%)'}, linewidths=0.5)
    
    plt.title('各分组细胞类型比例热图', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('compare/proportion_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('compare/proportion_heatmap.pdf', bbox_inches='tight')
    plt.close()

def plot_group_comparison(proportion_table, stats_df):
    """绘制组间比较点图"""
    if len(proportion_table) < 2:
        print("分组数量不足，跳过组间比较图")
        return
    
    # 选择有显著差异的细胞类型（p < 0.1）
    significant_celltypes = stats_df[stats_df['p_adjusted'] < 0.1]['celltype']
    
    if len(significant_celltypes) > 0:
        plt.figure(figsize=(10, 6))
        
        for i, celltype in enumerate(significant_celltypes):
            plt.subplot(1, len(significant_celltypes), i+1)
            proportions = proportion_table[celltype]
            
            # 绘制点图
            plt.scatter(range(len(proportions)), proportions, s=100, alpha=0.7)
            plt.plot(range(len(proportions)), proportions, 'r--', alpha=0.5)
            
            plt.title(f'{celltype}\n(p={stats_df[stats_df["celltype"]==celltype]["p_adjusted"].values[0]:.3f})')
            plt.xlabel('分组')
            plt.ylabel('比例 (%)')
            plt.xticks(range(len(proportions)), proportions.index, rotation=45)
        
        plt.tight_layout()
        plt.savefig('compare/group_comparison_dotplot.png', dpi=300, bbox_inches='tight')
        plt.savefig('compare/group_comparison_dotplot.pdf', bbox_inches='tight')
        plt.close()

def generate_summary_report(adata, count_table, proportion_table, stats_df):
    """生成分析总结报告"""
    with open('compare/analysis_summary.txt', 'w') as f:
        f.write("细胞类型比例差异分析报告\n")
        f.write("="*50 + "\n\n")
        
        f.write("1. 数据概览\n")
        f.write(f"   总细胞数: {adata.n_obs}\n")
        f.write(f"   分组数量: {len(adata.obs['group'].unique())}\n")
        f.write(f"   细胞类型数量: {len(adata.obs['celltype'].unique())}\n\n")
        
        f.write("2. 各分组细胞数量\n")
        group_counts = adata.obs['group'].value_counts()
        for group, count in group_counts.items():
            f.write(f"   {group}: {count} 个细胞\n")
        f.write("\n")
        
        f.write("3. 细胞类型比例统计\n")
        f.write(proportion_table.to_string())
        f.write("\n\n")
        
        f.write("4. 统计学显著性分析\n")
        f.write(stats_df.to_string())
        f.write("\n\n")
        
        # 找出显著差异的细胞类型
        significant = stats_df[stats_df['p_adjusted'] < 0.05]
        if len(significant) > 0:
            f.write("5. 显著差异的细胞类型 (p < 0.05)\n")
            for _, row in significant.iterrows():
                f.write(f"   {row['celltype']}: p = {row['p_adjusted']:.4f}\n")
        else:
            f.write("5. 未发现显著差异的细胞类型 (p < 0.05)\n")

def main(h5ad_file):
    """主函数"""
    print("开始细胞类型比例差异分析...")
    
    # 创建目录
    create_compare_directory()
    
    # 加载数据
    adata = load_and_validate_data(h5ad_file)
    if adata is None:
        return
    
    try:
        # 计算比例
        count_table, proportion_table = calculate_celltype_proportions(adata)
        
        # 统计学分析
        stats_df = statistical_analysis(adata, proportion_table)
        
        # 绘制图形
        plot_stacked_bar(proportion_table, adata)
        plot_heatmap(proportion_table)
        plot_group_comparison(proportion_table, stats_df)
        
        # 生成报告
        generate_summary_report(adata, count_table, proportion_table, stats_df)
        
        print("分析完成！结果保存在 compare 目录中")
        print(f"生成的文件:")
        for file in os.listdir('compare'):
            print(f"  - {file}")
            
    except Exception as e:
        print(f"分析过程中出现错误: {e}")

if __name__ == "__main__":
    # 替换为您的h5ad文件路径
    h5ad_file_path = "merged_single_cell_data_modified.h5ad"  # 请修改为实际文件路径
    
    if h5ad_file_path == "your_data.h5ad":
        print("请将脚本中的h5ad_file_path变量修改为您的实际文件路径")
    else:
        main(h5ad_file_path)
