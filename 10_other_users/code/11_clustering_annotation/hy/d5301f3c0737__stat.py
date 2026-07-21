import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import chi2_contingency
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图风格
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
sns.set_style("whitegrid")

def calculate_celltype_proportions(adata, sample_col='sample', celltype_col='celltype', group_col='group'):
    """
    计算细胞类型比例
    """
    # 检查必需的列是否存在
    required_cols = [sample_col, celltype_col, group_col]
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"列 '{col}' 不存在于 adata.obs 中")
    
    print(f"细胞数量: {adata.n_obs}")
    print(f"样本数量: {adata.obs[sample_col].nunique()}")
    print(f"细胞类型数量: {adata.obs[celltype_col].nunique()}")
    print(f"组别数量: {adata.obs[group_col].nunique()}")
    
    # 1. 样本级别的细胞类型比例统计
    sample_celltype_count = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col])
    sample_celltype_prop = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col], normalize='index') * 100
    
    # 2. 组级别的细胞类型比例统计
    group_celltype_count = pd.crosstab(adata.obs[group_col], adata.obs[celltype_col])
    group_celltype_prop = pd.crosstab(adata.obs[group_col], adata.obs[celltype_col], normalize='index') * 100
    
    # 3. 添加组信息到样本级别数据
    sample_to_group = adata.obs[[sample_col, group_col]].drop_duplicates().set_index(sample_col)[group_col]
    sample_celltype_prop['group'] = sample_celltype_prop.index.map(sample_to_group)
    sample_celltype_count['group'] = sample_celltype_count.index.map(sample_to_group)
    
    return {
        'sample_count': sample_celltype_count,
        'sample_proportion': sample_celltype_prop,
        'group_count': group_celltype_count,
        'group_proportion': group_celltype_prop
    }

def plot_stacked_barchart(proportion_data, title, figsize=(12, 8)):
    """
    绘制堆叠柱状图
    """
    plt.figure(figsize=figsize)
    
    # 准备数据
    plot_data = proportion_data.drop('group', axis=1, errors='ignore')
    celltypes = plot_data.columns
    
    # 绘制堆叠柱状图
    bottom = np.zeros(len(plot_data))
    colors = plt.cm.Set3(np.linspace(0, 1, len(celltypes)))
    
    for i, celltype in enumerate(celltypes):
        proportions = plot_data[celltype].values
        plt.bar(plot_data.index, proportions, bottom=bottom, 
                label=celltype, color=colors[i], alpha=0.8)
        bottom += proportions
    
    plt.ylabel('细胞类型比例 (%)', fontsize=12)
    plt.xlabel('')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    return plt.gcf()

def plot_group_comparison(proportion_data, level='sample', figsize=(15, 10)):
    """
    绘制组间比较图
    """
    if level == 'sample':
        data = proportion_data['sample_proportion'].reset_index().melt(
            id_vars=['sample', 'group'], var_name='celltype', value_name='proportion'
        )
        x_var = 'celltype'
        hue_var = 'group'
    else:
        data = proportion_data['group_proportion'].reset_index().melt(
            id_vars=['group'], var_name='celltype', value_name='proportion'
        )
        x_var = 'celltype'
        hue_var = 'group'
    
    # 创建子图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    
    # 箱线图
    sns.boxplot(data=data, x=x_var, y='proportion', hue=hue_var, ax=ax1)
    ax1.set_title(f'{level}级别细胞类型比例分布', fontsize=14, fontweight='bold')
    ax1.set_ylabel('比例 (%)')
    ax1.set_xlabel('细胞类型')
    ax1.tick_params(axis='x', rotation=45)
    
    # 小提琴图
    sns.violinplot(data=data, x=x_var, y='proportion', hue=hue_var, ax=ax2, split=True)
    ax2.set_title(f'{level}级别细胞类型比例分布（小提琴图）', fontsize=14, fontweight='bold')
    ax2.set_ylabel('比例 (%)')
    ax2.set_xlabel('细胞类型')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    return fig

def plot_bubble_chart(proportion_data, count_data, title_prefix, figsize=(12, 8)):
    """
    绘制气泡图显示数量和比例
    """
    # 准备数据
    prop_long = proportion_data.reset_index().melt(
        id_vars=proportion_data.index.names, var_name='celltype', value_name='proportion'
    )
    count_long = count_data.reset_index().melt(
        id_vars=count_data.index.names, var_name='celltype', value_name='count'
    )
    
    # 合并数据
    plot_data = pd.merge(prop_long, count_long, on=proportion_data.index.names + ['celltype'])
    
    plt.figure(figsize=figsize)
    
    # 创建气泡图
    scatter = plt.scatter(plot_data[plot_data.columns[0]], plot_data['celltype'], 
                        s=plot_data['count']/10,  # 调整气泡大小
                        c=plot_data['proportion'], 
                        cmap='viridis', alpha=0.7)
    
    plt.colorbar(scatter, label='比例 (%)')
    plt.xlabel(plot_data.columns[0])
    plt.ylabel('细胞类型')
    plt.title(f'{title_prefix} - 气泡图（大小=数量，颜色=比例）', fontsize=12, fontweight='bold')
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    return plt.gcf()

def statistical_test(adata, celltype_col='celltype', group_col='group'):
    """
    进行统计检验
    """
    # 卡方检验
    contingency_table = pd.crosstab(adata.obs[group_col], adata.obs[celltype_col])
    chi2, p_value, dof, expected = chi2_contingency(contingency_table)
    
    print("=" * 50)
    print("统计检验结果")
    print("=" * 50)
    print(f"卡方值: {chi2:.4f}")
    print(f"P值: {p_value:.4f}")
    print(f"自由度: {dof}")
    
    if p_value < 0.05:
        print("组间细胞组成存在显著差异 (p < 0.05)")
    else:
        print("组间细胞组成无显著差异 (p >= 0.05)")
    
    return chi2, p_value, dof

def create_comprehensive_report(adata, output_dir='celltype_proportion_analysis'):
    """
    创建完整的分析报告
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("开始细胞类型比例分析...")
    
    # 计算比例
    results = calculate_celltype_proportions(adata)
    
    # 保存统计结果
    with pd.ExcelWriter(os.path.join(output_dir, 'celltype_proportion_analysis.xlsx')) as writer:
        results['sample_count'].to_excel(writer, sheet_name='样本细胞数量')
        results['sample_proportion'].to_excel(writer, sheet_name='样本细胞比例')
        results['group_count'].to_excel(writer, sheet_name='组细胞数量')
        results['group_proportion'].to_excel(writer, sheet_name='组细胞比例')
    
    # 可视化
    print("生成可视化图表...")
    
    # 1. 样本级别堆叠柱状图
    fig1 = plot_stacked_barchart(
        results['sample_proportion'].drop('group', axis=1), 
        '各样本细胞类型比例分布'
    )
    fig1.savefig(os.path.join(output_dir, 'sample_stacked_barchart.png'), dpi=300, bbox_inches='tight')
    plt.close(fig1)
    
    # 2. 组级别堆叠柱状图
    fig2 = plot_stacked_barchart(
        results['group_proportion'], 
        '各组细胞类型比例分布'
    )
    fig2.savefig(os.path.join(output_dir, 'group_stacked_barchart.png'), dpi=300, bbox_inches='tight')
    plt.close(fig2)
    
    # 3. 组间比较图
    fig3 = plot_group_comparison(results, level='sample')
    fig3.savefig(os.path.join(output_dir, 'sample_group_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close(fig3)
    
    fig4 = plot_group_comparison(results, level='group')
    fig4.savefig(os.path.join(output_dir, 'group_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close(fig4)
    
    # 4. 气泡图
    fig5 = plot_bubble_chart(
        results['sample_proportion'].drop('group', axis=1),
        results['sample_count'].drop('group', axis=1),
        '样本级别'
    )
    fig5.savefig(os.path.join(output_dir, 'sample_bubble_chart.png'), dpi=300, bbox_inches='tight')
    plt.close(fig5)
    
    # 5. 统计检验
    chi2, p_value, dof = statistical_test(adata)
    
    # 保存统计检验结果
    stats_df = pd.DataFrame({
        'chi2': [chi2],
        'p_value': [p_value], 
        'dof': [dof],
        'significance': ['显著' if p_value < 0.05 else '不显著']
    })
    stats_df.to_excel(os.path.join(output_dir, 'statistical_test.xlsx'), index=False)
    
    # 6. 创建详细的比例表格（用于进一步分析）
    create_detailed_analysis_table(adata, results, output_dir)
    
    print(f"\n分析完成！结果保存在: {output_dir}/")
    print(f"统计表格: celltype_proportion_analysis.xlsx")
    print(f"可视化图表: 多个PNG文件")
    print(f"统计检验结果: statistical_test.xlsx")

def create_detailed_analysis_table(adata, results, output_dir):
    """
    创建详细分析表格
    """
    # 样本级别详细数据
    sample_detailed = results['sample_proportion'].reset_index()
    sample_detailed = sample_detailed.melt(id_vars=['sample', 'group'], 
                                          var_name='celltype', value_name='proportion')
    
    # 添加绝对数量
    sample_count_melted = results['sample_count'].reset_index().melt(
        id_vars=['sample', 'group'], var_name='celltype', value_name='count'
    )
    sample_detailed = pd.merge(sample_detailed, sample_count_melted, on=['sample', 'group', 'celltype'])
    
    # 计算每个细胞类型在各组中的统计量
    celltype_stats = sample_detailed.groupby(['celltype', 'group'])['proportion'].agg([
        'mean', 'std', 'min', 'max', 'count'
    ]).round(3)
    
    # 保存详细表格
    with pd.ExcelWriter(os.path.join(output_dir, 'detailed_analysis.xlsx')) as writer:
        sample_detailed.to_excel(writer, sheet_name='样本详细数据', index=False)
        celltype_stats.to_excel(writer, sheet_name='细胞类型统计')
        
        # 各样本总细胞数
        sample_totals = adata.obs['sample'].value_counts().reset_index()
        sample_totals.columns = ['sample', 'total_cells']
        sample_totals.to_excel(writer, sheet_name='样本总细胞数', index=False)

# 主执行函数
def main():
    """
    主函数：执行完整的细胞类型比例分析
    """
    # 读取h5ad文件（请修改为你的文件路径）
    h5ad_file_path = 'merged_single_cell_data.h5ad'  # 修改为你的文件路径
    
    try:
        print(f"读取数据文件: {h5ad_file_path}")
        adata = sc.read_h5ad(h5ad_file_path)
        
        # 检查数据基本信息
        print(f"数据形状: {adata.shape}")
        print(f"观察值(细胞)数量: {adata.n_obs}")
        print(f"变量(基因)数量: {adata.n_vars}")
        print(f"元数据列: {list(adata.obs.columns)}")
        
        # 执行分析
        create_comprehensive_report(adata, output_dir='celltype_proportion_analysis')
        
    except FileNotFoundError:
        print(f"错误: 文件 {h5ad_file_path} 未找到")
        print("请修改脚本中的文件路径为你的h5ad文件路径")
    except Exception as e:
        print(f"错误: {str(e)}")
        print("请检查数据文件格式和列名")

# 如果是直接运行脚本，则执行主函数
if __name__ == "__main__":
    main()

# 使用示例（取消注释以下代码来测试）
"""
# 示例1：基本使用
# adata = sc.read_h5ad('your_data.h5ad')
# create_comprehensive_report(adata)

# 示例2：自定义列名
# results = calculate_celltype_proportions(adata, sample_col='sample_id', celltype_col='cell_type', group_col='treatment')
"""
