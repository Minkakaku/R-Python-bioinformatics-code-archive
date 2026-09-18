import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def create_directories(base_dir):
    """创建输出目录"""
    dirs = ['g1', 'g2']
    for dir_name in dirs:
        os.makedirs(os.path.join(base_dir, dir_name), exist_ok=True)
    print(f"创建目录: {dirs}")

def calculate_sample_celltype_proportions(adata):
    """计算每个样本中细胞类型的比例"""
    # 确保必要的列存在
    required_cols = ['sample', 'celltype', 'g1', 'g2']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"缺少必要的列: {missing_cols}")
    
    # 计算每个样本的细胞类型数量
    sample_celltype_counts = pd.crosstab(adata.obs['celltype'], adata.obs['sample'])
    
    # 计算每个样本的细胞类型比例
    sample_celltype_proportions = sample_celltype_counts.div(sample_celltype_counts.sum(axis=0), axis=1) * 100
    
    return sample_celltype_counts, sample_celltype_proportions

def calculate_group_proportions(adata, sample_proportions, group_col):
    """根据分组计算细胞类型比例"""
    # 获取样本到分组的映射
    sample_to_group = dict(zip(adata.obs['sample'], adata.obs[group_col]))
    sample_to_group = {k: v for k, v in sample_to_group.items() if pd.notna(v)}
    
    # 筛选属于某个分组的样本
    groups = adata.obs[group_col].dropna().unique()
    
    # 创建结果字典
    group_results = {}
    
    for group in groups:
        # 获取该组的所有样本
        group_samples = [s for s, g in sample_to_group.items() if g == group]
        
        if group_samples:
            # 筛选属于该组的样本
            group_proportions = sample_proportions[group_samples]
            
            # 计算该组中每个样本的细胞类型比例
            group_results[group] = group_proportions
    
    return group_results

def plot_group_stacked_bar(group_results, group_name, output_dir):
    """绘制分组堆叠条形图"""
    for group, proportions_df in group_results.items():
        if proportions_df.empty:
            continue
            
        plt.figure(figsize=(max(10, len(proportions_df.columns) * 0.5), 8))
        
        # 转置以便于绘图
        plot_df = proportions_df.T
        
        # 创建颜色映射
        n_celltypes = len(plot_df.columns)
        colors = plt.cm.Set3(np.linspace(0, 1, n_celltypes))
        
        # 绘制堆叠条形图
        ax = plot_df.plot(kind='bar', stacked=True, color=colors, 
                         edgecolor='black', linewidth=0.5, figsize=(max(12, len(plot_df) * 0.6), 8))
        
        plt.title(f'{group_name}={group}-Sample proportion', fontsize=16, fontweight='bold')
        plt.xlabel('Sample', fontsize=12)
        plt.ylabel('Proportion (%)', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # 在柱子上添加百分比标签（只显示大于5%的）
        for container in ax.containers:
            ax.bar_label(container, labels=[f'{v.get_height():.1f}%' if v.get_height() > 5 else '' 
                                          for v in container], fontsize=8)
        
        plt.tight_layout()
        
        # 保存图形
        safe_group_name = str(group).replace('/', '_').replace('\\', '_')
        plt.savefig(os.path.join(output_dir, f'{group_name}_{safe_group_name}_stacked_bar.png'), 
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'{group_name}_{safe_group_name}_stacked_bar.pdf'), 
                   bbox_inches='tight')
        plt.close()
        
        # 保存数据
        plot_df.to_csv(os.path.join(output_dir, f'{group_name}_{safe_group_name}_proportions.csv'))

def plot_group_heatmap(group_results, group_name, output_dir):
    """绘制分组热图"""
    for group, proportions_df in group_results.items():
        if proportions_df.empty or proportions_df.shape[1] < 2:
            continue
            
        plt.figure(figsize=(max(10, proportions_df.shape[1] * 0.5), 
                           max(8, proportions_df.shape[0] * 0.3)))
        
        # 创建热图
        sns.heatmap(proportions_df, annot=True, fmt='.1f', cmap='YlOrRd', 
                   cbar_kws={'label': 'Cell proportion(%)'}, linewidths=0.5)
        
        plt.title(f'{group_name}={group} - Heatmap', fontsize=16, fontweight='bold')
        plt.xlabel('Sample', fontsize=12)
        plt.ylabel('Cell type', fontsize=12)
        
        plt.tight_layout()
        
        # 保存图形
        safe_group_name = str(group).replace('/', '_').replace('\\', '_')
        plt.savefig(os.path.join(output_dir, f'{group_name}_{safe_group_name}_heatmap.png'), 
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(output_dir, f'{group_name}_{safe_group_name}_heatmap.pdf'), 
                   bbox_inches='tight')
        plt.close()

def calculate_and_plot_group_comparisons(adata, sample_proportions, group_col, output_dir):
    """计算并绘制组间比较"""
    # 获取所有分组
    groups = adata.obs[group_col].dropna().unique()
    
    if len(groups) < 2:
        print(f"警告: {group_col}只有{len(groups)}个分组，无法进行组间比较")
        return
    
    # 计算每个分组的平均比例
    group_means = {}
    for group in groups:
        group_samples = adata.obs[adata.obs[group_col] == group]['sample'].unique()
        group_samples = [s for s in group_samples if s in sample_proportions.columns]
        
        if group_samples:
            group_data = sample_proportions[group_samples]
            group_means[group] = group_data.mean(axis=1)
    
    # 转换为DataFrame
    mean_df = pd.DataFrame(group_means)
    
    # 保存平均比例
    mean_df.to_csv(os.path.join(output_dir, f'{group_col}_mean_proportions.csv'))
    
    # 绘制分组平均比例堆叠条形图
    plt.figure(figsize=(10, 8))
    
    n_celltypes = len(mean_df)
    colors = plt.cm.Set3(np.linspace(0, 1, n_celltypes))
    
    ax = mean_df.T.plot(kind='bar', stacked=True, color=colors, 
                       edgecolor='black', linewidth=0.5, figsize=(10, 8))
    
    plt.title(f'{group_col} - Average cell proportion', fontsize=16, fontweight='bold')
    plt.xlabel(f'{group_col}Group', fontsize=12)
    plt.ylabel('Average proportion (%)', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Cell type', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 添加百分比标签
    for container in ax.containers:
        ax.bar_label(container, labels=[f'{v.get_height():.1f}%' if v.get_height() > 5 else '' 
                                      for v in container], fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{group_col}_group_means_stacked_bar.png'), 
               dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'{group_col}_group_means_stacked_bar.pdf'), 
               bbox_inches='tight')
    plt.close()
    
    # 绘制分组平均比例热图
    plt.figure(figsize=(8, max(6, len(mean_df) * 0.3)))
    sns.heatmap(mean_df, annot=True, fmt='.1f', cmap='YlOrRd', 
               cbar_kws={'label': 'Average proportion (%)'}, linewidths=0.5)
    plt.title(f'{group_col} - Average proportion for groups', fontsize=16, fontweight='bold')
    plt.xlabel(f'{group_col}Group', fontsize=12)
    plt.ylabel('Cell types', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{group_col}_group_means_heatmap.png'), 
               dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'{group_col}_group_means_heatmap.pdf'), 
               bbox_inches='tight')
    plt.close()
    
    # 进行统计检验（针对每个细胞类型）
    stats_results = []
    celltypes = mean_df.index
    
    for celltype in celltypes:
        # 准备数据
        group_data = []
        for group in groups:
            group_samples = adata.obs[adata.obs[group_col] == group]['sample'].unique()
            group_samples = [s for s in group_samples if s in sample_proportions.columns]
            
            if group_samples:
                celltype_proportions = sample_proportions.loc[celltype, group_samples].dropna()
                if len(celltype_proportions) > 0:
                    group_data.append(celltype_proportions.values)
        
        if len(group_data) >= 2:
            # 使用Kruskal-Wallis检验（非参数，不假设正态分布）
            try:
                h_stat, p_value = stats.kruskal(*group_data)
                stats_results.append({
                    'celltype': celltype,
                    'test': 'Kruskal-Wallis',
                    'statistic': h_stat,
                    'p_value': p_value
                })
            except:
                pass
    
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        # 校正p值
        _, stats_df['p_adjusted'], _, _ = multipletests(stats_df['p_value'], method='fdr_bh')
        stats_df.to_csv(os.path.join(output_dir, f'{group_col}_statistical_tests.csv'), index=False)
        
        # 标记显著差异的细胞类型
        significant = stats_df[stats_df['p_adjusted'] < 0.05]
        if not significant.empty:
            print(f"\n{group_col}中显著差异的细胞类型 (p < 0.05):")
            for _, row in significant.iterrows():
                print(f"  {row['celltype']}: p = {row['p_adjusted']:.4f}")

def generate_summary_report(adata, sample_counts, sample_proportions, output_dir):
    """生成分析总结报告"""
    with open(os.path.join(output_dir, 'analysis_summary.txt'), 'w') as f:
        f.write("细胞类型比例分析总结报告\n")
        f.write("="*50 + "\n\n")
        
        f.write("1. 数据概览\n")
        f.write(f"   总细胞数: {adata.n_obs}\n")
        f.write(f"   样本数量: {adata.obs['sample'].nunique()}\n")
        f.write(f"   细胞类型数量: {adata.obs['celltype'].nunique()}\n\n")
        
        f.write("2. 分组信息\n")
        f.write(f"   g1分组数量: {adata.obs['g1'].nunique()}\n")
        f.write(f"   g1分组: {', '.join([str(x) for x in adata.obs['g1'].dropna().unique()])}\n")
        f.write(f"   g2分组数量: {adata.obs['g2'].nunique()}\n")
        f.write(f"   g2分组: {', '.join([str(x) for x in adata.obs['g2'].dropna().unique()])}\n\n")
        
        f.write("3. 样本细胞数量统计\n")
        sample_stats = adata.obs['sample'].value_counts()
        f.write(sample_stats.to_string())
        f.write("\n\n")
        
        f.write("4. 细胞类型数量统计\n")
        celltype_stats = adata.obs['celltype'].value_counts()
        f.write(celltype_stats.to_string())

def main(h5ad_file, output_base_dir='celltype_analysis'):
    """主函数"""
    print("开始细胞类型比例差异分析...")
    
    # 创建输出目录
    os.makedirs(output_base_dir, exist_ok=True)
    create_directories(output_base_dir)
    
    # 加载数据
    print("加载h5ad数据...")
    adata = sc.read_h5ad(h5ad_file)
    
    # 检查必要的列
    required_cols = ['sample', 'celltype', 'g1', 'g2']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"数据缺少必要的列: {missing_cols}")
    
    print(f"数据加载成功: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    print(f"样本数量: {adata.obs['sample'].nunique()}")
    print(f"细胞类型数量: {adata.obs['celltype'].nunique()}")
    print(f"g1分组: {adata.obs['g1'].unique()}")
    print(f"g2分组: {adata.obs['g2'].unique()}")
    
    try:
        # 计算样本的细胞类型比例
        print("\n计算样本细胞类型比例...")
        sample_counts, sample_proportions = calculate_sample_celltype_proportions(adata)
        
        # 保存样本级别的比例数据
        sample_counts.to_csv(os.path.join(output_base_dir, 'sample_celltype_counts.csv'))
        sample_proportions.to_csv(os.path.join(output_base_dir, 'sample_celltype_proportions.csv'))
        
        # 为g1和g2分别进行分析
        for group_col, output_subdir in [('g1', 'g1'), ('g2', 'g2')]:
            print(f"\n分析 {group_col} 分组...")
            group_output_dir = os.path.join(output_base_dir, output_subdir)
            
            # 计算分组比例
            group_results = calculate_group_proportions(adata, sample_proportions, group_col)
            
            # 绘制分组堆叠条形图
            print(f"  绘制 {group_col} 分组堆叠条形图...")
            plot_group_stacked_bar(group_results, group_col, group_output_dir)
            
            # 绘制分组热图
            print(f"  绘制 {group_col} 分组热图...")
            plot_group_heatmap(group_results, group_col, group_output_dir)
            
            # 计算并绘制组间比较
            print(f"  计算 {group_col} 组间比较...")
            calculate_and_plot_group_comparisons(adata, sample_proportions, group_col, group_output_dir)
        
        # 生成总结报告
        print("\n生成分析总结报告...")
        generate_summary_report(adata, sample_counts, sample_proportions, output_base_dir)
        
        print("\n分析完成！")
        print(f"结果保存在目录: {output_base_dir}")
        print("\n生成的文件结构:")
        for root, dirs, files in os.walk(output_base_dir):
            level = root.replace(output_base_dir, '').count(os.sep)
            indent = ' ' * 2 * level
            print(f"{indent}{os.path.basename(root)}/")
            subindent = ' ' * 2 * (level + 1)
            for file in files:
                print(f"{subindent}{file}")
                
    except Exception as e:
        print(f"分析过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # 替换为您的h5ad文件路径
    h5ad_file_path = "epi_with_g1_g2.h5ad"  # 请修改为实际文件路径
    
    if h5ad_file_path == "your_file.h5ad":
        print("请将脚本中的h5ad_file_path变量修改为您的实际文件路径")
    else:
        main(h5ad_file_path)
