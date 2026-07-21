import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def check_and_fix_index_issues(adata):
    """
    检查并修复AnnData对象的索引问题
    """
    print("正在检查索引唯一性...")
    
    # 检查obs索引是否唯一
    if not adata.obs.index.is_unique:
        print(f"警告: 发现 {adata.obs.index.duplicated().sum()} 个重复的obs索引")
        print("重复的索引示例:", adata.obs.index[adata.obs.index.duplicated()].tolist()[:5])
        
        # 修复重复索引
        adata.obs_names_make_unique()
        print("已修复重复的obs索引")
    
    # 检查var索引是否唯一
    if not adata.var.index.is_unique:
        print(f"警告: 发现 {adata.var.index.duplicated().sum()} 个重复的var索引")
        adata.var_names_make_unique()
        print("已修复重复的var索引")
    
    return adata

def analyze_cell_composition(h5ad_path, output_dir="cell_composition_analysis"):
    """
    分析特定样本组与其他样本Adenoma阶段的细胞组成
    """
    # 创建输出目录
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # 1. 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 检查并修复索引问题
    adata = check_and_fix_index_issues(adata)
    
    # 检查必要的列
    required_columns = ['sample', 'celltype', 'stage']
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    if missing_columns:
        raise ValueError(f"缺少必要的列: {missing_columns}。可用的列有: {list(adata.obs.columns)}")
    
    print(f"数据概览:")
    print(f"总细胞数: {adata.n_obs}")
    print(f"样本列表: {list(adata.obs['sample'].unique())}")
    print(f"细胞类型: {list(adata.obs['celltype'].unique())}")
    print(f"Stage类型: {list(adata.obs['stage'].unique())}")
    
    # 2. 定义分析组
    # 组1: SR202500002C-1和SR202500002C-2样本的所有细胞
    target_samples = ['SR202500002C-1', 'SR202500002C-2']
    
    # 组2: 其他样本中stage为Adenoma的细胞
    # 修复索引问题：确保使用唯一索引
    other_samples_adenoma_mask = (
        (~adata.obs['sample'].isin(target_samples)) & 
        (adata.obs['stage'] == 'Adenoma')
    )
    
    # 使用安全的索引方式
    other_samples_adenoma_indices = adata.obs[other_samples_adenoma_mask].index
    print(f"其他样本Adenoma阶段细胞索引数量: {len(other_samples_adenoma_indices)}")
    
    # 检查索引唯一性
    if not other_samples_adenoma_indices.is_unique:
        print("警告: 其他样本Adenoma阶段的索引不唯一，正在修复...")
        # 去除重复索引
        other_samples_adenoma_indices = other_samples_adenoma_indices[~other_samples_adenoma_indices.duplicated()]
        print(f"修复后索引数量: {len(other_samples_adenoma_indices)}")
    
    if len(other_samples_adenoma_indices) == 0:
        print("警告: 未找到其他样本中的Adenoma阶段细胞")
        print(f"可用的stage值: {list(adata.obs['stage'].unique())}")
        return
    
    print(f"\n=== 分析组统计 ===")
    print(f"目标样本组 ({', '.join(target_samples)}) 细胞数: {len(adata.obs[adata.obs['sample'].isin(target_samples)])}")
    print(f"其他样本Adenoma阶段细胞数: {len(other_samples_adenoma_indices)}")
    
    # 3. 计算细胞比例
    def calculate_proportions(adata_subset, group_name):
        """计算指定细胞组的细胞类型比例"""
        if len(adata_subset) == 0:
            return pd.DataFrame()
        
        # 计算各细胞类型的数量
        celltype_counts = adata_subset.obs['celltype'].value_counts()
        
        # 计算比例
        total_cells = len(adata_subset)
        proportions = (celltype_counts / total_cells * 100).round(2)
        
        # 创建结果DataFrame
        result_df = pd.DataFrame({
            'celltype': proportions.index,
            'proportion': proportions.values,
            'count': celltype_counts.values,
            'total_cells': total_cells,
            'group': group_name
        })
        
        return result_df
    
    # 计算各组比例 - 使用安全的索引方式
    group1_mask = adata.obs['sample'].isin(target_samples)
    group1_indices = adata.obs[group1_mask].index
    if not group1_indices.is_unique:
        group1_indices = group1_indices[~group1_indices.duplicated()]
    
    group1_data = adata[group1_indices, :]
    group1_df = calculate_proportions(group1_data, 'Adenoma Polyps')
    
    group2_data = adata[other_samples_adenoma_indices, :]
    group2_df = calculate_proportions(group2_data, 'Adenoma Adjacent')
    
    # 合并结果
    all_results = pd.concat([group1_df, group2_df], ignore_index=True)
    
    # 4. 保存比例结果
    results_file = output_dir / "cell_composition_results.csv"
    all_results.to_csv(results_file, index=False)
    print(f"\n比例结果已保存: {results_file}")
    
    # 5. 可视化
    print("\n正在生成可视化图表...")
    
    # 设置图形样式
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    sns.set_style("whitegrid")
    
    # 创建颜色映射
    celltypes = all_results['celltype'].unique()
    palette = sns.color_palette('Set3', len(celltypes))
    color_dict = dict(zip(celltypes, palette))
    
    # 图1: 堆积条形图比较两组细胞组成
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 左图: 比例堆积条形图
    if not all_results.empty:
        pivot_df = all_results.pivot(index='group', columns='celltype', values='proportion').fillna(0)
        
        # 确保两组都有所有细胞类型（填充0）
        for celltype in celltypes:
            if celltype not in pivot_df.columns:
                pivot_df[celltype] = 0
        
        # 按比例排序细胞类型（按第一组）
        if not pivot_df.empty and len(pivot_df) > 0:
            celltype_order = pivot_df.iloc[0].sort_values(ascending=False).index
            pivot_df = pivot_df[celltype_order]
        
        # 绘制堆积条形图
        bottom = np.zeros(len(pivot_df))
        for i, celltype in enumerate(pivot_df.columns):
            proportions = pivot_df[celltype].values
            ax1.bar(pivot_df.index, proportions, bottom=bottom, 
                    label=celltype, color=color_dict[celltype], alpha=0.8)
            bottom += proportions
        
        ax1.set_ylabel('Proportion')
        ax1.set_title('Cell types')
        ax1.tick_params(axis='x', rotation=45)
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Celltypes')
    
    # 右图: 分组条形图显示主要细胞类型
    if not all_results.empty:
        # 选择比例最高的几种细胞类型进行详细比较
        top_celltypes = all_results.groupby('celltype')['proportion'].max().nlargest(min(10, len(celltypes))).index
        
        top_data = all_results[all_results['celltype'].isin(top_celltypes)]
        
        sns.barplot(data=top_data, x='celltype', y='proportion', hue='group', 
                    palette='viridis', ax=ax2)
        ax2.set_ylabel('Proportion')
        ax2.set_xlabel('Cell types')
        ax2.set_title('Main celltypes proportions')
        ax2.tick_params(axis='x', rotation=45)
        ax2.legend(title='Group')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'cell_composition_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 6. 打印统计摘要
    print("\n=== 统计分析结果 ===")
    print("各组细胞数量统计:")
    group_summary = all_results.groupby('group').agg({
        'total_cells': 'first',
        'celltype': 'count'
    }).rename(columns={'celltype': '细胞类型数量'})
    print(group_summary)
    
    if not all_results.empty:
        print("\n比例最高的5种细胞类型:")
        for group in all_results['group'].unique():
            group_data = all_results[all_results['group'] == group]
            top5 = group_data.nlargest(5, 'proportion')[['celltype', 'proportion']]
            print(f"\n{group}:")
            for _, row in top5.iterrows():
                print(f"  {row['celltype']}: {row['proportion']}%")
    
    # 7. 保存详细报告
    report_file = output_dir / "analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write("细胞组成分析报告\n")
        f.write("="*50 + "\n")
        f.write(f"分析时间: {pd.Timestamp.now()}\n")
        f.write(f"数据文件: {h5ad_path}\n")
        f.write(f"总细胞数: {adata.n_obs}\n")
        f.write(f"目标样本: {', '.join(target_samples)}\n")
        f.write(f"目标样本细胞数: {len(group1_data)}\n")
        f.write(f"其他样本Adenoma细胞数: {len(group2_data)}\n")
        f.write("\n详细比例:\n")
        f.write(all_results.to_string())
    
    print(f"\n分析完成！所有结果保存在: {output_dir}")
    return all_results

# 主程序
if __name__ == "__main__":
    # 配置文件路径
    h5ad_path = "adata_celltype_stage.h5ad"  # 替换为您的h5ad文件路径
    output_dir = "cell_composition_results"
    
    try:
        # 运行分析
        results = analyze_cell_composition(h5ad_path, output_dir)
        
    except Exception as e:
        print(f"错误: {str(e)}")
        import traceback
        traceback.print_exc()
