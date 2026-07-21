import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

def calculate_sample_celltype_proportions(adata, celltypes_of_interest, sample_col='sample'):
    """
    计算每个样本中指定细胞类型的比例
    """
    print("正在计算样本维度的细胞类型比例...")
    
    # 检查必要的列是否存在
    required_cols = [sample_col, 'MPCpercent', 'annotation']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"在adata.obs中未找到'{col}'列")
    
    # 提取所需数据
    df = adata.obs[[sample_col, 'MPCpercent', 'annotation']].copy()
    
    # 获取数据中实际存在的细胞类型
    actual_celltypes = adata.obs['annotation'].unique()
    
    # 筛选出我们感兴趣的且实际存在的细胞类型
    valid_celltypes = [ct for ct in celltypes_of_interest if ct in actual_celltypes]
    missing_celltypes = [ct for ct in celltypes_of_interest if ct not in actual_celltypes]
    
    if missing_celltypes:
        print(f"警告: 以下细胞类型在数据中不存在: {missing_celltypes}")
    
    if not valid_celltypes:
        raise ValueError("没有找到任何有效的细胞类型进行分析")
    
    print(f"将分析的细胞类型 ({len(valid_celltypes)}个): {valid_celltypes}")
    
    # 筛选只包含我们感兴趣细胞类型的数据
    df_filtered = df[df['annotation'].isin(valid_celltypes)].copy()
    
    # 计算每个样本中每种细胞类型的细胞数
    sample_cell_counts = pd.crosstab(df_filtered[sample_col], df_filtered['annotation'])
    
    # 确保所有感兴趣的细胞类型都有列（即使计数为0）
    for ct in valid_celltypes:
        if ct not in sample_cell_counts.columns:
            sample_cell_counts[ct] = 0
    
    # 重新排列列的顺序
    sample_cell_counts = sample_cell_counts[valid_celltypes]
    
    # 计算每个样本的总细胞数（针对这些细胞类型）
    sample_totals = sample_cell_counts.sum(axis=1)
    
    # 计算比例（百分比）
    sample_proportions = (sample_cell_counts.div(sample_totals, axis=0) * 100).round(4)
    
    # 添加分组信息
    sample_groups = df_filtered[[sample_col, 'MPCpercent']].drop_duplicates().set_index(sample_col)
    sample_proportions['MPCpercent'] = sample_proportions.index.map(sample_groups['MPCpercent'])
    
    # 重置索引，将样本名作为列
    sample_proportions = sample_proportions.reset_index()
    
    return sample_proportions, valid_celltypes

def kruskal_wallis_test_by_celltype(sample_proportions, celltypes, sample_col='sample'):
    """
    对每个细胞类型进行Kruskal-Wallis检验
    """
    print("\n进行Kruskal-Wallis检验...")
    
    results = {}
    
    for celltype in celltypes:
        # 准备各组数据
        groups = sorted([g for g in sample_proportions['MPCpercent'].unique() 
                         if pd.notna(g)])
        group_data = []
        group_labels = []
        
        for group in groups:
            group_values = sample_proportions[sample_proportions['MPCpercent'] == group][celltype].values
            # 过滤掉NaN值
            group_values = group_values[~np.isnan(group_values)]
            if len(group_values) > 0:
                group_data.append(group_values)
                group_labels.append(group)
        
        # 如果至少有2个组有数据，进行Kruskal-Wallis检验
        if len(group_data) >= 2:
            try:
                stat, p_value = stats.kruskal(*group_data)
                test_used = 'Kruskal-Wallis'
            except Exception as e:
                print(f"  细胞类型 {celltype}: 检验失败 - {e}")
                p_value = np.nan
                test_used = '检验失败'
        else:
            p_value = np.nan
            test_used = '数据不足'
        
        # 判断是否显著 (p < 0.05)
        significant = p_value < 0.05 if not np.isnan(p_value) else False
        
        results[celltype] = {
            'test': test_used,
            'p_value': p_value,
            'significant': significant
        }
    
    return results

def create_trend_summary_plot(sample_proportions, celltypes, results,
                              sample_col='sample', output_dir='MPC_trend_summary'):
    """
    创建所有指定细胞类型的趋势汇总图
    """
    print(f"\n创建趋势汇总图...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 确定分组顺序 (T1到T5)
    mpc_groups = sorted([g for g in sample_proportions['MPCpercent'].unique() 
                         if isinstance(g, str) and g.startswith('T')])
    
    if not mpc_groups:
        mpc_groups = sorted([g for g in sample_proportions['MPCpercent'].unique() 
                             if pd.notna(g)])
    
    print(f"分组顺序: {mpc_groups}")
    
    # 为每个分组创建数值映射
    group_to_num = {group: i for i, group in enumerate(mpc_groups)}
    
    # 创建图形
    n_celltypes = len(celltypes)
    n_cols = 4
    n_rows = (n_celltypes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*6, n_rows*5))
    
    # 如果只有一行，确保axes是二维数组
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # 设置颜色调色板
    colors = plt.cm.Set2(np.linspace(0, 1, len(celltypes)))
    
    # 遍历每个细胞类型
    for idx, celltype in enumerate(celltypes):
        row = idx // n_cols
        col = idx % n_cols
        
        ax = axes[row, col]
        
        # 准备绘图数据
        plot_data = []
        for group in mpc_groups:
            group_data = sample_proportions[sample_proportions['MPCpercent'] == group]
            values = group_data[celltype].dropna().values
            
            for value in values:
                plot_data.append({
                    'MPC_group_num': group_to_num[group],
                    'MPC_group': group,
                    'Proportion': value,
                    'CellType': celltype
                })
        
        plot_df = pd.DataFrame(plot_data)
        
        if plot_df.empty:
            ax.text(0.5, 0.5, f'No data for\n{celltype}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
            # 根据显著性添加*号
            title_suffix = "*" if results.get(celltype, {}).get('significant', False) else ""
            ax.set_title(f'{celltype}{title_suffix}', fontsize=12, fontweight='bold')
            continue
        
        # 计算每组的统计量
        stats_by_group = plot_df.groupby('MPC_group')['Proportion'].agg([
            'median', 'quantile', 'std', 'count'
        ]).rename(columns={'quantile': 'q1'})
        
        # 添加Q1和Q3
        q1 = plot_df.groupby('MPC_group')['Proportion'].quantile(0.25)
        q3 = plot_df.groupby('MPC_group')['Proportion'].quantile(0.75)
        stats_by_group['q1'] = q1
        stats_by_group['q3'] = q3
        
        # 重置索引
        stats_by_group = stats_by_group.reset_index()
        
        # 确保按照MPC分组顺序排序
        stats_by_group['group_num'] = stats_by_group['MPC_group'].map(group_to_num)
        stats_by_group = stats_by_group.sort_values('group_num')
        
        # 绘制简化箱线图 - 只显示上下线
        for i, row_data in stats_by_group.iterrows():
            group_num = row_data['group_num']
            q1_val = row_data['q1']
            q3_val = row_data['q3']
            
            # 绘制上下四分位线
            ax.hlines(y=q1_val, xmin=group_num-0.3, xmax=group_num+0.3, 
                     color='black', linewidth=1.5)
            ax.hlines(y=q3_val, xmin=group_num-0.3, xmax=group_num+0.3, 
                     color='black', linewidth=1.5)
            
            # 绘制中间的竖线连接上下线
            ax.vlines(x=group_num, ymin=q1_val, ymax=q3_val, 
                     color='black', linewidth=1.5)
        
        # 根据显著性设置折线样式
        is_significant = results.get(celltype, {}).get('significant', False)
        
        if is_significant:
            # 显著差异：红色实线
            line_style = '-'
            line_color = 'red'
            line_width = 3
        else:
            # 不显著：灰色虚线
            line_style = '--'
            line_color = 'gray'
            line_width = 2
        
        # 绘制中位数趋势线
        ax.plot(stats_by_group['group_num'], stats_by_group['median'], 
                line_style, linewidth=line_width, color=line_color, alpha=0.8, label='Median Trend')
        
        # 根据显著性添加*号
        title_suffix = "*" if is_significant else ""
        
        # 设置标题
        ax.set_title(f'{celltype}{title_suffix}', fontsize=12, fontweight='bold')
        
        # 设置x轴
        ax.set_xlabel('MPC Group', fontsize=10)
        ax.set_xticks(range(len(mpc_groups)))
        ax.set_xticklabels(mpc_groups)
        
        # 设置y轴标签
        if col == 0:
            ax.set_ylabel('Proportion (%)', fontsize=10)
        else:
            ax.set_ylabel('')
        
        # 添加网格
        ax.grid(True, axis='y', alpha=0.3, linestyle='--')
        
        # 统一设置y轴范围
        ax.set_ylim(0, 30)
    
    # 隐藏多余的子图
    for idx in range(len(celltypes), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)
    
    # 添加总标题
    plt.suptitle('Selected Cell Types by MPC Group (Simplified Boxplot with Median Trend)', 
                fontsize=16, fontweight='bold', y=1.02)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图形
    output_path = os.path.join(output_dir, "selected_celltypes_trend_summary.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"趋势汇总图已保存到: {output_path}")
    
    # 创建一个更简洁的版本，只包含趋势线
    create_trend_lines_only(sample_proportions, celltypes, mpc_groups, group_to_num, results, output_dir)

def create_trend_lines_only(sample_proportions, celltypes, mpc_groups, group_to_num, results, output_dir):
    """
    创建只包含趋势线的图形
    """
    print(f"创建趋势线图...")
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 遍历每个细胞类型绘制趋势线
    for idx, celltype in enumerate(celltypes):
        # 准备数据
        plot_data = []
        for group in mpc_groups:
            group_data = sample_proportions[sample_proportions['MPCpercent'] == group]
            values = group_data[celltype].dropna().values
            
            for value in values:
                plot_data.append({
                    'MPC_group_num': group_to_num[group],
                    'MPC_group': group,
                    'Proportion': value
                })
        
        plot_df = pd.DataFrame(plot_data)
        
        if plot_df.empty:
            continue
        
        # 计算每组的中位数
        median_by_group = plot_df.groupby('MPC_group_num')['Proportion'].median().sort_index()
        
        # 根据显著性设置折线样式
        is_significant = results.get(celltype, {}).get('significant', False)
        
        if is_significant:
            # 显著差异：红色实线
            line_style = '-'
            line_color = 'red'
            line_width = 2.5
            label = f"{celltype}*"
        else:
            # 不显著：灰色虚线
            line_style = '--'
            line_color = 'gray'
            line_width = 2
            label = celltype
        
        # 绘制趋势线
        ax.plot(median_by_group.index, median_by_group.values, 
                line_style, linewidth=line_width, color=line_color, 
                label=label, alpha=0.8)
    
    # 设置图表属性
    ax.set_xlabel('MPC Group', fontsize=12)
    ax.set_ylabel('Proportion (%)', fontsize=12)
    ax.set_title('Trend Lines of Selected Cell Types by MPC Group', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(mpc_groups)))
    ax.set_xticklabels(mpc_groups)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # 统一设置y轴范围
    ax.set_ylim(0, 30)
    
    # 添加图例
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图形
    output_path = os.path.join(output_dir, "selected_celltypes_trend_lines.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"趋势线图已保存到: {output_path}")

def main():
    # 配置参数
    h5ad_file = "adata_filtered_with_MPCpercent.h5ad"
    sample_col = 'sample'  # 如果样本列名不同，请修改此处
    
    # 感兴趣的细胞类型列表
    celltypes_of_interest = [
        "CD8+_Tem", "M2", "MKI67+_T", "B_follicular", "CD4+_Tex",
        "LYVE1+_TAM", "CD8+_Tcm", "NK", "IFN-T", "B_germinal_center",
        "EMT_M2", "Treg", "CD8+_Tex", "CD4+_Tem", "B_memory", "TAM"
    ]
    
    try:
        # 1. 读取数据
        print("="*60)
        print("读取h5ad数据...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        print(f"细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
        
        # 检查数据
        print(f"\n数据列: {adata.obs.columns.tolist()}")
        print(f"MPCpercent分组分布:")
        print(adata.obs['MPCpercent'].value_counts().sort_index())
        print(f"细胞类型数量: {adata.obs['annotation'].nunique()}")
        print(f"样本数量: {adata.obs[sample_col].nunique()}")
        
        # 2. 计算样本维度的细胞类型比例
        sample_proportions, valid_celltypes = calculate_sample_celltype_proportions(
            adata, celltypes_of_interest, sample_col
        )
        
        # 保存比例数据
        sample_proportions.to_csv("selected_celltype_proportions.csv", index=False, encoding='utf-8-sig')
        print(f"\n样本细胞类型比例数据已保存到: selected_celltype_proportions.csv")
        
        # 3. 进行Kruskal-Wallis检验
        results = kruskal_wallis_test_by_celltype(sample_proportions, valid_celltypes, sample_col)
        
        # 提取有显著差异的细胞类型
        significant_celltypes = [celltype for celltype, result in results.items() 
                                if result.get('significant', False)]
        
        # 保存统计结果
        results_df = pd.DataFrame.from_dict(results, orient='index')
        results_df.index.name = 'celltype'
        results_df.to_csv("selected_kruskal_wallis_results.csv", encoding='utf-8-sig')
        print(f"统计检验结果已保存到: selected_kruskal_wallis_results.csv")
        
        # 4. 打印结果摘要
        print("\n" + "="*60)
        print("分析结果摘要")
        print("="*60)
        print(f"总细胞类型数: {len(valid_celltypes)}")
        print(f"有显著差异的细胞类型数: {len(significant_celltypes)}")
        
        if significant_celltypes:
            print(f"\n有显著差异的细胞类型 (p < 0.05):")
            for celltype in significant_celltypes:
                p_val = results[celltype].get('p_value', np.nan)
                print(f"  {celltype}: p = {p_val:.4e}")
        else:
            print("\n没有发现显著差异的细胞类型 (p < 0.05)")
        
        # 5. 绘制趋势汇总图
        output_dir = "MPC_trend_summary_selected"
        create_trend_summary_plot(sample_proportions, valid_celltypes, results, sample_col, output_dir)
        
        # 6. 打印保存的文件列表
        print(f"\n生成的文件:")
        print(f"  1. selected_celltype_proportions.csv - 样本细胞类型比例数据")
        print(f"  2. selected_kruskal_wallis_results.csv - 统计检验结果")
        print(f"  3. {output_dir}/ - 包含趋势汇总图的目录")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
