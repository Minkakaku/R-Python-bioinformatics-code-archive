import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import warnings
import os
warnings.filterwarnings('ignore')

# 设置绘图风格
plt.style.use('default')
sns.set_palette("Set2")

def load_and_prepare_data(h5ad_file):
    """
    加载h5ad文件并检查必要的列
    """
    print("正在加载数据...")
    adata = sc.read_h5ad(h5ad_file)
    
    # 检查必要的列
    required_cols = ['sample', 'group', 'celltype']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"缺少必要的列: {missing_cols}")
    
    print(f"数据形状: {adata.shape}")
    print(f"样本数量: {adata.obs['sample'].nunique()}")
    print(f"分组数量: {adata.obs['group'].nunique()}")
    print(f"细胞类型数量: {adata.obs['celltype'].nunique()}")
    print(f"\n分组信息:")
    print(adata.obs['group'].value_counts())
    
    return adata

def calculate_celltype_proportions(adata):
    """
    计算每个样本中每种细胞类型的比例
    """
    print("\n正在计算细胞类型比例...")
    
    # 获取所有唯一的样本和细胞类型
    samples = adata.obs['sample'].unique()
    celltypes = adata.obs['celltype'].unique()
    groups = adata.obs['group'].unique()
    
    # 计算每个样本的总细胞数
    sample_total_cells = adata.obs['sample'].value_counts()
    
    # 创建样本-细胞类型的交叉表
    celltype_counts = pd.crosstab(adata.obs['sample'], adata.obs['celltype'])
    
    # 计算比例
    proportions = celltype_counts.div(sample_total_cells, axis=0)
    
    # 添加分组信息
    sample_groups = adata.obs[['sample', 'group']].drop_duplicates().set_index('sample')
    proportions['group'] = sample_groups['group']
    
    # 重置索引，将sample变为列
    proportions = proportions.reset_index()
    
    # 将宽格式转换为长格式，便于后续分析
    proportions_long = proportions.melt(id_vars=['sample', 'group'], 
                                        value_vars=celltypes,
                                        var_name='celltype',
                                        value_name='proportion')
    
    # 删除缺失值
    proportions_long = proportions_long.dropna()
    
    print(f"比例矩阵形状: {proportions.shape}")
    print(f"长格式数据形状: {proportions_long.shape}")
    
    return proportions, proportions_long, celltypes, groups

def statistical_analysis(proportions_long, celltypes, groups):
    """
    进行组间差异分析
    """
    print("\n正在进行组间差异分析...")
    
    results = {}
    
    for celltype in celltypes:
        # 提取该细胞类型的数据
        celltype_data = proportions_long[proportions_long['celltype'] == celltype]
        
        # 检查是否有足够的数据
        if celltype_data.empty:
            results[celltype] = {'test': '无数据', 'p_value': np.nan, 'pairwise': {}}
            continue
        
        # 按分组提取数据
        group_data = {}
        for group in groups:
            group_values = celltype_data[celltype_data['group'] == group]['proportion'].values
            if len(group_values) > 0:
                group_data[group] = group_values
        
        # 检查是否有至少两个组有数据
        if len(group_data) < 2:
            results[celltype] = {'test': '数据不足', 'p_value': np.nan, 'pairwise': {}}
            continue
        
        # 进行统计检验
        if len(groups) == 2:
            # 两组比较：Mann-Whitney U检验
            group1, group2 = groups
            if group1 in group_data and group2 in group_data:
                try:
                    stat, p_value = mannwhitneyu(group_data[group1], group_data[group2], 
                                                 alternative='two-sided')
                    test_used = 'Mann-Whitney U'
                except:
                    p_value = np.nan
                    test_used = '检验失败'
            else:
                p_value = np.nan
                test_used = '数据不足'
        else:
            # 多组比较：Kruskal-Wallis检验
            try:
                all_group_data = []
                for group in groups:
                    if group in group_data and len(group_data[group]) > 0:
                        all_group_data.append(group_data[group])
                
                if len(all_group_data) >= 2:
                    stat, p_value = kruskal(*all_group_data)
                    test_used = 'Kruskal-Wallis'
                else:
                    p_value = np.nan
                    test_used = '数据不足'
            except:
                p_value = np.nan
                test_used = '检验失败'
        
        # 两两比较
        pairwise_results = {}
        if len(groups) > 2:  # 无论整体是否显著，都进行两两比较
            for group1, group2 in combinations(groups, 2):
                if group1 in group_data and group2 in group_data:
                    try:
                        stat_pair, p_pair = mannwhitneyu(group_data[group1], group_data[group2], 
                                                        alternative='two-sided')
                        pairwise_results[f"{group1}_vs_{group2}"] = p_pair
                    except:
                        pairwise_results[f"{group1}_vs_{group2}"] = np.nan
                else:
                    pairwise_results[f"{group1}_vs_{group2}"] = np.nan
        
        # 存储结果
        results[celltype] = {
            'test': test_used,
            'p_value': p_value,
            'pairwise': pairwise_results,
            'significant': p_value < 0.05 if not np.isnan(p_value) else False
        }
    
    return results

def create_stacked_barplot(proportions, celltypes, output_dir="results"):
    """
    创建每个样本的细胞比例堆积条形图
    """
    print("\n正在生成样本细胞比例堆积条形图...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备数据
    plot_data = proportions.copy()
    
    # 设置图形大小
    n_samples = len(plot_data)
    fig_width = max(12, n_samples * 0.4)  # 根据样本数调整宽度
    fig_height = 8
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # 准备堆积条形图数据
    bar_data = plot_data[celltypes].values.T
    sample_labels = plot_data['sample']
    
    # 绘制堆积条形图
    bottom = np.zeros(n_samples)
    colors = plt.cm.Set3(np.linspace(0, 1, len(celltypes)))
    
    for i, celltype in enumerate(celltypes):
        ax.bar(range(n_samples), plot_data[celltype].values, 
               bottom=bottom, label=celltype, color=colors[i])
        bottom += plot_data[celltype].values
    
    # 设置x轴标签
    ax.set_xticks(range(n_samples))
    ax.set_xticklabels(sample_labels, rotation=45, ha='right')
    
    # 设置y轴
    ax.set_ylabel('Cell Proportion', fontsize=12)
    ax.set_ylim([0, 1])
    
    # 添加标题
    ax.set_title('Cell Type Proportions by Sample', fontsize=14, fontweight='bold')
    
    # 添加图例
    ax.legend(title='Cell Types', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # 保存图形
    plt.savefig(f'{output_dir}/sample_celltype_proportions_stacked.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/sample_celltype_proportions_stacked.pdf', 
                bbox_inches='tight')
    plt.close()
    
    print(f"已保存堆积条形图到: {output_dir}/sample_celltype_proportions_stacked.png/pdf")

def create_grouped_barchart(proportions_long, celltypes, groups, results, output_dir="results"):
    """
    创建多系列柱状图，显示每个细胞类型在不同分组中的平均比例
    """
    print("\n正在生成多系列柱状图...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 计算每个细胞类型在每个组中的平均比例和标准误差
    summary_data = []
    for celltype in celltypes:
        celltype_data = proportions_long[proportions_long['celltype'] == celltype]
        for group in groups:
            group_data = celltype_data[celltype_data['group'] == group]['proportion']
            if len(group_data) > 0:
                mean_prop = group_data.mean()
                std_prop = group_data.std()
                sem_prop = group_data.sem()  # 标准误差
                summary_data.append({
                    'celltype': celltype,
                    'group': group,
                    'mean': mean_prop,
                    'std': std_prop,
                    'sem': sem_prop,
                    'n': len(group_data)
                })
    
    summary_df = pd.DataFrame(summary_data)
    
    # 设置图形大小
    n_celltypes = len(celltypes)
    # 增加图形宽度
    fig_width = 18  # 固定宽度为18
    fig_height = 10
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # 设置柱状图的位置
    x_pos = np.arange(len(celltypes))
    # 增加细胞类型之间的间距
    x_pos = x_pos * 3  # 将间距扩大
    
    # 增加柱体宽度
    bar_width = 2 / len(groups)  # 增加柱体宽度
    
    # 为每个组选择颜色
    colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))
    
    # 存储每个柱体的y坐标，用于连线
    bar_y_positions = {}
    
    # 绘制柱状图
    for i, group in enumerate(groups):
        group_data = summary_df[summary_df['group'] == group]
        
        # 确保数据顺序与celltypes一致
        group_means = []
        group_sems = []
        for celltype in celltypes:
            cell_data = group_data[group_data['celltype'] == celltype]
            if len(cell_data) > 0:
                group_means.append(cell_data['mean'].values[0])
                group_sems.append(cell_data['sem'].values[0])
            else:
                group_means.append(0)
                group_sems.append(0)
        
        # 计算每个组的位置偏移
        x_offset = (i - (len(groups) - 1) / 2) * bar_width
        
        # 绘制柱状图，设置误差棒颜色与柱体相同
        bars = ax.bar(x_pos + x_offset, group_means, bar_width, 
                      label=group, color=colors[i], alpha=0.8,
                      yerr=group_sems, capsize=3, 
                      error_kw={'elinewidth': 1, 'ecolor': colors[i]})  # 设置误差棒颜色
        
        # 存储每个柱体的y坐标
        for j, bar in enumerate(bars):
            bar_y_positions[(j, i)] = bar.get_height()
    
    # 设置x轴标签
    ax.set_xticks(x_pos)
    ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=10)
    
    # 设置y轴
    ax.set_ylabel('Mean Proportion', fontsize=12)
    ax.set_ylim(bottom=0)
    
    # 添加标题
    ax.set_title('Mean Cell Type Proportions by Group', fontsize=14, fontweight='bold')
    
    # 添加图例
    ax.legend(title='Groups', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
    
    # 为每个细胞类型添加组间比较
    for i, celltype in enumerate(celltypes):
        celltype_results = results.get(celltype, {})
        pairwise_results = celltype_results.get('pairwise', {})
        
        if pairwise_results:
            # 获取这个细胞类型的所有组数据
            cell_data = summary_df[summary_df['celltype'] == celltype]
            group_means = {}
            for group in groups:
                group_data = cell_data[cell_data['group'] == group]
                if len(group_data) > 0:
                    group_means[group] = group_data['mean'].values[0]
                else:
                    group_means[group] = 0
            
            # 计算该细胞类型中最高的柱子高度
            max_height = max(group_means.values()) if group_means else 0
            
            # 获取所有显著差异的比较对
            significant_comparisons = []
            for comparison, p_val in pairwise_results.items():
                if not np.isnan(p_val) and p_val < 0.05:
                    significant_comparisons.append((comparison, p_val))
            
            # 按p值排序（从小到大）
            significant_comparisons.sort(key=lambda x: x[1])
            
            # 为每个比较对添加连线和p值，使用不同的高度
            for comp_idx, (comparison, p_val) in enumerate(significant_comparisons):
                # 解析比较对
                group1, group2 = comparison.split('_vs_')
                
                # 检查这两个组是否有数据
                if group1 not in group_means or group2 not in group_means:
                    continue
                
                # 获取组在x轴上的位置
                if group1 in groups and group2 in groups:
                    group1_idx = list(groups).index(group1)
                    group2_idx = list(groups).index(group2)
                    
                    # 计算每个组在柱状图中的x位置
                    x1 = x_pos[i] + (group1_idx - (len(groups) - 1) / 2) * bar_width
                    x2 = x_pos[i] + (group2_idx - (len(groups) - 1) / 2) * bar_width
                    
                    # 计算连线的y位置（基于最高柱体高度 + 0.08 + 0.03 * 比较索引）
                    y_line = max_height + 0.05 + (comp_idx * 0.02)
                    
                    # 绘制水平连线
                    ax.plot([x1, x2], [y_line, y_line], color='black', linewidth=0.5)
                    
                    # 在连线的末端向下延伸垂直的线
                    vertical_length = 0.005  # 垂直线的长度
                    ax.plot([x1, x1], [y_line, y_line - vertical_length], color='black', linewidth=0.5)
                    ax.plot([x2, x2], [y_line, y_line - vertical_length], color='black', linewidth=0.5)
                    
                    # 添加p值文本
                    if p_val < 0.001:
                        p_text = '***'
                    elif p_val < 0.01:
                        p_text = '**'
                    elif p_val < 0.05:
                        p_text = '*'
                    else:
                        p_text = f'p={p_val:.3f}'
                    
                    # 在连线上方添加p值
                    ax.text((x1 + x2) / 2, y_line + 0.005, p_text, 
                           ha='center', va='bottom', fontsize=6, fontweight='bold')
    
    # 调整y轴范围以容纳显著性标记
    y_max_plot = ax.get_ylim()[1]
    ax.set_ylim(top=y_max_plot + 0.2)  # 增加上端空间
    
    plt.tight_layout()
    
    # 保存图形
    plt.savefig(f'{output_dir}/celltype_grouped_barchart.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/celltype_grouped_barchart.pdf', 
                bbox_inches='tight')
    plt.close()
    
    print(f"已保存多系列柱状图到: {output_dir}/celltype_grouped_barchart.png/pdf")
    
    # 返回汇总数据，便于后续使用
    return summary_df

def save_results(proportions, proportions_long, results, summary_df, output_dir="results"):
    """
    保存分析结果到Excel文件
    """
    print("\n正在保存分析结果...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备描述性统计
    desc_stats = proportions.groupby('group')[proportions.columns[:-2]].describe()
    
    # 准备差异分析结果
    results_df = pd.DataFrame.from_dict(results, orient='index')
    results_df.index.name = 'celltype'
    
    # 提取两两比较结果
    pairwise_list = []
    for celltype, res in results.items():
        for comparison, p_val in res['pairwise'].items():
            pairwise_list.append({
                'celltype': celltype,
                'comparison': comparison,
                'p_value': p_val
            })
    
    if pairwise_list:
        pairwise_df = pd.DataFrame(pairwise_list)
    else:
        pairwise_df = pd.DataFrame(columns=['celltype', 'comparison', 'p_value'])
    
    # 保存到Excel
    with pd.ExcelWriter(f'{output_dir}/celltype_analysis_results.xlsx') as writer:
        proportions.to_excel(writer, sheet_name='样本比例矩阵', index=False)
        proportions_long.to_excel(writer, sheet_name='长格式数据', index=False)
        desc_stats.to_excel(writer, sheet_name='描述性统计')
        results_df.to_excel(writer, sheet_name='差异分析结果')
        summary_df.to_excel(writer, sheet_name='汇总数据', index=False)
        if not pairwise_df.empty:
            pairwise_df.to_excel(writer, sheet_name='两两比较结果', index=False)
    
    print(f"已保存结果到: {output_dir}/celltype_analysis_results.xlsx")

def main():
    """
    主函数
    """
    # 配置参数
    h5ad_file = "celltype_samples.h5ad"  # 修改为您的h5ad文件路径
    output_dir = "results"
    
    try:
        # 1. 加载数据
        adata = load_and_prepare_data(h5ad_file)
        
        # 2. 计算细胞类型比例
        proportions, proportions_long, celltypes, groups = calculate_celltype_proportions(adata)
        
        # 3. 组间差异分析
        results = statistical_analysis(proportions_long, celltypes, groups)
        
        # 4. 生成样本细胞比例堆积条形图
        create_stacked_barplot(proportions, celltypes, output_dir)
        
        # 5. 生成多系列柱状图
        summary_df = create_grouped_barchart(proportions_long, celltypes, groups, results, output_dir)
        
        # 6. 保存结果
        save_results(proportions, proportions_long, results, summary_df, output_dir)
        
        # 7. 打印摘要
        print("\n" + "="*60)
        print("分析完成摘要")
        print("="*60)
        
        # 计算显著细胞类型数量
        significant_count = sum(1 for r in results.values() if r.get('significant', False))
        
        print(f"分析的细胞类型数量: {len(celltypes)}")
        print(f"分组数量: {len(groups)}")
        print(f"具有显著组间差异的细胞类型数量: {significant_count}")
        
        if significant_count > 0:
            print("\n具有显著差异的细胞类型:")
            for celltype, res in results.items():
                if res.get('significant', False):
                    p_val = res.get('p_value', np.nan)
                    test = res.get('test', 'N/A')
                    print(f"  {celltype}: {test} p={p_val:.4f}")
        
        print(f"\n所有结果已保存到 '{output_dir}' 目录")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
