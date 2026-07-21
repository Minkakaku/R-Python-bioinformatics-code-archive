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

# 设置中文字体（如果需要显示中文标签）
# mpl.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
# mpl.rcParams['axes.unicode_minus'] = False

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
        if len(groups) > 2 and p_value < 0.05:  # 只在整体显著时进行两两比较
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

def create_boxplots_with_stats(proportions_long, celltypes, groups, results, output_dir="results"):
    """
    创建多系列箱线图，包含组间差异连线和P值
    """
    print("\n正在生成多系列箱线图...")
    
    # 确定子图的行列数
    n_celltypes = len(celltypes)
    
    # 如果细胞类型太多，分多页显示
    max_per_page = 12
    n_pages = (n_celltypes + max_per_page - 1) // max_per_page
    
    for page in range(n_pages):
        start_idx = page * max_per_page
        end_idx = min((page + 1) * max_per_page, n_celltypes)
        page_celltypes = celltypes[start_idx:end_idx]
        
        # 计算子图的行列
        n_subplots = len(page_celltypes)
        n_cols = min(3, n_subplots)
        n_rows = (n_subplots + n_cols - 1) // n_cols
        
        # 创建图形
        fig_width = n_cols * 6
        fig_height = n_rows * 5
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
        
        # 展平axes数组
        if n_rows * n_cols > 1:
            axes_flat = axes.ravel()
        else:
            axes_flat = [axes]
        
        for idx, celltype in enumerate(page_celltypes):
            ax = axes_flat[idx]
            
            # 提取该细胞类型的数据
            celltype_data = proportions_long[proportions_long['celltype'] == celltype]
            
            if celltype_data.empty:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
                ax.set_title(f'{celltype}\n(No data)')
                continue
            
            # 绘制箱线图
            plot_data = []
            group_labels = []
            
            for group in groups:
                group_values = celltype_data[celltype_data['group'] == group]['proportion'].values
                if len(group_values) > 0:
                    plot_data.append(group_values)
                    group_labels.append(group)
            
            if not plot_data:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
                ax.set_title(f'{celltype}\n(No data)')
                continue
            
            # 绘制箱线图
            bp = ax.boxplot(plot_data, labels=group_labels, patch_artist=True)
            
            # 设置箱线图颜色
            colors = plt.cm.Set2(np.linspace(0, 1, len(plot_data)))
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            # 添加散点显示个体数据点
            for i, group_vals in enumerate(plot_data):
                x_pos = np.random.normal(i + 1, 0.04, size=len(group_vals))
                ax.scatter(x_pos, group_vals, alpha=0.6, color='black', s=20, edgecolor='white', linewidth=0.5)
            
            # 设置标题
            p_value = results[celltype]['p_value']
            significance = "***" if results[celltype]['significant'] else "ns"
            title = f'{celltype}\n(p={p_value:.4f} {significance})'
            ax.set_title(title, fontsize=12, fontweight='bold')
            
            # 设置y轴标签
            ax.set_ylabel('Proportion', fontsize=10)
            ax.set_ylim(bottom=0)
            
            # 添加组间连线和p值
            if len(groups) > 1 and results[celltype]['significant']:
                # 获取两两比较结果
                pairwise_results = results[celltype]['pairwise']
                
                if pairwise_results:
                    # 计算y轴范围
                    y_min, y_max = ax.get_ylim()
                    y_range = y_max - y_min
                    
                    # 为每个比较对分配不同的高度
                    n_comparisons = len(pairwise_results)
                    y_positions = np.linspace(y_max + 0.02*y_range, 
                                            y_max + 0.15*y_range, 
                                            n_comparisons)
                    
                    for comp_idx, (comparison, p_val) in enumerate(pairwise_results.items()):
                        if np.isnan(p_val):
                            continue
                        
                        # 解析比较对
                        group1, group2 = comparison.split('_vs_')
                        
                        # 获取组在x轴上的位置
                        if group1 in group_labels and group2 in group_labels:
                            x1 = group_labels.index(group1) + 1
                            x2 = group_labels.index(group2) + 1
                            
                            # 绘制连线
                            y = y_positions[comp_idx]
                            ax.plot([x1, x2], [y, y], color='black', linewidth=1.5, alpha=0.7)
                            
                            # 在连线中间添加p值
                            x_mid = (x1 + x2) / 2
                            p_text = f'p={p_val:.3f}' if p_val >= 0.001 else 'p<0.001'
                            ax.text(x_mid, y + 0.005*y_range, p_text, 
                                   ha='center', va='bottom', fontsize=8,
                                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
                    
                    # 调整y轴范围以容纳连线和p值
                    ax.set_ylim([y_min, y_max + 0.2*y_range])
            
            # 添加网格
            ax.grid(True, alpha=0.3, axis='y')
        
        # 隐藏多余的子图
        for idx in range(len(page_celltypes), len(axes_flat)):
            axes_flat[idx].set_visible(False)
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图形
        page_suffix = f"_page{page+1}" if n_pages > 1 else ""
        plt.savefig(f'{output_dir}/celltype_boxplots{page_suffix}.png', 
                    dpi=300, bbox_inches='tight')
        plt.savefig(f'{output_dir}/celltype_boxplots{page_suffix}.pdf', 
                    bbox_inches='tight')
        plt.close()
    
    print(f"已保存箱线图到: {output_dir}/celltype_boxplots*.png/pdf")

def save_results(proportions, proportions_long, results, output_dir="results"):
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
        
        # 5. 生成多系列箱线图
        create_boxplots_with_stats(proportions_long, celltypes, groups, results, output_dir)
        
        # 6. 保存结果
        save_results(proportions, proportions_long, results, output_dir)
        
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
