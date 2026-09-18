import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import os

def analyze_celltype_sample_stats(adata, celltype_col='celltype', sample_col='sample'):
    """
    对h5ad数据中的celltype和sample进行统计分析
    
    Parameters:
    -----------
    adata : AnnData对象
        包含celltype和sample信息的单细胞数据
    celltype_col : str
        celltype所在的列名
    sample_col : str
        sample所在的列名
    """
    
    # 基础信息统计
    print("=== 数据基本信息 ===")
    print(f"总细胞数: {adata.n_obs}")
    print(f"总样本数: {adata.obs[sample_col].nunique()}")
    print(f"总细胞类型数: {adata.obs[celltype_col].nunique()}")
    
    # 样本和细胞类型的交叉统计
    print("\n=== 样本-细胞类型交叉统计 ===")
    cross_table = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col])
    print(cross_table)
    
    # 计算每个样本的细胞类型比例
    print("\n=== 各样本细胞类型比例 ===")
    sample_proportions = cross_table.div(cross_table.sum(axis=1), axis=0)
    print(sample_proportions.round(3))
    
    return cross_table, sample_proportions

def count_all_celltypes(adata, celltype_col='celltype'):
    """
    统计所有细胞类型的总数和详细信息
    
    Parameters:
    -----------
    adata : AnnData对象
        包含celltype信息的单细胞数据
    celltype_col : str
        celltype所在的列名
        
    Returns:
    --------
    celltype_stats : pandas DataFrame
        包含每种细胞类型统计信息的数据框
    """
    print("\n=== 所有细胞类型的总数统计 ===")
    
    # 统计每种细胞类型的数量
    celltype_counts = adata.obs[celltype_col].value_counts()
    
    # 创建详细统计信息的数据框
    celltype_stats = pd.DataFrame({
        '细胞类型': celltype_counts.index,
        '细胞数量': celltype_counts.values,
        '占比(%)': (celltype_counts.values / celltype_counts.sum() * 100).round(2)
    })
    
    # 添加排名
    celltype_stats['排名'] = range(1, len(celltype_stats) + 1)
    
    # 重新排列列
    celltype_stats = celltype_stats[['排名', '细胞类型', '细胞数量', '占比(%)']]
    
    # 打印详细统计
    print("细胞类型详细统计:")
    print(celltype_stats.to_string(index=False))
    
    # 打印汇总信息
    print(f"\n汇总信息:")
    print(f"1. 细胞类型总数: {len(celltype_counts)} 种")
    print(f"2. 总细胞数量: {celltype_counts.sum()} 个")
    print(f"3. 平均每种细胞类型的数量: {celltype_counts.mean():.1f} 个")
    print(f"4. 细胞数量最多的类型: {celltype_counts.idxmax()} ({celltype_counts.max()} 个, 占比{(celltype_counts.max()/celltype_counts.sum()*100):.1f}%)")
    print(f"5. 细胞数量最少的类型: {celltype_counts.idxmin()} ({celltype_counts.min()} 个, 占比{(celltype_counts.min()/celltype_counts.sum()*100):.1f}%)")
    
    return celltype_stats

def create_visualizations(adata, cross_table, sample_proportions, celltype_stats, 
                          celltype_col='celltype', sample_col='sample'):
    """
    创建统计可视化图表
    """
    # 设置图形样式
    plt.style.use('default')
    fig = plt.figure(figsize=(18, 15))
    
    # 1. 样本细胞数量分布
    ax1 = plt.subplot(3, 3, 1)
    sample_counts = adata.obs[sample_col].value_counts()
    colors1 = plt.cm.Set3(np.linspace(0, 1, len(sample_counts)))
    bars1 = ax1.bar(range(len(sample_counts)), sample_counts.values, color=colors1)
    ax1.set_title('各样本细胞数量分布', fontsize=12, fontweight='bold')
    ax1.set_xlabel('样本')
    ax1.set_ylabel('细胞数量')
    ax1.set_xticks(range(len(sample_counts)))
    ax1.set_xticklabels(sample_counts.index.astype(str), rotation=45, ha='right')
    
    # 在柱子上添加数值
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01*max(sample_counts.values),
                f'{int(height)}', ha='center', va='bottom', fontsize=9)
    
    # 2. 细胞类型分布（前20种，如果多于20种）
    ax2 = plt.subplot(3, 3, 2)
    celltype_counts = adata.obs[celltype_col].value_counts()
    top_n = min(20, len(celltype_counts))
    top_celltypes = celltype_counts.head(top_n)
    
    if len(celltype_counts) > 20:
        other_counts = celltype_counts[20:].sum()
        top_celltypes = pd.concat([top_celltypes, pd.Series([other_counts], index=['其他'])])
    
    colors2 = plt.cm.tab20c(np.linspace(0, 1, len(top_celltypes)))
    wedges, texts, autotexts = ax2.pie(top_celltypes.values, labels=top_celltypes.index, 
                                       autopct='%1.1f%%', colors=colors2, startangle=90)
    ax2.set_title(f'细胞类型分布 (前{top_n}种)', fontsize=12, fontweight='bold')
    
    # 3. 样本细胞类型组成热图
    ax3 = plt.subplot(3, 3, 3)
    # 限制热图显示的细胞类型数量，避免过于拥挤
    max_celltypes_heatmap = min(15, len(cross_table.columns))
    top_celltypes_heatmap = celltype_counts.head(max_celltypes_heatmap).index
    heatmap_data = sample_proportions[top_celltypes_heatmap].T
    
    im = ax3.imshow(heatmap_data, aspect='auto', cmap='YlOrRd')
    ax3.set_title(f'各样本细胞类型组成热图 (前{max_celltypes_heatmap}种)', fontsize=12, fontweight='bold')
    ax3.set_xlabel('样本')
    ax3.set_ylabel('细胞类型')
    ax3.set_xticks(range(len(heatmap_data.columns)))
    ax3.set_xticklabels(heatmap_data.columns, rotation=45, ha='right')
    ax3.set_yticks(range(len(heatmap_data.index)))
    ax3.set_yticklabels(heatmap_data.index)
    plt.colorbar(im, ax=ax3, label='比例')
    
    # 4. 堆叠柱状图显示细胞类型比例
    ax4 = plt.subplot(3, 3, 4)
    # 限制显示的细胞类型数量
    max_celltypes_bar = min(10, len(cross_table.columns))
    top_celltypes_bar = celltype_counts.head(max_celltypes_bar).index
    
    if len(cross_table.columns) > max_celltypes_bar:
        other_data = sample_proportions.drop(columns=top_celltypes_bar).sum(axis=1)
        plot_data = pd.concat([sample_proportions[top_celltypes_bar], 
                              other_data.rename('其他')], axis=1)
    else:
        plot_data = sample_proportions
    
    plot_data.plot(kind='bar', stacked=True, ax=ax4, colormap='tab20')
    ax4.set_title(f'各样本细胞类型比例 (前{max_celltypes_bar}种)', fontsize=12, fontweight='bold')
    ax4.set_xlabel('样本')
    ax4.set_ylabel('比例')
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax4.tick_params(axis='x', rotation=45)
    
    # 5. 细胞类型数量条形图（前20种）
    ax5 = plt.subplot(3, 3, 5)
    top_n_bar = min(20, len(celltype_counts))
    top_celltypes_bar = celltype_counts.head(top_n_bar)
    
    colors5 = plt.cm.viridis(np.linspace(0.2, 0.8, top_n_bar))
    bars5 = ax5.barh(range(top_n_bar), top_celltypes_bar.values[::-1], color=colors5)
    ax5.set_title(f'各细胞类型数量分布 (前{top_n_bar}种)', fontsize=12, fontweight='bold')
    ax5.set_xlabel('细胞数量')
    ax5.set_yticks(range(top_n_bar))
    ax5.set_yticklabels(top_celltypes_bar.index[::-1])
    ax5.invert_yaxis()  # 让最大的在顶部
    
    # 在条形上添加数值
    for i, bar in enumerate(bars5):
        width = bar.get_width()
        ax5.text(width + 0.01*max(top_celltypes_bar.values), bar.get_y() + bar.get_height()/2.,
                f'{int(width)}', ha='left', va='center', fontsize=9)
    
    # 6. 细胞类型占比饼图
    ax6 = plt.subplot(3, 3, 6)
    # 为饼图选择前8种
    top_n_pie = min(8, len(celltype_counts))
    top_for_pie = celltype_counts.head(top_n_pie)
    
    if len(celltype_counts) > top_n_pie:
        other_pie = celltype_counts[top_n_pie:].sum()
        top_for_pie = pd.concat([top_for_pie, pd.Series([other_pie], index=['其他'])])
    
    colors6 = plt.cm.Set3(np.linspace(0, 1, len(top_for_pie)))
    wedges, texts, autotexts = ax6.pie(top_for_pie.values, labels=top_for_pie.index, 
                                       autopct='%1.1f%%', colors=colors6, startangle=90)
    ax6.set_title(f'细胞类型占比 (前{top_n_pie}种)', fontsize=12, fontweight='bold')
    
    # 7. 细胞类型数量分布曲线
    ax7 = plt.subplot(3, 3, 7)
    sorted_counts = celltype_counts.sort_values(ascending=False).reset_index(drop=True)
    ax7.plot(range(1, len(sorted_counts)+1), sorted_counts.values, 'b-', linewidth=2, marker='o', markersize=4)
    ax7.set_title('细胞类型数量分布曲线', fontsize=12, fontweight='bold')
    ax7.set_xlabel('细胞类型排名')
    ax7.set_ylabel('细胞数量')
    ax7.set_xscale('log')
    ax7.set_yscale('log')
    ax7.grid(True, alpha=0.3)
    
    # 8. 样本数量 vs 细胞类型数量散点图
    ax8 = plt.subplot(3, 3, 8)
    sample_celltype_counts = cross_table.astype(bool).sum(axis=1)  # 每个样本有多少种细胞类型
    sample_total_counts = cross_table.sum(axis=1)  # 每个样本的总细胞数
    
    scatter = ax8.scatter(sample_total_counts, sample_celltype_counts, 
                         c=range(len(sample_total_counts)), cmap='viridis', s=100, alpha=0.7)
    ax8.set_title('样本细胞数量 vs 细胞类型数量', fontsize=12, fontweight='bold')
    ax8.set_xlabel('样本总细胞数')
    ax8.set_ylabel('样本细胞类型数')
    
    # 为每个点添加样本标签
    for i, (x, y) in enumerate(zip(sample_total_counts, sample_celltype_counts)):
        ax8.text(x, y, str(sample_total_counts.index[i]), fontsize=8, 
                ha='center', va='bottom')
    
    # 9. 细胞类型统计摘要
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    # 准备摘要文本
    summary_text = f"细胞类型统计摘要\n\n"
    summary_text += f"总细胞类型数: {len(celltype_counts)}\n"
    summary_text += f"总细胞数量: {celltype_counts.sum()}\n"
    summary_text += f"平均每类型: {celltype_counts.mean():.1f}\n"
    summary_text += f"中位数: {celltype_counts.median():.1f}\n"
    summary_text += f"标准差: {celltype_counts.std():.1f}\n\n"
    
    summary_text += f"最多的类型:\n"
    summary_text += f"  {celltype_counts.index[0]}: {celltype_counts.values[0]} ({celltype_stats.iloc[0, 3]}%)\n\n"
    
    summary_text += f"最少的类型:\n"
    summary_text += f"  {celltype_counts.index[-1]}: {celltype_counts.values[-1]} ({celltype_stats.iloc[-1, 3]}%)\n\n"
    
    summary_text += f"前5种类型占比:\n"
    top5_percent = celltype_counts.head(5).sum() / celltype_counts.sum() * 100
    summary_text += f"  {top5_percent:.1f}%"
    
    ax9.text(0.1, 0.95, summary_text, transform=ax9.transAxes, 
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.show()

def advanced_statistical_tests(cross_table):
    """
    进行高级统计检验
    """
    print("\n=== 高级统计分析 ===")
    
    # 卡方检验 - 检验样本和细胞类型是否独立
    chi2, p_value, dof, expected = chi2_contingency(cross_table)
    print(f"卡方检验结果:")
    print(f"卡方值: {chi2:.3f}")
    print(f"P值: {p_value:.3e}")
    print(f"自由度: {dof}")
    
    if p_value < 0.05:
        print("结论: 样本和细胞类型分布不独立 (p < 0.05)")
    else:
        print("结论: 样本和细胞类型分布独立 (p ≥ 0.05)")
    
    # 计算每个细胞类型的分布均匀性
    print("\n=== 细胞类型分布均匀性 ===")
    celltype_shannon = []
    for celltype in cross_table.columns:
        proportions = cross_table[celltype] / cross_table[celltype].sum()
        # 计算香农熵（均匀性度量）
        entropy = -np.sum(proportions * np.log(proportions + 1e-10))
        max_entropy = np.log(len(proportions))
        uniformity = entropy / max_entropy if max_entropy > 0 else 0
        celltype_shannon.append(uniformity)
        print(f"{celltype}: 均匀性 = {uniformity:.3f}")
    
    # 计算整体均匀性
    overall_uniformity = np.mean(celltype_shannon)
    print(f"\n整体均匀性: {overall_uniformity:.3f}")
    
    return chi2, p_value, celltype_shannon

def export_statistical_results(cross_table, sample_proportions, celltype_stats, output_prefix='celltype_sample_stats'):
    """
    导出统计结果到文件
    """
    # 创建输出目录
    output_dir = 'statistical_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存交叉表
    cross_table.to_csv(os.path.join(output_dir, f'{output_prefix}_cross_table.csv'))
    
    # 保存比例表
    sample_proportions.to_csv(os.path.join(output_dir, f'{output_prefix}_proportions.csv'))
    
    # 保存细胞类型统计
    celltype_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_celltype_stats.csv'), index=False)
    
    # 保存细胞类型详细数量
    celltype_counts = pd.Series(celltype_stats['细胞数量'].values, index=celltype_stats['细胞类型'])
    celltype_counts.to_csv(os.path.join(output_dir, f'{output_prefix}_celltype_counts.csv'), header=['细胞数量'])
    
    # 保存汇总统计
    summary_stats = pd.DataFrame({
        '指标': ['总细胞数', '总样本数', '总细胞类型数', '平均每样本细胞数', 
                '平均每细胞类型数', '最多细胞类型', '最多细胞类型数量', 
                '最少细胞类型', '最少细胞类型数量', '细胞类型均匀性'],
        '数值': [
            celltype_counts.sum(),
            len(cross_table),
            len(cross_table.columns),
            cross_table.sum().sum() / len(cross_table),
            celltype_counts.mean(),
            celltype_counts.idxmax(),
            celltype_counts.max(),
            celltype_counts.idxmin(),
            celltype_counts.min(),
            f"{celltype_counts.std()/celltype_counts.mean():.3f}"
        ]
    })
    summary_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_summary.csv'), index=False)
    
    # 保存为Excel文件（包含多个工作表）
    with pd.ExcelWriter(os.path.join(output_dir, f'{output_prefix}_all_results.xlsx')) as writer:
        cross_table.to_excel(writer, sheet_name='交叉表')
        sample_proportions.to_excel(writer, sheet_name='比例表')
        celltype_stats.to_excel(writer, sheet_name='细胞类型统计')
        summary_stats.to_excel(writer, sheet_name='汇总统计')
    
    print(f"\n统计结果已导出到目录: {output_dir}/")
    print(f"1. {output_prefix}_cross_table.csv - 样本-细胞类型交叉表")
    print(f"2. {output_prefix}_proportions.csv - 样本细胞类型比例表")
    print(f"3. {output_prefix}_celltype_stats.csv - 细胞类型统计表")
    print(f"4. {output_prefix}_celltype_counts.csv - 细胞类型数量表")
    print(f"5. {output_prefix}_summary.csv - 汇总统计表")
    print(f"6. {output_prefix}_all_results.xlsx - 所有结果的Excel文件")

# 主执行函数
def main():
    """
    主执行函数
    """
    # 读取h5ad文件
    # 请将下面的文件路径替换为你的实际文件路径
    h5ad_file = 'merged_single_cell_data.h5ad'  # 修改为你的文件路径
    
    try:
        print("正在读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"成功读取文件: {h5ad_file}")
    except FileNotFoundError:
        print(f"错误: 找不到文件 {h5ad_file}")
        print("请确保文件路径正确，并修改main()函数中的h5ad_file变量")
        return
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return
    
    # 检查必要的列是否存在
    required_cols = ['celltype', 'sample']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"错误: 缺少必要的列: {missing_cols}")
        print(f"可用的列: {list(adata.obs.columns)}")
        
        # 尝试猜测正确的列名
        celltype_candidates = [col for col in adata.obs.columns if 'cell' in col.lower() or 'type' in col.lower() or 'cluster' in col.lower()]
        sample_candidates = [col for col in adata.obs.columns if 'sample' in col.lower() or 'patient' in col.lower() or 'donor' in col.lower()]
        
        if celltype_candidates:
            print(f"\n可能的celltype列: {celltype_candidates}")
        if sample_candidates:
            print(f"可能的sample列: {sample_candidates}")
        
        # 如果找到候选列，询问用户是否使用
        if celltype_candidates and sample_candidates:
            use_default = input("\n是否使用第一个候选列? (y/n): ").lower()
            if use_default == 'y':
                celltype_col = celltype_candidates[0]
                sample_col = sample_candidates[0]
                print(f"使用: celltype_col='{celltype_col}', sample_col='{sample_col}'")
            else:
                return
        else:
            return
    else:
        celltype_col = 'celltype'
        sample_col = 'sample'
    
    # 检查数据质量
    print(f"\n=== 数据质量检查 ===")
    print(f"总细胞数: {adata.n_obs}")
    print(f"总基因数: {adata.n_vars}")
    
    # 检查是否有缺失值
    missing_celltype = adata.obs[celltype_col].isnull().sum()
    missing_sample = adata.obs[sample_col].isnull().sum()
    
    if missing_celltype > 0 or missing_sample > 0:
        print(f"警告: 发现缺失值")
        print(f"  - celltype列缺失: {missing_celltype} ({missing_celltype/adata.n_obs*100:.1f}%)")
        print(f"  - sample列缺失: {missing_sample} ({missing_sample/adata.n_obs*100:.1f}%)")
        
        # 移除缺失值
        original_size = adata.n_obs
        adata = adata[~adata.obs[celltype_col].isnull() & ~adata.obs[sample_col].isnull()].copy()
        print(f"已移除包含缺失值的细胞，剩余细胞数: {adata.n_obs} (移除了 {original_size - adata.n_obs} 个细胞)")
    
    # 执行统计分析
    cross_table, sample_proportions = analyze_celltype_sample_stats(adata, celltype_col, sample_col)
    
    # 统计所有细胞类型的总数
    celltype_stats = count_all_celltypes(adata, celltype_col)
    
    # 创建可视化
    print("\n正在生成可视化图表...")
    create_visualizations(adata, cross_table, sample_proportions, celltype_stats, celltype_col, sample_col)
    
    # 高级统计检验
    chi2, p_value, uniformity_scores = advanced_statistical_tests(cross_table)
    
    # 导出结果
    output_prefix = input("\n请输入输出文件的前缀 (直接回车使用默认前缀 'celltype_sample_stats'): ").strip()
    if not output_prefix:
        output_prefix = 'celltype_sample_stats'
    
    export_statistical_results(cross_table, sample_proportions, celltype_stats, output_prefix)
    
    print("\n" + "="*50)
    print("分析完成!")
    print("="*50)
    print(f"总计分析了 {adata.n_obs} 个细胞")
    print(f"来自 {adata.obs[sample_col].nunique()} 个样本")
    print(f"共 {adata.obs[celltype_col].nunique()} 种细胞类型")
    print("="*50)

# 如果直接运行此脚本，执行主函数
if __name__ == "__main__":
    main()
