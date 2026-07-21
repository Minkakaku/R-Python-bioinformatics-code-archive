import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import re
import matplotlib as mpl

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def analyze_celltype_sample_stats(adata, celltype_col='celltype', sample_col='sample'):
    """
    对h5ad数据中的celltype和sample进行统计分析
    """
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
    celltype_stats = celltype_stats[['排名', '细胞类型', '细胞数量', '占比(%)']]
    
    print("细胞类型详细统计:")
    print(celltype_stats.to_string(index=False))
    
    print(f"\n汇总信息:")
    print(f"细胞类型总数: {len(celltype_counts)} 种")
    print(f"总细胞数量: {celltype_counts.sum()} 个")
    print(f"平均每种细胞类型的数量: {celltype_counts.mean():.1f} 个")
    print(f"细胞数量最多的类型: {celltype_counts.idxmax()} ({celltype_counts.max()} 个, 占比{(celltype_counts.max()/celltype_counts.sum()*100):.1f}%)")
    print(f"细胞数量最少的类型: {celltype_counts.idxmin()} ({celltype_counts.min()} 个, 占比{(celltype_counts.min()/celltype_counts.sum()*100):.1f}%)")
    
    return celltype_stats

def analyze_tumor_cells(adata, celltype_col='celltype'):
    """
    统计细胞类型中包含'Tumor'字样的细胞数量并绘制饼图
    """
    print("\n=== Tumor相关细胞统计分析 ===")
    
    # 检测细胞类型名称中是否包含Tumor（不区分大小写）
    tumor_mask = adata.obs[celltype_col].str.contains('tumor', case=False, na=False)
    tumor_cell_count = tumor_mask.sum()
    non_tumor_cell_count = len(adata) - tumor_cell_count
    
    print(f"包含'Tumor'字样的细胞类型数量: {tumor_cell_count}")
    print(f"不包含'Tumor'字样的细胞类型数量: {non_tumor_cell_count}")
    print(f"Tumor细胞占比: {tumor_cell_count/len(adata)*100:.2f}%")
    
    # 找出具体的Tumor相关细胞类型
    tumor_celltypes = adata.obs.loc[tumor_mask, celltype_col].unique()
    print(f"Tumor相关细胞类型: {list(tumor_celltypes)}")
    
    # 统计每种Tumor细胞类型的数量
    tumor_celltype_counts = adata.obs.loc[tumor_mask, celltype_col].value_counts()
    print("\n各Tumor细胞类型数量分布:")
    for celltype, count in tumor_celltype_counts.items():
        percentage = count / tumor_cell_count * 100
        print(f"  {celltype}: {count} 个细胞 ({percentage:.1f}%)")
    
    return tumor_mask, tumor_cell_count, non_tumor_cell_count, tumor_celltypes

def create_pie_charts(adata, cross_table, sample_proportions, celltype_stats, 
                      tumor_cell_count, non_tumor_cell_count, tumor_celltypes,
                      celltype_col='celltype', sample_col='sample'):
    """
    创建饼图可视化
    """
    print("\n=== 生成饼图可视化 ===")
    
    # 设置图形样式
    plt.style.use('default')
    fig = plt.figure(figsize=(20, 16))
    
    # 1. 细胞类型分布饼图（前10种）
    ax1 = plt.subplot(2, 3, 1)
    celltype_counts = adata.obs[celltype_col].value_counts()
    top_n = min(10, len(celltype_counts))
    top_celltypes = celltype_counts.head(top_n)
    
    # 计算其他类型的总和
    if len(celltype_counts) > top_n:
        other_count = celltype_counts[top_n:].sum()
        pie_data = pd.concat([top_celltypes, pd.Series([other_count], index=['其他'])])
    else:
        pie_data = top_celltypes
    
    colors1 = plt.cm.Set3(np.linspace(0, 1, len(pie_data)))
    wedges1, texts1, autotexts1 = ax1.pie(pie_data.values, labels=pie_data.index, 
                                          autopct='%1.1f%%', colors=colors1, startangle=90,
                                          textprops={'fontsize': 10})
    
    # 设置百分比文本样式
    for autotext in autotexts1:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax1.set_title(f'细胞类型分布 (前{top_n}种)', fontsize=14, fontweight='bold', pad=20)
    
    # 2. Tumor vs Non-Tumor细胞分布饼图
    ax2 = plt.subplot(2, 3, 2)
    tumor_data = [tumor_cell_count, non_tumor_cell_count]
    tumor_labels = [f'Tumor细胞\n{tumor_cell_count}个', f'非Tumor细胞\n{non_tumor_cell_count}个']
    colors2 = ['#ff9999', '#66b3ff']
    
    wedges2, texts2, autotexts2 = ax2.pie(tumor_data, labels=tumor_labels, autopct='%1.1f%%',
                                         colors=colors2, startangle=90, explode=(0.1, 0))
    
    for autotext in autotexts2:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax2.set_title('Tumor vs 非Tumor细胞分布', fontsize=14, fontweight='bold', pad=20)
    
    # 3. 样本细胞数量分布饼图（前8个样本）
    ax3 = plt.subplot(2, 3, 3)
    sample_counts = adata.obs[sample_col].value_counts()
    top_samples = min(8, len(sample_counts))
    sample_data = sample_counts.head(top_samples)
    
    if len(sample_counts) > top_samples:
        other_sample_count = sample_counts[top_samples:].sum()
        sample_pie_data = pd.concat([sample_data, pd.Series([other_sample_count], index=['其他样本'])])
    else:
        sample_pie_data = sample_data
    
    colors3 = plt.cm.Pastel1(np.linspace(0, 1, len(sample_pie_data)))
    wedges3, texts3, autotexts3 = ax3.pie(sample_pie_data.values, labels=sample_pie_data.index,
                                         autopct='%1.1f%%', colors=colors3, startangle=90)
    
    for autotext in autotexts3:
        autotext.set_fontsize(9)
    
    ax3.set_title(f'样本细胞数量分布 (前{top_samples}个样本)', fontsize=14, fontweight='bold', pad=20)
    
    # 4. Tumor细胞类型详细分布
    ax4 = plt.subplot(2, 3, 4)
    if tumor_cell_count > 0:
        tumor_celltype_counts = adata.obs[adata.obs[celltype_col].str.contains('tumor', case=False, na=False)][celltype_col].value_counts()
        
        if len(tumor_celltype_counts) > 0:
            colors4 = plt.cm.Reds(np.linspace(0.4, 0.8, len(tumor_celltype_counts)))
            wedges4, texts4, autotexts4 = ax4.pie(tumor_celltype_counts.values, 
                                                  labels=tumor_celltype_counts.index,
                                                  autopct='%1.1f%%', colors=colors4,
                                                  startangle=90)
            
            for autotext in autotexts4:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
            
            ax4.set_title('Tumor细胞类型详细分布', fontsize=14, fontweight='bold', pad=20)
        else:
            ax4.text(0.5, 0.5, '无Tumor细胞类型数据', ha='center', va='center', 
                    transform=ax4.transAxes, fontsize=12)
            ax4.set_title('Tumor细胞类型详细分布', fontsize=14, fontweight='bold', pad=20)
    else:
        ax4.text(0.5, 0.5, '未检测到Tumor细胞', ha='center', va='center', 
                transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Tumor细胞类型详细分布', fontsize=14, fontweight='bold', pad=20)
    
    # 5. 细胞类型数量分布（环形图）
    ax5 = plt.subplot(2, 3, 5)
    top_celltypes_ring = min(8, len(celltype_counts))
    ring_data = celltype_counts.head(top_celltypes_ring)
    
    if len(celltype_counts) > top_celltypes_ring:
        other_ring_count = celltype_counts[top_celltypes_ring:].sum()
        ring_data = pd.concat([ring_data, pd.Series([other_ring_count], index=['其他'])])
    
    colors5 = plt.cm.viridis(np.linspace(0, 1, len(ring_data)))
    wedges5, texts5, autotexts5 = ax5.pie(ring_data.values, labels=ring_data.index,
                                        autopct='%1.1f%%', colors=colors5,
                                        wedgeprops=dict(width=0.3), startangle=90)
    
    ax5.set_title(f'细胞类型分布环形图 (前{top_celltypes_ring}种)', fontsize=14, fontweight='bold', pad=20)
    
    # 6. 细胞类型占比总结
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    # 准备摘要文本
    summary_text = "细胞类型统计摘要\n\n"
    summary_text += f"总细胞数: {len(adata):,}\n"
    summary_text += f"总细胞类型数: {len(celltype_counts)}\n"
    summary_text += f"Tumor细胞数: {tumor_cell_count:,} ({tumor_cell_count/len(adata)*100:.1f}%)\n"
    summary_text += f"非Tumor细胞数: {non_tumor_cell_count:,} ({non_tumor_cell_count/len(adata)*100:.1f}%)\n\n"
    
    if tumor_cell_count > 0:
        summary_text += f"Tumor细胞类型: {len(tumor_celltypes)} 种\n"
        summary_text += "主要Tumor类型:\n"
        tumor_top = adata.obs[adata.obs[celltype_col].str.contains('tumor', case=False, na=False)][celltype_col].value_counts().head(3)
        for i, (celltype, count) in enumerate(tumor_top.items()):
            summary_text += f"  {i+1}. {celltype}: {count} 个\n"
    
    summary_text += f"\n最多的细胞类型:\n"
    summary_text += f"  {celltype_counts.index[0]}: {celltype_counts.values[0]:,} 个 ({celltype_stats.iloc[0, 3]}%)\n"
    
    ax6.text(0.1, 0.95, summary_text, transform=ax6.transAxes, fontsize=11,
            verticalalignment='top', linespacing=1.5,
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    plt.tight_layout()
    plt.show()
    
    return fig

def advanced_statistical_tests(cross_table):
    """
    进行高级统计检验
    """
    print("\n=== 高级统计分析 ===")
    
    # 卡方检验
    chi2, p_value, dof, expected = chi2_contingency(cross_table)
    print(f"卡方检验结果:")
    print(f"卡方值: {chi2:.3f}")
    print(f"P值: {p_value:.3e}")
    print(f"自由度: {dof}")
    
    if p_value < 0.05:
        print("结论: 样本和细胞类型分布不独立 (p < 0.05)")
    else:
        print("结论: 样本和细胞类型分布独立 (p ≥ 0.05)")
    
    return chi2, p_value

def export_statistical_results(cross_table, sample_proportions, celltype_stats, 
                              tumor_cell_count, non_tumor_cell_count, tumor_celltypes,
                              output_prefix='celltype_sample_stats'):
    """
    导出统计结果到文件
    """
    import os
    output_dir = 'statistical_results'
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存交叉表
    cross_table.to_csv(os.path.join(output_dir, f'{output_prefix}_cross_table.csv'))
    
    # 保存比例表
    sample_proportions.to_csv(os.path.join(output_dir, f'{output_prefix}_proportions.csv'))
    
    # 保存细胞类型统计
    celltype_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_celltype_stats.csv'), index=False)
    
    # 保存Tumor细胞统计
    tumor_stats = pd.DataFrame({
        '类别': ['Tumor细胞', '非Tumor细胞'],
        '数量': [tumor_cell_count, non_tumor_cell_count],
        '占比(%)': [tumor_cell_count/(tumor_cell_count+non_tumor_cell_count)*100,
                  non_tumor_cell_count/(tumor_cell_count+non_tumor_cell_count)*100]
    })
    tumor_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_tumor_stats.csv'), index=False)
    
    # 保存汇总统计
    summary_stats = pd.DataFrame({
        '指标': ['总细胞数', '总样本数', '总细胞类型数', 'Tumor细胞数', 'Tumor细胞占比(%)',
                '非Tumor细胞数', '非Tumor细胞占比(%)', 'Tumor细胞类型数'],
        '数值': [
            cross_table.sum().sum(),
            len(cross_table),
            len(cross_table.columns),
            tumor_cell_count,
            tumor_cell_count/cross_table.sum().sum()*100,
            non_tumor_cell_count,
            non_tumor_cell_count/cross_table.sum().sum()*100,
            len(tumor_celltypes) if tumor_celltypes is not None else 0
        ]
    })
    summary_stats.to_csv(os.path.join(output_dir, f'{output_prefix}_summary.csv'), index=False)
    
    print(f"\n统计结果已导出到目录: {output_dir}/")
    print(f"1. {output_prefix}_cross_table.csv - 样本-细胞类型交叉表")
    print(f"2. {output_prefix}_proportions.csv - 样本细胞类型比例表")
    print(f"3. {output_prefix}_celltype_stats.csv - 细胞类型统计表")
    print(f"4. {output_prefix}_tumor_stats.csv - Tumor细胞统计表")
    print(f"5. {output_prefix}_summary.csv - 汇总统计表")

# 主执行函数
def main():
    """
    主执行函数
    """
    # 读取h5ad文件
    h5ad_file = 'merged_single_cell_data.h5ad'  # 请修改为你的文件路径
    
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
        print(f"celltype列缺失: {missing_celltype} ({missing_celltype/adata.n_obs*100:.1f}%)")
        print(f"sample列缺失: {missing_sample} ({missing_sample/adata.n_obs*100:.1f}%)")
        
        # 移除缺失值
        original_size = adata.n_obs
        adata = adata[~adata.obs[celltype_col].isnull() & ~adata.obs[sample_col].isnull()].copy()
        print(f"已移除包含缺失值的细胞，剩余细胞数: {adata.n_obs} (移除了 {original_size - adata.n_obs} 个细胞)")
    
    # 执行统计分析
    cross_table, sample_proportions = analyze_celltype_sample_stats(adata, celltype_col, sample_col)
    
    # 统计所有细胞类型的总数
    celltype_stats = count_all_celltypes(adata, celltype_col)
    
    # 统计Tumor相关细胞
    tumor_mask, tumor_cell_count, non_tumor_cell_count, tumor_celltypes = analyze_tumor_cells(adata, celltype_col)
    
    # 创建饼图可视化
    pie_chart_fig = create_pie_charts(adata, cross_table, sample_proportions, celltype_stats,
                                     tumor_cell_count, non_tumor_cell_count, tumor_celltypes,
                                     celltype_col, sample_col)
    
    # 高级统计检验
    chi2, p_value = advanced_statistical_tests(cross_table)
    
    # 导出结果
    output_prefix = input("\n请输入输出文件的前缀 (直接回车使用默认前缀 'celltype_sample_stats'): ").strip()
    if not output_prefix:
        output_prefix = 'celltype_sample_stats'
    
    export_statistical_results(cross_table, sample_proportions, celltype_stats,
                             tumor_cell_count, non_tumor_cell_count, tumor_celltypes,
                             output_prefix)
    
    print("\n" + "="*60)
    print("分析完成!")
    print("="*60)
    print(f"总计分析了 {adata.n_obs} 个细胞")
    print(f"来自 {adata.obs[sample_col].nunique()} 个样本")
    print(f"共 {adata.obs[celltype_col].nunique()} 种细胞类型")
    print(f"其中包含 {tumor_cell_count} 个Tumor相关细胞 ({tumor_cell_count/adata.n_obs*100:.1f}%)")
    print("="*60)

if __name__ == "__main__":
    main()
