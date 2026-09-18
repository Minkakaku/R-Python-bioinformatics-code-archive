import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import warnings
import os
warnings.filterwarnings('ignore')

# 设置绘图风格
plt.style.use('default')
sns.set_palette("Set3")

def load_and_prepare_data(h5ad_file):
    """
    加载h5ad文件并检查必要的列
    """
    print("正在加载数据...")
    adata = sc.read_h5ad(h5ad_file)
    
    # 检查必要的列
    if 'annotation' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'annotation'列")
    
    if 'patient' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'patient'列")
    
    print(f"数据形状: {adata.shape}")
    print(f"患者数量: {adata.obs['patient'].nunique()}")
    print(f"细胞类型数量: {adata.obs['annotation'].nunique()}")
    
    return adata

def get_top10_celltypes(adata):
    """
    统计数量排名前10的细胞类型
    """
    print("\n正在统计细胞类型数量...")
    
    # 统计每种细胞类型的数量
    celltype_counts = adata.obs['annotation'].value_counts()
    
    print(f"细胞类型数量统计 (前20):")
    for i, (celltype, count) in enumerate(celltype_counts.head(20).items()):
        print(f"  {i+1:2d}. {celltype}: {count}个细胞")
    
    # 获取前10细胞类型
    top10_celltypes = celltype_counts.head(10).index.tolist()
    other_celltypes = celltype_counts.index[10:].tolist()
    
    print(f"\n前10细胞类型: {top10_celltypes}")
    print(f"其他细胞类型数量: {len(other_celltypes)}")
    
    return top10_celltypes, other_celltypes, celltype_counts

def calculate_celltype_proportions(adata, top10_celltypes, other_celltypes):
    """
    计算每个患者中细胞类型的比例
    """
    print("\n正在计算细胞类型比例...")
    
    # 获取所有患者
    patients = adata.obs['patient'].unique()
    
    # 计算每个患者的总细胞数
    patient_total_cells = adata.obs['patient'].value_counts()
    
    # 创建患者-细胞类型的交叉表
    celltype_counts = pd.crosstab(adata.obs['patient'], adata.obs['annotation'])
    
    # 计算比例
    proportions = celltype_counts.div(patient_total_cells, axis=0)
    
    # 将其他细胞类型合并为"Others"
    if other_celltypes:
        proportions['Others'] = proportions[other_celltypes].sum(axis=1)
        proportions = proportions[top10_celltypes + ['Others']]
    else:
        proportions = proportions[top10_celltypes]
    
    # 重置索引，将patient变为列
    proportions = proportions.reset_index()
    proportions = proportions.rename(columns={'index': 'patient'})
    
    print(f"比例矩阵形状: {proportions.shape}")
    print(f"患者数量: {len(patients)}")
    
    return proportions, patients

def create_stacked_barplot(proportions, top10_celltypes, output_dir="results"):
    """
    创建细胞类型比例堆积柱状图
    """
    print("\n正在生成细胞类型比例堆积柱状图...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备数据
    plot_data = proportions.copy()
    
    # 设置图形大小
    n_patients = len(plot_data)
    fig_width = max(12, n_patients * 0.5)  # 根据患者数调整宽度
    fig_height = 8
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # 准备柱状图数据
    patient_labels = plot_data['patient']
    
    # 准备颜色
    n_categories = len(top10_celltypes) + 1  # 前10 + Others
    colors = plt.cm.Set3(np.linspace(0, 1, n_categories))
    
    # 绘制堆积柱状图
    bottom = np.zeros(n_patients)
    
    # 先绘制前10细胞类型
    for i, celltype in enumerate(top10_celltypes):
        ax.bar(range(n_patients), plot_data[celltype].values, 
               bottom=bottom, label=celltype, color=colors[i])
        bottom += plot_data[celltype].values
    
    # 最后绘制Others
    if 'Others' in plot_data.columns:
        ax.bar(range(n_patients), plot_data['Others'].values, 
               bottom=bottom, label='Others', color=colors[-1])
    
    # 设置x轴标签
    ax.set_xticks(range(n_patients))
    ax.set_xticklabels(patient_labels, rotation=45, ha='right', fontsize=10)
    
    # 设置y轴
    ax.set_ylabel('Cell Proportion', fontsize=12)
    ax.set_ylim([0, 1])
    
    # 添加标题
    ax.set_title('Top 10 Cell Type Proportions by Patient (with Others)', 
                fontsize=14, fontweight='bold')
    
    # 添加图例
    ax.legend(title='Cell Types', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # 保存图形
    plt.savefig(f'{output_dir}/top10_celltype_proportions_stacked.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/top10_celltype_proportions_stacked.pdf', 
                bbox_inches='tight')
    plt.close()
    
    print(f"已保存堆积柱状图到: {output_dir}/top10_celltype_proportions_stacked.png/pdf")

def create_pie_chart(adata, top10_celltypes, other_celltypes, output_dir="results"):
    """
    创建细胞类型分布的饼图
    """
    print("\n正在生成细胞类型分布饼图...")
    
    # 计算各类细胞的数量
    celltype_counts = adata.obs['annotation'].value_counts()
    
    # 计算前10细胞类型的数量
    top10_counts = [celltype_counts[celltype] for celltype in top10_celltypes]
    
    # 计算Others的数量
    others_count = celltype_counts[other_celltypes].sum() if other_celltypes else 0
    
    # 准备饼图数据
    sizes = top10_counts + [others_count]
    labels = top10_celltypes + ['Others']
    
    # 创建颜色
    n_categories = len(labels)
    colors = plt.cm.Set3(np.linspace(0, 1, n_categories))
    
    # 创建图形
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 第一个饼图：显示比例
    wedges, texts, autotexts = ax1.pie(sizes, labels=labels, colors=colors, 
                                      autopct='%1.1f%%', startangle=90)
    ax1.set_title('Cell Type Distribution (Percentage)', fontsize=14, fontweight='bold')
    
    # 第二个饼图：显示数量
    wedges2, texts2, autotexts2 = ax2.pie(sizes, labels=labels, colors=colors, 
                                         autopct=lambda p: f'{int(p * sum(sizes) / 100):,}', 
                                         startangle=90)
    ax2.set_title('Cell Type Distribution (Count)', fontsize=14, fontweight='bold')
    
    # 美化文本
    for autotext in autotexts + autotexts2:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    plt.tight_layout()
    
    # 保存图形
    plt.savefig(f'{output_dir}/celltype_distribution_pie.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/celltype_distribution_pie.pdf', 
                bbox_inches='tight')
    plt.close()

def create_summary_statistics(proportions, top10_celltypes, output_dir="results"):
    """
    创建汇总统计信息
    """
    print("\n正在生成汇总统计...")
    
    # 计算每个细胞类型的平均比例
    avg_proportions = proportions[top10_celltypes + (['Others'] if 'Others' in proportions.columns else [])].mean()
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 绘制水平条形图
    y_pos = np.arange(len(avg_proportions))
    bars = ax.barh(y_pos, avg_proportions.values)
    
    # 设置y轴标签
    ax.set_yticks(y_pos)
    ax.set_yticklabels(avg_proportions.index)
    
    # 设置x轴
    ax.set_xlabel('Average Proportion', fontsize=12)
    ax.set_xlim([0, 1])
    
    # 添加数值标签
    for i, (bar, proportion) in enumerate(zip(bars, avg_proportions.values)):
        width = bar.get_width()
        ax.text(width + 0.01, bar.get_y() + bar.get_height()/2, 
                f'{proportion:.3f}', ha='left', va='center', fontsize=10)
    
    # 添加标题
    ax.set_title('Average Cell Type Proportions Across All Patients', 
                fontsize=14, fontweight='bold')
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    # 保存图形
    plt.savefig(f'{output_dir}/celltype_average_proportions.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/celltype_average_proportions.pdf', 
                bbox_inches='tight')
    plt.close()

def save_results(proportions, celltype_counts, top10_celltypes, other_celltypes, output_dir="results"):
    """
    保存结果到Excel文件
    """
    print("\n正在保存结果到Excel文件...")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备细胞类型统计
    celltype_stats = []
    for celltype in top10_celltypes:
        celltype_stats.append({
            'CellType': celltype,
            'Count': celltype_counts[celltype],
            'Rank': f"Top {top10_celltypes.index(celltype) + 1}"
        })
    
    if other_celltypes:
        others_count = celltype_counts[other_celltypes].sum()
        celltype_stats.append({
            'CellType': 'Others',
            'Count': others_count,
            'Rank': 'Others'
        })
    
    celltype_stats_df = pd.DataFrame(celltype_stats)
    
    # 计算描述性统计
    desc_stats = proportions.describe()
    
    # 保存到Excel
    with pd.ExcelWriter(f'{output_dir}/celltype_analysis_results.xlsx') as writer:
        proportions.to_excel(writer, sheet_name='患者比例矩阵', index=False)
        celltype_stats_df.to_excel(writer, sheet_name='细胞类型统计', index=False)
        desc_stats.to_excel(writer, sheet_name='描述性统计')
    
    # 保存细胞类型数量到CSV
    celltype_counts_df = pd.DataFrame({
        'CellType': celltype_counts.index,
        'Count': celltype_counts.values
    })
    celltype_counts_df.to_csv(f'{output_dir}/celltype_counts.csv', index=False)
    
    print(f"已保存结果到: {output_dir}/celltype_analysis_results.xlsx")
    print(f"已保存细胞类型数量到: {output_dir}/celltype_counts.csv")

def main():
    """
    主函数
    """
    # 配置参数
    h5ad_file = "CN2_cells.h5ad"  # 修改为您的h5ad文件路径
    output_dir = "results"
    
    try:
        # 1. 加载数据
        adata = load_and_prepare_data(h5ad_file)
        
        # 2. 获取前10细胞类型
        top10_celltypes, other_celltypes, celltype_counts = get_top10_celltypes(adata)
        
        # 3. 计算细胞类型比例
        proportions, patients = calculate_celltype_proportions(adata, top10_celltypes, other_celltypes)
        
        # 4. 生成堆积柱状图
        create_stacked_barplot(proportions, top10_celltypes, output_dir)
        
        # 5. 生成饼图
        create_pie_chart(adata, top10_celltypes, other_celltypes, output_dir)
        
        # 6. 生成汇总统计
        create_summary_statistics(proportions, top10_celltypes, output_dir)
        
        # 7. 保存结果
        save_results(proportions, celltype_counts, top10_celltypes, other_celltypes, output_dir)
        
        # 8. 打印摘要
        print("\n" + "="*60)
        print("分析完成摘要")
        print("="*60)
        
        print(f"分析的细胞总数: {adata.shape[0]}")
        print(f"患者数量: {len(patients)}")
        print(f"前10细胞类型数量: {len(top10_celltypes)}")
        print(f"其他细胞类型数量: {len(other_celltypes) if other_celltypes else 0}")
        
        # 计算前10细胞类型占总细胞的比例
        top10_total = sum(celltype_counts[celltype] for celltype in top10_celltypes)
        top10_percentage = top10_total / adata.shape[0] * 100
        print(f"前10细胞类型占总细胞比例: {top10_percentage:.2f}%")
        
        if other_celltypes:
            others_total = celltype_counts[other_celltypes].sum()
            others_percentage = others_total / adata.shape[0] * 100
            print(f"其他细胞类型占总细胞比例: {others_percentage:.2f}%")
        
        print(f"\n所有结果已保存到 '{output_dir}' 目录")
        print(f"主要输出文件:")
        print(f"  - {output_dir}/top10_celltype_proportions_stacked.png/pdf: 堆积柱状图")
        print(f"  - {output_dir}/celltype_distribution_pie.png/pdf: 细胞类型分布饼图")
        print(f"  - {output_dir}/celltype_average_proportions.png/pdf: 平均比例条形图")
        print(f"  - {output_dir}/celltype_analysis_results.xlsx: 详细分析结果")
        print(f"  - {output_dir}/celltype_counts.csv: 细胞类型数量统计")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
