import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import warnings
warnings.filterwarnings('ignore')

# 设置字体和字号
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

# 创建输出目录
output_dir = "CNs15_celltype_analysis"
os.makedirs(output_dir, exist_ok=True)

# 1. 读取h5ad文件
print("Reading h5ad file...")
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15.h5ad")
print(f"Data shape: {adata.shape}")
print(f"obs columns: {adata.obs.columns.tolist()}")

# 检查必需列
required_cols = ['CNs_15', 'celltype']
missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    print(f"Error: Missing required columns: {missing_cols}")
    print(f"Available columns: {adata.obs.columns.tolist()}")
    exit(1)

# 创建数据副本
df = adata.obs[required_cols].copy()

# 检查数据
print(f"\nBasic information:")
print(f"Total cells: {len(df)}")
print(f"CNs_15 unique values: {df['CNs_15'].nunique()}")
print(f"celltype unique values: {df['celltype'].nunique()}")

# 检查CNs_15的取值
CNs_15_values = sorted(df['CNs_15'].dropna().unique().tolist())
print(f"CNs_15 groups: {CNs_15_values}")

# 2. 统计每个CNs_15中celltype的分布
print("\n" + "="*50)
print("Calculating cell type proportions in each CNs_15 cluster")
print("="*50)

# 2.1 计算原始计数
count_matrix = df.groupby(['CNs_15', 'celltype']).size().unstack(fill_value=0)

# 2.2 计算比例
proportion_matrix = count_matrix.div(count_matrix.sum(axis=1), axis=0) * 100

# 2.3 保存统计结果
count_matrix.to_csv(os.path.join(output_dir, "CNs15_celltype_counts.csv"), encoding='utf-8-sig')
proportion_matrix.to_csv(os.path.join(output_dir, "CNs15_celltype_proportions.csv"), encoding='utf-8-sig')

print(f"Count matrix saved to: {output_dir}/CNs15_celltype_counts.csv")
print(f"Proportion matrix saved to: {output_dir}/CNs15_celltype_proportions.csv")

# 3. 绘制热图
print("\n" + "="*50)
print("Generating heatmap")
print("="*50)

def plot_heatmap(proportion_df, output_dir, top_n=30):
    """
    绘制热图显示每个CNs_15中celltype的比例分布
    """
    # 如果celltype太多，只显示前top_n个最丰富的类型
    if proportion_df.shape[1] > top_n:
        # 计算每个celltype在所有CNs_15中的总比例
        celltype_totals = proportion_df.sum(axis=0)
        top_celltypes = celltype_totals.nlargest(top_n).index.tolist()
        heatmap_data = proportion_df[top_celltypes]
        print(f"Displaying top {top_n} cell types out of {proportion_df.shape[1]}")
    else:
        heatmap_data = proportion_df
        top_celltypes = proportion_df.columns.tolist()
    
    # 创建热图
    plt.figure(figsize=(max(12, len(heatmap_data.columns)*0.5), 
                       max(8, len(heatmap_data.index)*0.4)))
    
    # 使用颜色映射
    cmap = cm.get_cmap('YlOrRd')
    
    # 绘制热图
    sns.heatmap(heatmap_data, 
                cmap=cmap,
                annot=True, 
                fmt='.1f',
                annot_kws={'size': 8},
                linewidths=0.5,
                cbar_kws={'label': 'Proportion (%)'})
    
    plt.title('Cell Type Proportions in Each CNs_15 Cluster', fontsize=14, pad=20)
    plt.xlabel('Cell Type', fontsize=12)
    plt.ylabel('CNs_15 Cluster', fontsize=12)
    plt.tight_layout()
    
    # 保存热图
    output_file = os.path.join(output_dir, "CNs15_celltype_heatmap.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Heatmap saved to: {output_file}")
    
    return heatmap_data

# 绘制热图
heatmap_data = plot_heatmap(proportion_matrix, output_dir, top_n=30)

# 4. 为每个CNs_15绘制柱状图
print("\n" + "="*50)
print("Generating bar plots for each CNs_15 cluster")
print("="*50)

def plot_cn_bar_charts(proportion_df, count_df, output_dir, top_n=10):
    """
    为每个CNs_15绘制柱状图，显示前top_n个celltype和rest
    """
    all_cn_clusters = proportion_df.index.tolist()
    
    # 创建子目录保存柱状图
    barplot_dir = os.path.join(output_dir, "cn_bar_plots")
    os.makedirs(barplot_dir, exist_ok=True)
    
    # 存储每个CNs_15的前top_n celltypes
    cn_top_celltypes = {}
    
    for i, cn_cluster in enumerate(all_cn_clusters, 1):
        print(f"Processing CNs_15 cluster {i}/{len(all_cn_clusters)}: {cn_cluster}")
        
        # 获取该CNs_15的所有celltype比例
        cn_proportions = proportion_df.loc[cn_cluster]
        
        # 按比例排序
        cn_proportions_sorted = cn_proportions.sort_values(ascending=False)
        
        # 获取前top_n个celltype
        top_celltypes = cn_proportions_sorted.head(top_n)
        rest_celltypes = cn_proportions_sorted.iloc[top_n:]
        
        # 计算rest的总比例
        rest_proportion = rest_celltypes.sum()
        
        # 创建绘图数据
        plot_data = top_celltypes.copy()
        if rest_proportion > 0:
            plot_data['Rest'] = rest_proportion
        
        # 保存每个CNs_15的前top_n celltypes
        cn_top_celltypes[cn_cluster] = top_celltypes.index.tolist()
        
        # 绘制柱状图
        plt.figure(figsize=(max(10, len(plot_data)*0.6), 8))
        
        # 创建颜色数组
        colors = plt.cm.Set3(np.linspace(0, 1, len(plot_data)))
        
        # 绘制柱状图
        bars = plt.bar(range(len(plot_data)), plot_data.values, color=colors, edgecolor='black', linewidth=0.5)
        
        # 添加数值标签
        for idx, (bar, value) in enumerate(zip(bars, plot_data.values)):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    f'{value:.1f}%', ha='center', va='bottom', fontsize=9)
        
        # 设置X轴标签
        plt.xticks(range(len(plot_data)), plot_data.index.tolist(), rotation=45, ha='right', fontsize=10)
        
        # 设置标题和标签
        plt.title(f'Cell Type Composition in CNs_15: {cn_cluster}', fontsize=13, pad=20)
        plt.xlabel('Cell Type', fontsize=11)
        plt.ylabel('Proportion (%)', fontsize=11)
        
        # 添加网格
        plt.grid(True, axis='y', alpha=0.3, linestyle='--')
        
        # 设置Y轴范围
        plt.ylim(0, plot_data.max() * 1.2)
        
        # 添加总数标签
        total_cells = count_df.loc[cn_cluster].sum()
        plt.text(0.02, 0.98, f'Total cells: {int(total_cells)}', 
                transform=plt.gca().transAxes, fontsize=10, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        # 保存柱状图
        safe_cn = str(cn_cluster).replace('/', '_').replace('\\', '_').replace(':', '_')
        output_file = os.path.join(barplot_dir, f"CNs15_{safe_cn}_celltype_barplot.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    # 保存每个CNs_15的前top_n celltypes
    top_celltypes_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cn_top_celltypes.items()]))
    top_celltypes_df.to_csv(os.path.join(output_dir, "CNs15_top_celltypes_per_cluster.csv"), 
                           index=False, encoding='utf-8-sig')
    
    print(f"All bar plots saved to: {barplot_dir}/")
    print(f"Top cell types per cluster saved to: {output_dir}/CNs15_top_celltypes_per_cluster.csv")
    
    return cn_top_celltypes

# 绘制每个CNs_15的柱状图
cn_top_celltypes = plot_cn_bar_charts(proportion_matrix, count_matrix, output_dir, top_n=10)

# 5. 创建汇总柱状图
print("\n" + "="*50)
print("Generating summary bar plot")
print("="*50)

def plot_summary_bar_chart(proportion_df, output_dir, top_n=5):
    """
    创建汇总柱状图，显示每个CNs_15中top_n个主要的celltype
    """
    all_cn_clusters = proportion_df.index.tolist()
    
    # 收集每个CNs_15的主要celltype
    summary_data = []
    
    for cn_cluster in all_cn_clusters:
        cn_proportions = proportion_df.loc[cn_cluster]
        top_celltypes = cn_proportions.nlargest(top_n)
        
        for celltype, proportion in top_celltypes.items():
            summary_data.append({
                'CNs_15': cn_cluster,
                'celltype': celltype,
                'proportion': proportion,
                'rank': 1  # 这里可以改为实际排名
            })
    
    summary_df = pd.DataFrame(summary_data)
    
    # 创建汇总图
    fig, axes = plt.subplots(1, 1, figsize=(max(12, len(all_cn_clusters)*1.5), 8))
    
    # 准备数据用于分组柱状图
    plot_data = {}
    for cn_cluster in all_cn_clusters:
        cn_data = summary_df[summary_df['CNs_15'] == cn_cluster]
        if not cn_data.empty:
            # 只取前top_n个
            top_cn_celltypes = cn_data.nlargest(top_n, 'proportion')
            plot_data[cn_cluster] = top_cn_celltypes.set_index('celltype')['proportion'].to_dict()
    
    # 创建分组柱状图
    bar_width = 0.8 / top_n
    
    for i, cn_cluster in enumerate(all_cn_clusters):
        if cn_cluster in plot_data:
            celltypes_data = plot_data[cn_cluster]
            
            for j, (celltype, proportion) in enumerate(celltypes_data.items()):
                x_pos = i + j * bar_width - (top_n-1) * bar_width / 2
                axes.bar(x_pos, proportion, width=bar_width, 
                        label=celltype if i == 0 else "",  # 避免重复标签
                        alpha=0.7, edgecolor='black')
    
    axes.set_xlabel('CNs_15 Cluster', fontsize=12)
    axes.set_ylabel('Proportion (%)', fontsize=12)
    axes.set_title(f'Top {top_n} Cell Types in Each CNs_15 Cluster', fontsize=14, pad=20)
    axes.set_xticks(range(len(all_cn_clusters)))
    axes.set_xticklabels(all_cn_clusters, rotation=45, ha='right', fontsize=10)
    
    # 添加图例
    handles, labels = axes.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    axes.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    # 保存汇总图
    output_file = os.path.join(output_dir, "CNs15_summary_top_celltypes.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Summary bar plot saved to: {output_file}")

# 绘制汇总柱状图
plot_summary_bar_chart(proportion_matrix, output_dir, top_n=5)

# 6. 生成分析报告
print("\n" + "="*50)
print("Generating analysis report")
print("="*50)

with open(os.path.join(output_dir, "analysis_report.txt"), "w", encoding='utf-8') as f:
    f.write("Cell Type Composition Analysis in CNs_15 Clusters\n")
    f.write("="*60 + "\n\n")
    
    f.write("1. DATA OVERVIEW\n")
    f.write("-"*40 + "\n")
    f.write(f"Total cells: {len(df)}\n")
    f.write(f"Number of CNs_15 clusters: {len(CNs_15_values)}\n")
    f.write(f"Number of cell types: {df['celltype'].nunique()}\n")
    
    f.write("\n2. CNs_15 CLUSTER SIZES\n")
    f.write("-"*40 + "\n")
    cn_sizes = df['CNs_15'].value_counts()
    for cn_cluster, size in cn_sizes.items():
        proportion = (size / len(df) * 100)
        f.write(f"{cn_cluster}: {size} cells ({proportion:.2f}%)\n")
    
    f.write("\n3. CELL TYPE DISTRIBUTION ACROSS CLUSTERS\n")
    f.write("-"*40 + "\n")
    
    # 统计每个celltype出现在多少个CNs_15中
    celltype_presence = {}
    for celltype in df['celltype'].unique():
        presence = (count_matrix[celltype] > 0).sum()
        celltype_presence[celltype] = presence
    
    f.write("Cell types present in multiple CNs_15 clusters:\n")
    for celltype, presence in sorted(celltype_presence.items(), key=lambda x: x[1], reverse=True)[:20]:
        f.write(f"  {celltype}: {presence} clusters\n")
    
    f.write("\n4. TOP CELL TYPES IN EACH CLUSTER\n")
    f.write("-"*40 + "\n")
    
    for cn_cluster in proportion_matrix.index:
        top_celltypes = proportion_matrix.loc[cn_cluster].nlargest(5)
        f.write(f"\n{cn_cluster}:\n")
        for celltype, proportion in top_celltypes.items():
            f.write(f"  {celltype}: {proportion:.2f}%\n")
    
    f.write("\n5. OUTPUT FILES\n")
    f.write("-"*40 + "\n")
    f.write("Data files:\n")
    f.write(f"  {output_dir}/CNs15_celltype_counts.csv\n")
    f.write(f"  {output_dir}/CNs15_celltype_proportions.csv\n")
    f.write(f"  {output_dir}/CNs15_top_celltypes_per_cluster.csv\n")
    f.write("\nVisualization files:\n")
    f.write(f"  {output_dir}/CNs15_celltype_heatmap.png\n")
    f.write(f"  {output_dir}/cn_bar_plots/ (individual bar plots for each CNs_15)\n")
    f.write(f"  {output_dir}/CNs15_summary_top_celltypes.png\n")
    f.write(f"  {output_dir}/analysis_report.txt\n")

print(f"\nAnalysis complete! All outputs saved to: {output_dir}/")
print(f"Report saved to: {output_dir}/analysis_report.txt")
