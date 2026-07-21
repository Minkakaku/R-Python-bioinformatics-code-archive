import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 1. 加载数据
print("正在加载数据...")
adata = sc.read_h5ad('adata_cluster.h5ad')
print(f"数据加载完成，细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
print(f"聚类列: 'res_0.8'，可用聚类列: {[col for col in adata.obs.columns if 'res' in col or 'cluster' in col.lower()]}")

# 检查聚类列是否存在
if 'res_0.8' not in adata.obs.columns:
    print("错误: 未找到 'res_0.8' 列")
    print("可用的obs列:", list(adata.obs.columns))
    exit()

# 2. 统计每个cluster的细胞数量
print("\n" + "="*50)
print("Cluster统计信息")
print("="*50)

# 获取所有cluster
clusters = adata.obs['res_0.8'].unique()
clusters_sorted = sorted(clusters.astype(str))  # 排序便于阅读

# 统计每个cluster的细胞数量
cluster_counts = adata.obs['res_0.8'].value_counts()
cluster_counts_sorted = cluster_counts.sort_index(key=lambda x: x.astype(str))

# 创建统计表格
stats_df = pd.DataFrame({
    'Cluster': cluster_counts_sorted.index.astype(str),
    'Cell_Count': cluster_counts_sorted.values,
    'Percentage': (cluster_counts_sorted.values / len(adata) * 100).round(2)
})

# 3. 计算每个cluster的基因表达中位数
print("\n正在计算每个cluster的基因表达中位数...")

# 准备存储结果的字典
median_expressions = {}

# 对每个cluster计算基因表达中位数
for cluster in clusters_sorted:
    # 获取该cluster的细胞索引
    cluster_mask = adata.obs['res_0.8'] == cluster
    cluster_cells = adata[cluster_mask, :]
    
    # 获取表达矩阵
    if sparse.issparse(cluster_cells.X):
        # 对于稀疏矩阵
        expr_matrix = cluster_cells.X.toarray()
    else:
        expr_matrix = cluster_cells.X
    
    # 计算每个基因的中位数
    gene_medians = np.median(expr_matrix, axis=0)
    
    # 将结果存储在字典中
    median_expressions[cluster] = gene_medians
    
    # 计算该cluster的整体表达水平统计
    overall_median = np.median(expr_matrix)  # 所有基因所有细胞的中位数
    overall_mean = np.mean(expr_matrix)  # 所有基因所有细胞的均值
    
    # 添加到统计表格
    stats_df.loc[stats_df['Cluster'] == cluster, 'Overall_Median_Expression'] = overall_median
    stats_df.loc[stats_df['Cluster'] == cluster, 'Overall_Mean_Expression'] = overall_mean
    stats_df.loc[stats_df['Cluster'] == cluster, 'Zero_Expression_Percentage'] = (expr_matrix == 0).sum() / expr_matrix.size * 100

# 4. 创建基因表达中位数的DataFrame
# 获取基因名
gene_names = adata.var_names.tolist()

# 创建DataFrame，行为基因，列为cluster
median_df = pd.DataFrame(index=gene_names)
for cluster in clusters_sorted:
    median_df[cluster] = median_expressions[cluster]

# 5. 计算每个基因在所有cluster中的表达差异
median_df['Median_Across_Clusters'] = median_df.median(axis=1)  # 跨cluster的中位数
median_df['Std_Across_Clusters'] = median_df[clusters_sorted].std(axis=1)  # 跨cluster的标准差
median_df['Max_Min_Ratio'] = median_df[clusters_sorted].max(axis=1) / (median_df[clusters_sorted].min(axis=1) + 1e-10)  # 最大最小比

# 6. 输出统计结果
print("\n" + "="*50)
print("Cluster细胞数量统计")
print("="*50)
print(stats_df.to_string(index=False))

# 计算并显示总计
print(f"\n总计: {len(adata)} 个细胞, {len(clusters_sorted)} 个clusters")
print(f"平均每个cluster细胞数: {len(adata)/len(clusters_sorted):.1f}")
print(f"最大cluster: Cluster {stats_df.loc[stats_df['Cell_Count'].idxmax(), 'Cluster']} "
      f"({stats_df['Cell_Count'].max()} 个细胞, {stats_df['Cell_Count'].max()/len(adata)*100:.1f}%)")
print(f"最小cluster: Cluster {stats_df.loc[stats_df['Cell_Count'].idxmin(), 'Cluster']} "
      f"({stats_df['Cell_Count'].min()} 个细胞, {stats_df['Cell_Count'].min()/len(adata)*100:.1f}%)")

# 7. 可视化
print("\n正在生成可视化图表...")

# 创建输出目录
output_dir = Path("cluster_statistics")
output_dir.mkdir(exist_ok=True)

# 7.1 细胞数量条形图
plt.figure(figsize=(12, 6))
ax = sns.barplot(data=stats_df, x='Cluster', y='Cell_Count', palette='viridis')
plt.title(f'每个Cluster的细胞数量 (总计: {len(adata)} 个细胞)', fontsize=14, fontweight='bold')
plt.xlabel('Cluster', fontsize=12)
plt.ylabel('细胞数量', fontsize=12)

# 在柱子上添加数值
for i, (_, row) in enumerate(stats_df.iterrows()):
    ax.text(i, row['Cell_Count'] + max(stats_df['Cell_Count']) * 0.01, 
            f"{row['Cell_Count']}\n({row['Percentage']}%)", 
            ha='center', va='bottom', fontsize=9)

plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(output_dir / 'cluster_cell_counts.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'cluster_cell_counts.pdf', bbox_inches='tight')
plt.show()

# 7.2 表达水平热图（每个cluster整体表达中位数）
plt.figure(figsize=(10, 6))
cluster_order = stats_df.sort_values('Overall_Median_Expression', ascending=False)['Cluster'].tolist()
expr_matrix = stats_df.set_index('Cluster').loc[cluster_order][['Overall_Median_Expression', 'Overall_Mean_Expression']]

sns.heatmap(expr_matrix.T, annot=True, fmt='.3f', cmap='YlOrRd', 
            cbar_kws={'label': '表达水平'}, square=True)
plt.title('每个Cluster的整体表达水平', fontsize=14, fontweight='bold')
plt.xlabel('Cluster', fontsize=12)
plt.ylabel('表达指标', fontsize=12)
plt.tight_layout()
plt.savefig(output_dir / 'cluster_expression_levels.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'cluster_expression_levels.pdf', bbox_inches='tight')
plt.show()

# 7.3 表达分布箱线图
plt.figure(figsize=(14, 6))

# 采样部分基因以加快绘图速度
if len(gene_names) > 1000:
    sampled_genes = np.random.choice(gene_names, 1000, replace=False)
    plot_data = median_df.loc[sampled_genes, clusters_sorted]
else:
    plot_data = median_df[clusters_sorted]

# 将数据转换为长格式便于绘图
plot_data_long = plot_data.reset_index().melt(id_vars=['index'], 
                                                value_vars=clusters_sorted,
                                                var_name='Cluster', 
                                                value_name='Median_Expression')
plot_data_long = plot_data_long.rename(columns={'index': 'Gene'})

plt.subplot(1, 2, 1)
sns.boxplot(data=plot_data_long, x='Cluster', y='Median_Expression', palette='Set3')
plt.title('每个Cluster的基因表达中位数分布', fontsize=12)
plt.xlabel('Cluster', fontsize=10)
plt.ylabel('基因表达中位数', fontsize=10)
plt.xticks(rotation=45)

plt.subplot(1, 2, 2)
# 展示每个cluster的表达中位数的分布密度
for cluster in clusters_sorted[:10]:  # 限制显示前10个cluster避免过于拥挤
    cluster_data = plot_data_long[plot_data_long['Cluster'] == cluster]['Median_Expression']
    sns.kdeplot(data=cluster_data, label=f'Cluster {cluster}', alpha=0.7)
plt.title('基因表达中位数的分布密度', fontsize=12)
plt.xlabel('基因表达中位数', fontsize=10)
plt.ylabel('密度', fontsize=10)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.suptitle('基因表达分布分析', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(output_dir / 'gene_expression_distribution.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'gene_expression_distribution.pdf', bbox_inches='tight')
plt.show()

# 8. 保存结果到文件
print("\n正在保存结果到文件...")

# 保存cluster统计
stats_df.to_csv(output_dir / 'cluster_statistics.csv', index=False)
print(f"✓ Cluster统计已保存: {output_dir / 'cluster_statistics.csv'}")

# 保存详细的基因表达中位数
# 如果基因太多，可以保存为压缩格式或只保存部分
if len(gene_names) <= 10000:
    median_df.to_csv(output_dir / 'gene_median_expression_per_cluster.csv')
    print(f"✓ 基因表达中位数已保存: {output_dir / 'gene_median_expression_per_cluster.csv'}")
else:
    # 保存top变异基因
    top_variable_genes = median_df.nlargest(1000, 'Std_Across_Clusters')
    top_variable_genes.to_csv(output_dir / 'top_1000_variable_genes_per_cluster.csv')
    print(f"✓ Top 1000变异基因表达中位数已保存: {output_dir / 'top_1000_variable_genes_per_cluster.csv'}")
    
    # 保存所有基因的汇总统计
    summary_df = median_df[['Median_Across_Clusters', 'Std_Across_Clusters', 'Max_Min_Ratio']]
    summary_df.to_csv(output_dir / 'gene_expression_summary_statistics.csv')
    print(f"✓ 基因表达汇总统计已保存: {output_dir / 'gene_expression_summary_statistics.csv'}")

# 9. 生成汇总报告
print("\n" + "="*50)
print("汇总报告")
print("="*50)
print(f"1. 数据概览:")
print(f"   - 总细胞数: {len(adata)}")
print(f"   - 总基因数: {len(gene_names)}")
print(f"   - Cluster数: {len(clusters_sorted)}")

print(f"\n2. Cluster大小分布:")
print(f"   - 最大cluster: {stats_df['Cell_Count'].max()} 个细胞")
print(f"   - 最小cluster: {stats_df['Cell_Count'].min()} 个细胞")
print(f"   - 中位数cluster大小: {stats_df['Cell_Count'].median():.0f} 个细胞")

print(f"\n3. 表达水平统计:")
print(f"   - 所有细胞整体表达中位数: {stats_df['Overall_Median_Expression'].mean():.4f}")
print(f"   - 所有细胞整体表达均值: {stats_df['Overall_Mean_Expression'].mean():.4f}")
print(f"   - 平均零表达比例: {stats_df['Zero_Expression_Percentage'].mean():.2f}%")

print(f"\n4. 基因表达差异:")
if len(gene_names) <= 10000:
    high_var_genes = median_df.nlargest(5, 'Std_Across_Clusters')
    print(f"   - 变异最大的5个基因:")
    for idx, row in high_var_genes.iterrows():
        print(f"     {idx}: 标准差={row['Std_Across_Clusters']:.4f}, 最大最小比={row['Max_Min_Ratio']:.2f}")

print(f"\n5. 文件已保存到目录: {output_dir.resolve()}")
print("   包含: PNG/PDF格式图表、CSV格式数据表")

# 10. 显示聚类标签示例
print(f"\n" + "="*50)
print("聚类标签示例 (前10个细胞):")
print("="*50)
print(adata.obs[['res_0.8']].head(10))
