import scanpy as sc
import pandas as pd

# 1. 读取h5ad文件
adata = sc.read_h5ad('merged_all_celltypes.h5ad')  # 请将'your_file.h5ad'替换为您的实际文件路径

# 2. 从adata.obs中提取所需的列（sample和celltype）
# 确保列名与您的数据一致，这里假设列名分别为'sample'和'celltype'
df = adata.obs[['sample', 'celltype']].copy()

# 3. 计算每个sample-celltype组合的细胞数量
count_data = df.groupby(['sample', 'celltype']).size().reset_index(name='cell_count')

# 4. 计算每个样本的总细胞数
sample_total = df.groupby('sample').size().reset_index(name='total_cells')

# 5. 将总细胞数合并到count_data中，并计算百分比
result = count_data.merge(sample_total, on='sample')
result['percentage'] = (result['cell_count'] / result['total_cells'] * 100).round(2)  # 保留两位小数

# 6. 打印结果
print("每个样本中细胞类型的数量与百分比：")
print(result)

# 7. 可选：保存结果为CSV文件
result.to_csv('sample_celltype_stats.csv', index=False)
print("\n结果已保存到 sample_celltype_stats.csv")
