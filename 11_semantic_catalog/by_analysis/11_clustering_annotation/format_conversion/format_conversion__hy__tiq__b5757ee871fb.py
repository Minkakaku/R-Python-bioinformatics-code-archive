import scanpy as sc
import pandas as pd

# 读取 h5ad 文件
adata = sc.read_h5ad('adata_celltype.h5ad')  # 请将 'your_file.h5ad' 替换为您的实际文件路径

# 从 adata.obs 中提取 sample 和 group 列
# 使用双重括号 [[]] 来确保提取的结果是一个 DataFrame，而不是 Series
obs_data = adata.obs[['sample', 'group']]

# 将数据保存为 CSV 文件
# index=True 会将细胞的索引（通常是细胞条形码）也保存到 CSV 文件中
obs_data.to_csv('extracted_sample_group.csv', index=True)

print("数据提取完成，已保存为 extracted_sample_group.csv")
