import scanpy as sc
import numpy as np
from scipy.sparse import issparse

# 1. 读取h5ad文件
adata = sc.read_h5ad("Mixed_inflammatory_cells.h5ad")  # 替换为实际路径

# 2. 检查是否存在raw属性
if adata.raw is None:
    raise ValueError("AnnData对象未包含raw属性，请确保在预处理时已保存原始数据")

# 3. 获取raw矩阵（保留原始未过滤的基因）
raw_matrix = adata.raw.X  # 通常是稀疏矩阵（如csr_matrix）

# 4. 计算每个细胞的基因数（非零元素数量）
if issparse(raw_matrix):
    # 稀疏矩阵：按行统计非零值
    genes_per_cell = raw_matrix.getnnz(axis=1)  
else:
    # 密集矩阵：统计每行非零值
    genes_per_cell = np.count_nonzero(raw_matrix, axis=1)  

# 5. 计算基因数的中位数
median_genes = np.median(genes_per_cell)

# 6. 输出结果
print(f"所有细胞的基因数数组（前5个）: {genes_per_cell[:5]}")
print(f"基因数的中位数: {median_genes}")
