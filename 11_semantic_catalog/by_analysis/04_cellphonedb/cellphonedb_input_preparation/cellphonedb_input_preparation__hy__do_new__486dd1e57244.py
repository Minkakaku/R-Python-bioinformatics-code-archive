import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse

# 1. 读取 h5ad 文件
adata = sc.read_h5ad("integrated_data.h5ad")
#adata = adata1[adata1.obs["group"] == "normal"].copy()

# 2. 提取元数据 (meta.csv)
meta_df = adata.obs[["celltype_final"]].reset_index()
meta_df.columns = ["cell_id", "cell_annotation"]
meta_df.to_csv("meta.csv", index=False)

# 3. 提取原始计数矩阵
if adata.raw is not None:  # 优先从 raw 属性获取（原始计数）
    raw_counts = adata.raw.X.T
    gene_names = adata.raw.var_names
    cell_names = adata.obs_names  # 使用疾病组的细胞名
elif "counts" in adata.layers:  # 其次检查 layers
    raw_counts = adata.layers["counts"].T
    gene_names = adata.var_names
    cell_names = adata.obs_names
else:
    raise ValueError("原始计数矩阵未找到！请检查 raw.X 或 layers['counts']")

# 4. 标准化：仅进行文库大小归一化（CPM），去除log转换
# 4.1 计算每细胞的归一化因子（CPM）
cell_sums = raw_counts.sum(axis=0)  # 按列求和（每个细胞）
if issparse(cell_sums):
    cell_sums = cell_sums.A1  # 稀疏矩阵转一维数组
norm_factors = 1e6 / cell_sums

# 4.2 归一化（不进行log转换）
if issparse(raw_counts):
    # 稀疏矩阵优化：直接乘以缩放因子
    norm_data = raw_counts.multiply(norm_factors)
else:
    norm_data = raw_counts * norm_factors

# 5. 保存标准化矩阵 (count.csv)
count_df = pd.DataFrame(
    norm_data.toarray() if issparse(norm_data) else norm_data,
    index=gene_names,          # 行索引=基因名
    columns=cell_names         # 列索引=细胞ID
)
count_df.reset_index(inplace=True)
count_df.rename(columns={"index": "Gene_id"}, inplace=True)
count_df.to_csv("count.csv", index=False)  # 输出CPM矩阵（无量纲）
