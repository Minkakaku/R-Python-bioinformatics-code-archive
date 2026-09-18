import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse

# 1. 读取 h5ad 文件
adata1 = sc.read_h5ad("integrated_data.h5ad")
adata = adata1[adata1.obs["group"] == "disease"].copy()

# 2. 提取元数据 (meta.csv)
meta_df = adata.obs[["celltype_final"]].reset_index()  # 第一列为细胞ID
print(meta_df)
meta_df.columns = ["cell_id", "cell_annotation"]  # 重命名列
meta_df.to_csv("meta.csv", index=False)

# 3. 提取原始计数矩阵
if adata.raw is not None:  # 优先从 raw 属性获取
    raw_counts = adata.raw.X.T
    gene_names = adata.raw.var_names
    cell_names = adata.obs_names
    print(cell_names)  # 使用当前adata的细胞名
elif "counts" in adata.layers:  # 其次检查 layers
    raw_counts = adata.layers["counts"].T  # 需转置为基因×细胞
    gene_names = adata.var_names
    cell_names = adata.obs_names  # 使用当前adata的细胞名
else:
    raise ValueError("原始计数矩阵未找到！请检查 raw.X 或 layers['counts']")

# 4. 标准化：文库大小归一化 + log1p
# 4.1 计算每细胞的归一化因子（CPM）
cell_sums = raw_counts.sum(axis=0)  # 修正：按列求和（每个细胞）
if issparse(cell_sums):
    cell_sums = cell_sums.A1  # 稀疏矩阵转一维数组
norm_factors = 1e6 / cell_sums

# 4.2 归一化并取对数
if issparse(raw_counts):
    norm_data = raw_counts.multiply(norm_factors)  # 稀疏矩阵优化
else:
    norm_data = raw_counts * norm_factors
log_norm_data = np.log1p(norm_data)  # log(x+1) 转换

# 5. 保存标准化矩阵 (count.csv)
# 关键修正：确保列名为细胞ID
count_df = pd.DataFrame(
    log_norm_data.toarray() if issparse(log_norm_data) else log_norm_data,
    index=gene_names,          # 行索引=基因名
    columns=cell_names         # 列索引=细胞ID
)
count_df.reset_index(inplace=True)  # 将基因名从索引转为列
count_df.rename(columns={"index": "Gene_id"}, inplace=True)

# 修正列名：第一列名设为"Gene_id"，其他列名保持为细胞ID
count_df.to_csv("count.csv", index=False)  # 不保存行索引
