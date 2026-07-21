# Python 示例
import pandas as pd
import anndata
# 加载文件
meta = pd.read_csv("cell_phone/meta.tsv", sep="\t")
adata = anndata.read_h5ad("integrated_data.h5ad")  # 或读取文本矩阵
# 检查细胞数量
print("Meta细胞数:", meta.shape[0])
print("Counts细胞数:", adata.n_obs)  # 或 adata.shape[0]
# 检查ID一致性
meta_cells = set(meta["cell"])  # 假设列名为"Cell"
counts_cells = set(adata.obs_names)
print("ID差异:", meta_cells.symmetric_difference(counts_cells))
