import anndata as ad
import pandas as pd

# 读取 h5ad 文件
adata = ad.read_h5ad("adata_celltype.h5ad")

# 检查是否存在 'sample' 列
if 'sample' not in adata.obs.columns:
    raise KeyError("数据中没有 'sample' 列，请检查列名。")

# 统计每个样本的细胞数
sample_counts = adata.obs['sample'].value_counts().reset_index()
sample_counts.columns = ['sample', 'cell_count']

# 保存为 CSV
sample_counts.to_csv("sample_counts.csv", index=False)

print("统计完成，结果已保存至 sample_counts.csv")
print(sample_counts.head())
