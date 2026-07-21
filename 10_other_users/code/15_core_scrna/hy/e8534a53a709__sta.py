import scanpy as sc
import pandas as pd

# 读取并统计
adata = sc.read_h5ad("result.h5ad")
celltype_counts = adata.obs['cell_type'].value_counts(dropna=False)
total_cells = len(adata)

# 创建结果并保存
result = pd.DataFrame({
    'cell_type': celltype_counts.index,
    'count': celltype_counts.values,
    'percentage': (celltype_counts / total_cells * 100).round(2)
}).sort_values('count', ascending=False)

result.to_csv("flex.csv", index=False)
