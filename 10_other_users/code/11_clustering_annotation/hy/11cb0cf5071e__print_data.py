import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from matplotlib.colors import Normalize
import warnings


# 1. 读取h5ad数据并加载聚类信息
adata = sc.read_h5ad("cell_type_Epithelial_cells.h5ad")  # 替换为你的h5ad文件路径
print(adata)
print(adata.obs['group'])
print(adata.obs['n_genes_by_counts'])
print(adata.obs['pct_counts_mt'])
print(adata.obs['n_genes'])
print(adata.obs['pct_counts_ribo'])
#print(adata.raw.X)

# 定义要检查的列
columns_to_check = [
    'n_genes_by_counts',
    'pct_counts_mt',
    'n_genes',
    'pct_counts_ribo'
]

print("数值列的最小值和最大值:")
print("=" * 50)

for col in columns_to_check:
    if col in adata.obs.columns:
        col_data = adata.obs[col]
        
        # 尝试转换为数值类型
        try:
            numeric_data = pd.to_numeric(col_data, errors='coerce')
            if pd.api.types.is_numeric_dtype(col_data):
                print(f"{col}:")
                print(f"  最小值: {numeric_data.min():.4f}")
                print(f"  最大值: {numeric_data.max():.4f}")
            else:
                print(f"{col}: 非数值类型")
        except:
            print(f"{col}: 无法转换为数值类型")
    else:
        print(f"{col}: 列不存在")
