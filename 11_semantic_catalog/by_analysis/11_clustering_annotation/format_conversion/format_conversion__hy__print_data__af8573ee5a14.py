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
adata = sc.read_h5ad("merged_single_cell_data_modified.h5ad")  # 替换为你的h5ad文件路径
print(adata.obs['cell_id'])
print(adata.raw)
print(adata)
print(adata.obs['group'])
print(adata.obs['celltype'])
if 'celltype' in adata.obs.columns:
    celltype_series = adata.obs['celltype']
unique_celltypes = celltype_series.unique().tolist()
print("去重后的细胞类型列表 (共{}种):".format(len(unique_celltypes)))
for i, celltype in enumerate(unique_celltypes, 1):
    print(f"{i:2d}. {celltype}")
