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
adata = sc.read_h5ad("merged_data.h5ad")  # 替换为你的h5ad文件路径
print(adata)
print(adata.var)
print(adata.obs)
print(adata.raw)

#print(adata.obs['Sample'].unique().tolist())
#tient  统计每个样本的细胞数量（网页6中的实际应用）
#print(adata.obs['cell_type'].unique().tolist())
# 若样本信息存储在 'batch' 列
# 若不确定列名，先查看所有列名
print(adata.obs.columns)
print(adata.raw.X)
