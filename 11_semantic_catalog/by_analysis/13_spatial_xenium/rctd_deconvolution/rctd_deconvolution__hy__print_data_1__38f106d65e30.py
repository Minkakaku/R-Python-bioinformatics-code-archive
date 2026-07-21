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
adata = sc.read_h5ad("adata_celltype_RCTD_removed_with_patient.h5ad")  # 替换为你的h5ad文件路径
print(adata.obs)
print(adata)
#if 'sample' in adata.obs.columns:
#    celltype_series = adata.obs['sample']
#unique_celltypes = celltype_series.unique().tolist()
#print("去重后的细胞类型列表 (共{}种):".format(len(unique_celltypes)))
#for i, celltype in enumerate(unique_celltypes, 1):
#    print(f"{i:2d}. {celltype}")
#print(adata.obs['annotation'])
#if 'annotation' in adata.obs.columns:
#    celltype_series = adata.obs['annotation']
#unique_celltypes = celltype_series.unique().tolist()
#print("去重后的细胞类型列表 (共{}种):".format(len(unique_celltypes)))
#for i, celltype in enumerate(unique_celltypes, 1):
#    print(f"{i:2d}. {celltype}")
