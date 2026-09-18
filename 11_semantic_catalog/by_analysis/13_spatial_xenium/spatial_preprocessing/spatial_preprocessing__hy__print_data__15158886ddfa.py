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
adata = sc.read_h5ad("All_filter_with_region_CNs10.h5ad")  # 替换为你的h5ad文件路径
#print(adata.obs['sample']['POST_22'])
print(adata)
print(adata.obs['region'])
#print(adata.obs['celltype'])
#print(adata.obs['cell_type_source'])
print(adata.X)
