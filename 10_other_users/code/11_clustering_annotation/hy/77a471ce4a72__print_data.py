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
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15_T.h5ad")  


# 替换为你的h5ad文件路径
print(adata)
print(adata.obs['state'])
print(adata.obs['sample_new'])
print(adata.obs['T'])
print(adata.X)

