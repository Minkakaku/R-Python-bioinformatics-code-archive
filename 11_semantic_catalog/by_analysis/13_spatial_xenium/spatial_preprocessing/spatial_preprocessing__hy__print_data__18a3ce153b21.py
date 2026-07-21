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
adata = sc.read_h5ad("CN9_cells.h5ad")  # 替换为你的h5ad文件路径
print(adata.raw)
print(adata)
