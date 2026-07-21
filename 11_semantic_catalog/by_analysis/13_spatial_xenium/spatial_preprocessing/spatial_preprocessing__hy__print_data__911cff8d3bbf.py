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
adata = sc.read_h5ad("adata_filtered_with_MPCpercent.h5ad")  # 替换为你的h5ad文件路径
print(adata)
print(adata.obs['cell_id'])
print(adata.raw)
# 2. 从adata.obs中提取'sample'列
sample_data = adata.obs['sample']

# 3. 将'sample'列的数据写出到sample.txt文件
#sample_data.to_csv('sample.txt', sep='\t', index=True, header=False)


# 检查'sample'列是否存在
#if 'sample' not in adata.obs.columns:
#    print("错误：在adata.obs中未找到'sample'列")
#    print("可用的列有：", list(adata.obs.columns))
#else:
#    # 统计每个样本的细胞数量
#    cell_counts = adata.obs['sample'].value_counts().reset_index()
#    cell_counts.columns = ['sample', 'cell_count']
#    
#    # 按样本名排序（可选）
#    cell_counts = cell_counts.sort_values('sample')
#    
#    # 输出统计结果到屏幕
#    print("每个样本的细胞数量统计：")
#    print(cell_counts)
#    
#    # 将结果保存到sample_cell.txt文件
#    cell_counts.to_csv('sample_cell.txt', sep='\t', index=False)
#    
#    print("\n结果已成功导出到 sample_cell.txt")
