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
adata = sc.read_h5ad("MPCvsTNBC_all_sample_ct.h5ad")  # 替换为你的h5ad文件路径
print(adata)
print(adata.obs['group'])
print(adata.obs['annotation'])
# 2. 从adata.obs中提取'annotation'列
annotation_data = adata.obs['annotation']

# 3. 将'annotation'列的数据写出到annotation.txt文件
#annotation_data.to_csv('annotation.txt', sep='\t', index=True, header=False)
#
#
## 检查'annotation'列是否存在
#if 'annotation' not in adata.obs.columns:
#    print("错误：在adata.obs中未找到'annotation'列")
#    print("可用的列有：", list(adata.obs.columns))
#else:
#    # 统计每个样本的细胞数量
#    cell_counts = adata.obs['annotation'].value_counts().reset_index()
#    cell_counts.columns = ['annotation', 'cell_count']
#    
#    # 按样本名排序（可选）
#    cell_counts = cell_counts.sort_values('annotation')
#    
#    # 输出统计结果到屏幕
#    print("每个样本的细胞数量统计：")
#    print(cell_counts)
#    
#    # 将结果保存到annotation_cell.txt文件
#    cell_counts.to_csv('annotation_cell.txt', sep='\t', index=False)
#    
#    print("\n结果已成功导出到 annotation_cell.txt")
