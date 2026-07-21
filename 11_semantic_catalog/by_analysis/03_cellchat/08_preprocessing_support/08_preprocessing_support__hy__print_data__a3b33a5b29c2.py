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
adata = sc.read_h5ad("merged_annotated_final_group_reamin_new_with_CNs.h5ad")  # 替换为你的h5ad文件路径
print(adata.obs['singler_best'])
#print(adata.obs['scores.CD8+ Tem'])
print(adata.obs['group_new'])
print(adata.obs['CNs_10'])
print(adata.X)
print(adata)
#if 'group_new' in adata.obs.columns:
#    # 获取去重后的值
#    unique_groups = adata.obs['group_new'].unique()
#    
#    # 打印到控制台
#    print(f"Unique values in 'group_new' (共 {len(unique_groups)} 个):")
#    print(unique_groups)
#    print("\n" + "="*50)
#    
#    # 保存到 group.txt 文件
#    with open('group.txt', 'w') as f:
#        f.write(f"Unique values in 'group_new' (共 {len(unique_groups)} 个):\n")
#        f.write("="*50 + "\n")
#        for group in unique_groups:
#            f.write(f"{group}\n")
#    
#    print("✓ 已保存去重后的 group_new 值到 group.txt")
#    
#    # 额外：也可以保存为DataFrame格式（可选）
#    # 创建DataFrame并保存
#    df_unique = pd.DataFrame({
#        'group_new_unique': unique_groups
#    })
#    df_unique.to_csv('group_unique.csv', index=False)
#    print("✓ 同时保存了CSV格式到 group_unique.csv")
#    
#else:
#    print("⚠️ 错误: 'group_new' 列不存在于 adata.obs 中")
#    print("可用的列有:")
#    print(adata.obs.columns.tolist())
