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
#adata = sc.read_h5ad("merged_single_cell_data_with_raw.h5ad")  # 替换为你的h5ad文件路径
adata = sc.read_h5ad("merged_data_with_ct_ori.h5ad")
print(adata)
print(adata.var)
a = adata.obs['celltype'].to_list
print(adata.obs['cell_id'])
print(adata.obs['celltype'])
#print(adata.uns['pca'])
print(adata.raw)

#print("\n" + "="*50)
#print("开始检测重复的细胞ID")
#print("="*50)
#
## 确定用于检测重复的细胞标识列
#if 'cell_id' in adata.obs.columns:
#    cell_identifier = adata.obs['cell_id']
#    print("使用 'cell_id' 列进行重复检测")
#    identifier_source = "列"
#else:
#    cell_identifier = adata.obs_names
#    print("未找到 'cell_id' 列，使用细胞索引（obs_names）进行重复检测")
#    identifier_source = "索引"
#
# 检测重复值
#duplicates_mask = cell_identifier.duplicated(keep=False)  # 标记所有重复项
#duplicate_counts = cell_identifier.value_counts()
#true_duplicates = duplicate_counts[duplicate_counts > 1]  # 出现次数大于1的才是真正的重复项
#
#print(f"\n检测结果：")
#print(f"总细胞数量：{len(cell_identifier)}")
#print(f"唯一细胞{identifier_source}数量：{cell_identifier.nunique()}")
#print(f"重复细胞{identifier_source}数量：{len(true_duplicates)}")
#
#if len(true_duplicates) > 0:
#    print(f"\n重复的细胞{identifier_source}（出现次数 > 1）：")
#    for cell_id, count in true_duplicates.head(10).items():  # 只显示前10个重复项
#        print(f"  {cell_id}: 出现 {count} 次")
#    
#    if len(true_duplicates) > 10:
#        print(f"  ... 还有 {len(true_duplicates) - 10} 个重复项未显示")
#    
#    # 显示重复细胞的详细信息
#    print(f"\n重复细胞的详细注释信息：")
#    duplicate_cell_details = adata.obs[duplicates_mask]
#    if 'celltype' in duplicate_cell_details.columns:
#        print(duplicate_cell_details[['celltype']].head(10))  # 显示前10个重复细胞的类型
#else:
#    print(f"未发现重复的细胞{identifier_source}")
#
## 可选：处理重复项的方案建议
#print("\n" + "="*50)
#print("处理重复细胞的建议方案")
#print("="*50)
#
#if len(true_duplicates) > 0:
#    print("1. 调查重复原因：检查数据生成流程中是否出现细胞标识分配错误")
#    print("2. 数据清理方案：")
#    print("   - 保留第一个出现值：adata.obs = adata.obs.drop_duplicates(subset='cell_id', keep='first')")
#    print("   - 保留最后一个出现值：adata.obs = adata.obs.drop_duplicates(subset='cell_id', keep='last')")
#    print("   - 手动检查决定：比较重复细胞的注释信息，决定保留哪个")
#    print("3. 如果使用索引作为细胞ID，可执行：adata.obs_names_make_unique() 确保索引唯一")
#else:
#    print("细胞标识唯一，无需进一步处理。")
