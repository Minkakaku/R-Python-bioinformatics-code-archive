import scanpy as sc
import pandas as pd
import numpy as np

# 读取h5ad文件
adata = sc.read_h5ad('modify_tumor.h5ad')  # 替换为您的文件路径

# 检查anno_new列是否存在
if 'anno_new' not in adata.obs.columns:
    print("错误：adata.obs中没有找到'anno_new'列")
else:
    # 查看替换前的NaN值数量
    nan_count_before = adata.obs['anno_new'].isna().sum()
    print(f"替换前，anno_new列中的NaN值数量: {nan_count_before}")
    
    # 将所有NaN值替换为"LYVE1+_TAM"
    adata.obs['anno_new'] = adata.obs['anno_new'].fillna('LYVE1+_TAM')
    
    # 查看替换后的值分布
    nan_count_after = adata.obs['anno_new'].isna().sum()
    print(f"替换后，anno_new列中的NaN值数量: {nan_count_after}")
    
    # 查看替换后的值分布
    print("\nanno_new列替换后的值分布：")
    value_counts = adata.obs['anno_new'].value_counts()
    print(value_counts)
    
    # 特别查看"LYVE1+_TAM"的数量
    if 'LYVE1+_TAM' in adata.obs['anno_new'].values:
        lyve1_count = (adata.obs['anno_new'] == 'LYVE1+_TAM').sum()
        print(f"\nLYVE1+_TAM的数量: {lyve1_count}")
    
    # 保存修改后的数据
    adata.write('modify_tumor_new.h5ad')
    print("\n修改完成，已保存为新的h5ad文件")
