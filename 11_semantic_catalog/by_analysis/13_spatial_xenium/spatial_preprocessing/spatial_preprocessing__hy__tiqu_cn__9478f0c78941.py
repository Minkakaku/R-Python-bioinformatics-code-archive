#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import sys

def extract_cn7_cn8_simple():
    """简化版本，直接提取CN7和CN8细胞"""
    
    # 输入输出文件
    input_file = "adata_filtered_with_MPCpercent.h5ad"
    output_file = "CN9_cells.h5ad"
    
    print(f"加载数据: {input_file}")
    adata = sc.read_h5ad(input_file)
    
    print(f"原始数据形状: {adata.shape}")
    
    # 检查CNs_10列
    if 'CNs_10' not in adata.obs.columns:
        print("错误: 没有找到CNs_10列")
        return None
    
    # 查看CN分布
    print("\n原始CN分布:")
    cn_counts = adata.obs['CNs_10'].value_counts()
    for cn, count in cn_counts.items():
        print(f"  {cn}: {count}")
    
    # 提取CN7和CN8细胞
    cn_mask = adata.obs['CNs_10'].isin(['CN9'])
    adata_cn = adata[cn_mask].copy()
    
    print(f"\n提取后数据形状: {adata_cn.shape}")
    print(f"CN7细胞数: {cn_counts.get('CN9', 0)}")
    
    # 检查raw矩阵
    if adata.raw is not None:
        print(f"\n原始数据有raw矩阵")
        if adata_cn.raw is not None:
            print(f"提取后数据保留raw矩阵")
        else:
            print(f"警告: 提取后数据没有raw矩阵")
    else:
        print(f"\n原始数据没有raw矩阵")
    
    # 保存数据
    print(f"\n保存到: {output_file}")
    adata_cn.write_h5ad(output_file)
    
    return adata_cn

if __name__ == "__main__":
    adata_result = extract_cn7_cn8_simple()
    
    if adata_result is not None:
        print("\n处理完成!")
    else:
        print("\n处理失败!")
        sys.exit(1)
