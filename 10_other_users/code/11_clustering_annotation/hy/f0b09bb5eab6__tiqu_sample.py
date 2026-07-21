import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def extract_and_save_samples(h5ad_path, target_samples, matrix_type='X'):
    """
    提取指定样本并保存数据
    
    参数:
    h5ad_path: h5ad文件路径
    target_samples: 要提取的样本列表
    matrix_type: 要保存的矩阵类型，可以是'X', 'raw', 或layers中的键名
    """
    
    # 读取数据
    adata = sc.read_h5ad(h5ad_path)
    
    output_dir = "extracted_samples"
    Path(output_dir).mkdir(exist_ok=True)
    
    for sample in target_samples:
        if sample not in adata.obs['sample'].values:
            print(f"警告: 样本 {sample} 不存在")
            continue
            
        # 提取样本
        sample_mask = adata.obs['sample'] == sample
        adata_sample = adata[sample_mask, :].copy()
        
        # 根据矩阵类型提取表达矩阵
        if matrix_type == 'X':
            expr_matrix = adata_sample.X
        elif matrix_type == 'raw' and adata_sample.raw is not None:
            expr_matrix = adata_sample.raw.X
        elif matrix_type in adata_sample.layers:
            expr_matrix = adata_sample.layers[matrix_type]
        else:
            print(f"警告: 矩阵类型 '{matrix_type}' 不可用，使用默认的.X")
            expr_matrix = adata_sample.X
        
        # 保存h5ad
        h5ad_path = f"{output_dir}/{sample}.h5ad"
        adata_sample.write_h5ad(h5ad_path)
        
        # 保存表达矩阵
        if hasattr(expr_matrix, 'toarray'):
            # 稀疏矩阵
            sparse_matrix = expr_matrix
            dense_matrix = expr_matrix.toarray()
            
            np.save(f"{output_dir}/{sample}_{matrix_type}_dense.npy", dense_matrix)
            
            from scipy import sparse
            sparse.save_npz(f"{output_dir}/{sample}_{matrix_type}_sparse.npz", sparse_matrix)
        else:
            # 密集矩阵
            np.save(f"{output_dir}/{sample}_{matrix_type}.npy", expr_matrix)
        
        print(f"样本 {sample}: {adata_sample.shape[0]} 个细胞, 已保存")
    
    return output_dir

# 使用示例
if __name__ == "__main__":
    h5ad_path = "data_sp_ep.h5ad"
    target_samples = ["SR202500003B", "SR2025000016J", "SR2025000016K"]
    
    # 提取并保存，默认使用.X矩阵
    output_dir = extract_and_save_samples(h5ad_path, target_samples, matrix_type='X')
    
    # 如果数据有raw矩阵，可以这样提取:
    # output_dir = extract_and_save_samples(h5ad_path, target_samples, matrix_type='raw')
    
    print(f"\n处理完成! 文件保存在: {output_dir}")
