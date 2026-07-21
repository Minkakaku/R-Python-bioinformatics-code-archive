import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

def merge_h5ad_files(file1, file2, output_file, batch_key='batch'):
    """
    合并两个h5ad文件，保留共同基因并添加批次信息
    
    参数:
    file1: 第一个h5ad文件路径
    file2: 第二个h5ad文件路径
    output_file: 合并后的输出文件路径
    batch_key: obs中用于存储批次信息的列名
    """
    # 1. 读取数据
    print(f"读取文件: {file1}")
    adata1 = ad.read_h5ad(file1)
    print(f"读取文件: {file2}")
    adata2 = ad.read_h5ad(file2)
    
    # 2. 检查并统一基因名
    # 确保var_names一致（大小写敏感）
    adata1.var_names = adata1.var_names.astype(str)
    adata2.var_names = adata2.var_names.astype(str)
    
    # 3. 找到共同基因
    common_genes = list(set(adata1.var_names) & set(adata2.var_names))
    print(f"共同基因数量: {len(common_genes)}")
    
    if len(common_genes) == 0:
        raise ValueError("两个数据集没有共同基因，无法合并")
    
    # 4. 筛选共同基因
    adata1 = adata1[:, common_genes].copy()
    adata2 = adata2[:, common_genes].copy()
    
    # 5. 添加批次信息
    adata1.obs[batch_key] = 'batch1'
    adata2.obs[batch_key] = 'batch2'
    
    # 6. 合并数据集
    print("合并数据集...")
    merged_adata = ad.concat([adata1, adata2], 
                            join='outer', 
                            merge='same', 
                            index_unique=None)
    
    # 7. 保存合并后的数据
    print(f"保存合并数据到: {output_file}")
    merged_adata.write_h5ad(output_file)
    print("合并完成!")
    return merged_adata

# 使用示例
if __name__ == "__main__":
    # 替换为你的实际文件路径
    file1 = "adata_cluster.h5ad"
    file2 = "celltype_TNK.h5ad"
    output = "merged_data.h5ad"
    
    merged_data = merge_h5ad_files(file1, file2, output)
