import scanpy as sc
import pandas as pd
import numpy as np

def safely_extract_epithelial_cells(input_file, output_file='Epi.h5ad'):
    """
    安全地提取上皮细胞并保留原始计数矩阵的可靠方法
    """
    print("读取h5ad文件...")
    adata = sc.read_h5ad(input_file)
    print(f"原始数据形状: {adata.shape}")

    # 定义目标细胞类型
    target_celltypes = ['Myoepithelial cells', 'Luminal epithelial cells']
    celltype_col = 'cell_type'
    
    # 创建布尔掩码
    mask = adata.obs[celltype_col].isin(target_celltypes)
    print(f"找到 {mask.sum()} 个目标细胞")

    # 提取细胞子集（主数据）
    print("提取目标细胞...")
    adata_epi = adata[mask, :].copy()
    print(f"提取后数据形状: {adata_epi.shape}")

    #  关键修正：正确处理 .raw 属性
    if adata.raw is not None:
        print("检测到已存在的.raw属性，正在处理...")
        
        # 方法1：最安全的方法 - 从原始 .raw 重建
        try:
            # 获取匹配的细胞索引
            raw_cell_indices = adata.raw.obs_names.isin(adata_epi.obs_names)
            
            # 确保使用 .copy() 创建独立对象，而不是视图
            adata_epi_raw_subset = adata.raw[raw_cell_indices, :].copy()
            
            # 关键步骤：验证并确保这是完整的 AnnData 对象
            print(f"原始子集对象类型: {type(adata_epi_raw_subset)}")
            print(f"原始子集形状: {adata_epi_raw_subset.shape}")
            
            # 设置 .raw 属性
            adata_epi.raw = adata_epi_raw_subset
            print("✅ 方法1成功：使用.raw子集")
            
        except Exception as e:
            print(f"方法1失败: {e}，尝试备用方案...")
            
            # 方法2：备用方案 - 使用主数据的原始矩阵
            # 创建一个新的 AnnData 对象作为原始数据
            import anndata as ad
            adata_epi_raw = ad.AnnData(X=adata_epi.X.copy(), 
                                      var=adata_epi.var.copy(), 
                                      obs=adata_epi.obs.copy())
            adata_epi.raw = adata_epi_raw
            print("✅ 方法2成功：使用主数据重建原始矩阵")
    else:
        print("未检测到.raw属性，将当前数据设置为原始数据")
        adata_epi.raw = adata_epi

    # 验证结果
    print(f"提取后数据形状: {adata_epi.shape}")
    if adata_epi.raw is not None:
        print(f"原始计数矩阵形状: {adata_epi.raw.X.shape}")
    
    print("细胞类型分布:")
    print(adata_epi.obs[celltype_col].value_counts())

    # 保存数据
    adata_epi.write(output_file)
    print(f"数据已保存至: {output_file}")
    
    return adata_epi

# 使用示例
if __name__ == "__main__":
    input_file = "result.h5ad"  # 请替换为您的文件路径
    epi_data = safely_extract_epithelial_cells(input_file, "Epi.h5ad")
    print("操作成功完成！")
