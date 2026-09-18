import scanpy as sc
import pandas as pd
import numpy as np

def unify_epithelial_naming(adata, original_names=None, target_name="Epithelials"):
    """
    统一上皮细胞类型的命名
    
    参数:
    adata: AnnData对象
    original_names: 需要修改的原始名称列表，如果为None则使用默认值
    target_name: 统一后的目标名称
    """
    if original_names is None:
        original_names = [f'Epithelials-{i}' for i in range(1, 5)]  # Epithelials-1 到 Epithelials-4
    
    # 检查celltype列是否存在
    if 'celltype' not in adata.obs.columns:
        raise ValueError("数据中未找到'celltype'列，请检查列名是否正确")
    
    print("修改前的细胞类型分布:")
    print(adata.obs['celltype'].value_counts())
    
    # 创建副本以避免修改原始数据
    adata_modified = adata.copy()
    
    # 使用replace方法批量替换细胞类型名称[7,8](@ref)
    for original_name in original_names:
        if original_name in adata_modified.obs['celltype'].values:
            adata_modified.obs['celltype'] = adata_modified.obs['celltype'].replace(
                original_name, target_name
            )
            print(f"已将 {original_name} 替换为 {target_name}")
        else:
            print(f"未找到细胞类型: {original_name}，跳过")
    
    print("\n修改后的细胞类型分布:")
    print(adata_modified.obs['celltype'].value_counts())
    
    return adata_modified

# 主程序
if __name__ == "__main__":
    # 请将文件路径替换为您的实际h5ad文件路径
    file_path = "adata_celltype.h5ad"
    
    try:
        # 读取h5ad文件[1,4](@ref)
        print(f"正在读取文件: {file_path}")
        adata = sc.read_h5ad(file_path)
        
        print("数据基本信息:")
        print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
        print(f"细胞类型列包含的种类: {adata.obs['celltype'].unique()}")
        
        # 统一上皮细胞命名
        adata_modified = unify_epithelial_naming(adata)
        
        # 保存修改后的数据（可选）
        output_path = file_path.replace('.h5ad', '_modified.h5ad')
        adata_modified.write(output_path)
        print(f"\n修改后的数据已保存至: {output_path}")
        
        # 显示修改摘要
        original_epithelial_count = adata.obs['celltype'].str.contains('Epithelials').sum()
        new_epithelial_count = (adata_modified.obs['celltype'] == 'Epithelials').sum()
        
        print(f"\n修改摘要:")
        print(f"原上皮细胞总数: {original_epithelial_count}")
        print(f"现上皮细胞总数: {new_epithelial_count}")
        
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}，请检查文件路径是否正确")
    except Exception as e:
        print(f"处理过程中发生错误: {str(e)}")
