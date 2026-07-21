import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def simple_merge_h5ad_files(file_list, output_file="merged_data.h5ad"):
    """
    简单合并多个h5ad文件，完全不修改原始数据
    
    参数:
    file_list: h5ad文件路径列表
    output_file: 合并后的输出文件路径
    """
    
    print("开始简单合并h5ad文件...")
    print(f"需要合并的文件数量: {len(file_list)}")
    
    # 检查文件是否存在
    existing_files = [f for f in file_list if os.path.exists(f)]
    if len(existing_files) != len(file_list):
        missing_files = set(file_list) - set(existing_files)
        print(f"警告: 以下文件不存在: {missing_files}")
        file_list = existing_files
    
    if len(file_list) == 0:
        raise ValueError("没有找到可用的h5ad文件进行合并")
    
    # 存储所有AnnData对象
    adata_list = []
    
    # 读取所有文件，保持原样
    for file_path in file_list:
        try:
            print(f"读取文件: {file_path}")
            
            # 读取h5ad文件
            adata = sc.read_h5ad(file_path)
            
            # 完全不修改任何内容，包括obs_names
            # 不添加任何新列，不修改现有列
            
            adata_list.append(adata)
            print(f"成功读取: {adata.shape[0]} 细胞, {adata.shape[1]} 基因")
            
        except Exception as e:
            print(f"读取文件 {file_path} 时出错: {e}")
            continue
    
    if len(adata_list) == 0:
        raise ValueError("未能成功读取任何h5ad文件")
    
    print(f"\n成功读取 {len(adata_list)} 个文件")
    
    # 最简单的合并方式
    print("开始合并数据...")
    
    try:
        # 使用最简单的concat方式，不添加任何额外信息
        merged_adata = ad.concat(
            adata_list,
            axis=0,
            join='inner',  # 只保留共同基因
            index_unique=None  # 不处理重复索引
        )
        
        print("合并完成")
        
    except Exception as e:
        print(f"合并失败: {e}")
        print("尝试备用合并方法...")
        
        # 备用方法: 最基础的手动合并
        X_list = []
        obs_list = []
        var = adata_list[0].var
        
        for adata in adata_list:
            X_list.append(adata.X)
            obs_list.append(adata.obs)
        
        # 合并obs，保持原样
        merged_obs = pd.concat(obs_list, axis=0)
        
        # 合并X矩阵
        from scipy.sparse import vstack
        if hasattr(adata_list[0].X, 'shape'):
            merged_X = vstack(X_list) if len(X_list) > 1 else X_list[0]
        else:
            merged_X = np.vstack(X_list) if len(X_list) > 1 else X_list[0]
        
        # 创建新的AnnData对象
        merged_adata = ad.AnnData(X=merged_X, obs=merged_obs, var=var)
    
    # 保存合并结果
    print(f"保存合并结果到: {output_file}")
    merged_adata.write_h5ad(output_file)
    
    # 输出最简统计信息
    print("\n" + "="*50)
    print("合并完成!")
    print(f"最终数据维度: {merged_adata.shape[0]} 细胞, {merged_adata.shape[1]} 基因")
    print(f"输出文件: {output_file}")
    print("="*50)
    
    return merged_adata

def main():
    """主函数"""
    # 定义您的h5ad文件列表
    h5ad_files = [
        "celltype_B_cells.h5ad",
        "celltype_Endothelials.h5ad", 
        "celltype_GREM1+_CAFs.h5ad",
        "celltype_Macrophages.h5ad",
        "celltype_OLFML3+_CAFs.h5ad",
        "celltype_Plasmas.h5ad",
        "celltype_T_cells.h5ad",
        "celltype_CAFs-EMT.h5ad",
        "celltype_Epithelials.h5ad",
        "celltype_iCAFs.h5ad", 
        "celltype_Neutrophils.h5ad",
        "celltype_Pericytes.h5ad",
        "celltype_SMC.h5ad"
    ]
    
    # 检查文件是否存在
    existing_files = [f for f in h5ad_files if os.path.exists(f)]
    
    if len(existing_files) == 0:
        print("错误: 没有找到任何h5ad文件")
        return
    
    print("找到以下文件进行合并:")
    for f in existing_files:
        print(f"  {f}")
    
    # 设置输出文件名
    output_file = "merged_all_celltypes.h5ad"
    
    # 执行合并
    try:
        merged_adata = simple_merge_h5ad_files(existing_files, output_file)
        print("合并成功完成!")
        
    except Exception as e:
        print(f"合并过程中发生错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
