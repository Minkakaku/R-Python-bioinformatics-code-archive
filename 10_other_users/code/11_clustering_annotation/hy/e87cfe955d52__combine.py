import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
from scipy.sparse import csr_matrix, vstack

def merge_h5ad_files(file_list, output_file="merged_data.h5ad"):
    """
    合并多个h5ad文件，保留raw计数矩阵
    
    参数:
    file_list: h5ad文件路径列表
    output_file: 输出文件名
    """
    
    print("开始合并h5ad文件...")
    print(f"需要合并的文件数量: {len(file_list)}")
    
    adata_list = []
    sample_names = []
    
    # 1. 读取所有文件
    for i, file_path in enumerate(file_list):
        if not os.path.exists(file_path):
            print(f"警告: 文件 {file_path} 不存在，跳过")
            continue
            
        # 提取样本名（从文件名）
        sample_name = os.path.basename(file_path).replace('.h5ad', '')
        sample_names.append(sample_name)
        
        print(f"读取文件 {i+1}/{len(file_list)}: {file_path}")
        
        # 读取h5ad文件
        adata = sc.read_h5ad(file_path)
        
        # 确保有raw矩阵
        if adata.raw is None:
            print(f"警告: {sample_name} 没有raw矩阵，使用X矩阵")
            # 如果没有raw属性，将当前X矩阵设置为raw
            adata.raw = adata
        else:
            print(f"样本 {sample_name} 包含raw矩阵，形状: {adata.raw.X.shape}")
        
        # 为每个细胞添加样本来源标识
        adata.obs['batch'] = sample_name
        adata.obs['sample'] = sample_name
        
        adata_list.append(adata)
    
    if len(adata_list) == 0:
        raise ValueError("没有找到有效的h5ad文件进行合并")
    
    print(f"\n成功读取 {len(adata_list)} 个文件")
    
    # 2. 合并数据
    print("\n开始合并数据...")
    
    # 方法1: 使用concat函数（推荐）
    try:
        # 设置连接轴和合并方式
        merged_adata = ad.concat(
            adata_list,
            axis=0,
            join='outer',
            label='batch',
            keys=sample_names,
            index_unique='-',  # 使索引唯一
            fill_value=0       # 对于不同基因集的填充值
        )
        
        print("使用ad.concat合并成功")
        
    except Exception as e:
        print(f"方法1失败: {e}")
        print("尝试使用方法2: 手动合并...")
        
        # 方法2: 手动合并（备用方案）
        merged_adata = manual_merge_anndata(adata_list, sample_names)
    
    # 3. 处理raw矩阵
    print("\n处理raw计数矩阵...")
    
    # 检查是否所有样本都有raw矩阵
    has_raw_all = all(adata.raw is not None for adata in adata_list)
    
    if has_raw_all:
        print("所有样本都有raw矩阵，开始合并raw矩阵...")
        
        # 提取所有raw矩阵
        raw_matrices = []
        var_dfs = []
        
        for i, adata in enumerate(adata_list):
            raw_matrix = adata.raw.X
            raw_var = adata.raw.var
            
            # 确保矩阵为稀疏格式以节省内存
            if not hasattr(raw_matrix, 'shape'):
                raw_matrix = csr_matrix(raw_matrix)
            
            raw_matrices.append(raw_matrix)
            var_dfs.append(raw_var)
        
        # 合并var（基因信息）
        merged_raw_var = pd.concat(var_dfs).drop_duplicates()
        
        # 手动合并raw矩阵
        try:
            # 创建合并后的raw矩阵
            merged_raw_x = vstack(raw_matrices)
            
            # 创建raw的AnnData对象
            raw_merged = ad.AnnData(
                X=merged_raw_x,
                var=merged_raw_var,
                obs=merged_adata.obs
            )
            
            # 设置raw属性
            merged_adata.raw = raw_merged
            print("Raw矩阵合并成功")
            
        except Exception as e:
            print(f"Raw矩阵合并失败: {e}")
            print("将使用X矩阵作为raw矩阵")
            merged_adata.raw = merged_adata
    else:
        print("部分样本缺少raw矩阵，使用X矩阵作为raw矩阵")
        merged_adata.raw = merged_adata
    
    # 4. 保存合并后的数据
    print(f"\n保存合并后的数据到: {output_file}")
    merged_adata.write_h5ad(output_file)
    
    # 5. 输出统计信息
    print("\n" + "="*50)
    print("合并完成！统计信息:")
    print(f"总细胞数: {merged_adata.n_obs}")
    print(f"总基因数: {merged_adata.n_vars}")
    print(f"样本批次: {merged_adata.obs['batch'].unique().tolist()}")
    
    if merged_adata.raw is not None:
        print(f"Raw矩阵形状: {merged_adata.raw.X.shape}")
    
    print(f"合并后数据形状: {merged_adata.X.shape}")
    print("="*50)
    
    return merged_adata

def manual_merge_anndata(adata_list, sample_names):
    """手动合并AnnData对象（备用方法）"""
    
    # 获取所有基因的并集
    all_genes = set()
    for adata in adata_list:
        all_genes.update(adata.var_names)
    
    all_genes = sorted(list(all_genes))
    print(f"发现总基因数: {len(all_genes)}")
    
    # 创建新的var DataFrame
    merged_var = pd.DataFrame(index=all_genes)
    
    # 合并X矩阵
    merged_matrices = []
    merged_obs = []
    
    for i, (adata, sample_name) in enumerate(zip(adata_list, sample_names)):
        # 重新索引以匹配所有基因
        adata_reindexed = adata[:, all_genes].copy()
        
        # 处理缺失的基因（填充0）
        missing_genes = set(all_genes) - set(adata.var_names)
        if missing_genes:
            # 创建填充矩阵
            fill_matrix = csr_matrix((adata.n_obs, len(missing_genes)))
            
            # 水平拼接
            from scipy.sparse import hstack
            adata_reindexed.X = hstack([adata_reindexed.X, fill_matrix])
        
        merged_matrices.append(adata_reindexed.X)
        
        # 合并obs数据
        obs_data = adata.obs.copy()
        obs_data['original_batch'] = sample_name
        merged_obs.append(obs_data)
    
    # 垂直拼接所有矩阵
    merged_x = vstack(merged_matrices)
    merged_obs = pd.concat(merged_obs)
    
    # 创建合并后的AnnData对象
    merged_adata = ad.AnnData(
        X=merged_x,
        obs=merged_obs,
        var=merged_var
    )
    
    return merged_adata

# 使用示例
if __name__ == "__main__":
    # 你的文件列表
    file_list = [
        "B.h5ad",
        "Endo.h5ad", 
        "Epi.h5ad",
        "Macro.h5ad",
        "Neutrophils.h5ad",
        "Pericytes.h5ad",
        "Plasma.h5ad",
        "T.h5ad"
    ]
    
    # 检查文件是否存在
    existing_files = [f for f in file_list if os.path.exists(f)]
    missing_files = set(file_list) - set(existing_files)
    
    if missing_files:
        print("警告: 以下文件不存在:")
        for f in missing_files:
            print(f"  - {f}")
        print("将继续合并存在的文件...")
    
    if not existing_files:
        print("错误: 没有找到任何h5ad文件！")
        exit(1)
    
    # 合并文件
    merged_data = merge_h5ad_files(
        file_list=existing_files,
        output_file="merged_sc_data.h5ad"
    )
    
    print("合并完成！输出文件: merged_sc_data.h5ad")
