import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
from scipy.sparse import issparse, vstack

def merge_h5ad_files_preserve_sample(file_list, output_file="merged_data.h5ad"):
    """
    合并多个h5ad文件，保留原有的sample列，并添加文件名作为ct_ori列
    """
    
    print("开始合并h5ad文件...")
    print(f"需要合并的文件数量: {len(file_list)}")
    
    adata_list = []
    
    # 1. 读取所有文件
    for i, file_path in enumerate(file_list):
        if not os.path.exists(file_path):
            print(f"警告: 文件 {file_path} 不存在，跳过")
            continue
            
        # 提取文件名（不含扩展名）作为ct_ori的值
        file_basename = os.path.basename(file_path)
        ct_ori_value = file_basename.replace('.h5ad', '')
        
        print(f"读取文件 {i+1}/{len(file_list)}: {file_path}")
        print(f"CT来源标识: {ct_ori_value}")
        
        # 读取h5ad文件
        adata = sc.read_h5ad(file_path)
        
        # 检查是否已有sample列
        if 'sample' in adata.obs.columns:
            print(f"  检测到原有sample列，样本值示例: {adata.obs['sample'].iloc[0] if len(adata) > 0 else 'N/A'}")
        else:
            print(f"  警告: 未检测到sample列")
        
        # 添加ct_ori列（文件名）
        adata.obs['ct_ori'] = ct_ori_value
        
        # 确保细胞名称唯一
        adata.obs_names = [f"{ct_ori_value}_{name}" for name in adata.obs_names]
        
        adata_list.append(adata)
    
    if len(adata_list) == 0:
        raise ValueError("没有找到有效的h5ad文件进行合并")
    
    print(f"\n成功读取 {len(adata_list)} 个文件")
    
    # 2. 合并数据（使用outer join保留所有基因）[1,2](@ref)
    print("\n开始合并数据...")
    
    try:
        merged_adata = ad.concat(
            adata_list,
            axis=0,
            join='outer',  # 关键参数：保留所有基因[1](@ref)
            label='batch',
            keys=[os.path.basename(f).replace('.h5ad', '') for f in file_list if os.path.exists(f)],
            index_unique='-'
        )
        print("✅ 数据合并成功")
        
    except Exception as e:
        print(f"❌ 数据合并失败: {e}")
        # 备用合并方案
        merged_adata = manual_merge_with_ct_ori(adata_list, file_list)
    
    # 3. 检查合并结果
    print(f"\n合并后数据形状: {merged_adata.shape}")
    print(f"观察列(obs)包含: {list(merged_adata.obs.columns)}")
    
    # 检查sample列和ct_ori列
    if 'sample' in merged_adata.obs.columns:
        sample_counts = merged_adata.obs['sample'].value_counts()
        print(f"\n原有sample列的样本分布:")
        for sample, count in sample_counts.head().items():
            print(f"  {sample}: {count} cells")
    
    if 'ct_ori' in merged_adata.obs.columns:
        ct_ori_counts = merged_adata.obs['ct_ori'].value_counts()
        print(f"\nct_ori列（文件来源）分布:")
        for ct_ori, count in ct_ori_counts.items():
            print(f"  {ct_ori}: {count} cells")
    
    # 4. 保存结果
    print(f"\n保存结果到: {output_file}")
    merged_adata.write_h5ad(output_file)
    
    return merged_adata

def manual_merge_with_ct_ori(adata_list, file_list):
    """手动合并函数（备用方案）"""
    
    print("使用手动合并方案...")
    
    # 获取所有基因的并集[1](@ref)
    all_genes = set()
    for adata in adata_list:
        all_genes.update(adata.var_names)
    
    all_genes = sorted(list(all_genes))
    print(f"总基因数（并集）: {len(all_genes)}")
    
    # 合并矩阵
    merged_matrices = []
    merged_obs = []
    
    for i, (adata, file_path) in enumerate(zip(adata_list, file_list)):
        ct_ori_value = os.path.basename(file_path).replace('.h5ad', '')
        
        # 重新索引以匹配所有基因
        adata_reindexed = adata[:, all_genes].copy()
        
        merged_matrices.append(adata_reindexed.X)
        merged_obs.append(adata_reindexed.obs)
    
    # 创建合并后的对象
    merged_x = vstack(merged_matrices)
    merged_obs = pd.concat(merged_obs)
    
    merged_adata = ad.AnnData(
        X=merged_x,
        obs=merged_obs,
        var=pd.DataFrame(index=all_genes)
    )
    
    return merged_adata

# 使用示例
if __name__ == "__main__":
    # 您的文件列表
    file_list = [
        "B.h5ad",
        "Endo.h5ad", 
        "Epi.h5ad",
        "Macro.h5ad", 
        "Neutrophils.h5ad",
        "Pericytes.h5ad",
        "Plasma.h5ad", 
        "T.h5ad",
	"CAF.h5ad",
	"SMC.h5ad"
    ]   
 
    # 检查文件是否存在
    existing_files = [f for f in file_list if os.path.exists(f)]
    missing_files = set(file_list) - set(existing_files)
    
    if missing_files:
        print("⚠️ 以下文件不存在:")
        for f in missing_files:
            print(f"  - {f}")
    
    if not existing_files:
        print("❌ 错误: 没有找到任何h5ad文件！")
        exit(1)
    
    # 执行合并
    try:
        merged_data = merge_h5ad_files_preserve_sample(
            file_list=existing_files,
            output_file="merged_data_with_ct_ori.h5ad"
        )
        
        print("\n 合并完成！")
        print("列信息:")
        print(f"  - sample列: 保留原有样本信息")
        print(f"  - ct_ori列: 新增的文件来源标识")
        print(f"  - batch列: 合并时自动添加的批次信息")
        
    except Exception as e:
        print(f"\n 合并过程中出现错误: {e}")
