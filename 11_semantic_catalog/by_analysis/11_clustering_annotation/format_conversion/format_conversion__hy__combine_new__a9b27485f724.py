import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
from scipy.sparse import issparse, vstack

def merge_raw_matrices(file_list, output_file="merged_data_raw_preserved.h5ad"):
    """
    专门用于合并多个h5ad文件的raw计数矩阵，确保保留所有基因。
    
    参数:
    file_list: 需要合并的h5ad文件路径列表
    output_file: 输出文件名
    """
    
    print("开始合并h5ad文件的raw矩阵...")
    print(f"需要合并的文件数量: {len(file_list)}")
    
    raw_adata_list = []  # 用于存储从每个adata.raw提取的对象
    sample_names = []    # 样本标识
    all_genes_set = set() # 用于跟踪所有出现的基因
    
    # 1. 读取所有文件并检查raw矩阵
    print("\n步骤1: 检查每个文件的raw矩阵...")
    for i, file_path in enumerate(file_list):
        if not os.path.exists(file_path):
            print(f"警告: 文件 {file_path} 不存在，跳过")
            continue
            
        sample_name = os.path.basename(file_path).replace('.h5ad', '')
        sample_names.append(sample_name)
        print(f"检查文件 {i+1}/{len(file_list)}: {sample_name}")
        
        # 读取h5ad文件
        adata = sc.read_h5ad(file_path)
        
        # 检查是否存在raw矩阵
        if adata.raw is None:
            print(f"  ⚠️ 警告: {sample_name} 没有raw矩阵，将使用其X矩阵创建raw矩阵")
            # 如果没有raw，使用当前X创建
            adata.raw = adata
        
        # 提取raw矩阵对应的AnnData对象
        raw_adata = adata.raw.to_adata() if hasattr(adata.raw, 'to_adata') else adata.raw
        print(f"  ✅ Raw矩阵形状: {raw_adata.shape}")
        
        # 记录该样本的基因
        current_genes = set(raw_adata.var_names)
        all_genes_set.update(current_genes)
        print(f"  基因数: {len(current_genes)}")
        
        # 为细胞添加批次信息
        raw_adata.obs['batch'] = sample_name
        raw_adata.obs['sample'] = sample_name
        
        raw_adata_list.append(raw_adata)
    
    if len(raw_adata_list) == 0:
        raise ValueError("没有找到有效的h5ad文件进行合并")
    
    print(f"\n成功读取 {len(raw_adata_list)} 个文件的raw矩阵")
    print(f"所有文件中出现的总基因数（并集）: {len(all_genes_set)}")
    
    # 2. 合并raw矩阵（使用outer join保留所有基因）
    print("\n步骤2: 合并raw矩阵（保留所有基因）...")
    
    try:
        # 使用concat合并，join='outer'确保保留所有基因
        merged_raw = ad.concat(
            raw_adata_list,
            axis=0,              # 按行（细胞）合并
            join='outer',         # 关键参数：保留所有基因（并集）
            label='batch',        # 在obs中添加批次来源列
            keys=sample_names,    # 批次标识
            index_unique='-'     # 使细胞名称唯一
        )
        print("✅ Raw矩阵合并成功")
        
    except Exception as e:
        print(f"❌ Raw矩阵合并失败: {e}")
        raise
    
    # 3. 创建新的主AnnData对象（基于合并后的raw矩阵）
    print("\n步骤3: 创建最终的合并对象...")
    
    # 主X矩阵可以使用合并后raw矩阵的经过滤版本，或者直接使用raw矩阵
    # 这里我们创建一个与raw矩阵一致的主对象
    merged_adata = ad.AnnData(
        X=merged_raw.X.copy(),    # 表达矩阵
        obs=merged_raw.obs.copy(), # 细胞注释
        var=merged_raw.var.copy(), # 基因注释
        uns=merged_raw.uns.copy() if hasattr(merged_raw, 'uns') else {} # 非结构化数据
    )
    
    # 设置raw属性为合并后的raw矩阵（保持原始计数）
    merged_adata.raw = merged_raw
    print("✅ 最终对象创建完成")
    
    # 4. 保存结果
    print(f"\n步骤4: 保存结果到 {output_file}")
    merged_adata.write_h5ad(output_file)
    
    # 5. 输出详细统计
    print("\n" + "="*60)
    print("合并完成！详细统计信息:")
    print(f"总细胞数: {merged_adata.n_obs}")
    print(f"总基因数: {merged_adata.n_vars}")
    print(f"样本来源: {merged_adata.obs['batch'].unique().tolist()}")
    
    if merged_adata.raw is not None:
        print(f"Raw矩阵形状: {merged_adata.raw.X.shape}")
        print(f"Raw矩阵类型: {type(merged_adata.raw.X)}")
    
    # 检查每个批次的细胞数
    print("\n各样本细胞数量分布:")
    batch_counts = merged_adata.obs['batch'].value_counts()
    for batch, count in batch_counts.items():
        print(f"  {batch}: {count} 个细胞")
    
    print("="*60)
    
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
        "T.h5ad"
    ]
    
    # 检查文件是否存在
    existing_files = [f for f in file_list if os.path.exists(f)]
    missing_files = set(file_list) - set(existing_files)
    
    if missing_files:
        print("⚠️ 以下文件不存在，将合并存在的文件:")
        for f in missing_files:
            print(f"  - {f}")
    
    if not existing_files:
        print("❌ 错误: 没有找到任何h5ad文件！")
        exit(1)
    
    # 执行合并
    try:
        merged_data = merge_raw_matrices(
            file_list=existing_files,
            output_file="merged_data_with_raw.h5ad"
        )
        print("\n合并成功！输出文件: merged_data_with_raw.h5ad")
        
    except Exception as e:
        print(f"\n 合并过程中出现错误: {e}")
