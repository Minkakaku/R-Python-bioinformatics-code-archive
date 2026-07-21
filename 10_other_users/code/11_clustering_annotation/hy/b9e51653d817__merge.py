import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# 列出需要合并的h5ad文件
file_list = [
   "celltype_B_cells.h5ad",
   "celltype_DCs.h5ad",
   "celltype_Endothelials.h5ad",
   "celltype_Epi.h5ad",
   "celltype_Fibroblasts.h5ad",
   "celltype_Macrophages.h5ad",
   "celltype_Masts.h5ad",
   "celltype_Pericytes.h5ad",
   "celltype_Plasmas.h5ad",
   "celltype_T_cells.h5ad"
]

# 存储读取的AnnData对象和它们的raw矩阵
adatas = []
raw_matrices = []

print("开始读取和预处理数据...")

for i, file in enumerate(file_list):
    print(f"正在处理: {file}")
    
    # 读取数据
    adata = sc.read_h5ad(file)
    
    # 确保每个细胞的索引唯一，避免合并错误[6](@ref)
    adata.obs_names_make_unique()
    
    # 保存raw矩阵（如果存在）
    if adata.raw is not None:
        raw_matrices.append(adata.raw.to_adata())
        print(f"  - 检测到raw矩阵，形状: {adata.raw.X.shape}")
    else:
        raw_matrices.append(adata.copy())
        print(f"  - 未检测到raw矩阵，使用主矩阵")
    
    # 为每个数据集添加批次标识（可选，便于后续追踪）
    adata.obs['batch'] = f'batch_{i+1}'
    
    adatas.append(adata)
    print(f"  - 成功加载，细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")

print("\n开始合并数据...")

# 合并主AnnData对象[1](@ref)
# 使用merge='same'保留所有数据集共有的结构[2](@ref)
combined_adata = ad.concat(adatas, join='inner', merge='same')

print(f"合并后数据形状: {combined_adata.n_obs} 个细胞, {combined_adata.n_vars} 个基因")

# 处理raw矩阵的合并
print("\n处理raw矩阵合并...")
try:
    # 提取每个raw矩阵的AnnData对象
    raw_adatas = []
    for i, raw_adata in enumerate(raw_matrices):
        raw_adata.obs_names_make_unique()
        raw_adata.obs['batch'] = f'batch_{i+1}'
        raw_adatas.append(raw_adata)
    
    # 合并raw矩阵[1](@ref)
    combined_raw = ad.concat(raw_adatas, join='inner', merge='same')
    
    # 将合并后的raw矩阵关联到主对象
    combined_adata.raw = combined_raw
    print("成功合并raw矩阵")
    
except Exception as e:
    print(f"Raw矩阵合并出现警告: {e}")
    print("将使用主矩阵作为raw矩阵")
    combined_adata.raw = combined_adata

# 检查obs列是否完整保留
print(f"\n合并后的obs列: {list(combined_adata.obs.columns)}")

# 验证合并结果
print("\n验证合并结果:")
print(f"- 总细胞数: {combined_adata.n_obs}")
print(f"- 总基因数: {combined_adata.n_vars}") 
print(f"- 批次分布:")
print(combined_adata.obs['batch'].value_counts())

# 保存合并后的数据
output_file = 'Epithelials.h5ad'
combined_adata.write(output_file)
print(f"\n合并完成! 结果已保存为: {output_file}")

# 可选：显示详细摘要
print("\n数据摘要:")
print(combined_adata)
