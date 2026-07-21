import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack
import os

# 定义需要合并的h5ad文件列表
file_list = ["B_cells.celltype.h5ad","celltype_Plasmas.h5ad","Endothelials.celltype.h5ad","Epithelials.celltype.h5ad","Fibroblasts.celltype.h5ad","Macrophages.celltype.h5ad","Neutrophils.celltype.h5ad","T_cells.celltype.h5ad"]

print("开始合并多个h5ad文件的原始矩阵...")

# 步骤1: 获取完整的基因列表（所有数据集中出现的基因并集）
print("\n1. 收集所有数据集的基因列表...")
all_genes_set = set()

for file in file_list:
    if os.path.exists(file):
        adata = sc.read_h5ad(file)
        # 优先使用raw矩阵中的基因，如果没有raw则使用当前矩阵的基因
        if adata.raw is not None:
            genes = set(adata.raw.var_names)
        else:
            genes = set(adata.var_names)
        all_genes_set.update(genes)
        print(f"   {file}: 贡献了 {len(genes)} 个基因")

full_gene_list = sorted(list(all_genes_set))
print(f"   总基因数: {len(full_gene_list)}")

# 步骤2: 构建基因索引映射
gene_to_index = {gene: idx for idx, gene in enumerate(full_gene_list)}

# 步骤3: 逐个处理每个文件，重建完整的表达矩阵
print("\n2. 处理每个文件的原始矩阵...")
all_adatas = []
all_obs_dfs = []

for file in file_list:
    if not os.path.exists(file):
        print(f"   ⚠️ 文件不存在: {file}，跳过")
        continue
        
    print(f"   正在处理: {file}")
    adata = sc.read_h5ad(file)
    cell_type = file.replace('.celltype.h5ad', '')
    
    # 确定使用哪个矩阵作为原始计数
    if adata.raw is not None:
        source_adata = adata.raw.to_adata()
        matrix_source = "raw矩阵"
    else:
        source_adata = adata
        matrix_source = "主矩阵"
    
    print(f"      使用{matrix_source}, 细胞数: {source_adata.n_obs}, 基因数: {source_adata.n_vars}")
    
    # 重建完整的表达矩阵（包含所有基因）
    n_cells = source_adata.n_obs
    n_genes_total = len(full_gene_list)
    
    # 创建新的稀疏矩阵
    new_X = csr_matrix((n_cells, n_genes_total), dtype=source_adata.X.dtype)
    
    # 将现有基因的数据填充到新矩阵的正确位置
    for i, gene_name in enumerate(source_adata.var_names):
        if gene_name in gene_to_index:
            target_col_idx = gene_to_index[gene_name]
            # 复制该基因的表达数据到新矩阵
            new_X[:, target_col_idx] = source_adata.X[:, i]
        else:
            print(f"      警告: 基因 {gene_name} 不在总基因列表中")
    
    # 创建新的var（基因注释）
    new_var = pd.DataFrame(index=full_gene_list)
    
    # 保留原始obs信息，并添加来源标识
    new_obs = source_adata.obs.copy()
    new_obs['original_file'] = file
    new_obs['cell_type_source'] = cell_type
    
    # 为细胞ID添加前缀以避免重复
    new_obs_names = [f"{cell_type}_{cell_id}" for cell_id in source_adata.obs_names]
    
    # 创建新的AnnData对象
    new_adata = sc.AnnData(
        X=new_X,
        obs=new_obs,
        var=new_var
    )
    new_adata.obs_names = new_obs_names
    
    # 保留其他可能重要的信息
    if hasattr(source_adata, 'obsm') and source_adata.obsm:
        new_adata.obsm = source_adata.obsm.copy()
    if hasattr(source_adata, 'uns') and source_adata.uns:
        new_adata.uns = source_adata.uns.copy()
    
    all_adatas.append(new_adata)
    print(f"      完成处理 {file}")

# 步骤4: 合并所有重建后的AnnData对象
print("\n3. 合并所有数据集...")
if len(all_adatas) == 0:
    print("⚠️ 没有成功读取任何文件，请检查文件路径")
else:
    # 使用concat合并，join='outer'确保保留所有基因
    merged_adata = sc.concat(all_adatas, join='outer', label='batch', index_unique=None)
    
    print(f"合并完成！")
    print(f"总细胞数: {merged_adata.n_obs}")
    print(f"总基因数: {merged_adata.n_vars}")
    print(f"数据形状: {merged_adata.X.shape}")
    
    # 步骤5: 验证结果
    print("\n4. 验证合并结果...")
    print("obs列信息:", list(merged_adata.obs.columns))
    print("细胞来源分布:")
    print(merged_adata.obs['original_file'].value_counts())
    
    # 检查稀疏矩阵的非零元素比例
    nnz_ratio = merged_adata.X.nnz / (merged_adata.n_obs * merged_adata.n_vars)
    print(f"矩阵稀疏度: {nnz_ratio:.4%}")
    
    # 步骤6: 保存合并后的数据
    output_filename = "merged_all_celltypes_complete.h5ad"
    merged_adata.write_h5ad(output_filename)
    print(f"\n5. 合并后的数据已保存为: {output_filename}")
    
    # 最终统计信息
    print("\n=== 最终统计 ===")
    print(f"处理文件数: {len([f for f in file_list if os.path.exists(f)])}")
    print(f"成功合并文件数: {len(all_adatas)}")
    print(f"总基因数量: {merged_adata.n_vars}")
    print(f"总细胞数量: {merged_adata.n_obs}")
    
    # 检查是否保留了重要的obs信息
    if 'original_file' in merged_adata.obs.columns:
        print("✓ 成功保留了细胞来源信息")
    if 'cell_type_source' in merged_adata.obs.columns:
        print("✓ 成功保留了细胞类型来源信息")

print("\n合并流程完成！")
