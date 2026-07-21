import scanpy as sc
import os

# 1. 定义需要合并的h5ad文件列表
file_list = ["B_cells.celltype.h5ad","celltype_Plasmas.h5ad","Endothelials.celltype.h5ad","Epithelials.celltype.h5ad","Fibroblasts.celltype.h5ad","Macrophages.celltype.h5ad","Neutrophils.celltype.h5ad","T_cells.celltype.h5ad"]

# 2. 创建一个空列表来存储所有AnnData对象
adatas = []

# 3. 逐个读取h5ad文件，并进行关键预处理以保留原始矩阵
for file in file_list:
    print(f"正在读取并处理文件: {file}")
    
    # 提取细胞类型名称
    cell_type = file.replace('.celltype.h5ad', '')
    
    # 读取h5ad文件
    adata = sc.read_h5ad(file)
    
    # **关键步骤：在进行了任何可能修改X矩阵的操作之前，将当前状态（假定为原始计数）保存到 .raw 属性中[1,5,8](@ref)**
    # 这确保了后续分析始终能访问到原始计数矩阵[6](@ref)
    adata.raw = adata
    
    # 在obs中添加一列标识原始细胞类型
    adata.obs['original_celltype'] = cell_type
    
    # 为避免合并后细胞ID重复，为每个细胞的索引添加前缀
    adata.obs_names = [f"{cell_type}_{cell_id}" for cell_id in adata.obs_names]
    
    # 将处理后的AnnData对象添加到列表中
    adatas.append(adata)
    print(f"成功加载 {file}, 细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    print(f"  已为该数据集设置 .raw 属性。")

# 4. 合并所有AnnData对象
print("\n开始合并所有数据集...")
# 使用 `join='inner'` 确保只合并所有数据集共有的基因[4](@ref)
merged_adata = sc.concat(adatas, join='inner', label='source_file', index_unique=None)

print(f"合并完成！")
print(f"总细胞数: {merged_adata.n_obs}")
print(f"总基因数: {merged_adata.n_vars}")
print(f"合并后数据形状: {merged_adata.X.shape}")

# 5. **关键步骤：处理合并后对象的 .raw 属性**
# 合并操作不会自动处理每个单独对象的 .raw 属性[1](@ref)
# 因此，我们需要手动为合并后的对象创建一个新的 .raw 属性，其中包含所有细胞的原始计数。
print("\n正在处理合并后数据的原始矩阵（.raw 属性）...")

# 创建一个列表，用于存放从每个单独对象的 .raw 属性中提取的、经过相同过滤（基因交集）的原始数据
raw_adatas = []

for adata in adatas:
    # 获取当前数据集 .raw 中的 AnnData 视图
    raw_adata = adata.raw.to_adata() if adata.raw is not None else adata
    # 确保只取与合并后数据集中共同的基因[4](@ref)
    # merged_adata.var_names 包含了合并后保留的基因（交集）
    genes_to_keep = raw_adata.var_names.isin(merged_adata.var_names)
    raw_adata_filtered = raw_adata[:, genes_to_keep]
    raw_adatas.append(raw_adata_filtered)

# 现在，沿着观察轴（细胞）合并这些 .raw 数据，使用与主数据相同的连接方式和批次键[4](@ref)
merged_raw_adata = sc.concat(raw_adatas, join='inner', label='source_file', index_unique=None)

# 将合并后的原始数据赋值给最终合并对象的 .raw 属性
merged_adata.raw = merged_raw_adata
print("已成功设置合并后数据的 .raw 属性。")

# 6. 验证结果
print("\n=== 合并后数据验证 ===")
print(f"主数据矩阵 (merged_adata.X) 形状: {merged_adata.X.shape}")
if merged_adata.raw is not None:
    print(f"原始数据矩阵 (merged_adata.raw.X) 形状: {merged_adata.raw.X.shape}")
    # 检查维度是否匹配（细胞数应相同）
    if merged_adata.n_obs == merged_adata.raw.n_obs:
        print("✓ 验证通过：原始矩阵与主矩阵的细胞数量一致。")
    else:
        print("⚠️ 注意：原始矩阵与主矩阵的细胞数量不一致，请检查。")
else:
    print("⚠️ 警告：合并后对象的 .raw 属性未设置成功。")

print("\n合并后数据的obs列:", list(merged_adata.obs.columns))

# 7. 保存合并后的h5ad文件
output_filename = "merged_single_cell_data_with_raw.h5ad"
merged_adata.write_h5ad(output_filename)
print(f"\n合并后的数据（包含原始矩阵）已保存为: {output_filename}")
