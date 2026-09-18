import scanpy as sc
import pandas as pd
import numpy as np
import os

# 定义映射规则
mapping_rules = {
    "CXCL5+ Inf-stem-tumor": "CXCL5+ Inf stem tumor Epi",
    "Prol-Metab-Epi": "Prol Metab Epi", 
    "EMT-like Inv Ent": "EMT-like Inv Ent",
    "ESC-like Epi": "ESC-like Epi",
    "WNT+ ISC": "WNT+ ISC",
    "MUC5B+ CGC": "MUC5B+ CGC",
    "Hypo-Stem-tumor": "Hypo Stem tumor Epi",
    "Hypo-Metab": "Hypo Metab Epi",
    "Inf-Hyp Epi": "Inf Hyp Epi",
    "Early stress": "Early stress Epi",
    "Metab stress": "Metab stress Epi",
    "MUC5AC+ MGC": "MUC5AC+ MGC",
    "Mat Ent": "Mat Ent",
    "PCSK1+ AGC": "PCSK1+ AGC",
    "IGF2+ prolif-stem-tumor": "IGF2+ prolif stem tumor Epi"
}

# 1. 读取h5ad文件
print("读取h5ad文件...")
adata = sc.read_h5ad('/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/all_filter_modified.h5ad')
print(f"原始数据形状: {adata.shape}")
print(f"celltype列的数据类型: {type(adata.obs['celltype'])}")

# 获取原始的唯一值
original_celltypes = adata.obs['celltype'].unique()
print(f"celltype列的唯一值: {original_celltypes}")

# 2. 查看当前celltype的分布
print("\n当前celltype分布:")
celltype_counts = adata.obs['celltype'].value_counts()
print(celltype_counts)

# 3. 应用映射规则修改celltype
print("\n应用映射规则修改celltype...")

# 先确保celltype是字符串类型
if hasattr(adata.obs['celltype'], 'cat'):
    print("celltype是分类变量，转换为字符串进行修改...")
    adata.obs['celltype'] = adata.obs['celltype'].astype(str)

# 统计每个规则的修改数量
modification_counts = {}

# 创建反向映射，检查是否有重复映射
reverse_mapping = {}
for old_name, new_name in mapping_rules.items():
    if new_name in reverse_mapping:
        print(f"警告: 多个旧名称映射到同一个新名称 '{new_name}'")
    reverse_mapping[new_name] = old_name

# 应用映射规则
for old_name, new_name in mapping_rules.items():
    mask = adata.obs['celltype'] == old_name
    num_to_modify = mask.sum()
    
    if num_to_modify > 0:
        adata.obs.loc[mask, 'celltype'] = new_name
        modification_counts[old_name] = num_to_modify
        print(f"  将 '{old_name}' 改为 '{new_name}': 修改了 {num_to_modify} 个细胞")

# 4. 验证修改结果
print("\n验证修改结果...")
print(f"修改后celltype列的唯一值: {pd.unique(adata.obs['celltype'])}")

# 5. 查看修改后的分布
print("\n修改后celltype分布:")
celltype_counts_updated = adata.obs['celltype'].value_counts()
print(celltype_counts_updated)

# 6. 检查修改是否正确
print(f"\n修改结果统计:")
total_modified = sum(modification_counts.values())
print(f"  总共修改了 {total_modified} 个细胞")

# 检查是否有未映射的原始celltype
unmapped_original = set(original_celltypes) - set(mapping_rules.keys())
if unmapped_original:
    print(f"  未映射的原始celltype: {unmapped_original}")

# 检查是否有未使用的映射规则
unused_mappings = []
for old_name in mapping_rules.keys():
    if old_name not in original_celltypes:
        unused_mappings.append(old_name)

if unused_mappings:
    print(f"  未使用的映射规则: {unused_mappings}")


# 8. 保存修改后的数据
print("\n保存修改后的数据...")
output_path = '/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/all_filter_modified.h5ad'
adata.write_h5ad(output_path)
print(f"修改后的数据已保存到: {output_path}")

print("\n完成！")
