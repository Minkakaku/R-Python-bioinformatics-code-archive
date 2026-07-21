import scanpy as sc
import pandas as pd
import numpy as np
import os

# 1. 读取h5ad文件
print("读取h5ad文件...")
adata = sc.read_h5ad('/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/all_filter.h5ad')
print(f"原始数据形状: {adata.shape}")
print(f"celltype列的数据类型: {type(adata.obs['celltype'])}")
print(f"celltype列的唯一值: {adata.obs['celltype'].unique()}")

# 2. 查看当前celltype的分布
print("\n当前celltype分布:")
celltype_counts = adata.obs['celltype'].value_counts()
print(celltype_counts)

# 3. 修改Tumor Epi为ESC-like Epi
print("\n修改'Tumor Epi'为'ESC-like Epi'...")

# 方法1：将分类列转换为字符串，修改后再转回分类
if hasattr(adata.obs['celltype'], 'cat'):
    print("celltype是分类变量，转换为字符串进行修改...")
    adata.obs['celltype'] = adata.obs['celltype'].astype(str)

# 统计修改数量
mask = adata.obs['celltype'] == 'Tumor Epi'  # 注意这里是大写T
num_to_modify = mask.sum()
print(f"找到 {num_to_modify} 个'Tumor Epi'细胞需要修改")

# 执行修改
adata.obs.loc[mask, 'celltype'] = 'ESC-like Epi'

# 4. 验证修改结果
print("\n验证修改结果...")
print(f"修改后celltype列的唯一值: {pd.unique(adata.obs['celltype'])}")

# 5. 查看修改后的分布
print("\n修改后celltype分布:")
celltype_counts_updated = adata.obs['celltype'].value_counts()
print(celltype_counts_updated)

# 6. 检查修改是否正确
if num_to_modify > 0:
    # 确认修改后的数量
    esc_like_count = (adata.obs['celltype'] == 'ESC-like Epi').sum()
    tumor_epi_count = (adata.obs['celltype'] == 'Tumor Epi').sum()
    
    print(f"\n修改结果验证:")
    print(f"  ESC-like Epi细胞数: {esc_like_count}")
    print(f"  Tumor Epi细胞数: {tumor_epi_count}")
    
    if tumor_epi_count == 0 and esc_like_count == num_to_modify:
        print("✓ 修改成功！所有'Tumor Epi'已改为'ESC-like Epi'")
    else:
        print("⚠ 修改结果可能需要检查")

# 7. 保存修改后的数据
print("\n保存修改后的数据...")
output_path = '/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/all_filter_modified.h5ad'
adata.write_h5ad(output_path)
print(f"修改后的数据已保存到: {output_path}")

print("\n完成！")
