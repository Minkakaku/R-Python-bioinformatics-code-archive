import scanpy as sc

# 1. 加载h5ad文件
adata = sc.read_h5ad('raw_count.h5ad')  # 替换为你的文件路径

# 2. 定义要提取的细胞类型列表
target_cell_types = [
    "Stress-responsive_glandular_Epithelium",
    "Proliferative_Epithelium", 
    "Luminal_Epithelium",
    "FGFR2+_Basal_Epithelium",
    "Basal_Epithelium",
    "EMT_Epithelium",
    "Renal_Epithelium"
]

# 3. 筛选细胞
if 'annotation' not in adata.obs.columns:
    print("错误: 在adata.obs中未找到'annotation'列")
    print("可用的列有:", list(adata.obs.columns))
else:
    # 使用scanpy的subset函数，可以保留raw矩阵和其他属性
    subset_adata = adata[adata.obs['annotation'].isin(target_cell_types)].copy()
    
    # 输出结果
    print(f"原始数据: {adata.shape[0]} 个细胞, {adata.shape[1]} 个基因")
    print(f"提取后: {subset_adata.shape[0]} 个细胞, {subset_adata.shape[1]} 个基因")
    
    if subset_adata.raw is not None:
        print(f"Raw矩阵已保留: {subset_adata.raw.shape[0]} 个细胞, {subset_adata.raw.shape[1]} 个基因")
    
    # 显示细胞类型分布
    cell_type_counts = subset_adata.obs['annotation'].value_counts()
    print("\n提取的细胞类型分布:")
    for cell_type, count in cell_type_counts.items():
        print(f"{cell_type}: {count} 个细胞")
    
    print("\n提取完成！子集已保存到 subset_adata 变量中")
subset_adata.write('epithelial_raw.h5ad')
