import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# 定义输入文件列表（请确保这些文件在当前工作目录中）
files = [
    'B_cells.h5ad',
    'Endothelial_cells.h5ad',
    'Epithelial_cells.h5ad',
    'Fibroblasts.h5ad',
    'Macrophages.h5ad',
    'Plasma_cells.h5ad',
    'TNK_cells.h5ad'
]

# 初始化列表来存储主数据和raw数据
adatas = []
raw_adatas = []

print("开始加载单个.h5ad文件...")
for file in files:
    print(f"正在加载文件: {file}")
    adata = sc.read_h5ad(file)
    
    # 检查是否有raw矩阵，并提取raw数据[6](@ref)
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()  # 将raw转换回完整的AnnData对象
    else:
        # 如果没有raw矩阵，使用主数据作为后备
        raw_adata = adata.copy()
        print(f"警告: {file} 中没有找到raw矩阵，将使用主数据作为替代")
    
    adatas.append(adata)
    raw_adatas.append(raw_adata)

print("所有文件加载完成，开始合并数据...")

# 合并主数据[6](@ref)
# 使用join='inner'确保只保留所有样本共有的基因
# index_unique='-' 使细胞ID唯一，添加来源后缀
combined_adata = sc.concat(adatas, join='inner', axis=0)

# 合并raw数据
combined_raw_adata = sc.concat(raw_adatas, join='inner', axis=0)

# 将合并后的raw数据设置回主对象[6](@ref)
combined_adata.raw = combined_raw_adata

print("数据合并完成，正在保存结果...")

# 保存合并后的数据
combined_adata.write('combined_data.h5ad')

# 输出合并后的数据基本信息
print("=" * 50)
print("合并完成！数据已保存为 'combined_data.h5ad'")
print(f"合并后的细胞数量: {combined_adata.n_obs}")
print(f"合并后的基因数量: {combined_adata.n_vars}")
print(f"是否保留raw矩阵: {combined_adata.raw is not None}")

if combined_adata.raw is not None:
    print(f"Raw矩阵的细胞数量: {combined_adata.raw.n_obs}")
    print(f"Raw矩阵的基因数量: {combined_adata.raw.n_vars}")
