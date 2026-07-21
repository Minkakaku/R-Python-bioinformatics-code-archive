import scanpy as sc
import pandas as pd
import numpy as np

# 假设你的adata已经加载
# adata = sc.read_h5ad('your_data.h5ad')  # 如果数据需要从文件加载，取消注释此行

# 定义要提取的样本列表
sample_list = [
    'Q2220872-1', 'Q2236500-1', 'Q2243285-1', 'H2510098',
    'Q2406275-2', 'Q2324341-1', 'Z2101519-1', 'Z2370717-2', 'Z2279427-2'
]

# 创建一个字典来存储结果
results = []
adata = sc.read_h5ad('adata_raw.h5ad')

# 循环处理每个样本
for sample_name in sample_list:
    print(f"正在处理样本: {sample_name}")
    
    # 提取当前样本的细胞
    # 请确保 'sample' 是你的AnnData对象中存储样本名的列名，如果不是，请修改
    sample_mask = adata.obs['sample'].str.contains(sample_name)
    sample_adata = adata[sample_mask].copy()
    
    if sample_adata.n_obs == 0:
        print(f"  - 警告: 样本 {sample_name} 没有找到对应细胞，跳过处理。")
        results.append({
            'Sample': sample_name,
            'Total_Genes': 0,
            'Median_Genes_Per_Cell': 0,
            'Number_of_Cells': 0
        })
        continue
    
    # 统计基因总数 (在多少细胞中表达的基因数量)
    # 计算每个基因在多少细胞中表达 (计数 > 0 的细胞数)
    # 对于稀疏矩阵，adata.X > 0 会得到一个布尔矩阵，sum(axis=0) 对每列（基因）求和
    if hasattr(sample_adata.X, 'toarray'):
        # 如果是稀疏矩阵，先转换为稠密数组，但注意这可能消耗大量内存
        genes_expressed_per_gene = np.sum(sample_adata.X.toarray() > 0, axis=0)
    else:
        # 如果是稠密矩阵
        genes_expressed_per_gene = np.sum(sample_adata.X > 0, axis=0)
    gene_count = np.sum(genes_expressed_per_gene > 0) # 统计至少在一个细胞中表达的基因数量
    
    # 统计细胞中位数 (每个细胞表达的基因数量的中位数)
    # 计算每个细胞表达了多少个基因 (计数 > 0 的基因数)
    if hasattr(sample_adata.X, 'toarray'):
        genes_per_cell = np.sum(sample_adata.X.toarray() > 0, axis=1)
    else:
        genes_per_cell = np.sum(sample_adata.X > 0, axis=1)
    
    # 确保 genes_per_cell 是一个一维数组
    genes_per_cell = np.array(genes_per_cell).flatten()
    
    # 计算中位数
    median_genes_per_cell = np.median(genes_per_cell)
    
    # 记录结果
    results.append({
        'Sample': sample_name,
        'Total_Genes': int(gene_count), # 转换为整数
        'Median_Genes_Per_Cell': float(median_genes_per_cell), # 转换为浮点数
        'Number_of_Cells': int(sample_adata.n_obs)
    })
    
    print(f"  - 细胞数量: {sample_adata.n_obs}")
    print(f"  - 基因总数: {gene_count}")
    print(f"  - 细胞中位数: {median_genes_per_cell}")
    print()

# 将结果转换为DataFrame并保存
results_df = pd.DataFrame(results)
print("所有样本统计结果:")
print(results_df)

# 保存结果到CSV文件
results_df.to_csv('sample_statistics.csv', index=False)
print("结果已保存到 'sample_statistics.csv'")
