import scanpy as sc

# 读取H5AD文件
adata = sc.read_h5ad("integrated_data.h5ad")
print(adata.raw)
print(adata.obs) 
# 检查重复基因名
gene_names = adata.var_names
print(gene_names)
has_duplicates = gene_names.duplicated().any()
print("是否存在重复基因名:", has_duplicates)

# 统计重复数量
if has_duplicates:
    duplicates = gene_names[gene_names.duplicated(keep=False)]  # 保留所有重复项
    print(f"重复基因数量: {len(duplicates)}")
    print("重复基因示例:", duplicates.unique()[:5])  # 显示前5个重复基因
adata2 = adata.raw.to_adata()
sc.pp.normalize_total(adata2)
adata2.var_names_make_unique
adata2.obs_names_make_unique
adata2.write("new.h5ad")
print(adata2.obs["sub_CNs_10"])
