import scanpy as sc

file_list = ["B.h5ad", "Endo.h5ad", "Epi.h5ad", "Macro.h5ad", "Neutrophils.h5ad", "Pericytes.h5ad", "Plasma.h5ad", "T.h5ad"]
all_genes_sets = []

print("各文件基因数量及基因名前3个示例:")
for file in file_list:
    adata_temp = sc.read_h5ad(file)
    print(f"{file}: {adata_temp.n_vars} 个基因")
    print(f"    示例基因名: {list(adata_temp.var_names[:3])}")
    all_genes_sets.append(set(adata_temp.var_names))

# 计算所有文件基因名的交集和并集
from functools import reduce
common_genes = reduce(set.intersection, all_genes_sets)
all_genes_union = reduce(set.union, all_genes_sets)

print(f"\n所有文件共有的基因数（交集）: {len(common_genes)}")
print(f"所有文件出现的总基因数（并集）: {len(all_genes_union)}")

# 检查是否有基因名格式问题
for i, file in enumerate(file_list):
    adata_temp = sc.read_h5ad(file)
    unique_genes = set(adata_temp.var_names)
    only_in_this = unique_genes - all_genes_union.union() # 简化检查
    if only_in_this:
        print(f"\n警告：在 {file} 中发现疑似格式不一致的基因名: {list(only_in_this)[:3]} ...")
