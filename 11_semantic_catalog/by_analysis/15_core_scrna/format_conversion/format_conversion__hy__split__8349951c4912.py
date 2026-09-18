import scanpy as sc
import pandas as pd
import os
# 定义合并规则
celltype_mapping = {
    "Kerationocytes_disease_associated": "Kerationocytes_cell",
    "Kerationocytes": "Epithelial_cell",
    "T cells": "Immune_cell",
    "Macrophages": "Immune_cell"
}

def process_celltypes(adata_path, output_dir):
    # 读取h5ad数据
    adata = sc.read_h5ad(adata_path)
    
    # 创建新列celltype2
    adata.obs["celltype2"] = adata.obs["celltype"].copy()  # 假设原始列名为cell_type
    
    # 应用合并规则
    adata.obs["celltype2"] = adata.obs["celltype2"].replace(celltype_mapping)
    
    # 验证合并结果
    unique_types = adata.obs["celltype2"].unique()
    print(f"合并后的细胞类型: {unique_types}")
    
    # 分割数据
    os.makedirs(output_dir, exist_ok=True)
    for ct in unique_types:
        # 过滤当前细胞类型
        mask = adata.obs["celltype2"] == ct
        subset = adata[mask].copy()
        
        # 生成安全文件名
        safe_name = ct.replace(" ", "_").replace("/", "-")
        subset.write(os.path.join(output_dir, f"{safe_name}.h5ad"))
        print(f"保存 {ct}: {subset.n_obs}细胞 -> {safe_name}.h5ad")

# 使用示例
if __name__ == "__main__":
    process_celltypes(
        adata_path="adata_celltype.h5ad",  # 替换实际路径
        output_dir="./celltype2_split"
    )
