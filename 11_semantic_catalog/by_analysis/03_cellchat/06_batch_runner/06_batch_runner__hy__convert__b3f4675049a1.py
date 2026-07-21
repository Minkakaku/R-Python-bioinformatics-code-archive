import scanpy as sc
import scipy as sp
import pandas as pd
import os
import numpy as np

def h5ad_to_rds(h5ad_path, output_rds="output_seurat.rds"):
    # 1. 读取并校验数据
    adata = sc.read_h5ad(h5ad_path)
    if 'celltype_final' not in adata.obs:
        raise KeyError("adata.obs中必须包含'celltype_final'列")
    
    # 2. 修复obs重复问题
    if not adata.obs.index.is_unique:
        print("检测到重复的细胞ID，正在修复...")
        adata.obs_names_make_unique()  # 关键修复1[11](@ref)
    
    # 3. 处理矩阵格式
    os.makedirs("temp_export", exist_ok=True)
    matrix = adata.X
    
    # 检查并转换密集矩阵为稀疏格式
    if isinstance(matrix, np.ndarray):
        print("检测到密集矩阵，转换为CSR稀疏格式...")
        matrix = sp.sparse.csr_matrix(matrix)  # 关键修复2[1,8](@ref)
    
    # 4. 导出数据（转置为基因×细胞）
    sp.sparse.save_npz(os.path.join("temp_export", "matrix.npz"), matrix.T)
    adata.obs.to_csv(os.path.join("temp_export", "metadata.csv"))
    adata.var.to_csv(os.path.join("temp_export", "genes.csv"))

    # 5. 生成R重建脚本（修复编码问题）
    r_script = f"""
    library(Matrix)
    library(Seurat)
    
    # 读取数据
    matrix <- t(as(readMM("temp_export/matrix.npz"), "CsparseMatrix")) # 转置回细胞×基因
    metadata <- read.csv("temp_export/metadata.csv", row.names=1)
    genes <- read.csv("temp_export/genes.csv", row.names=1)
    
    # 设置行/列名
    rownames(matrix) <- rownames(genes)
    colnames(matrix) <- rownames(metadata)
    
    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(
        counts = matrix,
        meta.data = metadata,
        project = "Converted_h5ad"
    )
    
    # 显式添加细胞类型（防止元数据覆盖）
    seurat_obj$celltype_final <- metadata$celltype_final  
    
    # 保存RDS
    saveRDS(seurat_obj, file="{output_rds}")
    message(paste("RDS已保存至:", getwd(), "{output_rds}"))
    """
    
    # 6. 写入UTF-8编码的R脚本
    with open("reconstruct_seurat.R", "w", encoding="utf-8") as f:
        f.write(r_script)
    
    print("✅ 数据导出完成！请运行：\nRscript reconstruct_seurat.R")

# 使用示例
if __name__ == "__main__":
    h5ad_to_rds(
        h5ad_path="/home/hy/Project/FWSP20251095/Redo/subcluster/sub/integrated_data.h5ad", 
        output_rds="seurat_converted.rds"
    )
