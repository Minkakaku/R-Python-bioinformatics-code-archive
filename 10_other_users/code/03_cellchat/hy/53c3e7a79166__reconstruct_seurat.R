
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
    saveRDS(seurat_obj, file="seurat_converted.rds")
    message(paste("RDS已保存至:", getwd(), "seurat_converted.rds"))
    