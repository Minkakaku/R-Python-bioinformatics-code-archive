library(reticulate)
library(Seurat)
library(Matrix)

use_condaenv('/mnt/NexData_file/share/scRNA-rawdata/envs/scRNA_p')
anndata <- import("anndata", convert = FALSE)
scipy_sparse <- import("scipy.sparse", convert = FALSE)
np <- reticulate::import("numpy", convert = FALSE)

rds_files <- list.files(path = "/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou6/GSVA/",pattern = '.rds',full.names = T)
adata <- anndata$read_h5ad("/mnt/NexData_file/home/wmx/Project/SC2603132229-01/SC2603132229_report/03-Cell_filter/adata_counts.h5ad")
cell_names <- py_to_r(adata$obs_names$to_list())
gene_names <- py_to_r(adata$var_names$to_list())
x <- py_to_r(adata$X)
# 如果不是稀疏矩阵，可转为稀疏矩阵
if (!inherits(x, "sparseMatrix")) {
  x <- as(x, "dgCMatrix")
}
# AnnData 通常是 cell x gene
# Seurat 需要 gene x cell，所以要转置
x <- t(x)
x <- as(x, "TsparseMatrix")
x <- as(x, "dgCMatrix")
# 补上行列名
rownames(x) <- py_to_r(adata$var_names$to_list())
colnames(x) <- py_to_r(adata$obs_names$to_list())

for (files in rds_files){
  sce <- readRDS(files)
  # 对齐 Seurat 细胞和基因顺序
  common_cells <- intersect(colnames(sce), colnames(x))
  common_genes <- intersect(rownames(sce), rownames(x))
  sce <- subset(sce, cells = common_cells)
  counts <- x[common_genes, common_cells]
  sce[["RNA"]]$counts <- counts
  ff <- gsub(x = files,replacement = '_add_count.rds',pattern = '.rds',)
  saveRDS(sce, ff)
}
