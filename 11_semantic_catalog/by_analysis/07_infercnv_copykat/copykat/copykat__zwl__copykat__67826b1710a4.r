rm(list = ls())
gc()
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample8/')

#exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample8/all_celltype/raw.h5ad', from="anndata", to="seurat",
                      outFile='raw.rds')
obj <- readRDS(file = './raw.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)

Seurat::DimPlot(obj,reduction = 'umap',group.by = c('celltype','sample','group'))


#copykat############################################
library(copykat)
counts = as.matrix(obj@assays$RNA@counts)
sc_cnv = copykat(rawmat = counts,ngene.chr = 5,sam.name = 'OSCC',n.cores = 20)
saveRDS(file = 'OSCC_copykat.rds',sc_cnv)




