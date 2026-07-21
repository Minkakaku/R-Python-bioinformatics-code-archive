library(Seurat)

all_data <- readRDS("all_data_merge.rds")

sceasy::convertFormat(
	all_data, 
	from="seurat",
	to="anndata",
	outFile='all_data.h5ad')


# SeuratDisk
library(Seurat)
library(SeuratDisk)
SaveH5Seurat(all_data, filename = "all_data.h5Seurat")
Convert("all_data.h5Seurat", dest = "h5ad")
