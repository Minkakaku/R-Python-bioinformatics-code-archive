library(Seurat)

setwd('C:/Users/ztt/Desktop/FLEX_CRC')

cell_ann <- read.csv("Single_cell_data_cell_annotation.csv", header = TRUE)
head(cell_ann)

mat <- Read10X("CRC")
mat <- mat[,cell_ann$barcode]









