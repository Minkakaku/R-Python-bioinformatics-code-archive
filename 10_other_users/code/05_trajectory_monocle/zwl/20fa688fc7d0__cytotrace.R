setwd('/mnt/NexData_file/home/zwl/Project/FWAE20250109/require2/')

# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)

sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWAE20250109/require1/HSCs/HSPC_celltype.h5ad',
                      from="anndata", to="seurat",
                      outFile='HSPC_celltype.rds')
Seuset <- readRDS(file = './HSPC_celltype.rds')
Seuset <- FindVariableFeatures(Seuset,nfeatures = 2000)
Seurat::DimPlot(Seuset,reduction = 'umap',group.by = c('type','age','cell_type','phase'))
FeaturePlot(Seuset,features = c('G2M_score','S_score'))


# load the CytoTRACE 2 package
library(CytoTRACE) 

# extract expression data
expression_data <- Seurat$expression_data

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
cytotrace2_result <- cytotrace2(expression_data)

# extract annotation data
annotation <- data$annotation

# generate prediction and phenotype association plots with plotData function
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation,
                  expression_data = expression_data
)
