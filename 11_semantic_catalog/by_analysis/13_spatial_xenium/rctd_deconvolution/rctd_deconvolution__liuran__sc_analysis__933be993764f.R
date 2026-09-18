
library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)

scRNA_raw <- readRDS('scRNA/GSE131907_Lung_Cancer_raw_UMI_matrix.rds')
cell_annotation <- 
  read.table(
    "scRNA/GSE131907_Lung_Cancer_cell_annotation.txt", 
    sep = "\t",
    quote = "", 
    header = T,
    row.names = 1
  )  

scRNA <- CreateSeuratObject(
  counts = scRNA_raw,
  meta.data = cell_annotation
)

scRNA_plot <- scRNA

scRNA_plot <- NormalizeData(scRNA_plot,normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_plot <- FindVariableFeatures(object = scRNA_plot,selection.method = "vst",nfeatures = 2000)
scRNA_plot <- ScaleData(object = scRNA_plot,vars.to.regress ="nCount_RNA") 
scRNA_plot <- RunPCA(object = scRNA_plot)
scRNA_plot <- RunUMAP(scRNA_plot, dims = 1:20,n.neighbors=30,min.dist=0.3)



DimPlot(scRNA_plot,group.by="Cell_type.refined",pt.size=0.1,reduction="umap",raster=FALSE)+
  theme_bw()+ #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y=element_text(vjust=0.5,size=10,face = "bold"),
    axis.text.x=element_text(vjust=0.5,size=10,face = "bold"),
    axis.title.y=element_text(size=15,face = "bold"),
    axis.title.x=element_text(size=15,face = "bold"),
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title=element_text(size=20,face="bold",colour="black",hjust=0.5,lineheight = 0.5)
  )+
  guides(color = guide_legend(override.aes = list(size = 8),ncol=1))+
  labs(title="Cell Type")+ # 标题
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边




all_data_sc <- scRNA

all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$Cell_subtype 

Idents(all_data_sc) <- "annotation"

all_data_sc_subset<- all_data_sc[,all_data_sc@meta.data$annotation %in% NA]
cell_na <- colnames(all_data_sc_subset)

all_data_sc_subset <- all_data_sc[, !colnames(all_data_sc) %in% cell_na]

# saveRDS(all_data_sc, "scRNA/all_data_sc.rds")
saveRDS(all_data_sc_subset, "RCTD/all_data_sc_subset.rds")


