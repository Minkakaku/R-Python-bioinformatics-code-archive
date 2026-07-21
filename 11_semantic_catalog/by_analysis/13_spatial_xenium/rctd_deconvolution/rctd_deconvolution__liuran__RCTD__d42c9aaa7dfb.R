
library(tidyverse)
library(Seurat)
library(dplyr)
library(stringr)

# # only used fibroblast
# all_data_rctd <- readRDS("sc_stma.rds") %>%
#   .[, .$annotation %in% c(
#     'S3 fibroblasts',
#     'Myofibroblasts',
#     # Pericytes                                  
#     # Endothelial cells
#     'S1 fibroblasts',
#     'Immediate early-response (IER) fibroblasts',
#     'Peri-cryptal S2 fibroblasts'
#   )]

scRNA_raw <- readRDS('scRNA/GSE131907_Lung_Cancer_raw_UMI_matrix.rds')
cell_annotation <- read.table("scRNA/GSE131907_Lung_Cancer_cell_annotation.txt", sep = "\t", quote = "", header = T)

scRNA <- CreateSeuratObject(counts = scRNA_raw,
                            meta.data = cell_annotation)


all_data_sc <- scRNA

all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$Cell_subtype
# all_data_sc$annotation <- all_data_sc$Cell_type.refined
Idents(all_data_sc) <- "annotation"

all_data_sc_subset<- all_data_sc[,all_data_sc@meta.data$annotation %in% NA]
cell_na <- colnames(all_data_sc_subset)
all_data_sc_subset <- all_data_sc[, !colnames(all_data_sc) %in% cell_na]
all_data_sc_subset<- all_data_sc_subset[,all_data_sc_subset@meta.data$Sample_Origin %in% c("nLung","tLung","tL/B")]

all_data_sc <- all_data_sc_subset
dim(all_data_sc)

saveRDS(all_data_sc, "RCTD/all_data_sc_for_rctd.rds")


# 查看 annotation 列的原始值
unique(all_data_sc$annotation)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation <- gsub("/", "_", all_data_sc$annotation)
# all_data_sc$annotation <- gsub("CD8+/CD4+ Mixed Th", "CD8+_CD4+ Mixed Th", all_data_sc$annotation)

# 查看替换后的值
unique(all_data_sc$annotation)
table(all_data_sc$annotation)


# RCTD -----------------------
library(spacexr)
rctd_out_dir <- 'RCTD' 
dir.create(rctd_out_dir)

#2. single data
all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
# counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
counts_sc  <- all_data_sc %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$nCount_RNA

sc_cell_filter_idx <-  cluster_sc != "Undetermined"


# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)


##spatial--------
library(tidyverse)
library(Seurat)
library(patchwork)

# 3. create the RCTD query object
all_data_sp_rctd <- readRDS("RCTD/xenium_lung_cluster.rds")
# all_data_sp_rctd_ls <- all_data_sp_rctd  %>% SplitObject('sample')
# 
# rctd_query_ls <-
#   all_data_sp_rctd_ls %>%
#   imap(function(xenium, name){ # browser()
#     counts_hd <- xenium %>% LayerData(layer = "counts")
#     coords    <- 
#       cells_position_ls[[name]] %>% 
#       dplyr::select(
#         cell_id,
#         x = x_centroid, 
#         y = y_centroid
#       ) %>%
#       column_to_rownames('cell_id')
#     
#     # create the RCTD query object
#     query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
#     query
#   })

assay <- "Xenium"
coords_1 <- GetTissueCoordinates(all_data_sp_rctd,image = "PSC-I")
coords_2 <- GetTissueCoordinates(all_data_sp_rctd,image = "PSC-III")
coords <- rbind(coords_1, coords_2)
common_genes <- intersect(rownames(all_data_sp_rctd), rownames(all_data_sc))
rownames(coords) <- coords$cell
coords$cell <- NULL
### Create SpatialRNA object
puck <- SpatialRNA(coords = coords,
                   counts = GetAssayData(all_data_sp_rctd, assay = assay, layer = "counts")[common_genes, ],
                   require_int = TRUE)

myRCTD <- create.RCTD(puck, reference, max_cores = 10, UMI_min = 0,counts_MIN = 0)

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet') #doublet_mode = 'doublet'

saveRDS(myRCTD, "RCTD/myRCTD.rds")

############ result
results <- myRCTD@results
saveRDS(results, "RCTD/RctdResult.rds")


# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] 
spatialRNA <- myRCTD@spatialRNA



plot_weights(cell_type_names, spatialRNA, rctd_out_dir, norm_weights) 
# Plots all weights for each cell type as in full_mode. (saved as 
# 'results/cell_type_weights.pdf')
plot_weights_unthreshold(cell_type_names, spatialRNA, rctd_out_dir, norm_weights) 
# Plots the weights for each cell type as in doublet_mode. (saved as 
# 'results/cell_type_weights_doublets.pdf')
plot_weights_doublet(cell_type_names, spatialRNA, rctd_out_dir, results$weights_doublet, 
                     results$results_df) 


# Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
# 'results/cell_type_occur.pdf')
plot_cond_occur(cell_type_names, rctd_out_dir, norm_weights, spatialRNA)


myRCTD <- readRDS("RCTD/myRCTD.rds")
norm_weights = sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/') 
spaup <- "RCTD"
data = as.data.frame(myRCTD@results[["weights"]])
data$Max_CellType <- colnames(data)[apply(data, 1, which.max)]
results_df = myRCTD@results$results_df
results_df_data <- cbind(data,results_df)
row_indices <- match(rownames(results_df_data), rownames(all_data_sp_rctd@meta.data))
indices <- results_df_data[row_indices, ]
new_meta <- cbind(all_data_sp_rctd@meta.data,indices)
all_data_sp_rctd@meta.data <- new_meta
saveRDS(all_data_sp_rctd,"RCTD/RctdResult_final.rds")
save(all_data_sp_rctd, file = "RCTD/all_data_sp_rctd.rData")
all_data_sp_rctd <- readRDS("RCTD/RctdResult_final.rds")


seg_type <- "xenium"
class_label <- c("Max_CellType","first_type","second_type")
if(seg_type == "visiumHD"){
  for (i in class_label){
    p <- SpatialDimPlot(all_data_sp_rctd,images = fov,group.by = i,
                        image.alpha = 1.5, pt.size.factor = 12,ncol = 1) + 
      theme_bw()+ xlab("")+ylab("")+ 
      theme(axis.text = element_blank(),axis.ticks=element_blank())
    ggsave(p,filename = paste0(outdir,"/RCTD_Plots/",spe,cell_type,i,"_RCTD.png"),width = 15,height = 15,units = "in",dpi = 300)
  }
}else if (seg_type == "xenium"){
  for (i in class_label){
    Idents(all_data_sp_rctd) <- i
    DefaultBoundary(all_data_sp_rctd@images[[fov]] ) <- "segmentation"
    p3 <- ImageDimPlot(all_data_sp_rctd, fov = fov, axes = FALSE, border.color = "#003366", border.size = 0,cols = cols,
                       coord.fixed = FALSE,nmols = 10000)
    ggsave(p3,filename = paste0(rctd_out_dir,"/RCTD_Plots/",spe,cell_type,i,"_RCTD.png"),width = 15,height = 15,units = "in",dpi = 300)
  }
}


plan <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#ED5485",
           "#D9534F","#FFAD60","#57D1C9","#FFEEAD","#7F95D1","#FFE869",
           "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
           "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
           "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
           "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
           "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
           "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
           "#F17D80","#737495")



fov = "PSC-I"
spe = "PSC-I"

fov = "PSC-III"
spe = "PSC-III"

cell_type = "annotation"
for (i in class_label){
  Idents(all_data_sp_rctd) <- i
  DefaultBoundary(all_data_sp_rctd@images[[fov]] ) <- "segmentation"
  p3 <- ImageDimPlot(all_data_sp_rctd, fov = fov, axes = FALSE, border.color = "#003366", border.size = 0,cols = plan,
                     coord.fixed = FALSE,nmols = 10000)
  ggsave(p3,filename = paste0(rctd_out_dir,"/RCTD_Plots/",spe,"_",cell_type,"_",i,"_RCTD.png"),width = 26.5,height = 15,units = "in",dpi = 300)
  ggsave(p3,filename = paste0(rctd_out_dir,"/RCTD_Plots/",spe,"_",cell_type,"_",i,"_RCTD.pdf"),width = 26.5,height = 15)
}


# colnames(all_data_sp_rctd@meta.data) <- gsub("\\+","",colnames(all_data_sp_rctd@meta.data))
# fea <- colnames(data) %>% .[. != "Max_CellType"]
# fea <- gsub("\\+","",fea)
# fea <- gsub(" ","_",fea)
# if(seg_type == "visiumHD"){
#   for (i in fea){
#     p <- SpatialFeaturePlot(spatial,images = fov,features = i,
#                             image.alpha = 1.5, pt.size.factor = 12,ncol = 1) + 
#       theme_bw()+ xlab("")+ylab("")+ 
#       theme(axis.text = element_blank(),axis.ticks=element_blank())
#     ggsave(p,filename = paste0(outdir,"/RCTD_Plots/",spe,"_",cell_type,"_",i,"_RCTD.png"),width = 15,height = 15,units = "in",dpi = 300)
#   }
# }else if (seg_type == "xenium"){
#   for (i in fea){
#     DefaultBoundary(spatial@images[[fov]] ) <- "segmentation"
#     p2 <- ImageFeaturePlot(spatial, fov = fov,axes = FALSE, border.color = "#003366", border.size = 0,
#                            features = i,max.cutoff = "q90",coord.fixed = FALSE,nmols = 10000,
#                            cols = c('#003366', '#FFD700'))
#     ggsave(p2,filename = paste0(outdir,"/RCTD_Plots/",spe,"_",cell_type,"_",i,"_RCTD.png"),width = 15,height = 15,units = "in",dpi = 300)
#   }
# }

data = as.data.frame(myRCTD@results[["weights"]])
data$Max_CellType <- colnames(data)[apply(data, 1, which.max)]

colnames(all_data_sp_rctd@meta.data) <- gsub("\\+","",colnames(all_data_sp_rctd@meta.data))
colnames(all_data_sp_rctd@meta.data) <- gsub(" ","_",colnames(all_data_sp_rctd@meta.data))
fea <- colnames(data) %>% .[. != "Max_CellType"]
fea <- gsub("\\+","",fea)
fea <- gsub(" ","_",fea)


all_data_sp_rctd@meta.data$Max_CellType <- gsub("\\+","",all_data_sp_rctd@meta.data$Max_CellType)
all_data_sp_rctd@meta.data$Max_CellType <- gsub(" ","_",all_data_sp_rctd@meta.data$Max_CellType)


for (i in fea){
  DefaultBoundary(all_data_sp_rctd@images[[fov]] ) <- "segmentation"
  p2 <- ImageFeaturePlot(all_data_sp_rctd, fov = fov,axes = FALSE, border.color = "#003366", border.size = 0,
                         features = i,max.cutoff = "q90",coord.fixed = FALSE,nmols = 10000,
                         cols = c('#003366', '#FFD700'))
  ggsave(p2,filename = paste0(rctd_out_dir,"/RCTD_Plots/",spe,"_",cell_type,"_",i,"_RCTD.png"),width = 26.5,height = 15,units = "in",dpi = 300)
  ggsave(p2,filename = paste0(rctd_out_dir,"/RCTD_Plots/",spe,"_",cell_type,"_",i,"_RCTD.pdf"),width = 26.5,height = 15)
}






