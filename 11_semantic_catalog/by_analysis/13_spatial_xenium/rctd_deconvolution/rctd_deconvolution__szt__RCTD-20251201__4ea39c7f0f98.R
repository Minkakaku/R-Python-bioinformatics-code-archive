
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)


## spatial--------
# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 10Gb


all_data <- readRDS("data_for_analysis/all_data_projected_rctd_20250902.rds")

# RCTD -----------------------
library(spacexr)


# 2. single data --------
all_data_sc <- readRDS('FWSC20250984.25.6.9.rds')
all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$cell_type_new

# 查看 annotation 列的原始值
unique(all_data_sc$annotation)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation <- gsub("/", "_", all_data_sc$annotation)
all_data_sc$annotation <- gsub("fibroblast_mesenchymal stem cell", "Synovial fibroblasts", all_data_sc$annotation)
all_data_sc$annotation <- gsub("chondrocyte", "Chondrocytes", all_data_sc$annotation)
all_data_sc$annotation <- gsub("synovial intimal fibroblast", "Synovial intimal fibroblasts", all_data_sc$annotation)
all_data_sc$annotation <- gsub("Granulocyte", "Granulocytes", all_data_sc$annotation)
all_data_sc$annotation <- gsub("Osteoblast", "Osteoblasts", all_data_sc$annotation)
all_data_sc$annotation <- gsub("Osteoclast", "Osteoclasts", all_data_sc$annotation)


# 查看替换后的值
unique(all_data_sc$annotation)

out_gene_umi=as.matrix(all_data_sc@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20251201/all_data_sc_meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)



# subcluster数据处理
library(tidyverse)
library(reticulate)

use_python("/home/szt/miniforge3/envs/rpy/bin/python")
reticulate::py_config()

sc <- import("scanpy")
adata <- sc$read("RCTD_20250814/result.h5ad")
counts_x <- adata$X %>% t %>% as.matrix()
colnames(counts_x) <- rownames(adata$obs)

scRNA2 <- CreateSeuratObject(
  counts_x,
  meta.data = adata$obs
)

out_gene_umi=as.matrix(scRNA2@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file="RCTD_20251201/scRNA2_meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)


# 添加新的注释信息
annotation_new <- read.table(file = "RCTD_20251201/new_annotation.txt", header = T, sep = "\t")
annot_vec <- setNames(annotation_new$annotation_new, annotation_new$Cell)
all_data_sc_barcodes <- rownames(all_data_sc@meta.data)
new_annotations <- annot_vec[all_data_sc_barcodes]
all_data_sc@meta.data$annotation_new <- new_annotations



# 查看 annotation 列的原始值
unique(all_data_sc$annotation_new)
# 使用 gsub 函数替换 / 为 _
all_data_sc$annotation_new <- gsub("/", "_", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("fibroblast_mesenchymal stem cell", "Synovial fibroblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub0", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub3", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub5", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub7", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub8", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("chondrocyte-sub9", "Chondrocytes", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Osteoblast-sub1", "Osteoblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Osteoclasts-18", "Osteoclasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("synovial intimal fibroblast", "Synovial intimal fibroblasts", all_data_sc$annotation_new)
all_data_sc$annotation_new <- gsub("Granulocyte-sub4", "Granulocytes-sub4", all_data_sc$annotation_new)


# 查看替换后的值
unique(all_data_sc$annotation_new)

saveRDS(all_data_sc, "RCTD_20250902/all_data_sc_250902.rds")



# remove Erythrocyte --------------
rctd_outdir_remove_Erythrocyte <- "RCTD_20251201/remove_Erythrocyte"
dir.create("RCTD_20251201/remove_Erythrocyte", recursive = T)


all_data <- readRDS("data_for_analysis/all_data_projected_rctd_20250902.rds")

library(spacexr)

# 2. single data
all_data_sc <- readRDS("RCTD_20250902/all_data_sc_250902.rds")

# 单细胞数据处理
Idents(all_data_sc) <- "annotation_new"
all_data_sc <- all_data_sc[, !all_data_sc@meta.data$annotation_new %in% c("Erythrocytes-9", "Granulocytes-sub4", "Erythrocytes-17")]
unique(all_data_sc$annotation_new)

# all_data_sc <- subset(all_data_sc, 
#                       subset = annotation_new != "Erythrocytes-9" &
#                         annotation_new != "Erythrocytes-17")

all_data_sc$annotation_new <- case_when(
  all_data_sc$annotation_new == "Dendritic cells-21" ~ "Dendritic cells",
  all_data_sc$annotation_new == "Endothelial cells-19" ~ "Endothelial cells",
  all_data_sc$annotation_new == "Megakaryocytes-16" ~ "Megakaryocytes-16",
  all_data_sc$annotation_new == "Monocyte macrophages-5" ~ "Monocyte macrophages",
  all_data_sc$annotation_new == "Reticular cells-23" ~ "Reticular cells",
  all_data_sc$annotation_new == "Smooth muscle cells-25" ~ "Smooth muscle cells",
  all_data_sc$annotation_new == "Synovial fibroblasts-sub2" ~ "Synovial fibroblasts",
  all_data_sc$annotation_new == "Synovial intimal fibroblasts-sub10" ~ "Synovial intimal fibroblasts",
  all_data_sc$annotation_new == "T cells-12" ~ "T cells", 
  TRUE ~ as.character(all_data_sc$annotation_new)
)


table(all_data_sc$annotation_new)

count_sc_data <- readRDS('adata_sc_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")

count_sc_data <- subset(count_sc_data, cells = rownames(all_data_sc@meta.data))

counts_sc  <- count_sc_data %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation_new)
nUMI_sc    <- all_data_sc$total_counts

sc_cell_filter_idx <-  cluster_sc != "Unknown"
# sc_cell_filter_idx_new <-  !(cluster_sc %in% c("Unknown","Erythrocytes-9","Erythrocytes-17"))

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
visim_ls4 <- all_data %>% SplitObject(split.by = "sample")
rctd_query_ls <-
  visim_ls4 %>%
  map(function(visium){# browser()
    counts_hd <- visium %>% LayerData(layer = "counts")
    coords    <- visium %>% GetTissueCoordinates() %>% .[, 1:2]
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(
      query, reference, max_cores = 10,
      UMI_min = 50
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })


saveRDS(rctd_result_ls, 'RCTD_20251201/remove_Erythrocyte/rctd_result_ls.rds')


#5. addMetaData
visium_rctd_ls <- 
  visim_ls4 %>% 
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Other" # "Unknown"
    visium
  })

sample_list <- names(visium_rctd_ls)

for (i in sample_list) {
  visium_rctd_ls[[i]]@meta.data$first_type <- gsub("Megakaryocytes-16", "Megakaryocytes", visium_rctd_ls[[i]]@meta.data$first_type)
}




all_data_rctd <- merge(visium_rctd_ls[[1]], visium_rctd_ls[-1])
all_data_rctd$first_type <- gsub("Megakaryocytes-16", "Megakaryocytes", all_data_rctd$first_type)

saveRDS(all_data_rctd, file = "RCTD_20251201/remove_Erythrocyte/all_data_rctd.rds")


#6. spatial plot 
image_scale_ls <- 
  visium_rctd_ls %>% 
  # visium_roi_ls %>% 
  map(GetTissueCoordinates) %>% 
  map('x') %>% 
  map(range) %>% 
  map(diff) %>% 
  map(log) %>% 
  {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}

# RCTD Spatial plot
pdf(
  file.path(rctd_outdir_remove_Erythrocyte, 'spatial plot of RCTD.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 3 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv(na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
  })
dev.off()


# RCTD Spatial plot use custom_color

# plan <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#FFEEAD",
#           "#D9534F","#FFAD60","#57D1C9","#ED5485","#7F95D1","#FFE869",
#           "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
#           "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
#           "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
#           "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
#           "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
#           "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
#           "#F17D80","#737495")

plan1 <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#ED5485",
           "#D9534F","#FFAD60","#57D1C9","#FFEEAD","#7F95D1","#FFE869",
           "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
           "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
           "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
           "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
           "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
           "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
           "#F17D80","#737495")

plan2 <- c("#5050ff","#ce3d32","#749b58","#cfe685","#466983","#ba6338",
           "#ef5285","#802268","#6bd766","#d595a7","#ffc300","#96ceb4",
           "#c75127","#d58f5c","#7a65a5","#005995")

pdf(
  file.path(rctd_outdir_remove_Erythrocyte, 'spatial plot of RCTD - custom color.pdf'), 
  height = 15, 
  width = 15
)

visium_rctd_ls %>%
  imap(function(visium, name){ 
    # if(name == "HC") name_fix <- "NML" else 
    name_fix <- name
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 3 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = plan1, na.value = NA, labels = function(x) str_replace_na(x, "")) +
      guides(fill = guide_legend(override.aes = list(size = 5, color = NA)))
  })
dev.off()






### rctd remove erythrocyte other ------
meta_data_rctd_20251201 <- read.table("RCTD_20251201/remove_Erythrocyte/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(meta_data_rctd_20251201) <- gsub("first_type", "rctd_20251201", colnames(meta_data_rctd_20251201))
meta_data_rctd_20251201$rctd_20251201 <- gsub("Megakaryocytes-16", "Megakaryocytes", meta_data_rctd_20251201$rctd_20251201)

meta_data_rctd_20251201 <- meta_data_rctd_20251201 %>% 
  column_to_rownames("Cell")

head(meta_data_rctd_20251201)

all_data <- all_data %>% 
  AddMetaData(
    meta_data_rctd_20251201 %>% select(rctd_20251201)
  )
head(all_data@meta.data)
unique(all_data$rctd_20251201)

saveRDS(all_data, "data_for_analysis/all_data_projected_rctd_20251201.rds")







### 更改all_data image名
sample_mapping <- c(
  "X24M.16.008um" = "24M-1",
  "X24M.16.008um.2" = "24M-2",
  "X4M.7.008um" = "4M-1",
  "X4M.7.008um.4" = "4M-2"
)

all_data_renameImage <- all_data

current_img_name <- names(all_data_renameImage@images)
new_img_name <- sample_mapping[current_img_name]
unmapped <- current_img_name[is.na(new_img_name)]

names(all_data_renameImage@images) <- new_img_name
print(names(all_data_renameImage@images))

mapping_df <- data.frame(
  Original = current_img_name,
  New = new_img_name
)
write.csv(mapping_df, "image_name_mapping.csv", row.names = FALSE)

saveRDS(all_data_renameImage, file = "data_for_analysis/all_data_renameImage.rds")




# RCTD 空间图
all_data_renameImage$rctd_20251201_NA <- all_data_renameImage$rctd_20251201
all_data_renameImage$rctd_20251201_NA[all_data_renameImage$rctd_20251201_NA == "Other"] <- NA

for (i in unique(all_data_renameImage$sample)){
  p1 <- SpatialDimPlot(all_data_renameImage, images = i, image.scale = "hires", group.by = "rctd_20251201_NA", image.alpha = 2, stroke = 0, pt.size.factor = 3 * image_scale_ls[[i]], shape = 22)
  p2 <- p1 + theme_bw() + xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank()) + 
    scale_fill_manual(values = plan1, na.value = NA, labels = function(x) str_replace_na(x, ""))
  
  ggsave(p2, filename = paste0("RCTD_20251201/remove_Erythrocyte/", i, "_rctd_spatial.png"), width = 20, height = 20, dpi = 300)
  ggsave(p2, filename = paste0("RCTD_20251201/remove_Erythrocyte/", i, "_rctd_spatial.pdf"), width = 20, height = 20)
}





# 剔除骨髓腔中的滑膜细胞-------------
roi_ls <- list.files(
  'ROI_3',
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  # str_remove("\\..*")
  str_remove("_cells_to_remove.csv")


roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    dat <- read.csv(path)
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    dat<-dat %>% mutate(cell = barcode) %>% 
      column_to_rownames('barcode')
    dat
  }) 

roi_coords_ls <- map2(roi_coords_ls, names(roi_coords_ls), function(df, sample_name) {
  df %>% 
    mutate(cell_id = paste(sample_name, cell, sep = "_"), after = cell)
})

library(purrr)

cell_id <- roi_coords_ls %>% 
  map("cell_id") %>% 
  unlist(use.names = FALSE)

all_data <- readRDS("data_for_analysis/all_data_renameImage.rds")

all_data_filt <- subset(all_data, cells = setdiff(colnames(all_data), cell_id))
all_data_filt$rctd_20251201 <- gsub("Basophils-20", "Basophils", all_data_filt$rctd_20251201)

dim(all_data)
dim(all_data_filt)

names(all_data_filt@images)

original_images <- c("24M-1", "24M-2", "4M-1", "4M-2")
images_to_keep <- names(all_data_filt@images)[names(all_data_filt@images) %in% original_images]

# 如果原始图像都存在，直接子集化
if(length(images_to_keep) == 4) {
  all_data_filt@images <- all_data_filt@images[images_to_keep]
  cat("修复完成！当前图像名称:", paste(names(all_data_filt@images), collapse = ", "), "\n")
} else {
  # 如果原始图像不全，使用更安全的方法
  cat("采取备用方案...\n")
  # 方法：基于图像数据内容去重（更安全）
  image_data_list <- list()
  for(img_name in names(all_data_filt@images)) {
    img <- all_data_filt@images[[img_name]]
    # 使用关键标识来去重
    img_key <- paste(img@assay, img@key, sep = "_")
    if(!img_key %in% names(image_data_list)) {
      image_data_list[[img_key]] <- img
    }
  }
  all_data_filt@images <- image_data_list
}

# ### 更改all_data image名
# sample_mapping <- c(
#   "X24M.16.008um" = "24M-1",
#   "X24M.16.008um.2" = "24M-2",
#   "X4M.7.008um" = "4M-1",
#   "X4M.7.008um.4" = "4M-2"
# )
# 
# # all_data_renameImage <- all_data
# 
# current_img_name <- names(all_data_filt@images)
# new_img_name <- sample_mapping[current_img_name]
# unmapped <- current_img_name[is.na(new_img_name)]
# 
# names(all_data_filt@images) <- new_img_name
# print(names(all_data_filt@images))

# RCTD 空间图
all_data_filt$rctd_20251201_NA <- all_data_filt$rctd_20251201
all_data_filt$rctd_20251201_NA[all_data_filt$rctd_20251201_NA == "Other"] <- NA

for (i in unique(all_data_filt$sample)){
  p1 <- SpatialDimPlot(all_data_filt, images = i, image.scale = "hires", group.by = "rctd_20251201_NA", image.alpha = 2, stroke = 0, pt.size.factor = 3 * image_scale_ls[[i]], shape = 22)
  p2 <- p1 + theme_bw() + xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank()) + 
    scale_fill_manual(values = plan1, na.value = NA, labels = function(x) str_replace_na(x, ""))
  
  ggsave(p2, filename = paste0("RCTD_20251201/remove_Erythrocyte/", i, "_rctd_filt_spatial.png"), width = 20, height = 20, dpi = 300)
  ggsave(p2, filename = paste0("RCTD_20251201/remove_Erythrocyte/", i, "_rctd_filt_spatial.pdf"), width = 20, height = 20)
}


saveRDS(all_data_filt, file = "data_for_analysis/all_data_filt.rds")



# rctd remove other erythrocyte 注释信息导出 -------
rctd_1201 <- 
  all_data_filt[[]] %>% 
  rownames_to_column("Barcode") %>% 
  select(Barcode, rctd_20251201, sample)

project_groups <- split(rctd_1201, rctd_1201$sample)

# 为每个 project 创建 CSV 文件
walk(names(project_groups), function(proj) {
  # 提取当前 project 的数据
  proj_data <- project_groups[[proj]] %>%
    select(Barcode, CellType = rctd_20251201) %>% 
    mutate(Barcode = str_remove(Barcode, paste0(proj, "_")))# 重命名列
  
  # 创建文件名
  filename <- paste0("RCTD_20251201/remove_Erythrocyte/", proj, "_cell_annotations_filt.csv")
  
  # 保存为 CSV
  write.csv(proj_data, file = filename, row.names = FALSE)
  
  # 输出信息
  cat("已保存:", filename, "-", nrow(proj_data), "个细胞\n")
})





### 细胞比例 --------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

dir.create("RCTD_20251201/remove_Erythrocyte/Ratio", recursive = T)

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


# out_gene_umi=as.matrix(all_data@meta.data)
# row_name=c("Cell",rownames(out_gene_umi))
# out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
# out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
# write.table(out_gene_umi,file="RCTD_20251201/remove_Erythrocyte/Ratio/meta_data.txt",sep="\t",quote=F,row.names=F,col.names=F)
# 
# data <- read.table("RCTD_20251201/remove_Erythrocyte/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)

data <- all_data_filt[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 11, 10, 29)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data_filt <- data[data$Annotation != "Other"]
data_filt <- data %>% filter(Annotation != "Other")
data <- data_filt

mydata<-data %>% count(Sample,Annotation)
mymtrix<-dcast(mydata,Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20251201/remove_Erythrocyte/Ratio/Cell_type_stat.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Sample"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Sample",colnames(percent)),DATA)
write.table(DATA,file="RCTD_20251201/remove_Erythrocyte/Ratio/Cell_type_stat_percent.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot1)<-c("Sample","Celltype","Abund")

# samplemap <- read.table("sample_map_group.txt", header=T, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="RCTD_20251201/remove_Erythrocyte/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata<-final_data %>% count(Group,Annotation)
mymtrix<-dcast(mydata,Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out),file="RCTD_20251201/remove_Erythrocyte/Ratio/Cell_type_stat_group.txt",quote=F,row.names=T,col.names=F,sep="\t")

rownames(mymtrix)<-mymtrix[,"Group"]
mymtrix<-mymtrix[,-1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind,lapply(1:nrow(mymtrix),function(i){
  mymtrix[i,]<-round(mymtrix[i,]/rowSums(mymtrix)[i]*100,2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent),percent)
DATA<-rbind(c("Group",colnames(percent)),DATA)
write.table(DATA,file="RCTD_20251201/remove_Erythrocyte/Ratio/Cell_type_stat_percent_group.txt",quote=F,row.names=F,col.names=F,sep="\t")
colnames(data4plot2)<-c("Group","Celltype","Abund")

plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1,aes(x=Sample,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Sample",y="Percent(%)",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2,aes(x=Group,y=Abund,fill=Celltype)) +
  geom_bar(stat="identity",size=0.5,width=0.7)+
  labs(x="Group",y="",title="")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=15,colour="black",angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size=15,colour="black"),
    legend.text = element_text(size=15,colour="black"),
    legend.title = element_text(size=15,colour="black",face="bold"),
    axis.title = element_text(size=15,colour="black",face="bold"),
    plot.title = element_text(size=20,colour="black",face="bold",hjust=0.5,lineheight = 0.5),
    strip.background=element_rect(fill="white",size=2))

p <- wrap_plots(plot_list,ncol=2)+
  plot_layout(widths = c(2, 1)) # 两张图各自的宽度
ggsave(p,file="RCTD_20251201/remove_Erythrocyte/Ratio/ratio_barplot.pdf",width=12, height=6)
ggsave(p,file="RCTD_20251201/remove_Erythrocyte/Ratio/ratio_barplot.png",width=12, height=6,dpi=300)



# 箱线图（带显著性标记）
library(ggpubr)

# 普通版
data<-read.table(file = "RCTD_20251201/remove_Erythrocyte/Ratio/Cell_type_stat_percent.txt", header = T, sep="\t", comment.char = "", na.strings="",quote = "",check.names = F)
colnames(data)[1]<-"Gene"

data_filt <- data %>% filter(Gene != "Other")
data <- data_filt

data[,1]<-factor(data[,1])
data<-melt(data,value.name = "value",variable.name="Sample")

rename<-read.table(file="sample_map.txt", header = T, sep = "\t",quote = "",comment.char = "")
plot_data<-left_join(data,rename,by=c("Sample" = "COND1"))

sample<-read.table(file="sample_map_group.txt", header = T, sep = "\t",quote = "",comment.char = "")
colnames(sample) <- c("SampleID","Group")

plot_data<-left_join(plot_data,sample,by="SampleID")

groups<-unique(plot_data[,"Group"])

Y_Max<-max(plot_data[,"value"])*1.5
Y_Max_2<-max(plot_data[,"value"])*1.001

p <-ggplot(plot_data,aes(x=Gene,y=value,fill=Group))+
  geom_boxplot(color="black")+
  geom_jitter(position = position_dodge(0.75),size=1)+
  stat_compare_means(label = "p.format",size=3,
                     aes(group = Group))+
  theme_bw()+
  theme(axis.title.y=element_text(angle=90,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        axis.text.x=element_text(angle=30,colour="black",vjust=1,hjust = 1,
                                 size=12),
        axis.text.y=element_text(angle=0,colour="black",
                                 # hjust=0.5,vjust=0.5,
                                 size=12),
        legend.position="right",
        legend.title=element_text(angle=0,colour="black",face="bold",
                                  size=15,lineheight = 0.5),
        legend.text=element_text(angle=0,colour="black",size=12),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        plot.margin=unit(c(6,6,6,25),"mm"))+
  labs(x="",y="Percent of Celltype (%)")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) #边框粗细
ggsave(p,file="RCTD_20251201/remove_Erythrocyte/Ratio/diff_boxplot.pdf",width=20,height=5)
ggsave(p,file="RCTD_20251201/remove_Erythrocyte/Ratio/diff_boxplot.png",width=20,height=5)











# 对每种细胞类型进行正态性检验
shapiro_results <- plot_data %>%
  group_by(Gene) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(value)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(value)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "RCTD_20251201/remove_Erythrocyte/Ratio/shapiro_test_results.csv", row.names = FALSE)



plot_data_Chondrocytes  <- subset(plot_data, Gene == "Chondrocytes")
print(plot_data_Chondrocytes)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Chondrocytes, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 11),
    breaks = seq(0, 11, by = 2)
  ) +
  # 标签和标题
  labs(
    title = "Chondrocytes",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Chondrocytes_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Chondrocytes_plot.pdf", p_detailed, width = 10, height = 8)



plot_data_Osteoblasts <- subset(plot_data, Gene == "Osteoblasts")
p_detailed <- ggplot(plot_data_Osteoblasts, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(7, 11),
    breaks = seq(7, 11, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Osteoblasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Osteoblasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Osteoblasts_plot.pdf", p_detailed, width = 10, height = 8)





plot_data_Osteoclasts <- subset(plot_data, Gene == "Osteoclasts")
p_detailed <- ggplot(plot_data_Osteoclasts, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.1),
    breaks = seq(0, 2.1, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Osteoclasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Osteoclasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Osteoclasts_plot.pdf", p_detailed, width = 10, height = 8)





plot_data_Synovial_fibroblasts_sub2 <- subset(plot_data, Gene == "Synovial fibroblasts")
print(plot_data_Synovial_fibroblasts_sub2)

p_detailed <- ggplot(plot_data_Synovial_fibroblasts_sub2, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.35),
    breaks = seq(0.1, 0.35, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Synovial fibroblasts",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Synovial_fibroblasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Synovial_fibroblasts_plot.pdf", p_detailed, width = 10, height = 8)






# plot_data_Synovial_intimal_fibroblasts_sub10 <- subset(plot_data, Gene == "Synovial intimal fibroblasts")
# p_detailed <- ggplot(plot_data_Synovial_intimal_fibroblasts_sub10, aes(x = Group, y = value, fill = Group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.4, 0.65),
#     breaks = seq(0.4, 0.65, by = 0.05)
#   ) +
#   # 标签和标题
#   labs(
#     title = "Synovial intimal fibroblasts",
#     x = "",
#     y = "Expression Level",
#     fill = "Group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Synovial_intimal_fibroblasts_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Synovial_intimal_fibroblasts_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_B_cells_24  <- subset(plot_data, Gene == "B cells-24")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_24, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "B cells-24",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/B cells-24_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/B cells-24_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_B_cells_4  <- subset(plot_data, Gene == "B cells-4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_4, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(1, 2.8),
    breaks = seq(1, 2.8, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "B cells-4",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/B cells-4_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/B cells-4_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_B_cells_8  <- subset(plot_data, Gene == "B cells-8")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_B_cells_8, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 8),
    breaks = seq(0, 8, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "B cells-8",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/B cells-8_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/B cells-8_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Basophils_20  <- subset(plot_data, Gene == "Basophils")
print(plot_data_Basophils_20)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Basophils_20, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Basophils-20",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Basophils_20_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Basophils_20_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Dendritic_cells_21  <- subset(plot_data, Gene == "Dendritic cells")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Dendritic_cells_21, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.6, 2),
    breaks = seq(0.6, 2, by = 0.3)
  ) +
  # 标签和标题
  labs(
    title = "Dendritic cells",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Dendritic_cells_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Dendritic_cells_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Endothelial_cells_19  <- subset(plot_data, Gene == "Endothelial cells")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Endothelial_cells_19, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(2.4, 3),
    breaks = seq(2.4, 3, by = 0.1)
  ) +
  # 标签和标题
  labs(
    title = "Endothelial cells",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Endothelial_cells_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Endothelial_cells_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Granulocytes_0  <- subset(plot_data, Gene == "Granulocytes-0")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_0, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(18, 25),
    breaks = seq(18, 25, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-0",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_0_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_0_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_14  <- subset(plot_data, Gene == "Granulocytes-14")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_14, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(8, 11),
    breaks = seq(8, 11, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-14",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_14_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_14_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_15  <- subset(plot_data, Gene == "Granulocytes-15")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_15, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.6, 2),
    breaks = seq(0.6, 2, by = 0.4)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-15",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_15_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_15_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_2  <- subset(plot_data, Gene == "Granulocytes-2")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_2, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.3),
    breaks = seq(0, 0.3, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-2",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_2_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_2_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_3  <- subset(plot_data, Gene == "Granulocytes-3")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_3, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-3",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_3_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_3_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_Granulocytes_6  <- subset(plot_data, Gene == "Granulocytes-6")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_6, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 2.5),
    breaks = seq(0, 2.5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-6",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_6_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_6_plot.pdf", p_detailed, width = 10, height = 8)







plot_data_Granulocytes_7  <- subset(plot_data, Gene == "Granulocytes-7")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Granulocytes_7, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(15, 18),
    breaks = seq(15, 18, by = 1)
  ) +
  # 标签和标题
  labs(
    title = "Granulocytes-7",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_7_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Granulocytes_7_plot.pdf", p_detailed, width = 10, height = 8)






plot_data_Megakaryocytes_16  <- subset(plot_data, Gene == "Megakaryocytes")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Megakaryocytes_16, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(5, 14),
    breaks = seq(5, 14, by = 2)
  ) +
  # 标签和标题
  labs(
    title = "Megakaryocytes",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Megakaryocytes_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Megakaryocytes_plot.pdf", p_detailed, width = 10, height = 8)









plot_data_Monocyte_macrophages_5  <- subset(plot_data, Gene == "Monocyte macrophages")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Monocyte_macrophages_5, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Monocyte macrophages",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Monocyte_macrophages_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Monocyte_macrophages_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Reticular_cells_23  <- subset(plot_data, Gene == "Reticular cells")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Reticular_cells_23, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(3, 5),
    breaks = seq(3, 5, by = 0.5)
  ) +
  # 标签和标题
  labs(
    title = "Reticular cells",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Reticular_cells_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Reticular_cells_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_Smooth_muscle_cells_25  <- subset(plot_data, Gene == "Smooth muscle cells")
print(plot_data_Smooth_muscle_cells_25)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_Smooth_muscle_cells_25, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "Smooth muscle cells",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Smooth_muscle_cells_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/Smooth_muscle_cells_plot.pdf", p_detailed, width = 10, height = 8)








plot_data_T_cells_12  <- subset(plot_data, Gene == "T cells")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_T_cells_12, aes(x = Group, y = value, fill = Group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "T cells",
    x = "",
    y = "Expression Level",
    fill = "Group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/T_cells_plot.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("RCTD_20251201/remove_Erythrocyte/Ratio/T_cells_plot.pdf", p_detailed, width = 10, height = 8)




saveRDS(plot_data, file = "RCTD_20251201/remove_Erythrocyte/Ratio/plot_data.rds")










# rctd remove other erythrocyte 注释信息导出 -------
rctd_1201 <- 
  all_data[[]] %>% 
  rownames_to_column("Barcode") %>% 
  select(Barcode, rctd_20251201, sample)

project_groups <- split(rctd_1201, rctd_1201$sample)

# 为每个 project 创建 CSV 文件
walk(names(project_groups), function(proj) {
  # 提取当前 project 的数据
  proj_data <- project_groups[[proj]] %>%
    select(Barcode, CellType = rctd_20251201) %>% 
    mutate(Barcode = str_remove(Barcode, paste0(proj, "_")))# 重命名列
  
  # 创建文件名
  filename <- paste0("RCTD_20251201/remove_Erythrocyte/", proj, "_cell_annotations.csv")
  
  # 保存为 CSV
  write.csv(proj_data, file = filename, row.names = FALSE)
  
  # 输出信息
  cat("已保存:", filename, "-", nrow(proj_data), "个细胞\n")
})



