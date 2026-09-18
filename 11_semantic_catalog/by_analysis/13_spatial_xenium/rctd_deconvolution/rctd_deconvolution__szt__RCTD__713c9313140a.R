

library(tidyverse)
library(Seurat)

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

# all_data_sp_rctd <- readRDS("data_for_analysis/all_data_sp_rctd.rds")
cells_position_ls <- readRDS('data_for_analysis/cells_position_ls.rds')
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")

# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD - Fibroblast' 
dir.create(rctd_out_dir)

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

# 2. single data
all_data_sc <- readRDS("scRNA/scRNA_CAF.rds")
# all_data_rctd <- all_data
colnames(all_data_sc@meta.data)

all_data_sc$annotation <- all_data_sc$celltype

Idents(all_data_sc) <- "annotation"
counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$nCount_RNA

sc_cell_filter_idx <- ! cluster_sc == "other"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
sc_fib <- subset(all_data, annotation4 == 'Fibroblasts')



sc_fib@meta.data$sample %>% unique
sc_fib@meta.data$group %>% unique

all_data_sp_rctd_ls <- sc_fib %>% JoinLayers() %>% SplitObject('sample')

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords    <- 
      cells_position_ls[[name]] %>% 
      dplyr::select(
        cell_id,
        x = x_centroid, 
        y = y_centroid
      ) %>%
      column_to_rownames('cell_id')
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 10, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, 'RCTD - Fibroblast/rctd_result_ls.rds')
rctd_result_ls <- readRDS('RCTD - Fibroblast/rctd_result_ls.rds')

#5. addMetaData
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )


# 6. dimplot
meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

sc_fib <-
  sc_fib %>%
  AddMetaData(meta_rctd_result)

sc_fib$first_type <- as.character(sc_fib$first_type)
sc_fib$first_type[is.na(sc_fib$first_type)] <- "Unknown"

saveRDS(sc_fib, file.path(rctd_out_dir, "xenium_rctd.rds"))
# spatial plot --------------

sc_fib_meta <-
  sc_fib[[]] %>%
  rownames_to_column("cell_id")
cells_boundary_ls <- readRDS(file.path("R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))
# image_wh_r_ls<-readRDS(file.path("data_for_analysis", "image_wh_r_ls.rds"))

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

sp_p_ls1 <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(sc_fib_meta, by = "cell_id") %>% 
      dplyr::filter(!is.na(first_type))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = first_type, group = cell_id))  +
      scale_y_reverse() +
      coord_equal()+
      scale_fill_manual(values = plan) +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })


# # export image
# sp_dir <- "spatial plot" # dir.create(sp_dir)
# dir.create(sp_dir)

sp_p_ls1 %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(rctd_out_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })




# cell ratio --------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

rctd_ratio <- "RCTD - Fibroblast/Ratio"
dir.create(rctd_ratio)

# data <- read.table("Filter_Basicinfo_nGene_nUMI_mito.txt", header=T, sep="\t", check.names=F)
data <- sc_fib[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 21, 22, 47)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

data <- data %>% filter(Annotation != "Unknown")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(rctd_ratio, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(rctd_ratio, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(rctd_ratio, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(rctd_ratio, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

# orders <- read.table("order.txt", header = F, sep = "\t", quote = "", check.names = F)
# data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
# data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])


plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Sample", y = "Percent(%)", title = "") +
  scale_fill_manual(values = plan) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2)) + NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Group", y = "", title = "") +
  scale_fill_manual(values = plan) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(widths = c(3, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(rctd_ratio, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(rctd_ratio, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)






# rctd后回到all_data画umap
##映射------
library(tibble)
library(data.table)
meta.data.t <-
  sc_fib[['first_type']] %>%
  rownames_to_column('cells') %>%
  right_join(
    all_data[['annotation4']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_fbsub = ifelse((as.character(annotation4) == 'Fibroblasts'), as.character(first_type), as.character(annotation4))) %>%
  dplyr::select(cells, annotation_fbsub) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data, meta.data.t['annotation_fbsub'])

saveRDS(all_data, "all_data_fb_rctd_1223.rds")


# plan_rctd <- c("#5050FF","#CE3D32","#749B58","#F0E342","#5DB1DD","#924822","#9AFEFF","#802268","#6BD76B","#D595A7")
pdf(file.path(rctd_out_dir, str_glue("Dimplot annotation_fbsub.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_fbsub", label = TRUE, reduction = 'umap', cols = plan) + NoLegend()
dev.off()




# spatial plot --------------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_fbsub, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = plan) + 
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(rctd_out_dir, str_glue("spatial plot {grp} celltype.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })






## find markers -------
Idents(sc_fib) <- "annotation_fbsub"
DefaultAssay(sc_fib) <- "RNA"

all_markers <- FindAllMarkers(
  JoinLayers(sc_fib), 
  test.use = "wilcox", 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  only.pos = TRUE
)
write.csv(all_markers, file.path(rctd_out_dir, 'cluster diffmarker.csv'))


all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))
clean_top5 <- all_top5_gene[!duplicated(all_top5_gene)]

# top5 dotplot
p_dot_top5 <-
  DotPlot(
    sc_fib, 
    features = clean_top5, 
    dot.scale = 12, 
    col.min = 0, 
    # cols = c("white", "red")
    cols = c("lightgrey", "red")
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(rctd_out_dir, "top 5 gene dotplot.pdf"),
  height = 5,
  width  = 12
)





