# 1.RCTD 20260507 -------------

rctd_out_dir <- "20260509_RCTD"
dir.create(rctd_out_dir)


library(Seurat)
library(tidyverse)
library(reticulate)
library(dplyr)
library(spacexr)



### 1.scRNA ------
all_data_sc <- readRDS("scBPD.rds")
all_data_sc[["RNA"]] <- as(object = all_data_sc[["RNA"]], Class = "Assay5")
colnames(all_data_sc@meta.data)


all_data_sc@meta.data$annotation <- all_data_sc@meta.data$celltype

all_data_sc$annotation <- gsub("[ /]", "_", all_data_sc$annotation)



Idents(all_data_sc) <- "annotation"
counts_sc <- all_data_sc %>%
  JoinLayers() %>%
  LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc <- all_data_sc$nCount_RNA

sc_cell_filter_idx <- !cluster_sc == "other"


### 2.create the RCTD reference object ------
reference <- Reference(
  counts_sc[, sc_cell_filter_idx],
  droplevels(cluster_sc[sc_cell_filter_idx]),
  nUMI_sc[sc_cell_filter_idx]
)



### 3. create the RCTD query object ------
cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
all_data_sp_rctd <- readRDS("data_for_analysis/all_data_projected.rds")

all_data_sp_rctd@meta.data$project %>% unique()
all_data_sp_rctd@meta.data$group %>% unique()

all_data_sp_rctd_ls <- all_data_sp_rctd %>%
  JoinLayers() %>%
  SplitObject("sample")



rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    counts <- xenium %>% LayerData(layer = "counts")
    coords <-
      cells_position_ls[[name]] %>%
      dplyr::select(
        cell_id,
        x = x_centroid,
        y = y_centroid
      ) %>%
      column_to_rownames("cell_id")

    # create the RCTD query object
    query <- SpatialRNA(coords, counts, colSums(counts))
    query
  })





### 4.run RCTD ------------------
rctd_result_ls <-
  rctd_query_ls %>%
  imap(function(query, name) {
    # run RCTD
    RCTD <- create.RCTD(query, reference, 
                        max_cores     = 10, 
                        UMI_min       = 10, 
                        UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

    RCTD
  })

saveRDS(rctd_result_ls, file.path(rctd_out_dir, "rctd_result_ls.rds"))



### 5.addMetaData -----------
all_data_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)

    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

# all_data_anno <-
#   merge(
#     all_data_anno_ls[[1]],
#     all_data_anno_ls[-1]
#   )
# saveRDS(all_data_anno, file.path(rctd_out_dir, "all_data_anno.rds"))


meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_sp_rctd_anno <-
  all_data_sp_rctd %>%
  AddMetaData(meta_rctd_result)

saveRDS(all_data_sp_rctd_anno, file.path(rctd_out_dir, "all_data_sp_rctd_anno.rds"))





### spatial plot --------------
all_data_meta <-
  all_data_sp_rctd_anno[[]] %>%
  rownames_to_column("cell_id")

plan <- c(
  "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
  "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead",
  "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c",
  "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1",
  "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E",
  "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d",
  "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9",
  "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6",
  "#005995", "#fa625f", "#f17d80"
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project) { # browser()
    # image_wh_r <- image_wh_r_ls[[project]]

    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>%
      dplyr::filter(!is.na(first_type))

    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = first_type, group = cell_id)) +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan) +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })



sp_p_ls %>%
  iwalk(function(p, grp) {
    pdf(
      file.path(rctd_out_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })



### 注释信息加入到all_data中 --------
meta_data <- all_data_meta %>%
  mutate(RCTD = first_type) %>%
  dplyr::select(cell_id, RCTD) %>%
  column_to_rownames("cell_id")


all_data <- AddMetaData(
  object = all_data,
  metadata = meta_data
)


all_data@meta.data <- all_data@meta.data %>%
  mutate(group2 = gsub("Female|Male", "", group))

saveRDS(all_data, file = "data_for_analysis/all_data_rctd_260511.rds")





### cell ratio -------
cell_ratio_outdir <- file.path(rctd_out_dir, "Cell Ratio")
dir.create(cell_ratio_outdir)


library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

data <- all_data[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 6, 27, 26)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample ~ Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(cell_ratio_outdir, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i) {
  mymtrix[i, ] <- round(mymtrix[i, ] / rowSums(mymtrix)[i] * 100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(cell_ratio_outdir, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header = F, sep = "\t", check.names = F)
# colnames(samplemap)<-c("Sample", "Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data, file = "2. dimension reduction and cluster/res.0.5/Ratio/test.txt", quote = F, row.names = F, col.names = T, sep = "\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group ~ Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(cell_ratio_outdir, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i) {
  mymtrix[i, ] <- round(mymtrix[i, ] / rowSums(mymtrix)[i] * 100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(cell_ratio_outdir, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = T, sep = "\t", quote = "", check.names = F)
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
    strip.background = element_rect(fill = "white", size = 2)
  ) +
  NoLegend()

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
    strip.background = element_rect(fill = "white", size = 2)
  )

p <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(widths = c(4, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(cell_ratio_outdir, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(cell_ratio_outdir, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)





library(ggpubr)

data <- read.table(file = file.path(cell_ratio_outdir, "Cell_type_stat_percent.txt"), header = T, sep = "\t", comment.char = "", na.strings = "", quote = "", check.names = F)
colnames(data)[1] <- "Gene"
data[, 1] <- factor(data[, 1])
data <- melt(data, value.name = "value", variable.name = "Sample")

rename <- read.table(file = "sample_map_group.txt", header = T, sep = "\t", quote = "", comment.char = "")
plot_data <- left_join(data, rename, by = c("Sample" = "SampleID"))

sample <- read.table(file = "sample_map_group.txt", header = T, sep = "\t", quote = "", comment.char = "")
colnames(sample) <- c("SampleID", "Group")

plot_data <- left_join(plot_data, sample, by = c("Sample" = "SampleID"))

groups <- unique(plot_data[, "Group"])

Y_Max <- max(plot_data[, "value"]) * 1.5
Y_Max_2 <- max(plot_data[, "value"]) * 1.001

p <- ggplot(plot_data, aes(x = Gene, y = value, fill = Group)) +
  geom_boxplot(color = "black") +
  geom_jitter(position = position_dodge(0.75), size = 1) +
  stat_compare_means(
    label = "p.signif", size = 3,
    aes(group = Group)
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_text(angle = 90, colour = "black", face = "bold", size = 15, lineheight = 0.5),
    axis.text.x = element_text(angle = 30, colour = "black", vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(
      angle = 0, colour = "black",
      # hjust=0.5,vjust=0.5,
      size = 12
    ),
    legend.position = "right",
    legend.title = element_text(
      angle = 0, colour = "black", face = "bold",
      size = 15, lineheight = 0.5
    ),
    legend.text = element_text(angle = 0, colour = "black", size = 12),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = unit(c(6, 6, 6, 25), "mm")
  ) +
  labs(x = "", y = "Percent of Celltype (%)") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) + # 边框粗细
  scale_fill_manual(values = plan)
ggsave(p, file = file.path(cell_ratio_outdir, "diff_boxplot.pdf"), width = 15, height = 5)
ggsave(p, file = file.path(cell_ratio_outdir, "diff_boxplot.png"), width = 15, height = 5, dpi = 300)






# 2. CN analysis --------

library(Seurat)
library(tidyverse)
source("spatial analysis.R")

Date <- "20260511"
spatial_outdir <- str_glue("{Date} - Spatial analysis")
CN_out_dir <- file.path(spatial_outdir, "CN analysis")
dir.create(CN_out_dir, recursive = T)



all_data <- readRDS(file.path("data_for_analysis", "all_data_rctd_260511.rds"))
all_data@meta.data$annotation <- all_data@meta.data$RCTD
Idents(all_data) <- "annotation"

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

## CN analysis ---------
# options(warn = 0)
#
# Classification1 must be a factor
# add a column: obj_id


all_data_submeta_for_cn_analysis <-
  all_data[[c("sample", "group2", "annotation")]] %>%
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c("cell_id", "x_centroid", "y_centroid")], by = "cell_id") %>%
  dplyr::rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation = factor(annotation), # import!
    obj_id = cell_id # import!
  )

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x              = "X",
  y              = "Y",
  group.by       = "sample",
  classification = "annotation",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 7,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


### CN k15 ---------
res_cn_15 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x              = "X",
  y              = "Y",
  group.by       = "sample",
  classification = "annotation",
  k.nn           = 10,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 9,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()



### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, "CN10 spatial plot")
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")

All_type_meta_cn10 <-
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(CN_out_dir, "all_data_meta_cn10.rds"))

cn10_sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # image_wh_r <- image_wh_r_ls[[sample]]

    anno_data <- All_type_meta_cn10[c("sample", "group2")] %>% distinct()
    group_info <- anno_data$group2[anno_data$sample == sample]

    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")

    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id)) +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_"))
    p
  })


# export image
cn10_sp_p_ls %>%
  iwalk(function(p, name) {
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })



### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>%
  group_by(sample, group2, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))
All_type_cn10_count_df$CN <- factor(All_type_cn10_count_df$CN, levels = paste0("CN", 1:10))

p_cn10_props <-
  All_type_cn10_count_df %>%
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn10_props_group <-
  All_type_cn10_count_df %>%
  ggplot(aes(x = CN, y = freq, fill = group2)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print()
p_cn10_props_group %>% print()
dev.off()


### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, "CN15 spatial plot")
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")

All_type_meta_cn15 <-
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(CN_out_dir, "all_data_meta_cn15.rds"))

cn15_sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()

    anno_data <- All_type_meta_cn15[c("sample", "group2")] %>% distinct()
    group_info <- anno_data$group2[anno_data$sample == sample]

    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn15, by = "cell_id")

    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id)) +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_"))
    p
  })


# export image
cn15_sp_p_ls %>%
  iwalk(function(p, name) {
    pdf(
      file.path(cn15_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

### cell props ------
All_type_cn15_count_df <-
  All_type_meta_cn15 %>%
  group_by(sample, group2, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))
All_type_cn15_count_df$CN <- factor(All_type_cn15_count_df$CN, levels = paste0("CN", 1:15))

p_cn15_props <-
  All_type_cn15_count_df %>%
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn15_props_group <-
  All_type_cn15_count_df %>%
  ggplot(aes(x = CN, y = freq, fill = group2)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn15_props %>% print()
p_cn15_props_group %>% print()
dev.off()





