
# RCTD Myeloid -------------
library(tidyverse)
library(Seurat)
library(future)
plan("multisession", workers = 16, gc = TRUE)
options(future.globals.maxSize = 50 * 1024 ^ 3)

all_data <- readRDS("data_for_analysis/all_data_rctd_260716.rds")
all_data_sp_rctd <- all_data[, all_data@meta.data$RCTD == "Myeloid"]

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")

cells_boundary_ls_raw <- readRDS("data_for_analysis/cells_boundary_ls_raw.rds")
cells_position_ls_raw <- readRDS("data_for_analysis/cells_position_ls_raw.rds")

### RCTD ------------------
library(spacexr)
rctd_out_dir <- "260717-RCTD_Myeloid"
dir.create(rctd_out_dir)


### 2.single data ------
all_data_sc <- readRDS("scRNA annotation/Ps_Myeloid_annotation/all_data_sc.rds")
# all_data_rctd <- all_data
colnames(all_data_sc@meta.data)

all_data_sc$annotation <- all_data_sc$celltype

Idents(all_data_sc) <- "annotation"
counts_sc <- all_data_sc %>%
  JoinLayers() %>%
  LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc <- all_data_sc$nCount_RNA

sc_cell_filter_idx <- !cluster_sc == "other"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx],
  droplevels(cluster_sc[sc_cell_filter_idx]),
  nUMI_sc[sc_cell_filter_idx]
)

### 3.create the RCTD query object ------
all_data_sp_rctd@meta.data$sample %>% unique()
all_data_sp_rctd@meta.data$group %>% unique()
# all_data_sp_rctd$sample <- all_data_sp_rctd$project
all_data_sp_rctd_ls <- all_data_sp_rctd %>%
  JoinLayers() %>%
  SplitObject("sample")

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords <-
      cells_position_ls[[name]] %>%
      dplyr::select(
        cell_id,
        x = x_centroid,
        y = y_centroid
      ) %>%
      column_to_rownames("cell_id")
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

### 4.run RCTD ------
rctd_result_ls <-
  rctd_query_ls %>%
  imap(function(query, name) {
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 10, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, file.path(rctd_out_dir, "rctd_result_ls.rds"))


### 5.addMetaData ------
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })


### 6.dimplot ------
meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_sp_rctd_anno <-
  all_data_sp_rctd %>%
  AddMetaData(meta_rctd_result)

all_data_sp_rctd_anno$first_type <- as.character(all_data_sp_rctd_anno$first_type)
all_data_sp_rctd_anno$first_type[is.na(all_data_sp_rctd_anno$first_type)] <- "Unknown"

saveRDS(all_data_sp_rctd_anno, file.path(rctd_out_dir, "all_data_sp_rctd_anno.rds"))


### 7.spatial plot --------------
sp_outdir <- file.path(rctd_out_dir, "spatial plot")
dir.create(sp_outdir)

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

# celltype <- sort(unique(all_data_meta$first_type))
# plan2cell <- plan[1:length(celltype)]
# color_map <- setNames(plan2cell, celltype)
cell <- sort(unique(all_data_sp_rctd_anno$first_type))
writeLines(cell, file.path(rctd_out_dir, "cell.txt"))

color_map <- c(
  "ArtEC_01"        = "#e59572",
  "ArtEC_02"        = "#ea7070",
  "CapEC_01"        = "#600473",
  "CapEC_02"        = "#2694ab",
  "EndoMT"          = "#fdc4b6",
  "VenEC"           = "#96ceb4",
  "VenEC_activated" = "#7F95D1",
  "Unknown"         = "#D3D3D3"
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
      scale_fill_manual(values = color_map) +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls %>%
  iwalk(function(p, grp) {
    pdf(
      file.path(sp_outdir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })


sp_p_ls2 <-
  cells_boundary_ls_raw %>%
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
      scale_fill_manual(values = color_map) +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls2 %>%
  iwalk(function(p, grp) {
    pdf(
      file.path(sp_outdir, str_glue("spatial plot {grp}.pdf")),
      height = 20,
      width  = 15
    )
    print(p)
    dev.off()
  })



### 8.注释信息加入到all_data中 --------
library(tibble)
library(dplyr)

meta_data <- 
  all_data_sp_rctd_anno[["first_type"]] %>% 
  rownames_to_column("cell_id") %>% 
  right_join(
    all_data[["RCTD"]] %>% rownames_to_column("cell_id"),
    by = "cell_id"
  ) %>% 
  mutate(
    RCTD_FBsub = ifelse((
      as.character(RCTD) == "Myeloid"),
      as.character(first_type),
      as.character(RCTD)
    )
  ) %>% 
  dplyr::select(cell_id, RCTD_FBsub) %>% 
  column_to_rownames("cell_id")


all_data <- AddMetaData(
  object = all_data,
  metadata = meta_data )

saveRDS(all_data, file = "data_for_analysis/all_data_RCTD_FBsub.rds")



## sp plot2 ---------
sp_outdir <- file.path(rctd_out_dir, "spatial plot2")
dir.create(sp_outdir)

all_data_meta <-
  all_data[[]] %>%
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

# celltype <- sort(unique(all_data_meta$first_type))
# plan2cell <- plan[1:length(celltype)]
# color_map <- setNames(plan2cell, celltype)

color_map <- c(
  "ArtEC_01"        = "#ea7070",
  "ArtEC_02"        = "#fdc4b6",
  "Basal_KC"        = "#e59572",
  "CapEC_01"        = "#2694ab",
  "CapEC_02"        = "#96ceb4",
  "EndoMT"          = "#600473",
  "Fibroblast"      = "#7F95D1",
  "Hair_Follicle"   = "#dee2d1",
  "LVE"             = "#ef5285",
  "Mast"            = "#a5dff9",
  "Melanocyte"      = "#e8a0b8",
  "Myeloid"         = "#ffeead",
  "Schwann_cells"   = "#737495",
  "SMC"             = "#bccf3d",
  "Spino_KC"        = "#de4307",
  "T_cells"         = "#f29c2b",
  "VenEC"           = "#f6d04d",
  "VenEC_activated" = "#8bc24c",
  "Unknown"         = "#D3D3D3"
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project) { # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>%
      dplyr::filter(!is.na(RCTD_FBsub))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = RCTD_FBsub, group = cell_id)) +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = color_map) +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls %>%
  iwalk(function(p, grp) {
    pdf(
      file.path(sp_outdir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })


sp_p_ls2 <-
  cells_boundary_ls_raw %>%
  imap(function(cells_boundary, project) { # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>%
      dplyr::filter(!is.na(RCTD_FBsub))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = RCTD_FBsub, group = cell_id)) +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = color_map) +
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls2 %>%
  iwalk(function(p, grp) {
    pdf(
      file.path(sp_outdir, str_glue("spatial plot {grp}.pdf")),
      height = 20,
      width  = 15
    )
    print(p)
    dev.off()
  })




annotation_metadata <- all_data[[]] %>% 
  rownames_to_column("cell_id") %>% 
  select(cell_id, RCTD_FBsub, fov) # fov
project_groups <- split(annotation_metadata, annotation_metadata$fov)

walk(names(project_groups), function(proj) { # browser()
  proj_data <- project_groups[[proj]] %>% 
    select(cell_id, group = RCTD_FBsub) %>% 
    mutate(cell_id = str_remove(cell_id, paste0(proj, "_"))) %>% # 重命名列
    mutate(color = color_map[group])
  # mutate(cell_id = str_remove(cell_id, "^\\d+_")) # 重命名列
  
  filename <- paste0(rctd_out_dir, "/", proj, "_cell_annotations.csv")
  
  write.csv(proj_data, file = filename, row.names = FALSE)
})





### 9.cell ratio -------
cell_ratio_outdir <- file.path(rctd_out_dir, "Cell Ratio")
dir.create(cell_ratio_outdir)


library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

data <- all_data_sp_rctd_anno[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c("Cell", "sample", "group", "first_type")]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

data <- data %>% filter(Annotation != "Unknown")

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

orders <- read.table(file.path(rctd_out_dir, "order.txt"), header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])


plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7) +
  labs(x = "Sample", y = "Percent(%)", title = "") +
  scale_fill_manual(values = color_map) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
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
  scale_fill_manual(values = color_map) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2)
  )

p <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(widths = c(12, 1)) # 两张图各自的宽度
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
ggsave(p, file = file.path(cell_ratio_outdir, "diff_boxplot.pdf"), width = 12, height = 5)
ggsave(p, file = file.path(cell_ratio_outdir, "diff_boxplot.png"), width = 12, height = 5, dpi = 300)




library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

dat_props_immune <- all_data_sp_rctd_anno[[]] %>%
  filter(first_type != "Unknown") %>% 
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))

dat_props_immune$group <- factor(dat_props_immune$group, levels = c("Ps", "NC"))

dat_props_immune %>% head()
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153
# 2 Ctrl  Plasma cell  1515 0.162
# 3 Ctrl  T cell       4222 0.451
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127
# 6 Ctrl  Neutrophil    253 0.0270


p_props_immune_flow <-
  dat_props_immune %>%
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = group, y = freq, stratum = first_type, alluvium = first_type, fill = first_type, label = first_type)) +
  geom_stratum(width = 0.6, alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  # geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12)
  ) +
  scale_fill_manual(values = color_map) +
  # scale_x_discrete(expand = c(0.01, 0.01)) +
  # scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )

ggsave(file = file.path(cell_ratio_outdir, "props_flow.pdf"), p_props_immune_flow, width = 5, height = 6)
ggsave(file = file.path(cell_ratio_outdir, "props_flow.png"), p_props_immune_flow, width = 5, height = 6, dpi = 300)







# FindAllMarkers 计算每种RCTD细胞类型高表达基因 -------
Idents(all_data_sp_rctd_anno) <- "first_type"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data_sp_rctd_anno),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
all_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(rctd_out_dir, "marker.csv"), 
    row.names = FALSE
  )


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
    all_data_sp_rctd_anno, 
    features = clean_top5
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(rctd_out_dir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 15
)










