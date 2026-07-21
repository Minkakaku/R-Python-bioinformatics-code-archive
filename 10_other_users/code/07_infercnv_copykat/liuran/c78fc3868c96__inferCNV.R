# installed into linux server
# mamba install conda-forge::jags
# install.packages('rjags')
# BiocManager::install("infercnv")


# BiocManager::install("infercnv")
# https://sourceforge.net/projects/mcmc-jags/files
# If you are using R 4.1.3 or earlier then install JAGS-4.3.0.exe


# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024^3) # 10Gb


# 1.所有上皮细胞和down sample部分免疫细胞，跑inferCNV分析，
# 识别出正常肝细胞（肝母细胞/胆管上皮细胞？），恶性上皮细胞。
# 恶性上皮细胞跑cNMF分析，识别出上皮细胞程序。

library(infercnv)
library(Seurat)
library(tidyverse)

all_data <- readRDS("all_data_RCTD_subtype2.rds")
all_data <- subset(all_data, subset = RCTD_subtype2 != "Unknown")



# inferCNV PA1 ---------------
### used to infer epethial cells
### immunes cell used as control

infercnv_out_dir <- "inferCNV_PA1_outdir"
dir.create(infercnv_out_dir)


all_data_PA1 <- subset(all_data, subset = sample == "PA_1")

all_data_anno_meta <-
  all_data_PA1[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, annotation = RCTD_subtype2)

# all_data <- all_data %>% AddMetaData(all_data_anno_meta)
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.3)



# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file      = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep       = "\t"
  )

# undebug(CreateInfercnvObject)
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix         = JoinLayers(all_data_PA1)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
  annotations_file          = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file           = file.path("infercnv2", "gene_pos47.txt"),
  delim                     = "\t",
  # min_max_counts_per_cell = c(50,Inf),
  ref_group_names           = c("B cell", "Macrophage", "Mast cell", "T cell")
)


# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans() %>% { sum(. > 0.4) } # 291 <---
infercnv_obj@count.data %>% rowMeans() %>% { sum(. < 0.005) } # 1161 <---


infercnv_obj2 <- infercnv::run(
  infercnv_obj,
  cutoff            = 0.1,
  out_dir           = infercnv_out_dir,
  cluster_by_groups = TRUE,
  plot_steps        = FALSE,
  denoise           = TRUE,
  HMM               = FALSE,
  no_prelim_plot    = TRUE,
  png_res           = 300,
  # useRaster       = FALSE,
  output_format     = "pdf",
  no_plot           = T #  figure margin too large ----
)

infercnv::plot_cnv(
  infercnv_obj2,
  out_dir           = infercnv_out_dir,
  output_filename   = "infercnv",
  output_format     = "pdf",
  png_res           = 300,
  cluster_by_groups = TRUE,
  x.range           = "auto",
  x.center          = 1,
  title             = "infercnv",
  color_safe_pal    = FALSE
)


saveRDS(infercnv_obj2, file.path(infercnv_out_dir, "infercnv_obj2.rds"))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, "infercnv_obj.rds"))



infercnv_obj2_default_median_filtering <-
  infercnv_obj2 %>%
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           = infercnv_out_dir,
    output_filename   = "infercnv_median_filtering",
    output_format     = "pdf",
    cluster_by_groups = TRUE,
    x.range           = "auto",
    x.center          = 1,
    title             = "infercnv",
    color_safe_pal    = FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_out_dir, "infercnv_obj2_default_median_filtering.rds"))



### cnv score ------------
library(scales)
cnvScore <- function(data) {
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax(pmin(data, rng[2]), rng[1])
  
  data <- data %>%
    as.matrix() %>%
    t() %>%
    scale() %>%
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}



cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <- "CNV_score"
temp <-
  cnv_score %>%
  rownames_to_column("cell_id") %>%
  left_join(., as.data.frame(meta_data_infer_subset), "cell_id")

p_infer_cnv_score <-
  temp %>%
  ggplot(aes(annotation, log(CNV_score))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  # RotatedAxis() +
  geom_violin(aes(fill = annotation), color = "NA", drop = FALSE) +
  ggsci::scale_fill_igv()
  # scale_fill_manual(values = color_v) +
  

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height = 6,
  width = 12
)

cols_brewer <- c(
  "#f94144", "#f25c27", "#eb770a", "#f89f12", "#f9c74f", "#90be6d",
  "#43aa8b", "#4d908e", "#577590", "#277da1", "#967aa1",
  "#F2698D", "#F391AE", "#F3B9CF", "#AF8EB8", "#7B84B7", "#3F6FAF",
  "#035AA6", "#0473B4", "#8491B7", "#737AAB", "#5E5E97"
)

p_infer_cnv_score <-
  temp %>%
  ggplot(aes(annotation, log(CNV_score))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  geom_boxplot(aes(fill = annotation), color = "NA", outlier.colour = "black") +
  ggsci::scale_fill_igv() # +
  # scale_fill_manual(values = cols_brewer) # +
  # theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_boxplot.pdf"),
  height = 8,
  width = 10
)








# inferCNV PA2 ---------------
### used to infer epethial cells
### immunes cell used as control

infercnv_out_dir <- "inferCNV_PA2_outdir"
dir.create(infercnv_out_dir)


all_data_PA2 <- subset(all_data, subset = sample == "PA_2")

all_data_anno_meta <-
  all_data_PA2[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, annotation = RCTD_subtype2)

# all_data <- all_data %>% AddMetaData(all_data_anno_meta)
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.3)



# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )

# undebug(CreateInfercnvObject)
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data_PA2)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path("infercnv2", "gene_pos47.txt"),
  delim             = "\t",
  # min_max_counts_per_cell = c(50,Inf),
  ref_group_names   = c("B cell", "Macrophage", "Mast cell", "T cell")
)


# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans() %>% { sum(. > 0.4) } # 291 <---
infercnv_obj@count.data %>% rowMeans() %>% { sum(. < 0.005) } # 1161 <---


infercnv_obj2 <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = infercnv_out_dir,
  cluster_by_groups = TRUE,
  plot_steps = FALSE,
  denoise = TRUE,
  HMM = FALSE,
  no_prelim_plot = TRUE,
  png_res = 300,
  # useRaster   = FALSE,
  output_format = "pdf",
  no_plot = T #  figure margin too large ----
)

infercnv::plot_cnv(
  infercnv_obj2,
  out_dir = infercnv_out_dir,
  output_filename = "infercnv",
  output_format = "pdf",
  png_res = 300,
  cluster_by_groups = TRUE,
  x.range = "auto",
  x.center = 1,
  title = "infercnv",
  color_safe_pal = FALSE
)


saveRDS(infercnv_obj2, file.path(infercnv_out_dir, "infercnv_obj2.rds"))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, "infercnv_obj.rds"))



infercnv_obj2_default_median_filtering <-
  infercnv_obj2 %>%
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir = infercnv_out_dir,
    output_filename = "infercnv_median_filtering",
    output_format = "pdf",
    cluster_by_groups = TRUE,
    x.range = "auto",
    x.center = 1,
    title = "infercnv",
    color_safe_pal = FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_out_dir, "infercnv_obj2_default_median_filtering.rds"))



### cnv score ------------
library(scales)
cnvScore <- function(data) {
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax(pmin(data, rng[2]), rng[1])
  
  data <- data %>%
    as.matrix() %>%
    t() %>%
    scale() %>%
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}



cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <- "CNV_score"
temp <-
  cnv_score %>%
  rownames_to_column("cell_id") %>%
  left_join(., as.data.frame(meta_data_infer_subset), "cell_id")

p_infer_cnv_score <-
  temp %>%
  ggplot(aes(annotation, log(CNV_score))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  # RotatedAxis() +
  geom_violin(aes(fill = annotation), color = "NA", drop = FALSE) +
  ggsci::scale_fill_igv()
# scale_fill_manual(values = color_v) +


ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height = 6,
  width = 12
)

cols_brewer <- c(
  "#f94144", "#f25c27", "#eb770a", "#f89f12", "#f9c74f", "#90be6d",
  "#43aa8b", "#4d908e", "#577590", "#277da1", "#967aa1",
  "#F2698D", "#F391AE", "#F3B9CF", "#AF8EB8", "#7B84B7", "#3F6FAF",
  "#035AA6", "#0473B4", "#8491B7", "#737AAB", "#5E5E97"
)

p_infer_cnv_score <-
  temp %>%
  ggplot(aes(annotation, log(CNV_score))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  geom_boxplot(aes(fill = annotation), color = "NA", outlier.colour = "black") +
  ggsci::scale_fill_igv() # +
# scale_fill_manual(values = cols_brewer) # +
# theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_boxplot.pdf"),
  height = 8,
  width = 10
)










# inferCNV PC ---------------
### used to infer epethial cells
### immunes cell used as control

infercnv_out_dir <- "inferCNV_PC_outdir"
dir.create(infercnv_out_dir)


all_data_PC <- subset(all_data, subset = sample == "PC")

all_data_anno_meta <-
  all_data_PC[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, annotation = RCTD_subtype2)

# all_data <- all_data %>% AddMetaData(all_data_anno_meta)
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.3)



# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )

# undebug(CreateInfercnvObject)
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data_PC)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path("infercnv2", "gene_pos47.txt"),
  delim             = "\t",
  # min_max_counts_per_cell = c(50,Inf),
  ref_group_names   = c("B cell", "Macrophage", "Mast cell", "T cell")
)


# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans() %>% { sum(. > 0.4) } # 291 <---
infercnv_obj@count.data %>% rowMeans() %>% { sum(. < 0.005) } # 1161 <---


infercnv_obj2 <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = infercnv_out_dir,
  cluster_by_groups = TRUE,
  plot_steps = FALSE,
  denoise = TRUE,
  HMM = FALSE,
  no_prelim_plot = TRUE,
  png_res = 300,
  # useRaster   = FALSE,
  output_format = "pdf",
  no_plot = T #  figure margin too large ----
)

infercnv::plot_cnv(
  infercnv_obj2,
  out_dir = infercnv_out_dir,
  output_filename = "infercnv",
  output_format = "pdf",
  png_res = 300,
  cluster_by_groups = TRUE,
  x.range = "auto",
  x.center = 1,
  title = "infercnv",
  color_safe_pal = FALSE
)


saveRDS(infercnv_obj2, file.path(infercnv_out_dir, "infercnv_obj2.rds"))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, "infercnv_obj.rds"))



infercnv_obj2_default_median_filtering <-
  infercnv_obj2 %>%
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir = infercnv_out_dir,
    output_filename = "infercnv_median_filtering",
    output_format = "pdf",
    cluster_by_groups = TRUE,
    x.range = "auto",
    x.center = 1,
    title = "infercnv",
    color_safe_pal = FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_out_dir, "infercnv_obj2_default_median_filtering .rds"))



### cnv score ------------
library(scales)
cnvScore <- function(data) {
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax(pmin(data, rng[2]), rng[1])
  
  data <- data %>%
    as.matrix() %>%
    t() %>%
    scale() %>%
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}



cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <- "CNV_score"
temp <-
  cnv_score %>%
  rownames_to_column("cell_id") %>%
  left_join(., as.data.frame(meta_data_infer_subset), "cell_id")

p_infer_cnv_score <-
  temp %>%
  ggplot(aes(annotation, log(CNV_score))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  # RotatedAxis() +
  geom_violin(aes(fill = annotation), color = "NA", drop = FALSE) +
  ggsci::scale_fill_igv()
# scale_fill_manual(values = color_v) +


ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height = 6,
  width = 12
)

cols_brewer <- c(
  "#f94144", "#f25c27", "#eb770a", "#f89f12", "#f9c74f", "#90be6d",
  "#43aa8b", "#4d908e", "#577590", "#277da1", "#967aa1",
  "#F2698D", "#F391AE", "#F3B9CF", "#AF8EB8", "#7B84B7", "#3F6FAF",
  "#035AA6", "#0473B4", "#8491B7", "#737AAB", "#5E5E97"
)

p_infer_cnv_score <-
  temp %>%
  ggplot(aes(annotation, log(CNV_score))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  geom_boxplot(aes(fill = annotation), color = "NA", outlier.colour = "black") +
  ggsci::scale_fill_igv() # +
# scale_fill_manual(values = cols_brewer) # +
# theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_boxplot.pdf"),
  height = 8,
  width = 10
)














