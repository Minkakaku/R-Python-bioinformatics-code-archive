# 读入各个cluster或者细胞类型或者ECType的特征基因都在文件的cluster列区分基因
# 读入Pseudobulk化的矩阵
# GSVA分析各个cluster的样品打分
# GSVA相关性绘图
rm(list =ls())
Endo_sg <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Degs/Endo_all_cluster_marker.xls"
Epith_sg <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Degs/Epith_all_cluster_marker.xls"
Fibro_sg <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Degs/Fibro_all_cluster_marker.xls"
EC_sg <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Degs/ECType_cluster_diff_statistics.txt"
Pseudobulk_data <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/pseudobulk_logCPM_by_sample_for_GSVA.csv"

library(dplyr)
library(tibble)
library(readr)
library(GSVA)
library(pheatmap)
library(RColorBrewer)


outdir <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/GSVA_cluster_score"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# Function
# =========================
read_degs <- function(file) {
  
  message("[READ] ", file)
  
  df <- read.delim(file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  df <- as.data.frame(df)
  
  message("[COLNAMES] ", paste(colnames(df), collapse = ", "))
  
  if (!"cluster" %in% colnames(df)) {
    stop("文件中找不到 cluster 列：", file)
  }
  
  gene_col_candidates <- c(
    "gene", "Gene", "genes", "Genes",
    "symbol", "Symbol", "SYMBOL",
    "names", "Names",
    "gene_name", "GeneName",
    "feature", "Feature"
  )
  
  gene_col <- intersect(gene_col_candidates, colnames(df))
  
  if (length(gene_col) == 0) {
    stop("文件中找不到基因列，请检查列名：", file)
  }
  
  gene_col <- gene_col[1]
  
  df_use <- df %>%
    dplyr::select(
      cluster,
      gene = all_of(gene_col)
    ) %>%
    mutate(
      cluster = as.character(cluster),
      gene = as.character(gene),
      geneset = cluster
    ) %>%
    filter(
      !is.na(cluster),
      !is.na(gene),
      cluster != "",
      gene != ""
    ) %>%
    distinct(geneset, gene)
  
  return(df_use)
}

make_gene_sets <- function(df) {
  
  gene_sets <- split(df$gene, df$geneset)
  gene_sets <- lapply(gene_sets, unique)
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  
  return(gene_sets)
}

run_gsva_auto <- function(expr_mat, gene_sets) {
  
  gsva_res <- NULL
  
  if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    
    message("[GSVA] use new gsvaParam API")
    
    param <- GSVA::gsvaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      kcdf = "Gaussian"
    )
    
    gsva_res <- GSVA::gsva(param)
    
  } else {
    
    message("[GSVA] use old gsva API")
    
    gsva_res <- GSVA::gsva(
      expr = expr_mat,
      gset.idx.list = gene_sets,
      method = "gsva",
      kcdf = "Gaussian",
      verbose = TRUE
    )
  }
  
  return(gsva_res)
}

# =========================
# 1. Read signature genes
# =========================
Endo_df  <- read_degs(Endo_sg)
Epith_df <- read_degs(Epith_sg)
Fibro_df <- read_degs(Fibro_sg)
EC_df    <- read_degs(EC_sg)

sg_df <- bind_rows(
  Endo_df,
  Epith_df,
  Fibro_df,
  EC_df
)

write.csv(
  sg_df,
  file = file.path(outdir, "all_signature_genes_long_format.csv"),
  row.names = FALSE
)

gene_sets <- make_gene_sets(sg_df)

gene_set_size <- data.frame(
  geneset = names(gene_sets),
  n_gene_raw = lengths(gene_sets),
  stringsAsFactors = FALSE
)

# write.csv(
#   gene_set_size,
#   file = file.path(outdir, "geneset_size_raw.csv"),
#   row.names = FALSE
# )

message("[INFO] Gene set number: ", length(gene_sets))

# =========================
# 2. Read pseudobulk logCPM matrix
# rows = genes, cols = samples
# =========================
expr_df <- read.csv(
  Pseudobulk_data,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_mat <- as.data.frame(expr_df)

expr_mat[] <- lapply(expr_mat, function(x) as.numeric(as.character(x)))

expr_mat <- as.matrix(expr_mat)

expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

message("[INFO] Expression matrix dim:")
print(dim(expr_mat))

# =========================
# 3. Filter genesets by genes in expression matrix
# =========================
gene_sets_filter <- lapply(gene_sets, function(x) {
  intersect(unique(x), rownames(expr_mat))
})

gene_sets_filter <- gene_sets_filter[lengths(gene_sets_filter) >= 5]

gene_set_size_filter <- data.frame(
  geneset = names(gene_sets_filter),
  n_gene_used = lengths(gene_sets_filter),
  stringsAsFactors = FALSE
)

write.csv(
  gene_set_size_filter,
  file = file.path(outdir, "geneset_size_used_in_GSVA.csv"),
  row.names = FALSE
)

message("[INFO] Gene set number after filtering n >= 5: ", length(gene_sets_filter))

if (length(gene_sets_filter) == 0) {
  stop("过滤后没有可用基因集，请检查基因名是否和 pseudobulk 矩阵 rownames 一致")
}

# =========================
# 4. Run GSVA
# =========================
gsva_res <- run_gsva_auto(
  expr_mat = expr_mat,
  gene_sets = gene_sets_filter
)

gsva_res <- as.matrix(gsva_res)

write.csv(
  gsva_res,
  file = file.path(outdir, "GSVA_score_geneset_by_sample.csv")
)

# =========================
# 5. Correlation:
# EC genesets vs Endo / Epith / Fibro genesets
# =========================
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cor_method <- "spearman"   # 可改为 "pearson"

# ---------- 提取不同类型 geneset ----------
ec_sets    <- unique(EC_df$geneset)
epith_sets  <- unique(Epith_df$geneset)
endo_sets <- unique(Endo_df$geneset)
fibro_sets <- unique(Fibro_df$geneset)

target_sets <- c(endo_sets, epith_sets, fibro_sets)

cat("[INFO] EC genesets:", length(ec_sets), "\n")
cat("[INFO] Endo genesets:", length(endo_sets), "\n")
cat("[INFO] Epith genesets:", length(epith_sets), "\n")
cat("[INFO] Fibro genesets:", length(fibro_sets), "\n")

# ---------- 计算相关性 ----------
cor_one_pair <- function(x, y, method = "pearson") {
  
  idx <- complete.cases(x, y)
  x <- x[idx]
  y <- y[idx]
  
  if (length(x) < 3) {
    return(data.frame(
      cor = NA_real_,
      p_value = NA_real_,
      n = length(x)
    ))
  }
  
  if (sd(x, na.rm = TRUE) == 0 | sd(y, na.rm = TRUE) == 0) {
    return(data.frame(
      cor = NA_real_,
      p_value = NA_real_,
      n = length(x)
    ))
  }
  
  res <- cor.test(x, y, method = method)
  
  data.frame(
    cor = as.numeric(res$estimate),
    p_value = res$p.value,
    n = length(x)
  )
}

cor_df <- expand.grid(
  EC_geneset = ec_sets,
  Target_geneset = target_sets,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  do({
    x <- as.numeric(gsva_res[.$EC_geneset, ])
    y <- as.numeric(gsva_res[.$Target_geneset, ])
    
    tmp <- cor_one_pair(x, y, method = cor_method)
    
    data.frame(
      EC_geneset = .$EC_geneset,
      Target_geneset = .$Target_geneset,
      cor = tmp$cor,
      p_value = tmp$p_value,
      n = tmp$n,
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup() %>%
  mutate(
    Target_type = case_when(
      grepl("EC", Target_geneset) ~ "Endo",
      grepl("Cluster", Target_geneset) ~ "Epith",
      grepl("CAF", Target_geneset) ~ "Fibro",
      TRUE ~ "Other"
    ),
    EC_cluster = sub("^EC_", "", EC_geneset),
    Target_cluster = sub("^[^_]+_", "", Target_geneset),
    p_adj = p.adjust(p_value, method = "BH"),
    sig = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

write.csv(
  cor_df,
  file = file.path(outdir, paste0("Correlation_EC_vs_Endo_Epith_Fibro_", cor_method, ".csv")),
  row.names = FALSE
)

# =========================
# 6. Separate heatmap for Endo / Epith / Fibro
# =========================
plot_cor_heatmap_by_type <- function(type_use) {
  
  df_sub <- cor_df %>%
    filter(Target_type == type_use)
  
  if (nrow(df_sub) == 0) return(NULL)
  
  mat_sub <- df_sub %>%
    dplyr::select(EC_geneset, Target_geneset, cor) %>%
    tidyr::pivot_wider(
      names_from = Target_geneset,
      values_from = cor
    ) %>%
    as.data.frame()
  
  rownames(mat_sub) <- mat_sub$EC_geneset
  mat_sub <- mat_sub[, -1, drop = FALSE]
  mat_sub <- as.matrix(mat_sub)
  
  sig_sub <- df_sub %>%
    dplyr::select(EC_geneset, Target_geneset, sig) %>%
    tidyr::pivot_wider(
      names_from = Target_geneset,
      values_from = sig
    ) %>%
    as.data.frame()
  
  rownames(sig_sub) <- sig_sub$EC_geneset
  sig_sub <- sig_sub[, colnames(mat_sub), drop = FALSE]
  sig_sub <- as.matrix(sig_sub)
  
  pdf(
    file = file.path(outdir, paste0("Correlation_heatmap_EC_vs_", type_use, "_", cor_method, ".pdf")),
    width = max(6, ncol(mat_sub) * 0.45),
    height = max(5, nrow(mat_sub) * 0.35)
  )
  
  pheatmap::pheatmap(
    mat_sub,
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = sig_sub,
    number_color = "black",
    fontsize = 8,
    border_color = NA,
    main = paste0("Correlation: EC vs ", type_use, " GSVA score")
  )
  
  dev.off()
  
  png(
    file = file.path(outdir, paste0("Correlation_heatmap_EC_vs_", type_use, "_", cor_method, ".png")),
    width = max(1800, ncol(mat_sub) * 140),
    height = max(1800, nrow(mat_sub) * 120),
    res = 300
  )
  
  pheatmap::pheatmap(
    mat_sub,
    color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = sig_sub,
    number_color = "black",
    fontsize = 8,
    border_color = NA,
    main = paste0("Correlation: EC vs ", type_use, " GSVA score")
  )
  
  dev.off()
}

plot_cor_heatmap_by_type("Endo")
plot_cor_heatmap_by_type("Epith")
plot_cor_heatmap_by_type("Fibro")


# =========================
# 7. Per EC geneset correlation dotplot
# EC vs Endo / Epith / Fibro separately
# =========================
library(dplyr)
library(ggplot2)

outdir_ec_dot <- file.path(outdir, "Per_EC_correlation_dotplot")
dir.create(outdir_ec_dot, showWarnings = FALSE, recursive = TRUE)

plot_one_ec_one_type <- function(cor_df, ec_use, target_type_use, outdir) {
  
  df_plot <- cor_df %>%
    filter(
      EC_geneset == ec_use,
      Target_type == target_type_use
    ) %>%
    mutate(
      sub_cell_type = Target_geneset,
      pvalue = p_value
    ) %>%
    filter(!is.na(cor))
  
  if (nrow(df_plot) == 0) {
    message("[SKIP] ", ec_use, " vs ", target_type_use, ": no data")
    return(NULL)
  }
  
  df_plot <- df_plot %>%
    arrange(cor) %>%
    mutate(
      sub_cell_type = factor(sub_cell_type, levels = sub_cell_type),
      label = paste0(
        "r=", sprintf("%.3f", cor),
        ", p=", sprintf("%.3f", pvalue)
      ),
      neg_log10_p = -log10(pvalue + 1e-300)
    )
  
  p <- ggplot(
    df_plot,
    aes(x = cor, y = sub_cell_type)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_point(
      aes(size = neg_log10_p, color = cor)
    ) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey80"
    ) +
    labs(
      title = paste0(ec_use, " vs ", target_type_use),
      x = paste0("Pearson correlation with ", ec_use, " GSVA score"),
      y = NULL,
      size = "-log10(p)",
      color = "r"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(ec_use, "_vs_", target_type_use, "_correlation_dotplot.pdf")
    ),
    plot = p,
    width = 6,
    height = max(4, length(unique(df_plot$sub_cell_type)) * 0.35)
  )
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(ec_use, "_vs_", target_type_use, "_correlation_dotplot.png")
    ),
    plot = p,
    width = 6,
    height = max(4, length(unique(df_plot$sub_cell_type)) * 0.35),
    dpi = 300
  )
  
  return(p)
}

ec_list <- unique(cor_df$EC_geneset)
target_type_list <- c("Endo", "Epith", "Fibro")

for (ec_use in ec_list) {
  for (target_type_use in target_type_list) {
    
    message("[PLOT] ", ec_use, " vs ", target_type_use)
    
    plot_one_ec_one_type(
      cor_df = cor_df,
      ec_use = ec_use,
      target_type_use = target_type_use,
      outdir = outdir_ec_dot
    )
  }
}

message("[DONE] Per EC correlation dotplots saved to: ", outdir_ec_dot)




# =========================
# gsva_res: geneset x sample
# =========================
gsva_df <- as.data.frame(gsva_res) %>%
  rownames_to_column("geneset")

# =========================
# 转成长表
# =========================
gsva_long <- gsva_df %>%
  pivot_longer(
    cols = -geneset,
    names_to = "sample",
    values_to = "GSVA_score"
  ) %>%
  mutate(
    group = sub("^[0-9]+_", "", sample),
    category = sub("_.*$", "", geneset),
    cluster = sub("^[^_]+_", "", geneset)
  )

# =========================
# 每组每个 geneset 的平均值和中位数
# =========================
gsva_group_summary <- gsva_long %>%
  group_by(group, geneset, category, cluster) %>%
  summarise(
    n_sample = n(),
    mean_score = mean(GSVA_score, na.rm = TRUE),
    median_score = median(GSVA_score, na.rm = TRUE),
    sd_score = sd(GSVA_score, na.rm = TRUE),
    .groups = "drop"
  )

