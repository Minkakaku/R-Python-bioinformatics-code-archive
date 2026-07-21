suppressPackageStartupMessages({
  required_pkgs <- c("qs", "Seurat", "dplyr", "tidyr", "tibble", "ggplot2", "scales", "Matrix", "GSVA")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
  }
  library(qs)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(Matrix)
  library(GSVA)
})

set.seed(20260615)

work_dir <- "/mnt/NexData_file/home/lyh/sc_RNAseq_aftersale/Quanming/20260615_Epi"
obj_path <- file.path(work_dir, "Epi_recluster_4subgroups", "Epi_recluster.qs")

out_dir <- file.path(work_dir, "candidate_Epi_DEG_GSVA")
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

# Set to c("Epi_C1", "Epi_C3") to restrict analysis to selected candidate subclusters.
candidate_subclusters <- NULL
min_pseudobulk_cells <- 10
deg_min_pct <- 0.10
deg_logfc_threshold <- 0
deg_sig_logfc <- 0.25
deg_sig_padj <- 0.05
top_deg_n <- 20

write_df <- function(x, filename, row.names = FALSE) {
  csv_file <- file.path(table_dir, paste0(filename, ".csv"))
  utils::write.csv(x, csv_file, row.names = row.names, quote = TRUE)
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    xlsx_file <- file.path(table_dir, paste0(filename, ".xlsx"))
    openxlsx::write.xlsx(x, xlsx_file, rowNames = row.names)
  }
}

save_plot <- function(plot, filename, width, height) {
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

get_assay_data <- function(object, assay = "RNA", layer = "data") {
  mat <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) GetAssayData(object, assay = assay, slot = layer)
  )
  mat
}

average_by_group <- function(expr, groups) {
  groups <- droplevels(factor(groups))
  design <- Matrix::sparse.model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  counts <- Matrix::colSums(design)
  avg <- expr %*% design
  avg <- as.matrix(avg)
  avg <- sweep(avg, 2, counts, "/")
  counts_out <- as.integer(counts)
  names(counts_out) <- names(counts)
  list(expr = avg, counts = counts_out)
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

run_ssgsea <- function(expr_mat, gene_sets) {
  param <- GSVA::ssgseaParam(
    exprData = expr_mat,
    geneSets = gene_sets,
    minSize = 3,
    maxSize = 500,
    normalize = TRUE,
    verbose = FALSE
  )
  as.matrix(GSVA::gsva(param, verbose = FALSE))
}

pathway_gene_sets <- list(
  HER2_ERBB2_SIGNALING = c(
    "ERBB2", "ERBB3", "ERBB4", "EGFR", "GRB7", "STARD3", "PGAP3", "MIEN1",
    "NRG1", "NRG2", "EREG", "AREG", "HBEGF", "SHC1", "GRB2", "SOS1",
    "PIK3CA", "AKT1", "MAPK1", "MAPK3"
  ),
  RTK_SIGNALING = c(
    "EGFR", "ERBB2", "ERBB3", "MET", "FGFR1", "FGFR2", "FGFR3", "IGF1R",
    "PDGFRA", "PDGFRB", "KIT", "KDR", "FLT1", "AXL", "MERTK", "RET",
    "NTRK1", "NTRK2", "NTRK3", "IRS1", "IRS2", "SHC1", "GRB2", "SOS1", "PTPN11"
  ),
  PI3K_AKT_MTOR = c(
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PTEN", "AKT1",
    "AKT2", "AKT3", "MTOR", "RICTOR", "RPTOR", "TSC1", "TSC2", "PDPK1",
    "FOXO1", "FOXO3", "GSK3B", "RPS6KB1", "EIF4EBP1", "BAD", "BCL2L1", "MCL1"
  ),
  MAPK_SIGNALING = c(
    "KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1",
    "MAPK3", "MAPK8", "MAPK9", "MAPK14", "DUSP4", "DUSP6", "ELK1",
    "FOS", "JUN", "MYC", "ETS1", "SPRY2", "SPRY4"
  ),
  PROLIFERATION_CELL_CYCLE = c(
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
    "MCM7", "CDK1", "CCNB1", "CCNB2", "CCNA2", "CDC20", "CDC45",
    "UBE2C", "BIRC5", "AURKA", "AURKB", "CENPF", "TYMS", "STMN1", "HMGB2"
  ),
  MIGRATION_EMT = c(
    "VIM", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "CDH2",
    "ITGA5", "ITGB1", "ITGB6", "MMP2", "MMP7", "MMP9", "MMP14",
    "COL1A1", "COL1A2", "SPARC", "TGFBI", "SERPINE1", "CXCR4", "LAMC2", "LAMB3"
  ),
  STEMNESS = c(
    "PROM1", "ALDH1A1", "ALDH1A3", "SOX2", "SOX9", "NANOG", "POU5F1",
    "KLF4", "LGR5", "OLFM4", "BMI1", "EPCAM", "CD44", "ITGA6", "MYC",
    "NOTCH1", "NOTCH2", "HES1", "JAG1", "YAP1", "WWTR1"
  ),
  APOPTOSIS = c(
    "CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "BCL2L1",
    "MCL1", "BID", "BAD", "BBC3", "PMAIP1", "FAS", "FASLG", "TNFRSF10A",
    "TNFRSF10B", "TNFSF10", "TP53", "PARP1", "XIAP", "BIRC2", "BIRC3"
  ),
  ANTIGEN_PRESENTATION_MHC_I = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8",
    "PSMB9", "PSMB10", "NLRC5", "CALR", "CANX", "ERAP1", "ERAP2"
  ),
  ANTIGEN_PRESENTATION_MHC_II = c(
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD74", "CIITA", "HLA-DMA", "HLA-DMB", "CTSS", "CTSL", "LAMP2"
  ),
  CYTOTOXIC_IFN_RESPONSE = c(
    "GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "GNLY", "IFNG",
    "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "ISG15", "IFI6",
    "IFI27", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "GBP1", "CXCL16"
  )
)

key_genes <- unique(c(
  "ERBB2", "ERBB3", "EGFR", "ERBB4", "GRB7", "PIK3CA", "PIK3CB", "PIK3R1",
  "AKT1", "AKT2", "MTOR", "PTEN", "KRAS", "NRAS", "BRAF", "RAF1", "MAP2K1",
  "MAP2K2", "MAPK1", "MAPK3", "MKI67", "TOP2A", "PCNA", "VIM", "SNAI1",
  "SNAI2", "TWIST1", "ZEB1", "CDH2", "PROM1", "ALDH1A1", "SOX2", "LGR5",
  "BIRC5", "CASP3", "BAX", "BCL2", "HLA-A", "B2M", "HLA-DRA", "TAP1",
  "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "IRF1"
))

message("Reading upstream object: ", obj_path)
obj <- qread(obj_path)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[[DefaultAssay(obj)]], "Assay5") && exists("JoinLayers", mode = "function")) {
  obj <- JoinLayers(obj, assay = DefaultAssay(obj))
}

if (!"Epi_subcluster" %in% colnames(obj@meta.data)) {
  stop("Epi_subcluster was not found in metadata. Run recluster_Epi_4subgroups.R first.")
}
if (!"Group" %in% colnames(obj@meta.data)) {
  if ("group" %in% colnames(obj@meta.data)) {
    obj$Group <- obj$group
  } else {
    obj$Group <- "Unknown"
  }
}
if (!"Sample" %in% colnames(obj@meta.data)) {
  sample_col <- if ("sample" %in% colnames(obj@meta.data)) "sample" else "orig.ident"
  obj$Sample <- obj@meta.data[[sample_col]]
}

all_subclusters <- levels(factor(obj$Epi_subcluster))
if (is.null(candidate_subclusters)) {
  candidate_subclusters <- all_subclusters
}
candidate_subclusters <- intersect(candidate_subclusters, all_subclusters)
if (length(candidate_subclusters) < 2) {
  stop("At least two candidate subclusters are required for DEG analysis.")
}
obj <- subset(obj, subset = Epi_subcluster %in% candidate_subclusters)
obj$Epi_subcluster <- factor(as.character(obj$Epi_subcluster), levels = candidate_subclusters)
obj$Sample <- as.character(obj$Sample)
obj$Group <- factor(as.character(obj$Group), levels = unique(as.character(obj$Group)))
Idents(obj) <- obj$Epi_subcluster

cell_counts <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  count(Sample, Group, Epi_subcluster, name = "n_cells") %>%
  arrange(Epi_subcluster, Group, Sample)
write_df(cell_counts, "candidate_subcluster_cell_counts")

message("Preparing expression matrix")
expr <- get_assay_data(obj, assay = DefaultAssay(obj), layer = "data")
if (ncol(expr) == 0) {
  obj <- NormalizeData(obj, verbose = FALSE)
  expr <- get_assay_data(obj, assay = DefaultAssay(obj), layer = "data")
}

present_gene_sets <- lapply(pathway_gene_sets, function(gs) intersect(unique(gs), rownames(expr)))
gene_set_info <- tibble(
  pathway = names(pathway_gene_sets),
  n_input_genes = vapply(pathway_gene_sets, length, integer(1)),
  n_present_genes = vapply(present_gene_sets, length, integer(1)),
  present_genes = vapply(present_gene_sets, paste, character(1), collapse = ";")
)
write_df(gene_set_info, "candidate_pathway_gene_sets_used")
present_gene_sets <- present_gene_sets[vapply(present_gene_sets, length, integer(1)) >= 3]
if (length(present_gene_sets) == 0) {
  stop("No pathway gene set has at least 3 genes present in the object.")
}

message("Running DEG: each candidate subcluster vs the other candidate subclusters")
deg_list <- lapply(candidate_subclusters, function(cl) {
  other <- setdiff(candidate_subclusters, cl)
  res <- FindMarkers(
    object = obj,
    ident.1 = cl,
    ident.2 = other,
    assay = DefaultAssay(obj),
    only.pos = FALSE,
    test.use = "wilcox",
    min.pct = deg_min_pct,
    logfc.threshold = deg_logfc_threshold,
    verbose = FALSE
  )
  res %>%
    rownames_to_column("gene") %>%
    mutate(cluster = cl, comparison = paste0(cl, "_vs_other_candidates"), .before = 1)
})
deg_all <- bind_rows(deg_list)
fc_col <- if ("avg_log2FC" %in% colnames(deg_all)) "avg_log2FC" else "avg_logFC"
deg_all <- deg_all %>%
  mutate(
    direction = case_when(
      .data[[fc_col]] > 0 ~ "up_in_candidate",
      .data[[fc_col]] < 0 ~ "down_in_candidate",
      TRUE ~ "no_change"
    )
  ) %>%
  arrange(cluster, p_val_adj, desc(abs(.data[[fc_col]])))
write_df(deg_all, "candidate_subcluster_DEG_all")

deg_sig <- deg_all %>%
  filter(p_val_adj < deg_sig_padj, abs(.data[[fc_col]]) >= deg_sig_logfc) %>%
  arrange(cluster, p_val_adj, desc(abs(.data[[fc_col]])))
write_df(deg_sig, "candidate_subcluster_DEG_sig_padj0.05_logFC0.25")

deg_top <- bind_rows(
  deg_sig %>%
    filter(.data[[fc_col]] > 0) %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = top_deg_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(top_direction = "up"),
  deg_sig %>%
    filter(.data[[fc_col]] < 0) %>%
    group_by(cluster) %>%
    slice_min(order_by = .data[[fc_col]], n = top_deg_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(top_direction = "down")
) %>%
  arrange(cluster, top_direction, p_val_adj)
write_df(deg_top, paste0("candidate_subcluster_DEG_top", top_deg_n, "_up_down"))

message("Running GSVA/ssGSEA on subcluster pseudo-bulk expression")
sub_avg <- average_by_group(expr, obj$Epi_subcluster)
gsva_sub <- run_ssgsea(sub_avg$expr, present_gene_sets)
gsva_sub <- gsva_sub[names(present_gene_sets), candidate_subclusters, drop = FALSE]
gsva_sub_z <- row_zscore(gsva_sub)
write_df(as.data.frame(gsva_sub) %>% rownames_to_column("pathway"), "candidate_GSVA_scores_subcluster")
write_df(as.data.frame(gsva_sub_z) %>% rownames_to_column("pathway"), "candidate_GSVA_scores_subcluster_rowZ")

gsva_sub_long <- as.data.frame(gsva_sub) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "Epi_subcluster", values_to = "gsva_score") %>%
  left_join(
    as.data.frame(gsva_sub_z) %>%
      rownames_to_column("pathway") %>%
      pivot_longer(-pathway, names_to = "Epi_subcluster", values_to = "gsva_row_z"),
    by = c("pathway", "Epi_subcluster")
  ) %>%
  mutate(
    pathway = factor(pathway, levels = rev(names(present_gene_sets))),
    Epi_subcluster = factor(Epi_subcluster, levels = candidate_subclusters)
  )
write_df(gsva_sub_long, "candidate_GSVA_scores_subcluster_long")

message("Running GSVA/ssGSEA on sample-subcluster pseudo-bulk expression")
sample_groups <- paste(obj$Sample, obj$Group, obj$Epi_subcluster, sep = "__")
sample_group_counts <- table(sample_groups)
keep_cells <- sample_groups %in% names(sample_group_counts)[sample_group_counts >= min_pseudobulk_cells]
sample_avg <- average_by_group(expr[, keep_cells, drop = FALSE], sample_groups[keep_cells])
gsva_sample <- run_ssgsea(sample_avg$expr, present_gene_sets)
sample_meta <- tibble(sample_group = colnames(gsva_sample)) %>%
  separate(sample_group, into = c("Sample", "Group", "Epi_subcluster"), sep = "__", remove = FALSE) %>%
  left_join(
    tibble(sample_group = names(sample_avg$counts), pseudo_bulk_cell_count = as.integer(sample_avg$counts)),
    by = "sample_group"
  )
gsva_sample_long <- as.data.frame(gsva_sample) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "sample_group", values_to = "gsva_score") %>%
  left_join(sample_meta, by = "sample_group") %>%
  group_by(pathway) %>%
  mutate(gsva_row_z = as.numeric(scale(gsva_score))) %>%
  ungroup() %>%
  mutate(
    gsva_row_z = ifelse(is.na(gsva_row_z), 0, gsva_row_z),
    pathway = factor(pathway, levels = names(present_gene_sets)),
    Epi_subcluster = factor(Epi_subcluster, levels = candidate_subclusters)
  ) %>%
  arrange(pathway, Epi_subcluster, Sample)
write_df(gsva_sample_long, "candidate_GSVA_scores_sample_subcluster_long")
write_df(as.data.frame(gsva_sample) %>% rownames_to_column("pathway"), "candidate_GSVA_scores_sample_subcluster_matrix")

message("Saving GSVA function heatmap")
p_heatmap <- ggplot(gsva_sub_long, aes(x = Epi_subcluster, y = pathway, fill = gsva_row_z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", gsva_score)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Row Z") +
  labs(x = NULL, y = NULL, title = "Candidate Epi subcluster GSVA function heatmap") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
save_plot(p_heatmap, "candidate_function_GSVA_heatmap", width = 7.2, height = 5.8)

message("Saving pathway score plots")
p_dot <- ggplot(gsva_sub_long, aes(x = Epi_subcluster, y = pathway)) +
  geom_point(aes(size = abs(gsva_row_z), color = gsva_row_z)) +
  scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0, name = "Row Z") +
  scale_size_continuous(range = c(2, 8), name = "|Row Z|") +
  labs(x = NULL, y = NULL, title = "Candidate Epi pathway scores") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
save_plot(p_dot, "candidate_pathway_score_dotplot", width = 7.2, height = 5.6)

p_box <- ggplot(gsva_sample_long, aes(x = Epi_subcluster, y = gsva_score, fill = Epi_subcluster)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, color = "grey25", linewidth = 0.25) +
  geom_point(aes(shape = Group), position = position_jitter(width = 0.12, height = 0), size = 1.8, color = "black", alpha = 0.75) +
  facet_wrap(~ pathway, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = setNames(scales::hue_pal()(length(candidate_subclusters)), candidate_subclusters)) +
  labs(x = NULL, y = "GSVA/ssGSEA score", title = "Sample-level pathway scores by candidate Epi subcluster") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70")
  )
save_plot(p_box, "candidate_pathway_score_boxplot_sample_subcluster", width = 11, height = 8)

message("Saving key gene expression plots")
key_genes_present <- intersect(key_genes, rownames(expr))
write_df(tibble(gene = key_genes, present_in_object = key_genes %in% rownames(expr)), "candidate_key_genes_presence")
if (length(key_genes_present) > 0) {
  key_width <- max(9, min(18, length(key_genes_present) * 0.28))
  p_key_dot <- DotPlot(obj, features = key_genes_present, group.by = "Epi_subcluster") +
    scale_color_gradientn(colors = c("#2C7BB6", "#FFFFBF", "#D7191C")) +
    labs(x = NULL, y = NULL, title = "Key pathway genes in candidate Epi subclusters") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  save_plot(p_key_dot, "candidate_key_pathway_gene_dotplot", width = key_width, height = 4.8)
}

message("Saving top DEG average-expression heatmap")
top_heatmap_genes <- deg_sig %>%
  filter(.data[[fc_col]] > 0) %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()
top_heatmap_genes <- intersect(top_heatmap_genes, rownames(sub_avg$expr))
if (length(top_heatmap_genes) >= 2) {
  deg_avg <- sub_avg$expr[top_heatmap_genes, candidate_subclusters, drop = FALSE]
  deg_avg_z <- row_zscore(deg_avg)
  deg_avg_long <- as.data.frame(deg_avg_z) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "Epi_subcluster", values_to = "expr_row_z") %>%
    mutate(Epi_subcluster = factor(Epi_subcluster, levels = candidate_subclusters))
  p_deg_heat <- ggplot(deg_avg_long, aes(x = Epi_subcluster, y = gene, fill = expr_row_z)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Row Z") +
    labs(x = NULL, y = NULL, title = "Top up-regulated DEG average expression") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))
  heat_height <- max(5.5, min(12, length(top_heatmap_genes) * 0.24))
  save_plot(p_deg_heat, "candidate_top_DEG_average_expression_heatmap", width = 7.2, height = heat_height)
}

sink(file.path(out_dir, "sessionInfo.txt"))
cat("Input object:", obj_path, "\n")
cat("Candidate subclusters:", paste(candidate_subclusters, collapse = ", "), "\n")
cat("DEG comparison: each candidate subcluster vs other candidate subclusters\n")
cat("GSVA mode: ssGSEA on pseudo-bulk mean normalized expression\n")
cat("Minimum sample-subcluster cells:", min_pseudobulk_cells, "\n")
print(sessionInfo())
sink()

message("Done. Output directory: ", out_dir)
