suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(cellAdmix)
})

# ------------------------------- #
# User parameters
# ------------------------------- #

transcript_path <- "PATH/TO/transcripts.csv.gz"
cell_meta_path <- "PATH/TO/cell_metadata.csv"

# Optional: if you already have a coarse-clustered Seurat object, load it and
# export a two-column table: cell + coarse_cluster.
coarse_cluster_path <- "PATH/TO/coarse_cluster_mapping.csv"

out_dir <- "celladmix_cluster_clean"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# CellAdmix settings. Start conservatively.
k_nmf <- 8
h_nmf <- 20
nmol_dsamp <- 200000
num_nn_crf <- 10
same_label_ratio <- 5
knn_bridge <- 20
ncores <- 4

# Keep only reasonably sized coarse clusters for the first pass.
min_cells_per_coarse_cluster <- 30

# Bridge test threshold. Smaller = more conservative.
p_bridge <- 0.01

# Optional pre-filtering.
remove_unassigned <- TRUE
remove_control_features <- TRUE


# ------------------------------- #
# Helper functions
# ------------------------------- #

read_xenium_inputs <- function(transcript_path, cell_meta_path, coarse_cluster_path = NULL) {
  tx <- fread(transcript_path)
  meta <- fread(cell_meta_path)

  if (!is.null(coarse_cluster_path) && file.exists(coarse_cluster_path)) {
    coarse_map <- fread(coarse_cluster_path)
    meta <- meta %>% left_join(coarse_map, by = "cell")
  }

  list(tx = tx, meta = meta)
}


prepare_celladmix_df <- function(tx, meta) {
  tx <- as.data.frame(tx)
  meta <- as.data.frame(meta)

  # Adjust these column mappings to your Xenium export.
  # Expected output columns:
  # x, y, z, gene, cell, celltype, mol_id
  tx2 <- tx %>%
    transmute(
      x = x_location,
      y = y_location,
      z = if ("z_location" %in% colnames(tx)) z_location else 0,
      gene = feature_name,
      cell = as.character(cell_id),
      mol_id = if ("transcript_id" %in% colnames(tx)) {
        as.character(transcript_id)
      } else {
        paste0("mol_", seq_len(n()))
      }
    )

  meta2 <- meta %>%
    transmute(
      cell = as.character(cell),
      celltype = as.character(coarse_cluster)
    )

  df <- tx2 %>%
    left_join(meta2, by = "cell")

  if (remove_unassigned) {
    df <- df %>%
      filter(!is.na(cell), cell != "", cell != "UNASSIGNED", cell != "0")
  }

  if (remove_control_features) {
    df <- df %>%
      filter(
        !grepl("^NegControl", gene),
        !grepl("^BLANK", gene, ignore.case = TRUE),
        !grepl("^Unassigned", gene, ignore.case = TRUE)
      )
  }

  df <- df %>%
    filter(!is.na(celltype), celltype != "")

  df$mol_id <- make.unique(as.character(df$mol_id))
  rownames(df) <- df$mol_id
  df
}


filter_small_clusters <- function(df, min_cells_per_coarse_cluster = 30) {
  keep_types <- df %>%
    distinct(cell, celltype) %>%
    count(celltype, name = "n_cells") %>%
    filter(n_cells >= min_cells_per_coarse_cluster) %>%
    pull(celltype)

  df %>% filter(celltype %in% keep_types)
}


factor_fraction_per_cell <- function(df, crf_res) {
  tmp <- df %>%
    transmute(cell = cell, factor = as.character(crf_res[rownames(df), 1]))

  frac <- prop.table(table(tmp$cell, tmp$factor), margin = 1)
  frac <- as.matrix(frac)
  frac
}


removal_list_to_vector <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    return(character(0))
  }

  out <- c()
  for (f_nm in names(x)) {
    if (is.null(x[[f_nm]]) || length(x[[f_nm]]) == 0 || all(is.na(x[[f_nm]]))) {
      next
    }
    f_id <- sub("^f_", "", f_nm)
    out <- c(out, paste0(f_id, "_", x[[f_nm]]))
  }

  unique(out)
}


remove_factor_molecules <- function(df, crf_res, f_rm) {
  if (length(f_rm) == 0) {
    message("No factor-celltype pairs selected for removal.")
    return(df)
  }

  parsed <- lapply(f_rm, function(x) {
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    list(
      factor = parts[[1]],
      celltype = paste(parts[-1], collapse = "_")
    )
  })

  rm_tbl <- bind_rows(lapply(parsed, as.data.frame))
  colnames(rm_tbl) <- c("factor", "celltype")

  df2 <- df
  df2$factor <- as.character(crf_res[rownames(df2), 1])

  keep <- !(
    paste(df2$factor, df2$celltype, sep = "_") %in%
      paste(rm_tbl$factor, rm_tbl$celltype, sep = "_")
  )

  df2[keep, setdiff(colnames(df2), "factor"), drop = FALSE]
}


build_seurat_from_molecules <- function(df) {
  so <- get_counts_meta_seurat(df)
  Seurat::NormalizeData(so) %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA() %>%
    Seurat::FindNeighbors(dims = 1:30) %>%
    Seurat::FindClusters(resolution = 0.3) %>%
    Seurat::RunUMAP(dims = 1:30)
}


save_qc_summary <- function(df_before, df_after, file) {
  summary_tbl <- tibble(
    metric = c(
      "molecules_before",
      "molecules_after",
      "fraction_removed",
      "cells_before",
      "cells_after",
      "celltypes_before",
      "celltypes_after"
    ),
    value = c(
      nrow(df_before),
      nrow(df_after),
      1 - nrow(df_after) / nrow(df_before),
      dplyr::n_distinct(df_before$cell),
      dplyr::n_distinct(df_after$cell),
      dplyr::n_distinct(df_before$celltype),
      dplyr::n_distinct(df_after$celltype)
    )
  )

  fwrite(summary_tbl, file)
  summary_tbl
}


# ------------------------------- #
# 1. Read and prepare Xenium molecules
# ------------------------------- #

inputs <- read_xenium_inputs(
  transcript_path = transcript_path,
  cell_meta_path = cell_meta_path,
  coarse_cluster_path = coarse_cluster_path
)

df_raw <- prepare_celladmix_df(inputs$tx, inputs$meta)
df_use <- filter_small_clusters(
  df_raw,
  min_cells_per_coarse_cluster = min_cells_per_coarse_cluster
)

write.csv(
  df_use[, c("x", "y", "z", "gene", "cell", "celltype", "mol_id")],
  file.path(out_dir, "celladmix_input_molecules.csv"),
  row.names = FALSE
)


# ------------------------------- #
# 2. Learn spatial factors with NMF
# ------------------------------- #

message("Running NMF...")
res_nmf <- run_knn_nmf(
  df = df_use,
  k = k_nmf,
  h = h_nmf,
  nmol.dsamp = nmol_dsamp,
  n.cores = ncores
)

pdf(file.path(out_dir, "nmf_loadings.pdf"), width = 12, height = 8)
print(plot_nmf_loadings(res_nmf))
dev.off()


# ------------------------------- #
# 3. Smooth molecule-factor labels with CRF
# ------------------------------- #

message("Running CRF...")
crf_res <- run_crf_all(
  df = df_use,
  res = res_nmf,
  num.nn = num_nn_crf,
  same.label.ratio = same_label_ratio,
  n.cores = ncores
)

write.csv(
  data.frame(mol_id = rownames(crf_res), factor = crf_res[, 1]),
  file.path(out_dir, "crf_factor_labels.csv"),
  row.names = FALSE
)


# ------------------------------- #
# 4. Bridge test using coarse clusters as temporary cell types
# ------------------------------- #

# get_mol_knn() is not exported by the package, so use ::: explicitly.
message("Preparing molecule KNN matrices for bridge test...")
knn_obj <- cellAdmix:::get_mol_knn(df_use, knn_k = knn_bridge)
ncm_ct <- knn_obj[[1]]
ncm_cell <- knn_obj[[2]]

cell_annot <- df_use %>%
  distinct(cell, celltype) %>%
  as.data.frame()
rownames(cell_annot) <- cell_annot$cell

message("Running bridge test...")
bridge_combo <- run_bridge_test(
  df = df_use,
  crf.f.labs = crf_res,
  cell.annot = cell_annot,
  ncm.ct = ncm_ct,
  ncm.cell = ncm_cell,
  ncells.samp = 200,
  knn.k = knn_bridge,
  n.cores = ncores
)

bridge_res <- extract_bridge_res(
  combo.res = bridge_combo,
  all.ctypes = sort(unique(df_use$celltype)),
  n.factors = k_nmf,
  nmf.type = "joint",
  p.thresh = p_bridge,
  adj.pvals = TRUE
)

# bridge_res is a list:
# [[1]] factors to remove per cell type
# [[2]] likely source cell types
bridge_f_rm_raw <- removal_list_to_vector(bridge_res[[1]])
writeLines(bridge_f_rm_raw, file.path(out_dir, "bridge_factor_celltype_pairs_raw.txt"))


# ------------------------------- #
# 5. Optional refinement with check_f_rm()
# ------------------------------- #

cf_fracs_all <- factor_fraction_per_cell(df_use, crf_res)
cell_annot_vec <- setNames(cell_annot$celltype, cell_annot$cell)

bridge_checked <- check_f_rm(
  cell.annot = cell_annot_vec,
  f.ct.rm = bridge_res[[1]],
  f.from = bridge_res[[2]],
  cf.fracs.all = cf_fracs_all,
  ncm.ct.full = ncm_ct,
  do.clean = TRUE,
  median.thresh = 0.1
)

bridge_f_rm_final <- removal_list_to_vector(bridge_checked)
writeLines(bridge_f_rm_final, file.path(out_dir, "bridge_factor_celltype_pairs_final.txt"))


# ------------------------------- #
# 6. Remove flagged molecules
# ------------------------------- #

df_clean <- remove_factor_molecules(
  df = df_use,
  crf_res = crf_res,
  f_rm = bridge_f_rm_final
)

write.csv(
  df_clean[, c("x", "y", "z", "gene", "cell", "celltype", "mol_id")],
  file.path(out_dir, "molecules_cleaned.csv"),
  row.names = FALSE
)

qc_tbl <- save_qc_summary(
  df_before = df_use,
  df_after = df_clean,
  file = file.path(out_dir, "cleaning_summary.csv")
)
print(qc_tbl)


# ------------------------------- #
# 7. Rebuild Seurat and recluster
# ------------------------------- #

message("Rebuilding Seurat objects...")
so_before <- build_seurat_from_molecules(df_use)
so_after <- build_seurat_from_molecules(df_clean)

saveRDS(so_before, file.path(out_dir, "seurat_before_cleaning.rds"))
saveRDS(so_after, file.path(out_dir, "seurat_after_cleaning.rds"))

pdf(file.path(out_dir, "umap_before_after.pdf"), width = 12, height = 6)
print(DimPlot(so_before, group.by = "cell_type", label = TRUE) + ggtitle("Before cleaning"))
print(DimPlot(so_after, group.by = "cell_type", label = TRUE) + ggtitle("After cleaning"))
dev.off()


# ------------------------------- #
# 8. Optional: estimate contamination scores
# ------------------------------- #

# This step needs a reference scRNA-seq object with cell type labels.
# Example:
#
# so_rna <- readRDS("PATH/TO/reference_scRNA.rds")
# Idents(so_rna) <- "cell_type"
# Idents(so_before) <- "cell_type"
# Idents(so_after) <- "cell_type"
#
# spatial_cells <- so_before@meta.data %>%
#   tibble::rownames_to_column("cell") %>%
#   transmute(cell = cell, cell_type = cell_type)
#
# cell_pos <- df_use %>%
#   group_by(cell) %>%
#   summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
#   left_join(spatial_cells, by = "cell")
#
# cell_adj_df <- estimate_cell_adjacency(cell_pos, random.shift = 1e-2, n.cores = ncores)
# cell_type_adj_mat <- estimate_cell_type_adjacency(cell_adj_df)
#
# scores_before <- estimate_contamination_scores_seurat(
#   so.rna = so_rna,
#   so.spatial = so_before,
#   cell.type.adj.mat = cell_type_adj_mat,
#   idents = "cell_type"
# )
#
# scores_after <- estimate_contamination_scores_seurat(
#   so.rna = so_rna,
#   so.spatial = so_after,
#   cell.type.adj.mat = cell_type_adj_mat,
#   idents = "cell_type"
# )
#
# saveRDS(scores_before, file.path(out_dir, "scores_before.rds"))
# saveRDS(scores_after, file.path(out_dir, "scores_after.rds"))


message("Done. Review these files first:")
message("  - ", file.path(out_dir, "nmf_loadings.pdf"))
message("  - ", file.path(out_dir, "bridge_factor_celltype_pairs_final.txt"))
message("  - ", file.path(out_dir, "cleaning_summary.csv"))
message("  - ", file.path(out_dir, "umap_before_after.pdf"))

