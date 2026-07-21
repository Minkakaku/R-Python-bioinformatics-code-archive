setwd("/mnt/NexData_file/home/xjs/xutong_spatialEcoTyper/")
getwd()
rm(list = ls())

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

library(SpatialEcoTyper)

## Fused similarity matrix of spatial clusters across all samples
outdir = "integrated_SE/"
integrated = readRDS(file.path(outdir, "MultiSE_integrated_final.rds"))
dim(integrated)


## Single-cell metadata with an 'SE' column annotating the discovered SEs
# Column InitSE: Sample-level spatial clusters (before integration)
# Column SE: Conserved spatial ecotypes (after cross-sample integration)
finalmeta = readRDS(file.path(outdir, "MultiSE_metadata_final.rds"))
head(finalmeta[, c("Sample", "X", "Y", "CellType", "InitSE", "SE")])

# Prepare all-sample metadata for SE-specific state identification
# `finalmeta` comes from integrated_SE/MultiSE_metadata_final.rds.
stopifnot(all(c("CID", "Sample", "CellType", "SE") %in% colnames(finalmeta)))

scmeta <- finalmeta %>%
  mutate(
    CID = as.character(CID),
    Sample = as.character(Sample),
    CellType = as.character(CellType),
    SE = as.character(SE)
  )

scmeta$SE[is.na(scmeta$SE)] <- "NonSE"
sample_ids <- sort(unique(scmeta$Sample))
length(sample_ids)

# Read scdata for every sample in the all-sample integrated result
# Works whether the notebook is run from the project root or from spatial_ecotypes/.
scdata_root <- "./"
scdata_paths <- file.path(scdata_root, sample_ids, paste0(sample_ids, "_scdata.qs"))
names(scdata_paths) <- sample_ids

missing_scdata <- names(scdata_paths)[!file.exists(scdata_paths)]
if (length(missing_scdata) > 0) {
  stop(
    "Missing scdata files for samples: ",
    paste(missing_scdata, collapse = ", "),
    "
Expected paths:
",
    paste(scdata_paths[missing_scdata], collapse = "
")
  )
}

scdata_list <- lapply(scdata_paths, qs::qread)
length(scdata_list)

# Combine all sample expression matrices by shared genes and matched cells
# Each sample matrix is subset to cells present in the integrated metadata.
genes <- Reduce(intersect, lapply(scdata_list, rownames))
message("Shared genes across all samples: ", length(genes))

scdata_list <- lapply(names(scdata_list), function(sample_id) {
  sample_meta <- scmeta %>% filter(Sample == sample_id)
  sample_cells <- intersect(colnames(scdata_list[[sample_id]]), sample_meta$CID)
  
  if (length(sample_cells) == 0) {
    stop("No matched cells between scdata and scmeta for sample: ", sample_id)
  }
  
  scdata_list[[sample_id]][genes, sample_cells, drop = FALSE]
})
names(scdata_list) <- sample_ids

scdata <- do.call(cbind, scdata_list)
matched_cells <- intersect(colnames(scdata), scmeta$CID)
scdata <- scdata[, matched_cells, drop = FALSE]
scmeta <- scmeta[match(matched_cells, scmeta$CID), , drop = FALSE]
rownames(scmeta) <- scmeta$CID

stopifnot(identical(colnames(scdata), rownames(scmeta)))
dim(scdata)

# Leave-one-sample-out prediction on all samples
loocv_preds <- LoocvPredict(
  scdata = scdata,
  scmeta = scmeta,
  Sample = "Sample",
  SE = "SE",
  CellType = "CellType",
  repeats = 30,
  ncores = 20
)

saveRDS(loocv_preds, "loocv_preds.rds")

# Compute F1 matrix for all-sample SE-specific states
F1s <- ComputeMetrics(
  loocv_preds,
  SE = "SE",
  Pred = "cvPred",
  CellType = "CellType",
  Sample = "Sample",
  metric = "F1"
)
saveRDS(F1s, "F1s.rds")

# HeatmapView(F1s, c(0, 0.1, 0.3))
