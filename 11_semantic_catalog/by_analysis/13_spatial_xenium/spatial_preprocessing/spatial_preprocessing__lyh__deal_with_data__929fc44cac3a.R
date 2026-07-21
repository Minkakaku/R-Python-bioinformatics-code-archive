library(Seurat)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the input Seurat object [rdata/rds/qs]", metavar = "FILE"),
  make_option(c("-t", "--type"), type = "character", default = "rds",
              help = "Input Seurat object type [rdata/rds/qs] (default: rds)", metavar = "TYPE"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for split sample objects", metavar = "DIR"),
  make_option(c("-s", "--split_by"), type = "character", default = "sample",
              help = "Metadata column used to split the Seurat object (default: sample)", metavar = "COLNAME")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

required_args <- c("input", "output")
missing_args <- sapply(required_args, function(x) is.null(opt[[x]]))
if (any(missing_args)) {
  missing_names <- paste(names(which(missing_args)), collapse = ", ")
  stop(paste("Missing required arguments:", missing_names), call. = FALSE)
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (!file.exists(opt$input)) {
  stop(paste("Input file not found:", opt$input), call. = FALSE)
}

load_seurat_object <- function(file_path, data_type) {
  if (data_type == "rdata") {
    env <- new.env(parent = emptyenv())
    loaded_names <- load(file_path, envir = env)
    seurat_candidates <- loaded_names[sapply(loaded_names, function(name) inherits(env[[name]], "Seurat"))]
    if (length(seurat_candidates) == 0) {
      stop("No Seurat object found in the RData file.")
    }
    if (length(seurat_candidates) > 1) {
      warning("Multiple Seurat objects found in the RData file. Using the first one.")
    }
    return(env[[seurat_candidates[[1]]]])
  }

  if (data_type == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop("Package 'qs' is required to read .qs files. Please install it.")
    }
    return(qs::qread(file_path))
  }

  if (data_type == "rds") {
    return(readRDS(file_path))
  }

  stop("Unsupported data type. Please use 'rdata', 'rds', or 'qs'.")
}

cat("Loading Seurat object...\n")
seurat_obj <- load_seurat_object(opt$input, opt$type)
if (!inherits(seurat_obj, "Seurat")) {
  stop("The input file does not contain a Seurat object.")
}

if (!opt$split_by %in% colnames(seurat_obj@meta.data)) {
  stop(paste("Split column", opt$split_by, "not found in metadata"))
}

sample_ids <- unique(as.character(seurat_obj@meta.data[[opt$split_by]]))
cat("Found", length(sample_ids), "samples in column", opt$split_by, "\n")

for (sample_id in sample_ids) {
  cat("Processing sample:", sample_id, "\n")
  sample_dir <- file.path(opt$output, sample_id)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }

  sample_cells <- rownames(seurat_obj@meta.data)[as.character(seurat_obj@meta.data[[opt$split_by]]) == sample_id]
  sample_obj <- subset(seurat_obj, cells = sample_cells)
  sample_path <- file.path(sample_dir, "spatial_seurat.rds")
  saveRDS(sample_obj, file = sample_path)
  cat("Saved Seurat object to:", sample_path, "\n")
}

cat("Split finished successfully.\n")
