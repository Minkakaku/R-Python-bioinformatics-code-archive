###load library####
library(Seurat)
library(tidyverse)
library(openxlsx)
library(spacexr)
library(Matrix)
library(rhdf5)
library(optparse)

#
option_list <- list(
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="Path to single-cell reference data (required)", metavar="FILE"),
  make_option(c("-a", "--anno_column"), type="character", default=NULL, 
              help="Column name for cell type annotations (required)", metavar="COLNAME"),
  make_option(c("-t", "--type"), type="character", default="rds", 
              help="Reference data type [rdata/rds/qs] (default: rds)", metavar="TYPE"),
  make_option(c("-q", "--spatial_obj"), type="character", default=NULL,
              help="Path to spatial Seurat object (preferred input)", metavar="FILE"),
  make_option(c("-y", "--spatial_type"), type="character", default="rds",
              help="Spatial object type [rdata/rds/qs] (default: rds)", metavar="TYPE"),
  make_option(c("-s", "--spatial_ex"), type="character", default=NULL, 
              help="Path to spatial expression matrix (legacy HDF5 input)", metavar="FILE"),
  make_option(c("-c", "--coord"), type="character", default=NULL, 
              help="Path to spatial coordinates (legacy CSV input)", metavar="FILE"),
  make_option(c("-n", "--ncores"), type="integer", default=12, 
              help="Number of CPU cores to use (default: 12)", metavar="NUM"),
  make_option(c("-o", "--output"), type="character", default="./results", 
              help="Output directory path (default: ./results)", metavar="DIR")
)

# 更健壮的参数解析
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 强制检查必要参数
required_args <- c("ref", "anno_column")
missing_args <- sapply(required_args, function(x) is.null(opt[[x]]))
if (any(missing_args)) {
  missing_names <- paste(names(which(missing_args)), collapse=", ")
  stop(paste("Missing required arguments:", missing_names), call.=FALSE)
}

use_spatial_obj <- !is.null(opt$spatial_obj)
use_legacy_input <- !is.null(opt$spatial_ex) && !is.null(opt$coord)

if (!use_spatial_obj && !use_legacy_input) {
  stop("Provide either --spatial_obj or both --spatial_ex and --coord.", call. = FALSE)
}

if (use_spatial_obj && use_legacy_input) {
  warning("Both --spatial_obj and legacy matrix/coord inputs were provided. Using --spatial_obj.")
}

#setwd
work_dir <- opt$output
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}
setwd(work_dir)
cat("Working directory set to:", getwd(), "\n")
#load_data
cat("Loading reference data...\n")
load_seurat_object <- function(file_path, data_type) {
  if (!file.exists(file_path)) {
    stop(paste("Input file not found:", file_path))
  }

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
      stop("Package 'qs' needed for this function to work. Please install it.")
    }
    return(qs::qread(file_path))
  }

  if (data_type == "rds") {
    return(readRDS(file_path))
  }

  stop("Unsupported data type. Please use 'rdata', 'qs', or 'rds'.")
}
pbmc_filt <- load_seurat_object(opt$ref, opt$type)

#Ref
if (!opt$anno_column %in% names(pbmc_filt@meta.data)) {
  stop(paste("Annotation column", opt$anno_column, "not found in metadata"))
}
pbmc_filt$subtype  <- pbmc_filt@meta.data[, opt$anno_column]
pbmc_filt$subtype <- as.character(pbmc_filt$subtype)
subtype_counts <- table(pbmc_filt$subtype)
keep_subtypes <- names(subtype_counts)[subtype_counts > 25]
pbmc_filt <- subset(pbmc_filt, subset = subtype %in% keep_subtypes)
table(pbmc_filt$subtype)
cell_types <- as.factor(pbmc_filt$subtype)

scCounts <- as.matrix(GetAssayData(pbmc_filt, slot = "counts"))

#Spatial
cat("Loading spatial data...\n")
if (use_spatial_obj) {
  spatial_obj <- load_seurat_object(opt$spatial_obj, opt$spatial_type)

  if (!inherits(spatial_obj, "Seurat")) {
    stop("The spatial input does not contain a Seurat object.")
  }
  if (!"spatial" %in% names(spatial_obj@reductions)) {
    stop("Reduction 'spatial' not found in the spatial Seurat object.")
  }

  spatial_embeddings <- Embeddings(spatial_obj[["spatial"]])
  if (ncol(spatial_embeddings) < 2) {
    stop("Reduction 'spatial' must contain at least two dimensions for X/Y coordinates.")
  }

  coords <- spatial_embeddings[, 1:2, drop = FALSE]
  colnames(coords) <- c("X", "Y")
  coords <- as.data.frame(coords)

  spCounts <- as.matrix(GetAssayData(spatial_obj, slot = "counts"))
} else {
  if (!file.exists(opt$spatial_ex)) {
    stop(paste("Spatial expression file not found:", opt$spatial_ex))
  }
  Exp <- h5read(opt$spatial_ex, "mat") 
  #获取矩阵值
  spCounts <- Exp$block0_values
  #添加行名
  rownames(spCounts) <- Exp$axis0
  #添加列名
  colnames(spCounts) <- Exp$axis1
  #空间坐标
  coords <- read.csv(opt$coord, row.names = 1) 
}
#check
if (!all(colnames(spCounts) %in% rownames(coords))) {
  warning("Not all spots in expression matrix have coordinates. Subsetting...")
  common_spots <- intersect(colnames(spCounts), rownames(coords))
  if (length(common_spots) < 10) {
    stop("Less than 10 spots match between expression and coordinates")
  }
  spCounts <- spCounts[, common_spots]
  coords <- coords[common_spots, ]
}
coords <- coords[colnames(spCounts), , drop = FALSE]
#取交集
scCounts <- scCounts[rowSums(scCounts) > 3,]
spCounts <- spCounts[rowSums(spCounts) > 3,]
co_genes <- intersect(rownames(scCounts),rownames(spCounts))
scCounts <- scCounts[co_genes,]
spCounts <- spCounts[co_genes,]

##----------RCTD----------##
reference <- Reference(scCounts, cell_types)	#构建RCTD参考数据集
query <- SpatialRNA(coords, spCounts, colSums(spCounts))	#构建空转数据对象

RCTD <- create.RCTD(query, reference, max_cores = as.integer(opt$ncores),UMI_min = 30, UMI_min_sigma = 50)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD,file=paste0("myRCTD.rds"))
result <- RCTD@results$results_df
barcode <- rownames(result)
result <- cbind(barcode,result)
write.table(result,file=paste0("RCTD.txt"),sep="\t",col.names=T,row.names=F,quote=F)
