#!/usr/bin/env Rscript

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript spatial_ecotypes/run_spatialecotyper_all_samples.R [options]",
      "",
      "Purpose:",
      "  Run SpatialEcoTyper separately for each sample in a Seurat QS object.",
      "  Each sample is saved under <outdir>/<sample_safe>/ with scdata, scmeta, results, metadata, object and plots.",
      "",
      "Common options:",
      "  -h, --help                         Show this help message and exit.",
      "  --input-qs=FILE                    Input Seurat object saved by qs2::qs_save.",
      "                                     Default: /mnt/NexData_file/home/xjs/20260310_XE2510095462-01/spatialEcoTyper/all_data_4slides_sample_filtered_reUMAP.qs2",
      "  --outdir=DIR                       Output directory for per-sample folders.",
      "                                     Default: /mnt/NexData_file/home/xjs/20260310_XE2510095462-01/spatialEcoTyper",
      "  --sample-col=COL                   Metadata column containing sample IDs. Default: sample.",
      "  --celltype-col=COL                 Metadata column used as CellType. Default: singler_subtype.",
      "  --samples=S1,S2,S3                 Optional comma-separated subset of samples to run.",
      "  --assay=NAME                       Seurat assay name. Default: RNA.",
      "  --expr-layer=NAME                  Assay layer for expression data. Default: data.",
      "  --spatial-reduction=NAME           Reduction containing spatial coordinates. Default: spatial.",
      "",
      "SpatialEcoTyper options:",
      "  --radius=N                         Spatial neighborhood radius. Default: 50.",
      "  --ncores=N                         Parallel cores. Default: 15.",
      "",
      "Overwrite options:",
      "  --overwrite=true|false             Re-run samples with existing *_SpatialEcoTyper_results.rds.",
      "                                     Default: false.",
      "  --overwrite-plots=true|false       Re-save plots for samples with existing results.",
      "                                     Default: true.",
      "",
      "Plot options:",
      "  --plot-dpi=N                       Plot DPI. Default: 300.",
      "  --dimplot-width=N                  DimPlot width. Default: 7.",
      "  --dimplot-height=N                 DimPlot height. Default: 5.",
      "  --spatial-width=N                  SpatialView width. Default: 7.",
      "  --spatial-height=N                 SpatialView height. Default: 7.",
      "  --abundance-width=N                Cell type abundance plot width. Default: 14.",
      "  --abundance-height=N               Cell type abundance plot height. Default: 8.",
      "  --abundance-legend-ncol=N          Legend columns for abundance plot. Default: 4.",
      "  --max-se-id=N                      Preallocated SE color IDs, from SE0 to SE<N>. Default: 100.",
      "",
      "Examples:",
      "  Rscript spatial_ecotypes/run_spatialecotyper_all_samples.R --help",
      "  Rscript spatial_ecotypes/run_spatialecotyper_all_samples.R --samples=SS01,SS02 --ncores=15",
      "  Rscript spatial_ecotypes/run_spatialecotyper_all_samples.R --radius=50 --overwrite-plots=true",
      "",
      sep = "\n"
    )
  )
}

raw_args <- commandArgs(trailingOnly = TRUE)
if (any(raw_args %in% c("-h", "--help", "help"))) {
  print_usage()
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(SpatialEcoTyper))
suppressPackageStartupMessages(library(Polychrome))

parse_cli <- function(args) {
  opts <- list()
  flags <- character(0)

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }

    arg <- sub("^--", "", arg)
    if (grepl("=", arg, fixed = TRUE)) {
      key <- sub("=.*$", "", arg)
      value <- sub("^[^=]+=", "", arg)
      opts[[key]] <- value
    } else {
      flags <- c(flags, arg)
    }
  }

  opts$flags <- flags
  opts
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

opts <- parse_cli(raw_args)

input_qs <- "/mnt/NexData_file/home/xjs/20260310_XE2510095462-01/spatialEcoTyper/all_data_4slides_sample_filtered_reUMAP.qs2"
out_dir <- "/mnt/NexData_file/home/xjs/20260310_XE2510095462-01/spatialEcoTyper"

sample_col <- "sample"
celltype_col <- "annotation"
assay_name <- "RNA"
expr_layer <- "data"
spatial_reduction <- "spatial"

radius <- 50
ncores <- 15
overwrite <- FALSE
overwrite_plots <- TRUE

plot_dpi <- 300
dimplot_width <- 7
dimplot_height <- 5
spatial_width <- 7
spatial_height <- 7
abundance_width <- 14
abundance_height <- 8
abundance_legend_ncol <- 4
max_se_id <- 100

input_qs <- if (!is.null(opts[["input-qs"]])) opts[["input-qs"]] else input_qs
out_dir <- if (!is.null(opts$outdir)) opts$outdir else out_dir
sample_col <- if (!is.null(opts[["sample-col"]])) opts[["sample-col"]] else sample_col
celltype_col <- if (!is.null(opts[["celltype-col"]])) opts[["celltype-col"]] else celltype_col
assay_name <- if (!is.null(opts$assay)) opts$assay else assay_name
expr_layer <- if (!is.null(opts[["expr-layer"]])) opts[["expr-layer"]] else expr_layer
spatial_reduction <- if (!is.null(opts[["spatial-reduction"]])) opts[["spatial-reduction"]] else spatial_reduction
radius <- if (!is.null(opts$radius)) as.numeric(opts$radius) else radius
ncores <- if (!is.null(opts$ncores)) as.integer(opts$ncores) else ncores
overwrite <- as_bool(opts$overwrite, default = overwrite)
overwrite_plots <- as_bool(opts[["overwrite-plots"]], default = overwrite_plots)
plot_dpi <- if (!is.null(opts[["plot-dpi"]])) as.integer(opts[["plot-dpi"]]) else plot_dpi
dimplot_width <- if (!is.null(opts[["dimplot-width"]])) as.numeric(opts[["dimplot-width"]]) else dimplot_width
dimplot_height <- if (!is.null(opts[["dimplot-height"]])) as.numeric(opts[["dimplot-height"]]) else dimplot_height
spatial_width <- if (!is.null(opts[["spatial-width"]])) as.numeric(opts[["spatial-width"]]) else spatial_width
spatial_height <- if (!is.null(opts[["spatial-height"]])) as.numeric(opts[["spatial-height"]]) else spatial_height
abundance_width <- if (!is.null(opts[["abundance-width"]])) as.numeric(opts[["abundance-width"]]) else abundance_width
abundance_height <- if (!is.null(opts[["abundance-height"]])) as.numeric(opts[["abundance-height"]]) else abundance_height
abundance_legend_ncol <- if (!is.null(opts[["abundance-legend-ncol"]])) {
  as.integer(opts[["abundance-legend-ncol"]])
} else {
  abundance_legend_ncol
}
max_se_id <- if (!is.null(opts[["max-se-id"]])) as.integer(opts[["max-se-id"]]) else max_se_id

sample_include <- NULL
if (!is.null(opts$samples) && nzchar(opts$samples)) {
  sample_include <- strsplit(opts$samples, ",", fixed = TRUE)[[1]]
}

polychrome_palette <- function(n) {
  if (n == 0) {
    return(character(0))
  }

  set.seed(42)
  colors <- createPalette(
    n,
    seedcolors = c("#ff0000", "#00ff00", "#0000ff"),
    M = 50000
  )
  names(colors) <- NULL
  colors
}

save_spatialecotyper_plots <- function(obj, scmeta, sample_dir, sample_safe) {
  dimplot_gg <- DimPlot(obj, group.by = "SE") +
    scale_color_manual(values = global_se_colors)
  ggsave(
    filename = file.path(sample_dir, paste0(sample_safe, "_DimPlot_SE.png")),
    plot = dimplot_gg,
    width = dimplot_width,
    height = dimplot_height,
    dpi = plot_dpi
  )

  spatial_gg <- SpatialView(scmeta, by = "SE")
  spatial_gg <- spatial_gg +
    scale_color_manual(values = global_se_colors) +
    scale_fill_manual(values = global_se_colors)
  ggsave(
    filename = file.path(sample_dir, paste0(sample_safe, "_SpatialView_SE.png")),
    plot = spatial_gg,
    width = spatial_width,
    height = spatial_height,
    dpi = plot_dpi
  )

  abundance_df <- scmeta %>% filter(!is.na(SE)) %>% count(SE, CellType)
  abundance_gg <- ggplot(abundance_df, aes(SE, n, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = global_celltype_colors) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = grid::unit(0.28, "cm"),
      plot.margin = grid::unit(c(0.4, 0.4, 0.4, 0.4), "cm")
    ) +
    guides(fill = guide_legend(ncol = abundance_legend_ncol, byrow = TRUE)) +
    coord_flip() +
    labs(y = "Cell type abundance")
  ggsave(
    filename = file.path(sample_dir, paste0(sample_safe, "_SE_CellType_abundance.png")),
    plot = abundance_gg,
    width = abundance_width,
    height = abundance_height,
    dpi = plot_dpi
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(format(Sys.time()), " Read input: ", input_qs)
adata <- qs2::qs_read(input_qs)

needed_meta_cols <- c(sample_col, celltype_col)
missing_meta_cols <- setdiff(needed_meta_cols, colnames(adata@meta.data))
if (length(missing_meta_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}
if (!assay_name %in% names(adata@assays)) {
  stop("Missing assay: ", assay_name)
}
if (!spatial_reduction %in% names(adata@reductions)) {
  stop("Missing reduction: ", spatial_reduction)
}

global_celltype_levels <- sort(unique(as.character(adata@meta.data[[celltype_col]])))
global_celltype_levels <- global_celltype_levels[!is.na(global_celltype_levels) & global_celltype_levels != ""]
global_celltype_colors <- setNames(polychrome_palette(length(global_celltype_levels)), global_celltype_levels)

global_se_levels <- paste0("SE", seq.int(0, max_se_id))
global_se_colors <- setNames(polychrome_palette(length(global_se_levels)), global_se_levels)

samples <- sort(unique(as.character(adata@meta.data[[sample_col]])))
samples <- samples[!is.na(samples) & samples != ""]
if (!is.null(sample_include)) {
  missing_samples <- setdiff(sample_include, samples)
  if (length(missing_samples) > 0) {
    stop("Requested samples are not present in ", sample_col, ": ", paste(missing_samples, collapse = ", "))
  }
  samples <- sample_include
}

message(format(Sys.time()), " Samples to run: ", paste(samples, collapse = ", "))
message(format(Sys.time()), " Global CellType colors: ", length(global_celltype_colors))
orig_wd <- getwd()

for (sample_id in samples) {
  sample_safe <- gsub("[^A-Za-z0-9_.-]+", "_", sample_id)
  sample_dir <- file.path(out_dir, sample_safe)
  result_rds <- file.path(sample_dir, paste0(sample_safe, "_SpatialEcoTyper_results.rds"))
  metadata_qs <- file.path(sample_dir, paste0(sample_safe, "_SpatialEcoTyper_metadata.qs2"))
  obj_rds <- file.path(sample_dir, paste0(sample_safe, "_SpatialEcoTyper_obj.rds"))
  plot_files <- file.path(
    sample_dir,
    paste0(
      sample_safe,
      c("_DimPlot_SE.png", "_SpatialView_SE.png", "_SE_CellType_abundance.png")
    )
  )

  if (file.exists(result_rds) && !overwrite) {
    if (all(file.exists(plot_files)) && !overwrite_plots) {
      message(format(Sys.time()), " Skip existing sample with plots: ", sample_id, " -> ", result_rds)
      next
    }

    if (file.exists(metadata_qs) && file.exists(obj_rds)) {
      message(format(Sys.time()), " Existing result found, save missing plots: ", sample_id)
      obj <- readRDS(obj_rds)
      scmeta <- qs2::qs_read(metadata_qs)
      save_spatialecotyper_plots(obj, scmeta, sample_dir, sample_safe)
      rm(obj, scmeta)
      gc()
      next
    }

    message(format(Sys.time()), " Skip existing sample without cached metadata for plotting: ", sample_id, " -> ", result_rds)
    next
  }

  message(format(Sys.time()), " Start sample: ", sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  sample_cells <- rownames(adata@meta.data)[as.character(adata@meta.data[[sample_col]]) == sample_id]
  adata_sample <- subset(adata, cells = sample_cells)

  scdata <- GetAssayData(adata_sample, assay = assay_name, layer = expr_layer)

  spatial_coord <- Embeddings(adata_sample, spatial_reduction)
  scmeta <- adata_sample@meta.data
  scmeta$X <- spatial_coord[rownames(scmeta), 1]
  scmeta$Y <- spatial_coord[rownames(scmeta), 2]
  scmeta <- scmeta %>%
    select(X, Y, all_of(celltype_col)) %>%
    rename(CellType = all_of(celltype_col))

  keep_cells <- rownames(scmeta)[
    !is.na(scmeta$X) &
      !is.na(scmeta$Y) &
      !is.na(scmeta$CellType) &
      scmeta$CellType != ""
  ]
  if (length(keep_cells) == 0) {
    message(format(Sys.time()), " Skip sample with no valid cells: ", sample_id)
    next
  }

  scmeta <- scmeta[keep_cells, , drop = FALSE]
  scdata <- scdata[, keep_cells, drop = FALSE]

  qs2::qs_save(scdata, file.path(sample_dir, paste0(sample_safe, "_scdata.qs2")))
  qs2::qs_save(scmeta, file.path(sample_dir, paste0(sample_safe, "_scmeta.qs2")))

  setwd(sample_dir)
  se_results <- tryCatch(
    SpatialEcoTyper(
      scdata,
      scmeta,
      outprefix = sample_safe,
      radius = radius,
      ncores = ncores
    ),
    finally = setwd(orig_wd)
  )

  qs2::qs_save(se_results$metadata, metadata_qs)
  saveRDS(se_results$obj, obj_rds)

  obj <- se_results$obj
  scmeta <- se_results$metadata
  save_spatialecotyper_plots(obj, scmeta, sample_dir, sample_safe)

  message(format(Sys.time()), " Done sample: ", sample_id)
  rm(adata_sample, scdata, scmeta, spatial_coord, se_results, obj)
  gc()
}

setwd(orig_wd)
message(format(Sys.time()), " All done. DEG/marker analysis was not run.")
