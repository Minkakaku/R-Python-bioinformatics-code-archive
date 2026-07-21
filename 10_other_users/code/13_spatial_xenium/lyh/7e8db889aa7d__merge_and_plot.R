library(optparse)
library(Seurat)
library(qs)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(parallel)
library(plotthis)

safe_library <- function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

safe_library("ggpubr")

COLOR_SCRIPT <- "/mnt/NexData_file/home/lyh/Script/R/R_fuc/color.R"

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the original Xenium/Seurat object", metavar = "FILE"),
  make_option(c("-t", "--input_type"), type = "character", default = "rds",
              help = "Input object type [rdata/rds/qs] (default: rds)", metavar = "TYPE"),
  make_option(c("-r", "--result_dir"), type = "character", default = NULL,
              help = "Directory containing per-sample *_RCTD.txt files", metavar = "DIR"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory for merged object and plots", metavar = "DIR"),
  make_option(c("-q", "--qs_output"), type = "character", default = "Xenium_RCTD.qs",
              help = "Output filename for the merged Xenium object in qs format (default: Xenium_RCTD.qs)", metavar = "FILE"),
  make_option(c("-m", "--merged_txt"), type = "character", default = "RCTR.txt",
              help = "Merged RCTD table filename to write under output_dir (default: RCTR.txt)", metavar = "FILE"),
  make_option(c("-s", "--sample_column"), type = "character", default = "sample",
              help = "Sample column in metadata (default: sample)", metavar = "COL"),
  make_option(c("-g", "--group_column"), type = "character", default = "group",
              help = "Group column in metadata for group plots (default: group)", metavar = "COL"),
  make_option(c("-c", "--celltype_column"), type = "character", default = "first_type",
              help = "Metadata column name used to store merged RCTD label (default: first_type)", metavar = "COL"),
  make_option(c("-n", "--ncores"), type = "integer", default = 8,
              help = "Cores used for spatial plotting on Linux (default: 8)", metavar = "NUM")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("input", "result_dir", "output_dir")
missing_args <- sapply(required_args, function(x) is.null(opt[[x]]))
if (any(missing_args)) {
  missing_names <- paste(required_args[missing_args], collapse = ", ")
  stop(paste("Missing required arguments:", missing_names), call. = FALSE)
}

if (!dir.exists(opt$result_dir)) {
  stop(paste("Result directory not found:", opt$result_dir), call. = FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("Input object not found:", opt$input), call. = FALSE)
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$output_dir, "plot")
spatial_plot_dir <- file.path(plot_dir, "spatial_plot")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_plot_dir, recursive = TRUE, showWarnings = FALSE)

load_seurat_object <- function(file_path, data_type) {
  if (data_type == "rdata") {
    env <- new.env(parent = emptyenv())
    loaded_names <- load(file_path, envir = env)
    seurat_candidates <- loaded_names[sapply(loaded_names, function(name) inherits(env[[name]], "Seurat"))]
    if (length(seurat_candidates) == 0) {
      stop("No Seurat object found in the RData file.")
    }
    return(env[[seurat_candidates[[1]]]])
  }

  if (data_type == "qs") {
    return(qs::qread(file_path))
  }

  if (data_type == "rds") {
    return(readRDS(file_path))
  }

  stop("Unsupported data type. Please use 'rdata', 'rds', or 'qs'.")
}

celltype_palette <- function(labels) {
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels)]
  if (!exists("color3", inherits = TRUE)) {
    stop("color3 was not found after sourcing color.R", call. = FALSE)
  }
  if (length(labels) <= length(color3)) {
    palette_values <- color3[seq_along(labels)]
  } else {
    extra_n <- length(labels) - length(color3)
    extra_colors <- grDevices::hcl.colors(extra_n, palette = "Dynamic")
    palette_values <- c(color3, extra_colors)
  }
  setNames(palette_values[seq_along(labels)], labels)
}

group_palette <- function(labels) {
  labels <- unique(as.character(labels))
  labels <- labels[!is.na(labels)]
  if (!exists("color4", inherits = TRUE)) {
    stop("color4 was not found after sourcing color.R", call. = FALSE)
  }
  if (length(labels) <= length(color4)) {
    palette_values <- color4[seq_along(labels)]
  } else {
    extra_n <- length(labels) - length(color4)
    extra_colors <- grDevices::hcl.colors(extra_n, "Set 2")
    palette_values <- c(color4, extra_colors)
  }
  setNames(palette_values[seq_along(labels)], labels)
}

save_plot_dual <- function(plot_obj, pdf_path, png_path, width = 8, height = 6) {
  ggsave(pdf_path, plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plot_obj, width = width, height = height, dpi = 300)
}

read_result_tables <- function(result_dir) {
  result_files <- list.files(result_dir, pattern = "_RCTD\\.txt$", full.names = TRUE)
  if (length(result_files) == 0) {
    stop("No per-sample *_RCTD.txt files were found.", call. = FALSE)
  }

  result_list <- lapply(result_files, function(path) {
    read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  })

  bind_rows(result_list)
}

plot_spatial_maps <- function(obj, sample_column, celltype_column, colors, output_dir, ncores) {
  if (!sample_column %in% colnames(obj@meta.data)) {
    warning("Sample column not found. Skipping spatial plots.")
    return(invisible(NULL))
  }

  samples <- unique(as.character(obj@meta.data[[sample_column]]))
  samples <- samples[!is.na(samples)]
  runner <- if (.Platform$OS.type == "windows") lapply else function(X, FUN) parallel::mclapply(X, FUN, mc.cores = ncores)

  runner(samples, function(sample_id) {
    if (!sample_id %in% names(obj@images)) {
      warning(paste("Image slot not found for sample:", sample_id, "- skipping spatial image plot."))
      return(NULL)
    }

    image_obj <- obj@images[[sample_id]]
    if ("segmentation" %in% names(image_obj@boundaries)) {
      DefaultBoundary(image_obj) <- "segmentation"
      obj@images[[sample_id]] <- image_obj
    }

    p <- ImageDimPlot(
      obj,
      fov = sample_id,
      nmols = 20000,
      dark.background = FALSE,
      border.color = "black",
      border.size = 0.1,
      cols = colors,
      group.by = celltype_column
    ) + labs(title = paste0(sample_id, " RCTD Subtype"))

    pdf_path <- file.path(output_dir, paste0(sample_id, "_spatial.pdf"))
    png_path <- file.path(output_dir, paste0(sample_id, "_spatial.png"))
    save_plot_dual(p, pdf_path, png_path, width = 12, height = 12)
    NULL
  })
}

plot_umap <- function(obj, group_by, colors, split_by = NULL, output_prefix = "celltype_umap", output_dir = ".") {
  if (!"umap" %in% names(obj@reductions)) {
    warning("UMAP reduction not found. Skipping UMAP plots.")
    return(invisible(NULL))
  }

  if (!requireNamespace("scplotter", quietly = TRUE)) {
    warning("Package 'scplotter' is not installed. Skipping UMAP plots.")
    return(invisible(NULL))
  }

  p <- scplotter::CellDimPlot(
    obj,
    group_by = group_by,
    split_by = split_by,
    reduction = "umap",
    label = FALSE,
    palcolor = colors,
    pt_size = 0.1
  ) &
    theme_bw() &
    theme(
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )

  width <- if (is.null(split_by)) 9 else 19
  save_plot_dual(
    p,
    file.path(output_dir, paste0(output_prefix, ".pdf")),
    file.path(output_dir, paste0(output_prefix, ".png")),
    width = width,
    height = 6
  )
}

plot_group_stacked <- function(obj, group_column, celltype_column, colors, output_dir) {
  sum_tbl <- table(obj@meta.data[[group_column]], obj@meta.data[[celltype_column]])
  write.xlsx(sum_tbl, file.path(output_dir, "Celltype_sum_group.xlsx"), rowNames = TRUE)

  sum_pct <- as.data.frame(t(sum_tbl)) %>%
    group_by(Var2) %>%
    mutate(Freq_pct = Freq / sum(Freq) * 100) %>%
    dplyr::select(-Freq)

  sum_pct_wide <- sum_pct %>%
    pivot_wider(names_from = Var2, values_from = Freq_pct) %>%
    column_to_rownames("Var1") %>%
    t() %>%
    as.data.frame()
  write.xlsx(sum_pct_wide, file.path(output_dir, "Celltype_sum_group_pct.xlsx"), rowNames = TRUE)

  plot_df <- sum_pct %>%
    rename(celltype = Var1, group = Var2, percentage = Freq_pct)
  p <- TrendPlot(plot_df,x=group,y=percentage,group_by = celltype,palcolor=colors)+
    labs(x = "Group", y = "Cell Proportion (%)", fill = "Cell Type") 

  save_plot_dual(
    p,
    file.path(output_dir, "celltype_Stacked_group.pdf"),
    file.path(output_dir, "celltype_Stacked_group.png"),
    width = 9,
    height = 6
  )
}

plot_sample_proportion <- function(obj, sample_column, group_column, celltype_column, group_colors, output_dir) {
  metadata_df <- obj@meta.data %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::transmute(
      Cell = .data[["Cell"]],
      sample = .data[[sample_column]],
      group = .data[[group_column]],
      celltype = .data[[celltype_column]]
    )

  cell_num_df <- metadata_df %>%
    dplyr::count(sample, group, name = "Num")

  ratio_df <- metadata_df %>%
    dplyr::count(sample, group, celltype, name = "Sum") %>%
    tidyr::complete(
      tidyr::nesting(sample, group),
      celltype,
      fill = list(Sum = 0)
    ) %>%
    dplyr::left_join(cell_num_df, by = c("sample", "group")) %>%
    dplyr::mutate(Percentage = Sum / Num)

  wide_df <- ratio_df %>%
    dplyr::select(sample, group, celltype, Percentage) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = Percentage)

  openxlsx::write.xlsx(wide_df, file.path(output_dir, "Cell_Proportion_sample.xlsx"))

  p1 <- ggplot2::ggplot(ratio_df, ggplot2::aes(x = celltype, y = Percentage, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.25, position = ggplot2::position_dodge(width = 0.5), size = 0.1, colour = "#FFFFFF") +
    ggplot2::geom_point(pch = 21, position = ggplot2::position_jitterdodge(), colour = "white", size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggpubr::stat_compare_means(
      ggplot2::aes(group = group),
      method = "t.test",
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
      label = "p.signif",
      size = 3
    ) +
    ggplot2::labs(title = "Comparison of Cell Proportion", x = "Cell Type", y = "Cell Proportion") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title = ggplot2::element_text(size = 10),
      legend.key.size = grid::unit(0.4, "cm")
    )

  save_plot_dual(
    p1,
    file.path(output_dir, "Cell_Proportion_BoxPlot_sign1.pdf"),
    file.path(output_dir, "Cell_Proportion_BoxPlot_sign1.png"),
    width = 12,
    height = 5
  )

  dodge_width <- 0.8
  p2 <- ggplot2::ggplot(ratio_df, ggplot2::aes(x = celltype, y = Percentage, fill = group)) +
    ggplot2::geom_bar(stat = "summary", fun = "mean", position = ggplot2::position_dodge(width = dodge_width), width = 0.7, color = "black", size = 0.4) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      position = ggplot2::position_dodge(width = dodge_width),
      width = 0.25,
      size = 0.5
    ) +
    ggplot2::geom_point(
      pch = 21,
      position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = dodge_width),
      colour = "white",
      stroke = 0.2,
      size = 2,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    ggpubr::stat_compare_means(
      ggplot2::aes(group = group),
      method = "t.test",
      label = "p.signif",
      hide.ns = TRUE,
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),
      vjust = -0.5,
      size = 6
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(x = NULL, y = "Cell Proportion") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )

  save_plot_dual(
    p2,
    file.path(output_dir, "Cell_Proportion_BoxPlot_sign2.pdf"),
    file.path(output_dir, "Cell_Proportion_BoxPlot_sign2.png"),
    width = 12,
    height = 5
  )
}

plot_marker_dot <- function(obj, marker_df, output_dir, continuous_palette, celltype_column) {
  marker_df <- marker_df %>% distinct(gene, .keep_all = TRUE) %>% group_by(cluster) %>% filter(pct.1 >= 0.1) %>% slice_max(avg_log2FC, n = 3)
  marker_df$cluster <- as.factor(marker_df$cluster)

  p <- NULL
  if (requireNamespace("scRNAtoolVis", quietly = TRUE)) {
    p <- scRNAtoolVis::jjDotPlot(
      object = obj,
      markerGene = marker_df,
      anno = TRUE,
      ytree = FALSE,
      xtree = FALSE,
      id = celltype_column,
      textRot = 45,
      textSize = 14,
      lwd = 1.5,
      plot.margin = c(4, 1, 1, 1),
      same.pos.label = TRUE,
      segWidth = 0.8,
      midpoint = 0,
      col.min = -2.5,
      col.max = 2.5
    ) + scale_fill_gradientn(colours = continuous_palette)
  } else {
    features <- unique(marker_df$gene)
    p <- DotPlot(obj, features = features, group.by = celltype_column) +
      RotatedAxis() +
      scale_color_gradientn(colours = continuous_palette)
  }

  save_plot_dual(
    p,
    file.path(output_dir, "Marker_DotPlot.pdf"),
    file.path(output_dir, "Marker_DotPlot.png"),
    width = 16,
    height = 8
  )
}

plot_marker_volcano <- function(marker_df, output_dir, label_colors) {
  p <- NULL
  if (requireNamespace("scRNAtoolVis", quietly = TRUE)) {
    p <- scRNAtoolVis::markerVolcano(marker_df, labelCol = label_colors, topn = 3)
  } else {
    marker_df <- marker_df %>%
      mutate(log10_padj = -log10(p_val_adj + 1e-300)) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>%
      ungroup()

    p <- ggplot(marker_df, aes(x = avg_log2FC, y = log10_padj, color = cluster)) +
      geom_point(alpha = 0.7, size = 1.2) +
      geom_text(aes(label = gene), size = 2.5, check_overlap = TRUE, vjust = -0.3) +
      theme_classic() +
      labs(x = "avg_log2FC", y = "-log10(adj p value)")
  }

  save_plot_dual(
    p,
    file.path(output_dir, "Marker_Volcano.pdf"),
    file.path(output_dir, "Marker_Volcano.png"),
    width = 25,
    height = 8
  )
}

message("Loading original Xenium object...")
pbmc_filt <- load_seurat_object(opt$input, opt$input_type)
if (!inherits(pbmc_filt, "Seurat")) {
  stop("Input file does not contain a Seurat object.", call. = FALSE)
}

if (!file.exists(COLOR_SCRIPT)) {
  stop(paste("Color script not found:", COLOR_SCRIPT), call. = FALSE)
}
source(COLOR_SCRIPT)

message("Reading per-sample RCTD result tables...")
rctd_result <- read_result_tables(opt$result_dir)
if (!"barcode" %in% colnames(rctd_result) || !"first_type" %in% colnames(rctd_result)) {
  stop("Merged RCTD table must contain 'barcode' and 'first_type' columns.", call. = FALSE)
}

rctd_result <- rctd_result %>% distinct(barcode, .keep_all = TRUE)
write.table(
  rctd_result,
  file = file.path(opt$output_dir, opt$merged_txt),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message("Merging RCTD labels into Xenium metadata...")
meta_df <- pbmc_filt@meta.data %>%
  rownames_to_column("barcodes") %>%
  left_join(rctd_result %>% rename(barcodes = barcode), by = "barcodes") %>%
  column_to_rownames("barcodes")
pbmc_filt@meta.data <- meta_df

if (!opt$celltype_column %in% colnames(pbmc_filt@meta.data)) {
  if ("first_type" %in% colnames(pbmc_filt@meta.data)) {
    pbmc_filt@meta.data[[opt$celltype_column]] <- pbmc_filt@meta.data[["first_type"]]
  } else {
    stop("Merged metadata does not contain the requested celltype column.", call. = FALSE)
  }
}

qs_output_path <- file.path(opt$output_dir, opt$qs_output)
message("Saving merged Xenium object to: ", qs_output_path)
qs::qsave(pbmc_filt, qs_output_path)

celltype_values <- pbmc_filt@meta.data[[opt$celltype_column]]
celltype_values <- as.character(celltype_values)
celltype_values[is.na(celltype_values) | celltype_values == ""] <- "Unknown"
pbmc_filt@meta.data[[opt$celltype_column]] <- factor(celltype_values, levels = unique(celltype_values))
celltype_colors <- celltype_palette(pbmc_filt@meta.data[[opt$celltype_column]])
continuous_palette <- if (exists("color_s2", inherits = TRUE)) color_s2 else grDevices::hcl.colors(11, "Blue-Red 3")

message("Generating spatial plots...")
plot_spatial_maps(pbmc_filt, opt$sample_column, opt$celltype_column, celltype_colors, spatial_plot_dir, opt$ncores)

message("Generating UMAP plots if available...")
plot_umap(pbmc_filt, opt$celltype_column, celltype_colors, NULL, "celltype_umap", plot_dir)
if (opt$group_column %in% colnames(pbmc_filt@meta.data)) {
  pbmc_filt@meta.data[[opt$group_column]] <- factor(pbmc_filt@meta.data[[opt$group_column]])
  plot_umap(pbmc_filt, opt$celltype_column, celltype_colors, opt$group_column, "celltype_umap_split", plot_dir)
}

if (opt$group_column %in% colnames(pbmc_filt@meta.data)) {
  message("Generating group summary plots...")
  group_colors <- group_palette(pbmc_filt@meta.data[[opt$group_column]])
  plot_group_stacked(pbmc_filt, opt$group_column, opt$celltype_column, celltype_colors, plot_dir)
  plot_sample_proportion(pbmc_filt, opt$sample_column, opt$group_column, opt$celltype_column, group_colors, plot_dir)
} else {
  warning("Group column not found. Group-level plots were skipped.")
}

message("Running marker analysis...")
Idents(pbmc_filt) <- opt$celltype_column
obj_markers <- FindAllMarkers(pbmc_filt, only.pos = TRUE, return.thresh = 0.01, logfc.threshold = 0.26)
write.xlsx(obj_markers, file.path(plot_dir, "Find_all_markers.xlsx"))

plot_marker_dot(pbmc_filt, obj_markers, plot_dir, continuous_palette, opt$celltype_column)
plot_marker_volcano(obj_markers, plot_dir, celltype_colors)

message("Post-processing finished successfully.")
