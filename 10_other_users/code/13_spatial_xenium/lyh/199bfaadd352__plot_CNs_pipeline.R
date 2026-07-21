#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(scplotter)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(openxlsx)
  library(qs)
  library(jsonlite)
})

DEFAULT_COLORS <- c(
  "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", "#d9534f", "#ffad60",
  "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c",
  "#4695d6", "#fed95c", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879",
  "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83",
  "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d",
  "#e8a0b8", "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473",
  "#f17d80", "#737495"
)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- list(
    work_dir = NULL,
    rdata_path = NULL,
    data_type = "qs",
    cns_df = NULL,
    config_path = NULL,
    heatmap_script = NULL,
    group_col = "Group",
    sample_col = "sample",
    celltype_col = "predicted.celltype",
    cn_col = NULL,
    mode = NULL,
    object_name = NULL,
    output_object = "add_CNS.qs"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1]] else NULL
    if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript plot_CNs_pipeline.R --work_dir <dir> --rdata_path <file> --cns_df <csv> --cn_col <name> [options]\n",
        "Options:\n",
        "  --data_type       qs | rds | rdata (default: qs)\n",
        "  --config_path     JSON config exported by the Python stage\n",
        "  --heatmap_script  Path to heatmap.R used for CN heatmap plotting\n",
        "  --group_col       Metadata column used for stacked group plots (default: Group)\n",
        "  --sample_col      Metadata column used for spatial split plots (default: sample)\n",
        "  --celltype_col    Metadata column used for cell type spatial plots (default: predicted.celltype)\n",
        "  --cn_col          CN column name to plot\n",
        "  --mode            knn or radius. Used for output naming\n",
        "  --object_name     Object name to extract when loading .RData\n",
        "  --output_object   Output qs filename (default: add_CNS.qs)\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    }

    if (!startsWith(key, "--") || is.null(value)) {
      stop(sprintf("Unknown or incomplete parameter: %s", key))
    }

    name <- substring(key, 3)
    if (!name %in% names(parsed)) {
      stop(sprintf("Unknown parameter: %s", key))
    }

    parsed[[name]] <- value
    i <- i + 2
  }

  required <- c("work_dir", "rdata_path", "cns_df")
  missing <- required[vapply(parsed[required], is.null, logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required parameters: %s", paste(missing, collapse = ", ")))
  }

  parsed
}

load_seurat_object <- function(data_path, data_type, object_name = NULL) {
  if (data_type == "qs") {
    return(qs::qread(data_path))
  }
  if (data_type == "rds") {
    return(readRDS(data_path))
  }
  if (data_type == "rdata") {
    env <- new.env(parent = emptyenv())
    loaded_names <- load(data_path, envir = env)
    if (!is.null(object_name)) {
      if (!exists(object_name, envir = env, inherits = FALSE)) {
        stop(sprintf("Object '%s' not found in %s", object_name, data_path))
      }
      return(get(object_name, envir = env, inherits = FALSE))
    }
    seurat_candidates <- loaded_names[vapply(
      loaded_names,
      function(x) inherits(get(x, envir = env, inherits = FALSE), "Seurat"),
      logical(1)
    )]
    if (length(seurat_candidates) != 1) {
      stop("Unable to determine which Seurat object to use from the .RData file. Please provide --object_name.")
    }
    return(get(seurat_candidates[[1]], envir = env, inherits = FALSE))
  }
  stop("Unsupported data_type. Use 'qs', 'rds', or 'rdata'.")
}

named_palette <- function(values) {
  if (is.null(values)) {
    return(setNames(character(0), character(0)))
  }
  clean_values <- sort(unique(as.character(values)))
  if (length(clean_values) <= length(DEFAULT_COLORS)) {
    colors <- DEFAULT_COLORS[seq_along(clean_values)]
  } else {
    colors <- grDevices::hcl.colors(length(clean_values), palette = "Dynamic")
  }
  setNames(colors, clean_values)
}

save_named_palette <- function(palette, output_path) {
  palette_df <- data.frame(label = names(palette), color = unname(palette), stringsAsFactors = FALSE)
  write.csv(palette_df, output_path, row.names = FALSE)
}

merge_cns_metadata <- function(seurat_obj, cns_df, cn_col) {
  common_ids <- intersect(rownames(seurat_obj@meta.data), rownames(cns_df))
  if (length(common_ids) == 0) {
    stop("No overlapping cell IDs were found between the Seurat object and the CN metadata.")
  }
  if (!cn_col %in% colnames(cns_df)) {
    stop(sprintf("CN column '%s' was not found in CN metadata CSV.", cn_col))
  }

  aligned_meta <- seurat_obj@meta.data
  aligned_meta[[cn_col]] <- NA
  aligned_meta[common_ids, cn_col] <- cns_df[common_ids, cn_col]
  seurat_obj@meta.data <- aligned_meta
  seurat_obj
}

prefix_cn_values <- function(seurat_obj, cn_col) {
  if (!cn_col %in% colnames(seurat_obj@meta.data)) {
    return(seurat_obj)
  }
  seurat_obj@meta.data[[cn_col]] <- paste0("CNs_", seurat_obj@meta.data[[cn_col]])
  seurat_obj
}

save_heatmap_plot <- function(seurat_obj, celltype_col, cn_col, output_prefix) {
  if (!exists("cn_exp_heatmap", mode = "function")) {
    stop("cn_exp_heatmap was not loaded. Please check the heatmap.R script path.")
  }
  if (!celltype_col %in% colnames(seurat_obj@meta.data) || !cn_col %in% colnames(seurat_obj@meta.data)) {
    message(sprintf("Skip heatmap because %s or %s is missing.", celltype_col, cn_col))
    return(invisible(NULL))
  }

  plot_df <- seurat_obj@meta.data[, c(celltype_col, cn_col), drop = FALSE]
  plot_df[[celltype_col]] <- as.character(plot_df[[celltype_col]])
  plot_df[[cn_col]] <- as.character(plot_df[[cn_col]])
  p <- cn_exp_heatmap(plot_df, cluster_col = celltype_col, cn_col = cn_col)
  ggsave(paste0(output_prefix, ".pdf"), p, height = 5, width = 8)
  ggsave(paste0(output_prefix, ".png"), p, height = 5, width = 8, dpi = 300)
}

save_group_summary_plot <- function(seurat_obj, group_col, cn_col, palette, output_prefix) {
  if (!group_col %in% colnames(seurat_obj@meta.data) || !cn_col %in% colnames(seurat_obj@meta.data)) {
    message(sprintf("Skip group summary for %s because %s or %s is missing.", output_prefix, group_col, cn_col))
    return(invisible(NULL))
  }

  sum_table <- table(seurat_obj@meta.data[[group_col]], seurat_obj@meta.data[[cn_col]])
  sum_pct <- as.data.frame(t(sum_table)) %>%
    group_by(Var2) %>%
    mutate(Freq_pct = Freq / sum(Freq) * 100) %>%
    ungroup()
  sum_pct_matrix <- sum_pct %>%
    select(Var1, Var2, Freq_pct) %>%
    pivot_wider(names_from = "Var2", values_from = "Freq_pct") %>%
    column_to_rownames("Var1") %>%
    t() %>%
    as.data.frame()

  openxlsx::write.xlsx(as.data.frame.matrix(sum_table), paste0(output_prefix, ".xlsx"), rowNames = TRUE)
  openxlsx::write.xlsx(sum_pct_matrix, paste0(output_prefix, "_pct.xlsx"), rowNames = TRUE)

  p <- ggplot(sum_pct, aes(x = Var2, y = Freq_pct, fill = Var1)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(name = "CNs", values = palette) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black", angle = 45, hjust = 1)
    ) +
    labs(x = "", y = "Percentage", title = "")

  ggsave(paste0(output_prefix, ".pdf"), p, height = 6, width = 6)
  ggsave(paste0(output_prefix, ".png"), p, height = 6, width = 6, dpi = 300)
}

save_spatial_plots <- function(seurat_obj, sample_col, cn_col, celltype_col, cn_palette, celltype_palette) {
  if (!sample_col %in% colnames(seurat_obj@meta.data)) {
    message(sprintf("Skip spatial plots because %s is missing.", sample_col))
    return(invisible(NULL))
  }

  dir.create("./spatial_plot_CNs", showWarnings = FALSE)

  clean_theme <- theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

  for (sample_id in unique(seurat_obj@meta.data[[sample_col]])) {
    sample_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[sample_col]] == sample_id]
    tmp <- subset(seurat_obj, cells = sample_cells)

    p_list <- list()
    if (cn_col %in% colnames(tmp@meta.data)) {
      p_list[[length(p_list) + 1]] <- CellDimPlot(
        tmp,
        group_by = cn_col,
        reduction = "spatial",
        raster = FALSE,
        palcolor = cn_palette,
        pt_size = 0.03,
        theme = "theme_bw",
        aspect.ratio = NULL
      ) + clean_theme
    }
    if (celltype_col %in% colnames(tmp@meta.data)) {
      p_list[[length(p_list) + 1]] <- CellDimPlot(
        tmp,
        group_by = celltype_col,
        reduction = "spatial",
        raster = FALSE,
        palcolor = celltype_palette,
        pt_size = 0.03,
        theme = "theme_bw",
        aspect.ratio = NULL
      ) + clean_theme
    }

    if (length(p_list) == 0) {
      next
    }

    plot_width <- max(6, length(p_list) * 6)
    combined_plot <- wrap_plots(p_list, ncol = length(p_list))
    ggsave(paste0("./spatial_plot_CNs/", sample_id, "_CNs.pdf"), combined_plot, height = 6, width = plot_width, device = cairo_pdf)
    ggsave(paste0("./spatial_plot_CNs/", sample_id, "_CNs.png"), combined_plot, height = 6, width = plot_width, dpi = 300)
  }
}

main <- function() {
  args <- parse_args()
  work_dir <- normalizePath(args$work_dir, mustWork = FALSE)
  data_path <- normalizePath(args$rdata_path, mustWork = TRUE)
  cns_path <- normalizePath(args$cns_df, mustWork = TRUE)

  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(work_dir)
  message("Working directory: ", getwd())

  if (!is.null(args$heatmap_script)) {
    source(normalizePath(args$heatmap_script, mustWork = TRUE))
  }

  config <- NULL
  if (!is.null(args$config_path)) {
    config <- jsonlite::fromJSON(normalizePath(args$config_path, mustWork = TRUE))
  }

  cn_col <- if (!is.null(args$cn_col)) args$cn_col else config$cn_column
  mode <- if (!is.null(args$mode)) args$mode else config$mode
  if (is.null(cn_col) || !nzchar(cn_col)) {
    stop("CN column is required. Please provide --cn_col or a valid --config_path.")
  }
  if (is.null(mode) || !nzchar(mode)) {
    mode <- "cn"
  }

  seurat_obj <- load_seurat_object(data_path, tolower(args$data_type), args$object_name)
  cns_df <- read.csv(cns_path, row.names = 1, check.names = FALSE)

  seurat_obj <- merge_cns_metadata(seurat_obj, cns_df, cn_col = cn_col)
  seurat_obj <- prefix_cn_values(seurat_obj, cn_col = cn_col)

  cn_palette <- if (cn_col %in% colnames(seurat_obj@meta.data)) named_palette(seurat_obj@meta.data[[cn_col]]) else setNames(character(0), character(0))
  celltype_palette <- if (args$celltype_col %in% colnames(seurat_obj@meta.data)) named_palette(seurat_obj@meta.data[[args$celltype_col]]) else setNames(character(0), character(0))

  save_named_palette(cn_palette, paste0(mode, "_palette.csv"))
  save_named_palette(celltype_palette, "celltype_palette.csv")

  save_heatmap_plot(
    seurat_obj,
    celltype_col = args$celltype_col,
    cn_col = cn_col,
    output_prefix = paste0("CNs_heatmap_", mode)
  )
  save_group_summary_plot(
    seurat_obj,
    group_col = args$group_col,
    cn_col = cn_col,
    palette = cn_palette,
    output_prefix = paste0("CNs_Stacked_group_", mode)
  )
  save_spatial_plots(
    seurat_obj,
    sample_col = args$sample_col,
    cn_col = cn_col,
    celltype_col = args$celltype_col,
    cn_palette = cn_palette,
    celltype_palette = celltype_palette
  )

  qs::qsave(seurat_obj, args$output_object)
  message("Saved updated object to: ", file.path(work_dir, args$output_object))
}

main()
