# Infiltration zone / spatial proximity analysis
#
# Input: a data.frame with cell type, sample, region/group, X, and Y columns.
# Method: for each sample, calculate each cell's nearest distance to cells in the
# specified target region, then assign distance-band labels.
#
# Required for distance calculation:
# install.packages("FNN")
#
# Required only for plotting:
# install.packages("ggplot2")

require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required. Install it with install.packages('",
      pkg, "').",
      call. = FALSE
    )
  }
}

check_columns <- function(dat, cols) {
  missing_cols <- setdiff(cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

make_one_col_formula <- function(col) {
  safe_col <- gsub("`", "``", col, fixed = TRUE)
  stats::as.formula(paste0("~ `", safe_col, "`"))
}

format_break_value <- function(x) {
  if (isTRUE(all.equal(x, round(x)))) {
    as.character(as.integer(round(x)))
  } else {
    as.character(x)
  }
}

make_distance_breaks <- function(interval = 50, max_distance = 500, breaks = NULL) {
  if (is.null(breaks)) {
    if (!is.numeric(interval) || length(interval) != 1 || interval <= 0) {
      stop("'interval' must be one positive number.", call. = FALSE)
    }
    if (!is.numeric(max_distance) || length(max_distance) != 1 || max_distance <= 0) {
      stop("'max_distance' must be one positive number.", call. = FALSE)
    }
    breaks <- seq(0, max_distance, by = interval)
    if (tail(breaks, 1) < max_distance) {
      breaks <- c(breaks, max_distance)
    }
  }

  if (!is.numeric(breaks) || length(breaks) < 2) {
    stop("'breaks' must be a numeric vector with at least two values.", call. = FALSE)
  }
  breaks <- sort(unique(breaks))
  if (breaks[1] != 0) {
    breaks <- c(0, breaks)
  }
  if (any(diff(breaks) <= 0)) {
    stop("'breaks' must be strictly increasing.", call. = FALSE)
  }
  breaks
}

make_band_labels <- function(breaks) {
  paste0(
    vapply(head(breaks, -1), format_break_value, character(1)),
    "_to_",
    vapply(tail(breaks, -1), format_break_value, character(1))
  )
}

assign_distance_bands <- function(distance, breaks, band_labels, outside_label) {
  as.character(cut(
    distance,
    breaks = c(breaks, Inf),
    labels = c(band_labels, outside_label),
    right = FALSE,
    include.lowest = TRUE
  ))
}

make_zone_summary <- function(dat, sample_col, celltype_col, zone_col, sample_group_col = NULL) {
  group_cols <- character(0)
  if (!is.null(sample_group_col)) {
    check_columns(dat, sample_group_col)
    group_cols <- sample_group_col

    sample_group_map <- unique(dat[c(sample_col, sample_group_col)])
    sample_group_n <- aggregate(
      sample_group_map[[sample_group_col]],
      by = sample_group_map[sample_col],
      FUN = function(x) length(unique(x))
    )
    if (any(sample_group_n$x > 1)) {
      bad_samples <- sample_group_n[[sample_col]][sample_group_n$x > 1]
      stop(
        "Each sample must map to one group in '", sample_group_col,
        "'. Multiple groups found for sample(s): ",
        paste(bad_samples, collapse = ", "),
        call. = FALSE
      )
    }
  }

  tab <- as.data.frame(table(
    sample = dat[[sample_col]],
    zone = dat[[zone_col]],
    celltype = dat[[celltype_col]],
    useNA = "no"
  ))
  names(tab) <- c(sample_col, zone_col, celltype_col, "n")
  tab <- tab[tab$n > 0, , drop = FALSE]

  if (length(group_cols) > 0) {
    sample_group_map <- unique(dat[c(sample_col, group_cols)])
    tab <- merge(tab, sample_group_map, by = sample_col, all.x = TRUE)
  }

  zone_totals <- aggregate(
    tab$n,
    by = tab[c(group_cols, sample_col, zone_col)],
    FUN = sum
  )
  names(zone_totals)[names(zone_totals) == "x"] <- "zone_n"

  sample_totals <- aggregate(
    zone_totals$zone_n,
    by = zone_totals[c(group_cols, sample_col)],
    FUN = sum
  )
  names(sample_totals)[names(sample_totals) == "x"] <- "sample_n"

  zone_summary <- merge(zone_totals, sample_totals, by = c(group_cols, sample_col), all.x = TRUE)
  zone_summary$zone_prop_in_sample <- zone_summary$zone_n / zone_summary$sample_n

  celltype_summary <- merge(tab, zone_totals, by = c(group_cols, sample_col, zone_col), all.x = TRUE)
  celltype_summary$celltype_prop_in_zone <- celltype_summary$n / celltype_summary$zone_n

  list(
    zone_summary = zone_summary,
    celltype_summary = celltype_summary
  )
}

summarize_celltype_percent_by_sample <- function(
    celltype_summary,
    sample_col = "sample",
    fov_col = "fov",
    zone_col = "Distance_Group",
    celltype_col = "rctd_dominant_type",
    count_col = "n",
    zone_n_col = "zone_n",
    group_col = NULL,
    sample_group = NULL,
    percent_col = "celltype_percent_in_zone") {
  check_columns(celltype_summary, c(
    sample_col,
    fov_col,
    zone_col,
    celltype_col,
    count_col,
    zone_n_col
  ))

  dat <- as.data.frame(celltype_summary)
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(dat)) {
      if (is.null(sample_group)) {
        stop(
          "Column '", group_col, "' is not in celltype_summary. ",
          "Please provide sample_group with columns '", sample_col, "' and '",
          group_col, "'.",
          call. = FALSE
        )
      }
      check_columns(sample_group, c(sample_col, group_col))
      sample_group <- unique(as.data.frame(sample_group)[c(sample_col, group_col)])
      dat <- merge(dat, sample_group, by = sample_col, all.x = TRUE)
    }
    if (any(is.na(dat[[group_col]]))) {
      stop("Some samples have no group value after merging sample_group.", call. = FALSE)
    }
  }

  if (is.factor(celltype_summary[[zone_col]])) {
    zone_levels <- levels(celltype_summary[[zone_col]])
  } else {
    zone_values <- unique(as.character(dat[[zone_col]]))
    zone_start <- suppressWarnings(as.numeric(sub("_to_.*$", "", zone_values)))
    zone_levels <- zone_values[order(
      is.na(zone_start),
      zone_start,
      zone_values
    )]
  }
  dat[[zone_col]] <- factor(as.character(dat[[zone_col]]), levels = zone_levels)

  count_by <- c(sample_col, zone_col, celltype_col)
  if (!is.null(group_col)) {
    count_by <- c(group_col, count_by)
  }

  count_summary <- aggregate(
    dat[[count_col]],
    by = dat[count_by],
    FUN = sum
  )
  names(count_summary)[names(count_summary) == "x"] <- "n"

  zone_cols <- c(sample_col, fov_col, zone_col, zone_n_col)
  if (!is.null(group_col)) {
    zone_cols <- c(group_col, zone_cols)
  }
  zone_dat <- unique(dat[zone_cols])
  zone_by <- c(sample_col, zone_col)
  if (!is.null(group_col)) {
    zone_by <- c(group_col, zone_by)
  }
  zone_summary <- aggregate(
    zone_dat[[zone_n_col]],
    by = zone_dat[zone_by],
    FUN = sum
  )
  names(zone_summary)[names(zone_summary) == "x"] <- "zone_n"

  merge_by <- c(sample_col, zone_col)
  if (!is.null(group_col)) {
    merge_by <- c(group_col, merge_by)
  }
  out <- merge(count_summary, zone_summary, by = merge_by, all.x = TRUE)
  out[[percent_col]] <- out$n / out$zone_n * 100
  out[[zone_col]] <- factor(as.character(out[[zone_col]]), levels = zone_levels)

  order_cols <- c(sample_col, zone_col, celltype_col)
  if (!is.null(group_col)) {
    order_cols <- c(group_col, order_cols)
  }
  out <- out[do.call(order, out[order_cols]), , drop = FALSE]
  rownames(out) <- NULL

  out
}

# Main function ---------------------------------------------------------------
#
# dat:
#   data.frame containing cell coordinates and annotations.
# region_values:
#   value(s) in group_col used as the target region. Example: "TLS".
# interval:
#   distance-band width. Default is 50.
# max_distance:
#   maximum distance to split into bands. Cells at or beyond this value are
#   assigned to outside_label, for example ">500".
# breaks:
#   optional custom break vector. Example: c(0, 50, 100, 200, 300, 500).
#
# Boundary rule:
#   distance bands are left-closed and right-open: [0,50), [50,100), ...
#   target-region cells are always assigned to Inside_<region_label>.
infiltration_zone_analysis <- function(
    dat,
    celltype_col,
    sample_col,
    group_col,
    x_col = "X",
    y_col = "Y",
    region_values,
    region_label = NULL,
    sample_group_col = NULL,
    interval = 50,
    max_distance = 500,
    breaks = NULL,
    distance_col = "distance_to_region",
    zone_col = "Distance_Group",
    no_region_label = "No_target_region") {
  require_namespace("FNN")

  if (!is.data.frame(dat)) {
    stop("'dat' must be a data.frame.", call. = FALSE)
  }
  check_columns(dat, c(celltype_col, sample_col, group_col, x_col, y_col))
  if (!is.null(sample_group_col)) {
    check_columns(dat, sample_group_col)
  }

  if (missing(region_values) || length(region_values) == 0) {
    stop("'region_values' must contain at least one target region value.", call. = FALSE)
  }

  if (is.null(region_label)) {
    region_label <- if (length(region_values) == 1) {
      as.character(region_values[[1]])
    } else {
      "Region"
    }
  }
  region_label <- gsub("[^A-Za-z0-9_]+", "_", region_label)

  breaks <- make_distance_breaks(
    interval = interval,
    max_distance = max_distance,
    breaks = breaks
  )
  band_labels <- make_band_labels(breaks)
  outside_label <- paste0(">", format_break_value(tail(breaks, 1)))
  inside_label <- paste0("Inside_", region_label)
  zone_levels <- c(inside_label, band_labels, outside_label, no_region_label)

  out <- as.data.frame(dat)
  out[[distance_col]] <- NA_real_
  out[[zone_col]] <- NA_character_

  target_flag <- out[[group_col]] %in% region_values
  samples <- unique(out[[sample_col]])
  empty_samples <- character(0)

  for (sample_id in samples) {
    sample_idx <- which(out[[sample_col]] == sample_id)
    sample_target_idx <- sample_idx[target_flag[sample_idx]]

    if (length(sample_target_idx) == 0) {
      out[[zone_col]][sample_idx] <- no_region_label
      empty_samples <- c(empty_samples, as.character(sample_id))
      next
    }

    out[[distance_col]][sample_target_idx] <- 0
    out[[zone_col]][sample_target_idx] <- inside_label

    query_idx <- setdiff(sample_idx, sample_target_idx)
    if (length(query_idx) > 0) {
      target_xy <- as.matrix(out[sample_target_idx, c(x_col, y_col)])
      query_xy <- as.matrix(out[query_idx, c(x_col, y_col)])

      if (any(!stats::complete.cases(target_xy))) {
        stop("Target region contains NA coordinates in sample: ", sample_id, call. = FALSE)
      }
      if (any(!stats::complete.cases(query_xy))) {
        stop("Non-target cells contain NA coordinates in sample: ", sample_id, call. = FALSE)
      }

      nn <- FNN::get.knnx(data = target_xy, query = query_xy, k = 1)
      out[[distance_col]][query_idx] <- as.numeric(nn$nn.dist[, 1])
      out[[zone_col]][query_idx] <- assign_distance_bands(
        distance = out[[distance_col]][query_idx],
        breaks = breaks,
        band_labels = band_labels,
        outside_label = outside_label
      )
    }
  }

  out[[zone_col]] <- factor(out[[zone_col]], levels = zone_levels)

  if (length(empty_samples) > 0) {
    warning(
      "No target region cells found in sample(s): ",
      paste(unique(empty_samples), collapse = ", "),
      call. = FALSE
    )
  }

  summary <- make_zone_summary(
    dat = out,
    sample_col = sample_col,
    celltype_col = celltype_col,
    zone_col = zone_col,
    sample_group_col = sample_group_col
  )

  list(
    annotated = out,
    zone_summary = summary$zone_summary,
    celltype_summary = summary$celltype_summary,
    params = list(
      region_values = region_values,
      region_label = region_label,
      interval = interval,
      max_distance = max_distance,
      breaks = breaks,
      inside_label = inside_label,
      outside_label = outside_label,
      no_region_label = no_region_label,
      distance_col = distance_col,
      zone_col = zone_col,
      sample_group_col = sample_group_col
    )
  )
}

make_zone_palette <- function(zone_levels) {
  zone_levels <- as.character(zone_levels)
  palette <- rep("#D9D9D9", length(zone_levels))
  names(palette) <- zone_levels

  inside_idx <- grepl("^Inside_", zone_levels)
  outside_idx <- grepl("^>", zone_levels)
  no_region_idx <- grepl("^No_", zone_levels)
  band_idx <- grepl("_to_", zone_levels)

  palette[inside_idx] <- "#E41A1C"
  palette[outside_idx] <- "#D9D9D9"
  palette[no_region_idx] <- "#969696"

  if (any(band_idx)) {
    n_band <- sum(band_idx)
    palette[band_idx] <- grDevices::colorRampPalette(
      c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1")
    )(n_band)
  }

  palette
}

plot_infiltration_map <- function(
    annotated_df,
    sample_col,
    x_col = "X",
    y_col = "Y",
    zone_col = "Distance_Group",
    sample_id = NULL,
    point_size = 0.2,
    point_alpha = 0.8,
    reverse_y = FALSE,
    show_counts = TRUE,
    facet = TRUE) {
  require_namespace("ggplot2")
  check_columns(annotated_df, c(sample_col, x_col, y_col, zone_col))

  plot_dat <- annotated_df
  if (!is.null(sample_id)) {
    plot_dat <- plot_dat[plot_dat[[sample_col]] %in% sample_id, , drop = FALSE]
  }
  if (nrow(plot_dat) == 0) {
    stop("No cells left to plot after filtering sample_id.", call. = FALSE)
  }

  if (!is.factor(plot_dat[[zone_col]])) {
    plot_dat[[zone_col]] <- factor(plot_dat[[zone_col]])
  }
  zone_levels <- levels(plot_dat[[zone_col]])
  zone_levels <- zone_levels[zone_levels %in% as.character(plot_dat[[zone_col]])]
  palette <- make_zone_palette(zone_levels)

  legend_labels <- zone_levels
  if (show_counts) {
    counts <- table(plot_dat[[zone_col]])
    legend_labels <- paste0(zone_levels, "(", as.integer(counts[zone_levels]), ")")
  }
  names(legend_labels) <- zone_levels

  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[zone_col]]
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_manual(
      values = palette,
      breaks = zone_levels,
      labels = legend_labels,
      drop = FALSE,
      na.value = "#BDBDBD"
    ) +
    ggplot2::labs(
      title = paste0("N = ", nrow(plot_dat)),
      color = zone_col,
      x = x_col,
      y = y_col
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  if (isTRUE(facet) && is.null(sample_id) && length(unique(plot_dat[[sample_col]])) > 1) {
    p <- p + ggplot2::facet_wrap(make_one_col_formula(sample_col))
  }

  p
}

plot_celltype_by_zone <- function(
    celltype_summary,
    sample_col,
    celltype_col,
    zone_col = "Distance_Group",
    fill_by = c("celltype", "zone")) {
  require_namespace("ggplot2")
  fill_by <- match.arg(fill_by)
  check_columns(celltype_summary, c(
    sample_col,
    zone_col,
    celltype_col,
    "celltype_prop_in_zone"
  ))

  fill_col <- if (fill_by == "celltype") celltype_col else zone_col

  ggplot2::ggplot(
    celltype_summary,
    ggplot2::aes(
      x = .data[[zone_col]],
      y = .data[["celltype_prop_in_zone"]],
      fill = .data[[fill_col]]
    )
  ) +
    ggplot2::geom_col(position = "stack", width = 0.85) +
    ggplot2::facet_wrap(make_one_col_formula(sample_col)) +
    ggplot2::labs(
      x = zone_col,
      y = "Cell type proportion in zone",
      fill = fill_col
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

plot_celltype_spatial_distribution <- function(
    celltype_summary,
    sample_col,
    celltype_col,
    zone_col = "Distance_Group",
    group_col = "Group",
    prop_col = "celltype_prop_in_zone",
    sample_group = NULL,
    celltypes = NULL,
    zone_levels = NULL,
    group_colors = NULL,
    title = "Cell Type Spatial Distribution",
    x_label = "Distance from TLS boundary (um)",
    y_label = "Mean Percentage (%) +/- SE",
    ncol = 4,
    point_size = 1.4,
    line_width = 0.7,
    errorbar_width = 0.18,
    free_y = TRUE) {
  require_namespace("ggplot2")
  check_columns(celltype_summary, c(sample_col, zone_col, celltype_col, prop_col))

  plot_dat <- as.data.frame(celltype_summary)

  if (!group_col %in% colnames(plot_dat)) {
    if (is.null(sample_group)) {
      stop(
        "Column '", group_col, "' is not in celltype_summary. ",
        "Please provide sample_group with columns '", sample_col, "' and '",
        group_col, "'.",
        call. = FALSE
      )
    }
    check_columns(sample_group, c(sample_col, group_col))
    sample_group <- unique(as.data.frame(sample_group)[c(sample_col, group_col)])
    plot_dat <- merge(plot_dat, sample_group, by = sample_col, all.x = TRUE)
  }

  if (any(is.na(plot_dat[[group_col]]))) {
    stop("Some samples have no group value after merging sample_group.", call. = FALSE)
  }

  if (!is.null(celltypes)) {
    plot_dat <- plot_dat[plot_dat[[celltype_col]] %in% celltypes, , drop = FALSE]
  }
  if (nrow(plot_dat) == 0) {
    stop("No rows left for plotting after filtering celltypes.", call. = FALSE)
  }

  if (is.null(zone_levels)) {
    if (is.factor(plot_dat[[zone_col]])) {
      zone_levels <- levels(plot_dat[[zone_col]])
      zone_levels <- zone_levels[zone_levels %in% as.character(plot_dat[[zone_col]])]
    } else {
      zone_levels <- unique(as.character(plot_dat[[zone_col]]))
    }
  }
  plot_dat[[zone_col]] <- factor(as.character(plot_dat[[zone_col]]), levels = zone_levels)
  plot_dat[[group_col]] <- factor(plot_dat[[group_col]])
  plot_dat[[celltype_col]] <- factor(plot_dat[[celltype_col]])

  # Convert sample-level proportions to percentages before calculating mean and SE.
  plot_dat$.__percent <- as.numeric(plot_dat[[prop_col]]) * 100

  stat_dat <- aggregate(
    plot_dat$.__percent,
    by = plot_dat[c(group_col, zone_col, celltype_col)],
    FUN = function(x) {
      x <- x[is.finite(x)]
      n <- length(x)
      mean_x <- if (n == 0) NA_real_ else mean(x)
      se_x <- if (n <= 1) 0 else stats::sd(x) / sqrt(n)
      c(mean = mean_x, se = se_x, n = n)
    }
  )

  stat_values <- if (is.list(stat_dat$x)) {
    do.call(rbind, stat_dat$x)
  } else {
    as.matrix(stat_dat$x)
  }
  stat_dat$x <- NULL
  stat_dat$mean_percent <- stat_values[, "mean"]
  stat_dat$se <- stat_values[, "se"]
  stat_dat$n_sample <- as.integer(stat_values[, "n"])
  stat_dat$ymin <- pmax(0, stat_dat$mean_percent - stat_dat$se)
  stat_dat$ymax <- stat_dat$mean_percent + stat_dat$se

  if (is.null(group_colors)) {
    default_colors <- c("#FF6F61", "#00AFC4", "#2E7D32", "#8E44AD", "#F39C12")
    group_levels <- levels(plot_dat[[group_col]])
    group_colors <- rep(default_colors, length.out = length(group_levels))
    names(group_colors) <- group_levels
  }

  facet_scales <- if (isTRUE(free_y)) "free_y" else "fixed"

  ggplot2::ggplot(
    stat_dat,
    ggplot2::aes(
      x = .data[[zone_col]],
      y = .data[["mean_percent"]],
      color = .data[[group_col]],
      group = .data[[group_col]]
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["ymin"]], ymax = .data[["ymax"]]),
      width = errorbar_width,
      linewidth = line_width
    ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(make_one_col_formula(celltype_col), scales = facet_scales, ncol = ncol) +
    ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label,
      color = group_col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0, size = 13),
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Example --------------------------------------------------------------------
#
# res <- infiltration_zone_analysis(
#   dat = df,
#   celltype_col = "cell_type",
#   sample_col = "sample",
#   group_col = "region",
#   x_col = "X",
#   y_col = "Y",
#   region_values = "TLS",
#   region_label = "TLS",
#   sample_group_col = "Group",
#   interval = 50,
#   max_distance = 500
# )
#
# annotated_df <- res$annotated
# zone_summary <- res$zone_summary
# celltype_summary <- res$celltype_summary
#
# p1 <- plot_infiltration_map(
#   annotated_df,
#   sample_col = "sample",
#   x_col = "X",
#   y_col = "Y",
#   zone_col = "Distance_Group",
#   sample_id = "Sample_1",
#   point_size = 0.25
# )
# print(p1)
#
# p2 <- plot_celltype_by_zone(
#   celltype_summary,
#   sample_col = "sample",
#   celltype_col = "cell_type",
#   zone_col = "Distance_Group"
# )
# print(p2)
#
# p3 <- plot_celltype_spatial_distribution(
#   celltype_summary = celltype_summary,
#   sample_col = "sample",
#   celltype_col = "cell_type",
#   zone_col = "Distance_Group",
#   group_col = "Group",
#   sample_group = unique(df[c("sample", "Group")]),
#   zone_levels = c(
#     "Inside_TLS", "0_to_50", "50_to_100", "100_to_200",
#     "200_to_300", "300_to_400", "400_to_500", ">500"
#   ),
#   group_colors = c(T_NR = "#FF6F61", T_R = "#00AFC4"),
#   ncol = 4
# )
# print(p3)
#
# ggsave("celltype_spatial_distribution_relative_to_TLS.pdf", p3, width = 16, height = 18)
# ggsave("celltype_spatial_distribution_relative_to_TLS.png", p3, width = 16, height = 18, dpi = 300)
#
# ggsave("Sample_1_infiltration_zone.pdf", p1, width = 8, height = 6)
# write.csv(annotated_df, "infiltration_zone_annotated.csv", row.names = FALSE)
# write.csv(zone_summary, "infiltration_zone_summary.csv", row.names = FALSE)
# write.csv(celltype_summary, "infiltration_celltype_summary.csv", row.names = FALSE)
