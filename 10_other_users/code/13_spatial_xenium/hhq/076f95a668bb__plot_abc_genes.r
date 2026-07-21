library(Seurat)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(rstatix)
library(paletteer)

abc_genes <- c(
  "ABCA1", "ABCA3", "ABCA4", "ABCA7", "ABCB1", "ABCB4", "ABCB6",
  "ABCC1", "ABCC12", "ABCC2", "ABCC3", "ABCC4", "ABCC6", "ABCC8",
  "ABCC9", "ABCD1", "ABCD3", "ABCD4", "ABCF3", "ABCG1", "ABCG2"
)

meta_vars <- c("group", "singler", "sample", "singler_subtype")

plot_data <- FetchData(
  adata,
  vars = c(abc_genes, meta_vars)
)

head(plot_data)

safe_name <- function(x) {
  gsub("[\\\\/:*?\"<>|]", "_", x)
}

get_facet_layout <- function(n, max_ncol = 6) {
  if (n <= 1) {
    return(list(ncol = 1, nrow = 1))
  }

  ncol <- min(max_ncol, ceiling(sqrt(n)))
  nrow <- ceiling(n / ncol)

  list(ncol = ncol, nrow = nrow)
}

comparisons <- list(
  c("PreR", "PostR"),
  c("PreR", "PreNR"),
  c("PreNR", "PostNR"),
  c("PostR", "PostNR")
)

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(x)
}

get_signif_stats <- function(data_bar, comparisons) {
  data_for_test <- data_bar %>%
    filter(!is.na(CN), !is.na(group), is.finite(gene_expr))

  if (nrow(data_for_test) == 0) {
    return(data.frame())
  }

  stat_list <- lapply(split(data_for_test, data_for_test$CN), function(df) {
    df$group <- droplevels(df$group)

    valid_comparisons <- Filter(function(comp) {
      all(vapply(comp, function(group_name) {
        sum(df$group == group_name, na.rm = TRUE) > 0
      }, logical(1)))
    }, comparisons)

    if (length(valid_comparisons) == 0) {
      return(NULL)
    }

    stat_test <- tryCatch(
      {
        df %>%
          wilcox_test(
            gene_expr ~ group,
            comparisons = valid_comparisons
          ) %>%
          filter(!is.na(p.adj.signif), p.adj.signif != "ns")
      },
      error = function(e) {
        NULL
      }
    )

    if (is.null(stat_test) || nrow(stat_test) == 0) {
      return(NULL)
    }

    y_range <- range(df$gene_expr, na.rm = TRUE)
    y_step <- diff(y_range)
    if (!is.finite(y_step) || y_step <= 0) {
      y_step <- max(abs(y_range[2]) * 0.1, 0.1)
    }

    stat_test$CN <- df$CN[1]
    stat_test$y.position <- y_range[2] + seq_len(nrow(stat_test)) * y_step * 0.15
    stat_test
  })

  bind_rows(stat_list)
}

major_celltypes <- c(
  "Epithelial cells",
  "Fibroblasts",
  "CD4T",
  "CD8T",
  "Macrophages",
  "B cells"
)

celltype_height <- 15
celltype_width <- 30

for (gene in abc_genes) {

  gene_dir <- safe_name(gene)
  dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)

  # =========================
  # 1. Celltype boxplot
  # =========================

  data_bar <- plot_data %>%
    group_by(group, sample, singler) %>%
    summarise(
      gene_expr = mean_or_na(.data[[gene]]),
      .groups = "drop"
    )

  data_bar$CN <- as.factor(data_bar$singler)
  data_bar$group <- factor(
    data_bar$group,
    levels = c("PreR", "PostR", "PreNR", "PostNR")
  )

  g <- ggplot(data_bar, aes(group, gene_expr)) +
    geom_boxplot(aes(fill = CN), outlier.shape = NA, alpha = 0.8) +
    geom_jitter(size = 0.5, color = "gray") +
    xlab("") +
    ylab("Expression") +
    ggtitle(gene) +
    ggforce::facet_wrap_paginate(
      ~ CN,
      ncol = 6,
      nrow = 3,
      scales = "free"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 15),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(color = "black", size = 14),
      plot.title = element_text(face = "bold", size = 21),
      strip.background = element_blank(),
      strip.text.x = element_text(hjust = 0)
    ) +
    paletteer::scale_fill_paletteer_d(palette = "ggsci::category20_d3") +
    NoLegend()

  stat.test <- get_signif_stats(data_bar, comparisons)

  if (nrow(stat.test) > 0) {
    g <- g +
      stat_pvalue_manual(
        stat.test,
        label = "p.adj.signif",
        tip.length = 0,
        hide.ns = TRUE
      )
  }

  options(
    repr.plot.height = celltype_height,
    repr.plot.width = celltype_width,
    repr.plot.res = 200
  )

  print(g)

  ggsave(
    file.path(gene_dir, paste0(gene, "_celltype_boxplot.pdf")),
    g,
    dpi = 300,
    height = celltype_height,
    width = celltype_width
  )

  ggsave(
    file.path(gene_dir, paste0(gene, "_celltype_boxplot.png")),
    g,
    dpi = 300,
    height = celltype_height,
    width = celltype_width
  )

  # =========================
  # 2. Subtype boxplot
  #    Save by major celltype folder
  #    Dynamic ncol / nrow
  # =========================

  for (ct in major_celltypes) {

    ct_dir <- file.path(gene_dir, safe_name(ct))
    dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

    data_bar <- plot_data %>%
      filter(singler == ct) %>%
      group_by(group, sample, singler_subtype) %>%
      summarise(
        gene_expr = mean_or_na(.data[[gene]]),
        .groups = "drop"
      )

    if (nrow(data_bar) == 0) {
      next
    }

    data_bar$CN <- as.factor(data_bar$singler_subtype)
    data_bar$group <- factor(
      data_bar$group,
      levels = c("PreR", "PostR", "PreNR", "PostNR")
    )

    n_subtypes <- dplyr::n_distinct(data_bar$CN)
    facet_layout <- get_facet_layout(n_subtypes, max_ncol = 6)

    subtype_height <- max(5, facet_layout$nrow * 4)
    subtype_width <- max(8, facet_layout$ncol * 5)

    g <- ggplot(data_bar, aes(group, gene_expr)) +
      geom_boxplot(aes(fill = CN), outlier.shape = NA, alpha = 0.8) +
      geom_jitter(size = 0.5, color = "gray") +
      xlab("") +
      ylab("Expression") +
      ggtitle(paste0(gene, " - ", ct)) +
      ggforce::facet_wrap_paginate(
        ~ CN,
        ncol = facet_layout$ncol,
        nrow = facet_layout$nrow,
        scales = "free"
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 15),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(color = "black", size = 14),
        plot.title = element_text(face = "bold", size = 21),
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0)
      ) +
      paletteer::scale_fill_paletteer_d(palette = "ggsci::category20_d3") +
      NoLegend()

    stat.test <- get_signif_stats(data_bar, comparisons)

    if (nrow(stat.test) > 0) {
      g <- g +
        stat_pvalue_manual(
          stat.test,
          label = "p.adj.signif",
          tip.length = 0,
          hide.ns = TRUE
        )
    }

    options(
      repr.plot.height = subtype_height,
      repr.plot.width = subtype_width,
      repr.plot.res = 200
    )

    print(g)

    ggsave(
      file.path(
        ct_dir,
        paste0(gene, "_", safe_name(ct), "_subtype_boxplot.pdf")
      ),
      g,
      dpi = 300,
      height = subtype_height,
      width = subtype_width
    )

    ggsave(
      file.path(
        ct_dir,
        paste0(gene, "_", safe_name(ct), "_subtype_boxplot.png")
      ),
      g,
      dpi = 300,
      height = subtype_height,
      width = subtype_width
    )
  }
}
