library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(fgsea)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(plotly)
library(jsonlite)

##################### Parameters ###############################################

db.dir <- snakemake@params[["db_dir"]]
diff.data.file <- snakemake@input[["dep_results"]]
out.dir <- snakemake@params[["out_dir"]]

##################### Helper Functions #########################################

# Save empty placeholder (png + pdf + json)
save_empty_plot <- function(png_path, msg = "No data available") {
  for (dev_fn in list(
    list(fn = png, path = png_path, args = list(width = 800, height = 400)),
    list(fn = pdf, path = sub("\\.png$", ".pdf", png_path), args = list(width = 8, height = 4))
  )) {
    do.call(dev_fn$fn, c(list(dev_fn$path), dev_fn$args))
    plot.new()
    text(0.5, 0.5, msg, cex = 2, col = "grey50")
    dev.off()
  }
  writeLines("{}", sub("\\.png$", ".json", png_path))
}

# Fisher's exact test for enrichment (one-sided, greater)
enrich.test <- function(x) {
  a <- x[1]
  b <- x[2] - a
  c <- x[3] - a
  d <- x[4] - a - b - c
  fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
}

# Perform enrichment analysis for a single regulation group
enrich.single <- function(go.data, diff.list) {
  if (is.null(diff.list) || nrow(diff.list) == 0) {
    return(data.table())
  }
  if (!"Protein accession" %in% colnames(diff.list)) {
    return(data.table())
  }

  go.data.diff <- merge(go.data, diff.list[, c("Protein accession", "Gene name")], by = "Protein accession")
  if (nrow(go.data.diff) == 0) {
    return(data.table())
  }

  # Gene/Protein aggregation per term
  gene_protein <- go.data.diff %>%
    group_by(Terms) %>%
    summarise(
      `Gene names` = paste(unique(`Gene name`), collapse = ";"),
      Proteins = paste(unique(`Protein accession`), collapse = ";"),
      .groups = "drop"
    )

  # Enrichment statistics
  mapping <- as.data.frame(table(go.data.diff$Terms), stringsAsFactors = FALSE)
  background <- as.data.frame(table(go.data$Terms), stringsAsFactors = FALSE)
  stats <- merge(mapping, background, by = "Var1")
  colnames(stats) <- c("Terms", "mapping", "background")
  stats$`All mapping` <- length(unique(go.data.diff$`Protein accession`))
  stats$`All background` <- length(unique(go.data$`Protein accession`))
  stats$`Fold enrich` <- round(stats$mapping / stats$background / stats$`All mapping` * stats$`All background`, 2)

  stats <- stats[stats$mapping > 1 & stats$`Fold enrich` > 1, ]
  if (nrow(stats) == 0) {
    return(data.table())
  }

  stats$p.value <- signif(apply(stats[, c("mapping", "background", "All mapping", "All background")], 1, enrich.test), 3)
  stats$adj.P.Val <- signif(p.adjust(stats$p.value, method = "BH"), 3)
  stats <- stats[, c("Terms", "mapping", "Fold enrich", "p.value", "adj.P.Val")]

  merge(stats, gene_protein, by = "Terms")
}

# Combine enrichment results across regulation groups (up/down)
enrich.combine <- function(go.data, diff.list) {
  if (nrow(diff.list) == 0) {
    return(data.table())
  }
  groups <- split(diff.list, diff.list$Regulation)
  results <- lapply(groups, function(x) {
    tryCatch(enrich.single(go.data, x), error = function(e) {
      message("enrich.single error: ", e$message)
      data.table()
    })
  })
  results <- Filter(function(x) !is.null(x) && nrow(x) > 0, results)
  if (length(results) == 0) {
    return(data.table())
  }
  rbindlist(results, idcol = "Regulation", fill = TRUE)
}

# Calculate running enrichment score for GSEA plot
calc_ES <- function(geneList, pathway, p = 1) {
  hits <- names(geneList) %in% pathway
  Nh <- sum(hits)
  Nm <- length(geneList) - Nh
  Phit <- ifelse(hits, abs(geneList)^p / sum(abs(geneList[hits])^p), 0)
  Pmiss <- ifelse(!hits, 1 / Nm, 0)
  RES <- cumsum(Phit - Pmiss)
  ES <- ifelse(max(RES) > abs(min(RES)), max(RES), min(RES))
  list(ES = ES, DF = data.frame(rank = seq_along(geneList), score = geneList, hit = as.integer(hits), runningES = RES))
}

# Generate GSEA enrichment score plot
gsea.plot <- function(pathway_genes, geneList, gsea.res) {
  gsea.res <- gsea.res[1, ]
  pathway_name <- gsea.res$pathway
  es_result <- calc_ES(geneList, pathway_genes[[pathway_name]])
  max_es <- es_result$ES
  df <- es_result$DF
  yn <- round(max_es / 10, 2)

  # Annotation table position
  if (max_es > 0) {
    xmin <- 0.6 * max(df$rank)
    xmax <- 0.9 * max(df$rank)
  } else {
    xmin <- 0.1 * max(df$rank)
    xmax <- 0.4 * max(df$rank)
  }

  ggplot(df, aes(x = rank, y = runningES)) +
    geom_tile(aes(y = 0, fill = score), height = 2 * yn) +
    scale_fill_gradientn(colors = c("steelblue", "white", "red"), limits = c(-1, 1), oob = scales::squish) +
    geom_line(size = 1.1, color = "#1f78b4") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(data = df[df$hit == 1, ], aes(x = rank, xend = rank, y = -yn, yend = yn), color = "black", alpha = 0.6) +
    geom_point(data = df[which.max(abs(df$runningES)), ], aes(x = rank, y = max_es), size = 4, color = "firebrick") +
    geom_text_repel(
      data = df[which.max(abs(df$runningES)), ],
      aes(x = rank, y = max_es, label = paste0("ES: ", round(max_es, 3))),
      nudge_y = yn, size = 5, color = "firebrick", segment.color = "firebrick", segment.size = 0.5, show.legend = FALSE
    ) +
    annotation_custom(
      grob = tableGrob(
        gsea.res[, c("pval", "NES", "size")] %>% mutate(pval = signif(pval, 2), NES = round(NES, 3)),
        rows = NULL
      ),
      xmin = xmin, xmax = xmax, ymin = 0.5 * max_es, ymax = max_es
    ) +
    labs(x = "Rank in Ordered Gene List", y = "Running Enrichment Score", title = pathway_name) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(), axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), plot.title = element_text(size = 18, face = "bold"),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
}

# Build interactive JSON specification for bubble plot
create.json <- function(plot_data) {
  plot_data$Terms <- reorder(plot_data$Terms, plot_data$`Fold enrich`)
  plot_data <- plot_data[order(plot_data$Terms), ]

  axis_label_style <- list(
    show = TRUE, color = "#000", fontStyle = "normal",
    fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
  )
  empty_name_style <- list(color = "", fontStyle = "normal", fontWeight = "normal", fontFamily = "", fontSize = "")

  list(
    width = 1000, height = 800, outerBorder = TRUE,
    margin = list(top = 40, right = 60, bottom = 60, left = 60),
    grid = list(list(left = 0, width = 280), list(left = 350, width = 280)),
    dataset = list(list(source = list(
      terms = plot_data$Terms, log2FoldEnrich = log2(plot_data$`Fold enrich`),
      pvalue = plot_data$p.value, mapping = plot_data$mapping,
      regulation = plot_data$Regulation, sig = rep("sig", nrow(plot_data))
    ))),
    series = list(list(
      type = "facetedBubble", datasetIndex = 0,
      encode = list(
        x = "log2FoldEnrich", y = "terms", size = "mapping", color = "pvalue",
        facet = "regulation", tooltip = c("terms", "log2FoldEnrich", "pvalue", "mapping", "sig", "regulation"),
        annotation = "sig"
      ),
      border = list(show = TRUE)
    )),
    xAxis = list(list(
      type = "scaleLinear", name = "log2(`Fold enrich`)",
      nameTextStyle = empty_name_style, axisLabel = axis_label_style
    )),
    yAxis = list(list(
      type = "scaleBand", name = "Terms",
      nameTextStyle = empty_name_style, axisLabel = axis_label_style
    )),
    legend = list(show = TRUE, site = "right", orient = "vertical", top = 30, right = 0, type = "circle", value = list()),
    sizeLegend = list(
      show = TRUE, title = "mapping", orient = "vertical", count = 4, right = 10, top = 390,
      bordered = TRUE, color = "#000000", fontSize = 12, site = "right"
    ),
    colorbar = list(
      show = TRUE, title = "p.value", width = 20, height = 150, site = "right", right = 10, top = 220,
      fontSize = 12, tickCount = 5, useScientificNotation = TRUE
    ),
    axisTick = list(show = TRUE, labelType = "number", count = 3)
  )
}

# Merge enrichment results with annotation metadata, with error handling
safe_merge <- function(list_data, enrich_data, category_name) {
  tryCatch(
    {
      if (is.null(enrich_data) || nrow(enrich_data) == 0) {
        return(data.table())
      }
      result <- merge(list_data, enrich_data, by = "Terms")
      if (nrow(result) > 0) result$Category <- category_name
      result
    },
    error = function(e) {
      message("Enrichment merge skipped for ", category_name, ": ", e$message)
      data.table()
    }
  )
}

##################### Enrichment Pipeline ######################################

enrich.pipeline <- function(data.diff, out.dir) {
  # Run enrichment per annotation category
  enrich.out <- lapply(split(gmt.data, gmt.data$Category), function(x) enrich.combine(x, data.diff))

  # Merge with annotation metadata (apply GO level / term length filters)
  enrich.out.list <- list(
    GOBP = safe_merge(go.list %>% filter(nchar(as.character(Terms)) < 60, `GO levels` > 5, `GO levels` < 10), enrich.out$`Biological process`, "GOBP"),
    GOCC = safe_merge(go.list %>% filter(nchar(as.character(Terms)) < 60, `GO levels` > 3), enrich.out$`Cellular component`, "GOCC"),
    GOMF = safe_merge(go.list %>% filter(nchar(as.character(Terms)) < 60, `GO levels` > 3), enrich.out$`Molecular function`, "GOMF"),
    KEGG = safe_merge(kegg.list %>% filter(nchar(as.character(Terms)) < 60), enrich.out$KEGG, "KEGG"),
    PFAM = safe_merge(pfam.list %>% filter(nchar(as.character(Terms)) < 60), enrich.out$Pfam, "PFAM")
  )

  # Collect full (unfiltered) results for the combined table
  all_results <- list()

  # Process each category: generate plots + per-category tables
  for (cat_name in names(enrich.out.list)) {
    x <- enrich.out.list[[cat_name]]
    category_name <- if (nrow(x) > 0) x$Category[1] else cat_name

    if (nrow(x) == 0) {
      save_empty_plot(file.path(out.dir, paste0("Figure 4. ", category_name, " enrichment of ", compare.name, ".png")),
        msg = paste0("No ", category_name, " enrichment results")
      )
      fwrite(data.table(), file = file.path(out.dir, paste0("Table 4. ", category_name, " enrichment of ", compare.name, ".xls")), sep = "\t")
      next
    }

    # Save per-category table (all results, unfiltered)
    fwrite(x, file = file.path(out.dir, paste0("Table 4. ", category_name, " enrichment of ", compare.name, ".xls")), sep = "\t")

    # Collect for combined table (pre-filter: all results)
    all_results[[cat_name]] <- x

    # Filter for plot & JSON: preserve the legacy per-regulation dense-rank selection,
    # then strictly cap the combined result at 20 terms.
    x_sig <- x[p.value < 0.05]
    if (nrow(x_sig) == 0) {
      save_empty_plot(file.path(out.dir, paste0("Figure 4. ", category_name, " enrichment of ", compare.name, ".png")),
        msg = paste0("No significant ", category_name, " enrichment results (p < 0.05)")
      )
      next
    }
    N <- length(unique(x_sig$Regulation))
    df.plot <- x_sig[, .SD[frankv(p.value, ties.method = "dense") <= 20 / N], by = "Regulation"]
    if (nrow(df.plot) > 20) {
      plot_quota <- floor(20 / N)
      setorder(df.plot, Regulation, p.value, Terms)
      df.plot[, plot_rank := seq_len(.N), by = "Regulation"]
      df.plot.keep <- df.plot[plot_rank <= plot_quota]
      plot_slots <- 20 - nrow(df.plot.keep)
      if (plot_slots > 0) {
        df.plot.extra <- df.plot[plot_rank > plot_quota][order(p.value, Terms)]
        df.plot.keep <- rbind(df.plot.keep, head(df.plot.extra, plot_slots), use.names = TRUE)
      }
      df.plot <- df.plot.keep[, plot_rank := NULL]
    }

    if (nrow(df.plot) > 0) {
      # Save interactive JSON (same filtered data as plot)
      tryCatch(
        {
          write_json(create.json(df.plot), file.path(out.dir, paste0("Figure 4. ", category_name, " enrichment of ", compare.name, ".json")), auto_unbox = TRUE, pretty = TRUE)
        },
        error = function(e) message("JSON export skipped for ", category_name, ": ", e$message)
      )

      maxn <- max(nchar(as.character(df.plot$Terms)), na.rm = TRUE)
      if (is.infinite(maxn)) maxn <- 20
      N <- length(unique(df.plot$Regulation))
      height_num <- max(4, 2 + length(unique(df.plot$Terms)) * 0.2)
      width_num <- max(8, (0.08 - 0.02 * N) * maxn + 1 + N * 3)

      p <- ggplot(df.plot, aes(y = reorder(Terms, `Fold enrich`), x = log2(`Fold enrich`), fill = p.value)) +
        geom_point(aes(size = mapping), shape = 21, color = "grey30", stroke = 0.5) +
        theme_test(base_size = 14) +
        ylab("Terms") +
        facet_grid(. ~ Regulation) +
        scale_fill_material(reverse = TRUE, palette = "red") +
        scale_size(range = c(3, 6)) +
        guides(size = guide_legend(override.aes = list(shape = 16, color = "black", fill = "black"))) +
        scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
        scale_y_discrete(expand = expansion(add = 0.6)) +
        theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

      ggsave(plot = p, filename = file.path(out.dir, paste0("Figure 4. ", category_name, " enrichment of ", compare.name, ".pdf")), width = width_num, height = height_num)
      ggsave(plot = p, filename = file.path(out.dir, paste0("Figure 4. ", category_name, " enrichment of ", compare.name, ".png")), width = width_num, height = height_num, dpi = 300)
    }
  }

  # Return combined unfiltered results (Category, Terms, Regulation, p.value only)
  if (length(all_results) > 0) {
    enrich.result <- rbindlist(all_results, fill = TRUE)[, c("Category", "Terms", "Regulation", "p.value")]
    setorder(enrich.result, p.value)
    unique(enrich.result, by = c("Category", "Terms", "Regulation"))
  } else {
    data.table()
  }
}

##################### Main Program #############################################

compare.name <- sub("\\.xls$", "", sub(".*Table 2\\. DEPs list of ", "", diff.data.file))

# Load data
data.diff.all <- fread(diff.data.file)
go.data <- fread(file.path(db.dir, "GO_annotation.tsv"))
nog.data <- fread(file.path(db.dir, "NOG_annotation.tsv"))
pfam.data <- fread(file.path(db.dir, "Pfam_annotation.tsv"))
kegg.data <- fread(file.path(db.dir, "KEGG_annotation.tsv"))
annot.data <- fread(file.path(db.dir, "annotation.tsv"))

# Build annotation lookup tables
go.list <- distinct(go.data, ID, Category, Terms, `GO levels`)
kegg.list <- distinct(kegg.data, `Pathway ID`, `Pathway level 1`, `Pathway level 2`, Terms) %>% filter(`Pathway level 1` != "Human Diseases")
pfam.list <- distinct(pfam.data, Terms)

# Combine annotation data for enrichment
pfam.data$Category <- "Pfam"
kegg.data$Category <- "KEGG"
gmt.data <- rbindlist(list(
  go.data[, c("Protein accession", "Category", "Terms")],
  pfam.data[, c("Protein accession", "Category", "Terms")],
  kegg.data[, c("Protein accession", "Category", "Terms")] %>% filter(Terms %in% kegg.list$Terms)
)) %>%
  distinct_all() %>%
  filter(`Protein accession` %in% data.diff.all$`Protein accession`)

# --- Enrichment Analysis ---
if (!dir.exists(out.dir)) dir.create(out.dir)

# Up/Down-regulated separately
enrich.out.combine <- enrich.pipeline(data.diff.all %>% filter(Regulation != "unchanged"), out.dir)
fwrite(enrich.out.combine, file = file.path(out.dir, paste0("Table 4. Functional enrichment of ", compare.name, ".xls")), sep = "\t")

# All DEPs combined
out.dir.all <- file.path(out.dir, "ALL_DEPs_enrichment")
if (!dir.exists(out.dir.all)) dir.create(out.dir.all)
enrich.out.combine <- enrich.pipeline(data.diff.all %>% filter(Regulation != "unchanged") %>% mutate(Regulation = "All DEPs"), out.dir.all)
fwrite(enrich.out.combine, file = file.path(out.dir.all, paste0("Table 4. Functional enrichment of ", compare.name, ".xls")), sep = "\t")

# --- GSEA Analysis ---
gmt.list <- lapply(split(gmt.data, gmt.data$Category), function(x) split(x$`Protein accession`, x$Terms))
geneList <- sort(setNames(data.diff.all$log2FC, data.diff.all$`Protein accession`) %>% na.omit(), decreasing = TRUE)

set.seed(42)
kegg.gsea <- fgsea(pathways = gmt.list$KEGG, stats = geneList, minSize = 10, maxSize = 500)
kegg.gsea$Genes <- lapply(kegg.gsea$leadingEdge, function(x) {
  unique(data.diff.all$`Gene name`[data.diff.all$`Protein accession` %in% x])
})

out.dir.gsea <- file.path(out.dir, "GSEA")
if (!dir.exists(out.dir.gsea)) dir.create(out.dir.gsea)

# Pathway ID mapping for safe filenames
pathway_id_map <- setNames(kegg.list$`Pathway ID`, kegg.list$Terms)

fwrite(kegg.gsea, file = file.path(out.dir.gsea, paste0("Table 5. KEGG GSEA of ", compare.name, ".xls")), sep = "\t")

kegg.gsea.sig <- kegg.gsea %>%
  filter(pval < 0.05) %>%
  slice_min(pval, n = 20)

if (nrow(kegg.gsea.sig) > 0) {
  # GSEA bar plot
  maxn <- max(nchar(as.character(kegg.gsea.sig$pathway)), na.rm = TRUE)
  if (is.infinite(maxn)) maxn <- 20

  p <- ggplot(kegg.gsea.sig, aes(y = reorder(pathway, NES), x = NES, fill = -log10(pval))) +
    geom_bar(stat = "identity", width = 0.8, color = "black") +
    theme_classic(base_size = 16) +
    geom_vline(xintercept = 0) +
    scale_fill_material("red") +
    labs(y = "KEGG pathway", fill = "-log10(p-value)")

  ggsave(
    plot = p, filename = file.path(out.dir.gsea, paste0("Figure 5. KEGG GSEA of ", compare.name, ".pdf")),
    width = maxn * 0.06 + 8, height = 2 + nrow(kegg.gsea.sig) * 0.3
  )
  ggsave(
    plot = p, filename = file.path(out.dir.gsea, paste0("Figure 5. KEGG GSEA of ", compare.name, ".png")),
    width = maxn * 0.06 + 8, height = 2 + nrow(kegg.gsea.sig) * 0.3
  )

  # Individual GSEA enrichment plots
  out.dir.gsea.plots <- file.path(out.dir.gsea, paste("Enrich plot of", compare.name))
  if (!dir.exists(out.dir.gsea.plots)) dir.create(out.dir.gsea.plots)

  lapply(split(kegg.gsea.sig, kegg.gsea.sig$pathway), function(x) {
    tryCatch(
      {
        pid <- pathway_id_map[x$pathway[1]]
        if (is.na(pid)) pid <- gsub("[^A-Za-z0-9_-]", "_", x$pathway[1])
        file_label <- paste0(pid, "_", compare.name)
        p <- gsea.plot(pathway_genes = gmt.list$KEGG, geneList = geneList, gsea.res = x)
        ggsave(plot = p, filename = file.path(out.dir.gsea.plots, paste0("Figure 6. KEGG GSEA plot of ", file_label, ".pdf")), width = 8, height = 5)
        ggsave(plot = p, filename = file.path(out.dir.gsea.plots, paste0("Figure 6. KEGG GSEA plot of ", file_label, ".png")), width = 8, height = 5)
      },
      error = function(e) message("GSEA plot skipped for ", x$pathway[1], ": ", e$message)
    )
  })
} else {
  message("No significant KEGG GSEA pathways found (pval < 0.05). Generating placeholder.")
  save_empty_plot(file.path(out.dir.gsea, paste0("Figure 5. KEGG GSEA of ", compare.name, ".png")),
    msg = "No significant KEGG GSEA pathways (p < 0.05)"
  )
}

# --- Finalize ---
write.csv("", file.path(out.dir, "done"))
