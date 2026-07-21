if (length(snakemake@log) > 0) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log)
  sink(log, type = "message")
}

project_dir <- getwd()
to_abs <- function(path) {
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(project_dir, path)
}

options(stringsAsFactors = FALSE)
library(tidyverse)
library(igraph)
library(data.table)

contrast <- gsub("_diffexp.DEG.xls", "", basename(snakemake@input[[1]]))
top6_clusters_plot_png <- to_abs(snakemake@output[[1]])
top6_clusters_plot_pdf <- to_abs(snakemake@output[[2]])
string_link <- to_abs(snakemake@output[[3]])
work_dir <- dirname(top6_clusters_plot_png)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

stringdb.version <- sub("^v", "", as.character(snakemake@params[["STRINGdb_version"]]))
stringdb.dir <- to_abs(snakemake@params[["STRINGdb_dir"]])
version.dir <- file.path(stringdb.dir, paste0("v", stringdb.version))
if (dir.exists(version.dir)) {
  stringdb.dir <- version.dir
}
score.threshold <- as.numeric(snakemake@params[["STRINGdb_score_threshold"]])
if (is.na(score.threshold)) {
  score.threshold <- 700
}

write_empty_ppi <- function(text = "可用于PPI分析的差异蛋白少于2个，无法构建互作网络。") {
  library(ggplot2)
  library(showtext)
  showtext_auto()

  dir.create(dirname(string_link), recursive = TRUE, showWarnings = FALSE)
  writeLines("node1,node2,combined_score", string_link)

  p <- ggplot() +
    annotate("text", x = 4, y = 25, size = 6, label = text) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white"))
  ggsave(
    p,
    filename = top6_clusters_plot_pdf,
    device = "pdf",
    width = 10,
    height = 10
  )
  ggsave(
    p,
    filename = top6_clusters_plot_png,
    device = "png",
    width = 4,
    height = 4
  )
}

stringdb_file <- function(species.id, file.kind, ext) {
  file.path(
    stringdb.dir,
    paste0(species.id, ".protein.", file.kind, ".v", stringdb.version, ext)
  )
}

existing_stringdb_file <- function(species.id, file.kind, exts) {
  for (ext in exts) {
    path <- stringdb_file(species.id, file.kind, ext)
    if (file.exists(path) && file.info(path)$size > 0) {
      return(path)
    }
  }
  ""
}

infer_species_id <- function(string_ids) {
  species.count <- data.frame(
    species.id = sub("\\..*", "", string_ids),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(species.id), species.id != "") %>%
    count(species.id, name = "Freq") %>%
    arrange(desc(Freq))

  if (nrow(species.count) == 0) {
    return("")
  }

  local.species <- species.count$species.id[
    vapply(
      species.count$species.id,
      function(id) existing_stringdb_file(id, "links", c(".txt.gz", ".txt")) != "",
      logical(1)
    )
  ]
  if (length(local.species) > 0) {
    return(local.species[[1]])
  }
  species.count$species.id[[1]]
}

read_string_links <- function(link.file, hits) {
  hits.file <- tempfile()
  writeLines(hits, hits.file)
  on.exit(unlink(hits.file), add = TRUE)

  reader <- if (grepl("\\.gz$", link.file)) "gzip -dc" else "cat"
  awk.script <- sprintf(
    "BEGIN { while ((getline line < \"%s\") > 0) hit[line]=1 } NR==1 || (($1 in hit) && ($2 in hit))",
    hits.file
  )
  cmd <- paste(reader, shQuote(link.file), "| awk", shQuote(awk.script))
  fread(cmd = cmd)
}

plot_ppi_cluster_local <- function(edges_df, cluster_nodes, dat_map, cluster_id = 1) {
  sub_edges <- edges_df %>% filter(from %in% cluster_nodes & to %in% cluster_nodes)
  if (nrow(sub_edges) == 0) {
    plot.new()
    text(0.5, 0.5, paste0("Cluster ", cluster_id, "\nNo interactions"), cex = 1.2)
    return(invisible(NULL))
  }
  g <- graph_from_data_frame(sub_edges[, c("from", "to")], directed = FALSE)
  node_ids <- V(g)$name
  node_labels <- dat_map$gene_id[match(node_ids, dat_map$STRING_id)]
  node_labels[is.na(node_labels)] <- node_ids[is.na(node_labels)]
  edge_weights <- sub_edges$combined_score / max(sub_edges$combined_score, na.rm = TRUE)
  node_fc <- dat_map$log2FoldChange[match(node_ids, dat_map$STRING_id)]
  node_fc[is.na(node_fc)] <- 0
  max_abs_fc <- max(abs(node_fc), na.rm = TRUE)
  if (!is.finite(max_abs_fc) || max_abs_fc == 0) {
    max_abs_fc <- 1
  }
  color_scale_up <- colorRampPalette(c("#FFCCCC", "#CC0000"))(100)
  color_scale_down <- colorRampPalette(c("#CCCCFF", "#0000CC"))(100)
  color_index <- pmin(ceiling(abs(node_fc) / max_abs_fc * 99) + 1, 100)
  node_colors <- ifelse(node_fc > 0, color_scale_up[color_index], color_scale_down[color_index])
  n_nodes <- length(V(g))
  node_size <- ifelse(n_nodes > 50, 4, ifelse(n_nodes > 20, 6, 8))
  label_cex <- ifelse(n_nodes > 50, 0.4, ifelse(n_nodes > 20, 0.55, 0.7))
  layout <- layout_with_kk(g)
  plot(
    g,
    layout = layout,
    vertex.label = node_labels,
    vertex.label.cex = label_cex,
    vertex.label.color = "black",
    vertex.label.family = "sans",
    vertex.size = node_size,
    vertex.color = node_colors,
    vertex.frame.color = "grey60",
    edge.width = edge_weights * 0.5 + 0.5,
    edge.color = adjustcolor("grey50", alpha.f = 0.6),
    edge.curved = 0.1,
    main = paste0("Cluster ", cluster_id, " (", n_nodes, " nodes)")
  )
}

DEG <- read.table(to_abs(snakemake@input[[1]]), sep = "\t", header = TRUE, check.names = FALSE)
string_mapping <- read.table(
  to_abs(snakemake@input[["string_mapping"]]),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

DEG$gene_id <- DEG$gene
ppi_df <- DEG[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
ppi_df <- na.omit(ppi_df)
ppi_df <- merge(ppi_df, string_mapping, by = "gene_id")
ppi_df <- dplyr::distinct(ppi_df, string_id, .keep_all = TRUE)

if (nrow(ppi_df) < 2) {
  write_empty_ppi()
  quit(save = "no")
}

dat_map <- ppi_df %>%
  mutate(STRING_id = string_id)
hits <- unique(dat_map$STRING_id)
hits_2000 <- hits
if (length(hits_2000) > 2000) {
  hits_2000 <- hits_2000[1:2000]
}
spe <- infer_species_id(hits)
message("STRING species ID: ", spe)
message("STRINGdb directory: ", stringdb.dir)
message("PPI cluster input STRING IDs: ", length(hits_2000), " / ", length(hits))

link.file <- existing_stringdb_file(spe, "links", c(".txt.gz", ".txt"))
if (link.file == "") {
  write_empty_ppi(
    paste0("未找到该物种的 STRING v", stringdb.version, " 本地互作数据库，无法生成PPI网络。")
  )
  quit(save = "no")
}

tryCatch({
  links <- read_string_links(link.file, hits)
  required.columns <- c("protein1", "protein2", "combined_score")
  missing.columns <- setdiff(required.columns, colnames(links))
  if (length(missing.columns) > 0) {
    stop("STRING links file is missing columns: ", paste(missing.columns, collapse = ", "))
  }

  all_interactions <- links %>%
    rename(from = protein1, to = protein2) %>%
    mutate(combined_score = as.numeric(combined_score)) %>%
    filter(combined_score >= score.threshold) %>%
    mutate(
      pair_a = pmin(from, to),
      pair_b = pmax(from, to)
    ) %>%
    group_by(pair_a, pair_b) %>%
    slice_max(combined_score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(from = pair_a, to = pair_b, combined_score = combined_score)

  if (nrow(all_interactions) == 0) {
    write_empty_ppi(
      paste0("在 STRING v", stringdb.version, " 中未检索到差异蛋白之间的高置信互作。")
    )
    quit(save = "no")
  }

  cluster_interactions <- all_interactions %>%
    filter(from %in% hits_2000, to %in% hits_2000)

  if (nrow(cluster_interactions) == 0) {
    cluster_interactions <- all_interactions
  }

  g_cluster <- graph_from_data_frame(cluster_interactions[, c("from", "to")], directed = FALSE)
  g_cluster <- simplify(g_cluster)
  ppi.clusters <- cluster_fast_greedy(g_cluster)
  clustersList <- split(V(g_cluster)$name, membership(ppi.clusters))
  clustersList <- clustersList[order(lengths(clustersList), decreasing = TRUE)]
  n_cluster <- min(6, length(clustersList))

  png(
    filename = top6_clusters_plot_png,
    units = "in",
    width = 15,
    height = 10,
    res = 400
  )
  par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
  for (i in 1:n_cluster) {
    plot_ppi_cluster_local(cluster_interactions, clustersList[[i]], dat_map, cluster_id = i)
  }
  if (n_cluster < 6) {
    for (i in (n_cluster + 1):6) {
      plot.new()
    }
  }
  dev.off()

  pdf(file = top6_clusters_plot_pdf, width = 15, height = 10)
  par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
  for (i in 1:n_cluster) {
    plot_ppi_cluster_local(cluster_interactions, clustersList[[i]], dat_map, cluster_id = i)
  }
  if (n_cluster < 6) {
    for (i in (n_cluster + 1):6) {
      plot.new()
    }
  }
  dev.off()

  dat_link <- all_interactions %>%
    transmute(node1 = from, node2 = to, combined_score = combined_score)
  write.csv(dat_link, string_link, row.names = FALSE, quote = FALSE)
}, error = function(e) {
  message("STRINGdb PPI failed; writing placeholder outputs: ", conditionMessage(e))
  write_empty_ppi(
    paste0("STRING v", stringdb.version, " 本地互作数据库读取或解析失败，已生成占位结果。")
  )
})
