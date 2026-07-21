log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

work_dir <- getwd()

options(timeout = 1200)
options(stringsAsFactors = F)
library(tidyverse)
library(igraph)

# Input
## DEG results
DEG <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE)
## Species
spe <- 10116
## Local STRINGdb
STRINGdb_version <- snakemake@params[["STRINGdb_version"]]
STRINGdb_dir <- snakemake@params[["STRINGdb_dir"]]
STRINGdb_file_version <- STRINGdb_version
if (!grepl(".", STRINGdb_file_version, fixed=TRUE)) {
  STRINGdb_file_version <- paste0(STRINGdb_file_version, ".0")
}
if (STRINGdb_file_version == "11.0b") {
  STRINGdb_file_version <- "11.0"
}

protein_info_file <- file.path(STRINGdb_dir, paste0(spe, ".protein.info.v", STRINGdb_file_version, ".txt.gz"))
protein_aliases_file <- file.path(STRINGdb_dir, paste0(spe, ".protein.aliases.v", STRINGdb_file_version, ".txt.gz"))
protein_links_file <- file.path(STRINGdb_dir, paste0(spe, ".protein.links.v", STRINGdb_file_version, ".txt.gz"))
required_stringdb_files <- c(protein_info_file, protein_aliases_file, protein_links_file)
missing_stringdb_files <- required_stringdb_files[!file.exists(required_stringdb_files) | file.info(required_stringdb_files)$size == 0]
if (length(missing_stringdb_files) > 0) {
  stop(
    "Local STRINGdb files are missing or empty. PPI is running in offline mode and will not download files.\nMissing files:\n",
    paste(missing_stringdb_files, collapse = "\n")
  )
}

# Output
top6_clusters_plot_png <- snakemake@output[[1]]
top6_clusters_plot_pdf <- snakemake@output[[2]]
string_link <- snakemake@output[[3]]

# Extract DEG
ppi_df <- DEG[, c("log2FoldChange","pvalue","padj","gene_name")]
ppi_df <- na.omit(ppi_df)
ppi_df <- dplyr::distinct(ppi_df, gene_name, .keep_all=TRUE)
colnames(ppi_df) <- c('log2FC', 'pvalue', 'padj', 'gene')
gene_diff <- ppi_df

gene_diff <- gene_diff[order(gene_diff$pvalue),]
if (dim(gene_diff)[1] > 2000) {
  gene_diff <- head(gene_diff,2000)
}

# Load STRINGdb files locally and avoid online version checks.
load_stringdb_proteins <- function(path, file_version) {
  if (file_version %in% c("11.0", "11.0b")) {
    proteins <- read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE, quote="")
    proteins <- proteins[, c("protein_external_id", "preferred_name", "protein_size", "annotation")]
  } else {
    proteins <- read.table(path, sep="\t", skip=1, header=FALSE, stringsAsFactors=FALSE, fill=TRUE, quote="")
    colnames(proteins) <- c("protein_external_id", "preferred_name", "protein_size", "annotation")
  }
  proteins
}

map_to_string_ids <- function(my_data_frame, proteins_file, aliases_file, file_version) {
  proteins <- load_stringdb_proteins(proteins_file, file_version)
  aliases <- read.table(aliases_file, skip=1, sep="\t", header=FALSE, quote="", stringsAsFactors=FALSE, fill=TRUE)
  colnames(aliases) <- c("STRING_id", "alias", "sources")
  alias_df <- rbind(
    data.frame(STRING_id=proteins$protein_external_id, alias=proteins$preferred_name, stringsAsFactors=FALSE),
    data.frame(STRING_id=proteins$protein_external_id, alias=proteins$protein_external_id, stringsAsFactors=FALSE),
    aliases[, c("STRING_id", "alias")]
  )
  alias_df$alias_key <- toupper(iconv(alias_df$alias, "WINDOWS-1252", "UTF-8"))
  alias_df <- alias_df[!is.na(alias_df$alias_key) & alias_df$alias_key != "", ]
  alias_df <- dplyr::distinct(alias_df, alias_key, .keep_all=TRUE)
  mapped_df <- my_data_frame
  mapped_df$gene_key <- toupper(iconv(mapped_df$gene, "WINDOWS-1252", "UTF-8"))
  mapped_df$STRING_id <- alias_df$STRING_id[match(mapped_df$gene_key, alias_df$alias_key)]
  mapped_df <- mapped_df[!is.na(mapped_df$STRING_id), ]
  mapped_df$gene_key <- NULL
  mapped_df
}

load_stringdb_interactions <- function(path, string_ids, score_threshold=700) {
  links <- read.table(path, sep=" ", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
  links <- links[links$combined_score >= score_threshold, ]
  links <- links[links$protein1 %in% string_ids & links$protein2 %in% string_ids, ]
  data.frame(from=links$protein1, to=links$protein2, combined_score=links$combined_score, stringsAsFactors=FALSE)
}

get_stringdb_clusters <- function(edges_df, string_ids) {
  string_ids <- unique(string_ids[!is.na(string_ids)])
  edges_df <- edges_df[edges_df$from %in% string_ids & edges_df$to %in% string_ids, ]
  if (nrow(edges_df) == 0) {
    return(list())
  }
  graph <- graph_from_data_frame(edges_df[, c("from", "to")], directed=FALSE)
  if (ecount(graph) == 0 || vcount(graph) < 2) {
    return(list())
  }
  clusters <- cluster_fast_greedy(simplify(graph))
  cluster_list <- split(names(membership(clusters)), membership(clusters))
  cluster_list[order(lengths(cluster_list), decreasing=TRUE)]
}

# Get STRING id
dat_map <- map_to_string_ids(gene_diff, protein_info_file, protein_aliases_file, STRINGdb_file_version)
hits <- dat_map$STRING_id

# Limit the number of input genes to 2000 for STRINGdb::plot_network()
if (length(hits) > 2000) {
  hits_2000 <- hits[1:2000]
} else {
  hits_2000 <- hits
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
  node_labels <- dat_map$gene[match(node_ids, dat_map$STRING_id)]
  node_labels[is.na(node_labels)] <- node_ids[is.na(node_labels)]
  edge_weights <- sub_edges$combined_score / max(sub_edges$combined_score, na.rm = TRUE)
  node_fc <- dat_map$log2FC[match(node_ids, dat_map$STRING_id)]
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

if (length(hits) > 0) {
  # iGraph clustering
  ## algorithm: fastgreedy(default), walktrap, edge.betweenness
  all_interactions <- load_stringdb_interactions(protein_links_file, hits, score_threshold=700)
  clustersList <- get_stringdb_clusters(all_interactions, hits)
  if (length(clustersList) == 0) {
    write.csv(data.frame(node1=character(), node2=character(), combined_score=numeric()),
              string_link,
              row.names=FALSE,
              quote=FALSE)
    png(filename=top6_clusters_plot_png, units="in", width=4, height=4, res=400)
    plot.new()
    text(0.5, 0.5, "No PPI interactions", cex=1.2)
    dev.off()
    pdf(file=top6_clusters_plot_pdf, width=4, height=4)
    plot.new()
    text(0.5, 0.5, "No PPI interactions", cex=1.2)
    dev.off()
    quit(save="no", status=0)
  }
  # Plot first 6 clusters
  n_cluster <- 6
  if (length(clustersList) < 6) {
    n_cluster <- length(clustersList)
  } 
  png(filename = top6_clusters_plot_png,
      units = "in",
      width = 15,
      height = 10,
      res = 400)
  par(mfrow = c(2,3), mar = c(1, 1, 2, 1))
  for(i in 1:n_cluster){
    plot_ppi_cluster_local(all_interactions, clustersList[[i]], dat_map, cluster_id = i)
  }
  if (n_cluster < 6) {
    for (i in (n_cluster + 1):6) {
      plot.new()
    }
  }
  dev.off()
  pdf(file = top6_clusters_plot_pdf,
      width = 15,
      height = 10)
  par(mfrow = c(2,3), mar = c(1, 1, 2, 1))
  for(i in 1:n_cluster){
    plot_ppi_cluster_local(all_interactions, clustersList[[i]], dat_map, cluster_id = i)
  }
  if (n_cluster < 6) {
    for (i in (n_cluster + 1):6) {
      plot.new()
    }
  }
  dev.off()
  
  # Export interaction for downstream Cytoscape
  dat_link <- all_interactions
  ## Convert stringID to gene symbol
  dat_link$from <- dat_map[match(dat_link$from,dat_map$STRING_id),'gene']
  dat_link$to <- dat_map[match(dat_link$to,dat_map$STRING_id),'gene']  
  colnames(dat_link) <- c('node1','node2','combined_score')
  ## Remove duplicates
  dat_link <- dat_link %>% distinct(node1, node2, .keep_all = T)
  write.csv(dat_link,
            string_link,
            row.names = F,
            quote = F)
} else {
  # if no hits, generate empty table and plots
  library(glue)
  library(ggplot2)
  library(showtext)
  showtext_auto()
  
  # generate empty results table
  headers <- glue("node1\tnode2\tcombined_score\n")
  
  fileConn<-file(string_link)
  writeLines(headers, fileConn)
  close(fileConn)
  
  # generate empty plots
  text = "差异表达基因过少，无法进行富集分析！"
  p <- ggplot() + 
    annotate("text", x = 4, y = 25, size=8, label = text) + 
    theme_void() + 
    theme(panel.background = element_rect(fill = 'white'))
  ggsave(p, 
         filename = top6_clusters_plot_pdf, 
         device = "pdf",
         width = 10, 
         height = 10)
  ggsave(p, 
         filename = top6_clusters_plot_png, 
         device = "png",
         width = 4, 
         height = 4)
}
