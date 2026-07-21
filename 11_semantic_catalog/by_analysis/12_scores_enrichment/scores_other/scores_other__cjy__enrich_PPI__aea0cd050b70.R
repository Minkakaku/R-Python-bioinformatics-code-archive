log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

work_dir <- getwd()

options(timeout = 1200)
options(stringsAsFactors = F)
library(tidyverse)
library(STRINGdb)
library(igraph)

# Input
## DEG results
DEG <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE)
## Species
if (snakemake@params[["species"]] == "homo_sapiens") {
  spe <- 9606
} else if (snakemake@params[["species"]] == "mus_musculus") {
  spe <- 10090
} 
## Local STRINGdb
STRINGdb_version <- snakemake@params[["STRINGdb_version"]]
STRINGdb_dir <- snakemake@params[["STRINGdb_dir"]]

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

# Set up STRINGdb
print(1)
string_db <- STRINGdb$new(version = STRINGdb_version,
                          species = spe,   #human9606，mouse10090 
                          score_threshold = 700, #PPI score: default400, low150，high700, ultra high900
                          input_directory = STRINGdb_dir,
                          protocol="http") 
print(2)

# Get STRING id
dat_map <- string_db$map(my_data_frame = gene_diff, 
                         my_data_frame_id_col_names = "gene", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )
hits <- dat_map$STRING_id

print(3)

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


if (length(hits) > 1) {
  # iGraph clustering
  ## algorithm: fastgreedy(default), walktrap, edge.betweenness
  clustersList <- string_db$get_clusters(string_ids = hits_2000 ,
                                         algorithm  = "fastgreedy" ) 
  all_interactions <- string_db$get_interactions(hits)
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
