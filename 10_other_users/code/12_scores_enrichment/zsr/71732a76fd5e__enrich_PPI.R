log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

work_dir <- getwd()

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
string_db <- STRINGdb$new(version = STRINGdb_version,
                          species = spe,   #human9606，mouse10090 
                          score_threshold = 700, #PPI score: default400, low150，high700, ultra high900
                          input_directory = STRINGdb_dir) 

# Get STRING id
dat_map <- string_db$map(my_data_frame = gene_diff, 
                         my_data_frame_id_col_names = "gene", #使用gene symbol或ENTREZID都可
                         removeUnmappedRows = TRUE )
hits <- dat_map$STRING_id

# Limit the number of input genes to 2000 for STRINGdb::plot_network()
if (length(hits) > 2000) {
  hits_2000 <- hits[1:2000]
} else {
  hits_2000 <- hits
}


if (length(hits) > 1) {
  # iGraph clustering
  ## algorithm: fastgreedy(default), walktrap, edge.betweenness
  clustersList <- string_db$get_clusters(string_ids = hits_2000 ,
                                         algorithm  = "fastgreedy" ) 
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
  par(mfrow = c(2,3))
  for(i in 1:n_cluster){
    string_db$plot_network(clustersList[[i]])
  }
  dev.off()
  pdf(file = top6_clusters_plot_pdf,
      width = 15,
      height = 10)
  par(mfrow = c(2,3))
  for(i in 1:n_cluster){
    string_db$plot_network(clustersList[[i]])
  }
  dev.off()
  
  # Export interaction for downstream Cytoscape
  dat_link <- string_db$get_interactions(hits)
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

