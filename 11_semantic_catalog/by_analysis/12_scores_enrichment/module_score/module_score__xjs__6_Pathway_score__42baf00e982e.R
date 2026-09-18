
pathways <- readRDS("xenium5k_related_pathway.rds")

pathways_use <- pathways[c("HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_DNA_REPAIR", "HALLMARK_GLYCOLYSIS", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_NOTCH_SIGNALING", "HALLMARK_P53_PATHWAY", "HALLMARK_COMPLEMENT")]


all_data_score <- 
  AddModuleScore(all_data, 
                 features = pathways_use, 
                 name = names(pathways_use))
colnames(all_data_score@meta.data)[51:67] <- gsub("(HALLMARK_.+?)\\d+$", "\\1", colnames(all_data_score@meta.data)[51:67])
colnames(all_data_score@meta.data)[63] <- gsub("HALLMARK_MYC_TARGETS_V", "HALLMARK_MYC_TARGETS_V1", colnames(all_data_score@meta.data)[63])

dir.create("Pathways_score")


data_pathway_score <- FetchData(all_data_score, 
                                vars = c(
                                  "RCTD_subtype",
                                  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                                  "HALLMARK_HYPOXIA",
                                  "HALLMARK_APOPTOSIS",
                                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                  "HALLMARK_DNA_REPAIR",
                                  "HALLMARK_GLYCOLYSIS",
                                  "HALLMARK_TGF_BETA_SIGNALING", 
                                  "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
                                  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                                  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                                  "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                                  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                                  "HALLMARK_MYC_TARGETS_V1", 
                                  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
                                  "HALLMARK_NOTCH_SIGNALING", 
                                  "HALLMARK_P53_PATHWAY", 
                                  "HALLMARK_COMPLEMENT"))



pathway_cols <- colnames(data_pathway_score)[-1]
plot_list <- list()
for (i in seq_along(pathway_cols)) {
  pathway <- pathway_cols[i]
  
  p <- ggplot(data_pathway_score, aes_string(x = "RCTD_subtype", y = pathway, fill = "RCTD_subtype")) + 
    geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.8) + 
    scale_fill_manual(values = plan) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "none") + 
    labs(title = pathway, x = "CellType", y = "Pathway score")
  
  plot_list[[pathway]] <- p
}

pdf("Pathways_score/All_pathway_score_boxplot.pdf", width = 10, height = 6)
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
}
dev.off()








library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(gtools)


celltype <- sort(unique(all_data_score$RCTD_subtype))
celltype_raw <- unique(all_data_score$RCTD_subtype)
celltype <- mixedsort(as.character(celltype_raw))
print(celltype)

pdf("Pathways_score/All_celltype_pathway_scores.pdf", width = 15, height = 6)

for (celltype_id in celltype) {
  celltype_data <- all_data_score@meta.data %>% 
    filter(RCTD_subtype == celltype_id) %>% 
    tibble::rownames_to_column("cell_id") %>% 
    dplyr::select(cell_id, all_of(pathway_cols))
  
  cell_count <- nrow(celltype_data)
  
  if(cell_count > 0) {
    # 转换为长格式
    celltype_long <- celltype_data %>%
      pivot_longer(
        cols = all_of(pathway_cols),
        names_to = "Score_Type",
        values_to = "Score_Value"
      )
    
    # 创建图形
    p <- ggplot(celltype_long, aes(x = Score_Type, y = Score_Value, fill = Score_Type)) +
      geom_violin(alpha = 0.6) +
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
      labs(
        title = paste(celltype_id, "Pathway Score"),
        subtitle = paste("cell_count:", cell_count),
        x = "Pathway",
        y = "Pathway Score"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none"
      ) +
      scale_fill_viridis_d()
    
    # 打印图形到当前PDF页面
    print(p)
    cat("已添加 Cluster", celltype_id, "到PDF文件\n")
  } else {
    cat("Cluster", celltype_id, "没有细胞，跳过\n")
  }
}

dev.off()






# pathway_gene_available <- lapply(pathways_use, function(genes) {
#   intersect(genes, rownames(all_data_score))
# })
# 
# pathway_gene_available_gene_counts <- sapply(pathway_gene_available, length)
# pathways_use_gene_counts <- sapply(pathways_use, length)
# 
# pathway_gene_available_df <- data.frame(
#   Pathway = names(gene_counts),
#   Pathway_gene_available_gene_counts = pathway_gene_available_gene_counts, 
#   Pathways_use_gene_counts = pathways_use_gene_counts
# )
# print(pathway_gene_available_df)




lines <- lapply(names(pathways_use), function(pathway_name) {
  genes <- pathways_use[[pathway_name]]
  # 将通路名称和基因用制表符连接
  paste(c(pathway_name, genes), collapse = "\t")
})

writeLines(unlist(lines), "Pathways_score/pathways_output.txt")


