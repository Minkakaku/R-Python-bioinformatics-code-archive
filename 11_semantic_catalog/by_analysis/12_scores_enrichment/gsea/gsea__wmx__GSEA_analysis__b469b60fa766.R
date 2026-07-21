library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

# =========================
# 输入输出
# =========================
deg_files <- list.files("/mnt/NexData_file/home/wmx/Project/SC2603022014-01/shouhou6/DEG/",pattern = 'All_differentially_expressed_genes.xls',recursive = T,full.names = T)

output <- "/mnt/NexData_file/home/wmx/Project/SC2603022014-01/shouhou6/GSEA_GO_KEGG/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# =========================
# 读取差异基因
# =========================
for(deg_file in deg_files){
  deg <- read.csv(deg_file, header = TRUE, check.names = FALSE,sep = '\t')
  dirr = basename(dirname(deg_file))
  outdir = paste(output,dirr,sep = '/')
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  deg <- deg %>%
    dplyr::rename(
      gene = names,
      logFC = logfoldchanges
    ) %>%
    filter(!is.na(gene), !is.na(logFC))
  
  # 去掉重复基因，排序
  deg <- deg %>%
    arrange(desc(scores)) %>%
    distinct(gene, .keep_all = TRUE)
  # =========================
  # gene symbol 转 Entrez ID
  # =========================
  gene_df <- bitr(
    deg$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  deg2 <- deg %>%
    inner_join(gene_df, by = c("gene" = "SYMBOL"))
  
  # 如果多个 SYMBOL 对应同一个 ENTREZID，保留scores 最大的
  deg2 <- deg2 %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = scores, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  geneList <- deg2$scores
  names(geneList) <- deg2$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  cat("GSEA genes:", length(geneList), "\n")
  
  # =========================
  # GO GSEA
  # =========================
  gsea_go <- gseGO(
    geneList = geneList,
    OrgDb = org.Mm.eg.db,
    ont = "ALL",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE
  )
  
  go_res <- as.data.frame(gsea_go)
  
  # =========================
  # core_enrichment EntrezID 转 Symbol
  # =========================
  entrez_all <- unique(unlist(strsplit(go_res$core_enrichment, "/")))
  entrez_all <- entrez_all[!is.na(entrez_all) & entrez_all != ""]
  
  id2symbol <- bitr(
    entrez_all,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Mm.eg.db
  )
  
  id2symbol_vec <- setNames(id2symbol$SYMBOL, id2symbol$ENTREZID)
  
  go_res$core_enrichment_symbol <- sapply(go_res$core_enrichment, function(x) {
    ids <- unlist(strsplit(x, "/"))
    syms <- id2symbol_vec[ids]
    syms <- syms[!is.na(syms)]
    paste(syms, collapse = "/")
  })
  
  write.csv(
    go_res,
    file = file.path(outdir, "GSEA_GO_all.csv"),
    row.names = FALSE
  )
  
  go_sig <- go_res %>%
    filter(pvalue < 0.05)
  
  write.csv(
    go_sig,
    file = file.path(outdir, "GSEA_GO_pvalue0.05.csv"),
    row.names = FALSE
  )
  
  # =========================
  # KEGG GSEA
  # =========================
  gsea_kegg <- gseKEGG(
    geneList = geneList,
    organism = "mmu",
    keyType = "kegg",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE
  )
  
  kegg_res <- as.data.frame(gsea_kegg)
  
  # =========================
  # core_enrichment EntrezID 转 Symbol
  # =========================
  entrez_all <- unique(unlist(strsplit(kegg_res$core_enrichment, "/")))
  entrez_all <- entrez_all[!is.na(entrez_all) & entrez_all != ""]
  
  id2symbol <- bitr(
    entrez_all,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Mm.eg.db
  )
  
  id2symbol_vec <- setNames(id2symbol$SYMBOL, id2symbol$ENTREZID)
  
  kegg_res$core_enrichment_symbol <- sapply(kegg_res$core_enrichment, function(x) {
    ids <- unlist(strsplit(x, "/"))
    syms <- id2symbol_vec[ids]
    syms <- syms[!is.na(syms)]
    paste(syms, collapse = "/")
  })
  
  write.csv(
    kegg_res,
    file = file.path(outdir, "GSEA_KEGG_all.csv"),
    row.names = FALSE
  )
  
  kegg_sig <- kegg_res %>%
    filter(pvalue < 0.05)
  
  write.csv(
    kegg_sig,
    file = file.path(outdir, "GSEA_KEGG_pvalue0.05.csv"),
    row.names = FALSE
  )
  
  # =========================
  # 输出 GO pvalue 显著通路单通路图
  # =========================
  go_plot_dir <- file.path(outdir, "GO_single_pathway_pvalue0.05")
  dir.create(go_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(go_sig) > 0) {
    
    for (pathway in go_sig$ID) {
      
      pathway_name <- go_sig$Description[match(pathway, go_sig$ID)]
      pathway_title <- stringr::str_wrap(pathway_name, width = 45)
      
      p <- gseaplot2(
        gsea_go,
        geneSetID = pathway,
        title = pathway_title,
        base_size = 12
      )
      
      filename <- gsub("[^A-Za-z0-9_]", "_", paste0(pathway, "_", pathway_name))
      
      ggsave(
        file.path(go_plot_dir, paste0(filename, ".pdf")),
        p,
        width = 7,
        height = 5
      )
      
      ggsave(
        file.path(go_plot_dir, paste0(filename, ".png")),
        p,
        width = 7,
        height = 5,
        dpi = 300
      )
    }
  }
  
  # =========================
  # 输出 KEGG pvalue 显著通路单通路图
  # =========================
  kegg_plot_dir <- file.path(outdir, "KEGG_single_pathway_pvalue0.05")
  dir.create(kegg_plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(kegg_sig) > 0) {
    
    for (pathway in kegg_sig$ID) {
      
      pathway_name <- kegg_sig$Description[match(pathway, kegg_sig$ID)]
      pathway_name <- gsub(x = pathway_name,pattern = " - Mus musculus (house mouse)",replacement = "",fixed = TRUE)
      pathway_title <- stringr::str_wrap(pathway_name, width = 45)
      p <- gseaplot2(
        gsea_kegg,
        geneSetID = pathway,
        title = pathway_title,
        base_size = 12
      )
      
      filename <- gsub("[^A-Za-z0-9_]", "_", paste0(pathway, "_", pathway_name))
      
      ggsave(
        file.path(kegg_plot_dir, paste0(filename, ".pdf")),
        p,
        width = 7,
        height = 5
      )
      
      ggsave(
        file.path(kegg_plot_dir, paste0(filename, ".png")),
        p,
        width = 7,
        height = 5,
        dpi = 300
      )
    }
  }
  
  cat("GO pvalue < 0.05:", nrow(go_sig), "\n")
  cat("KEGG pvalue < 0.05:", nrow(kegg_sig), "\n")
  cat("Done:", outdir, "\n")
}