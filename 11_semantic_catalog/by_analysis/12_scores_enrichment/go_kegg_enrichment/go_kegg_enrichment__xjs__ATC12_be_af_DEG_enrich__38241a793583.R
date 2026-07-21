library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create('lubin')
out_dir <- 'lubin/ATC12_be_af_DEG_enrich'
dir.create(out_dir)
all_data <- readRDS("R_OUTPUT/all_data_lubin.rds")

all_data$sample %>% table()

ATC12 <- all_data[, all_data$sample %in% c('ATC12_af_T_1', 'ATC12_af_T_2', 'ATC12_be_T_1', 'ATC12_be_T_2')]
ATC12$group <- str_replace_all(ATC12$sample, c(
  'ATC12_af_T_1' = 'after', 
  'ATC12_af_T_2' = 'after', 
  'ATC12_be_T_1' = 'before',
  'ATC12_be_T_2' = 'before'
))
ATC12$group %>% table()


ATC12_markers <- 
  FindMarkers(
    JoinLayers(ATC12), 
    group.by = 'group',
    test.use = "wilcox", 
    ident.1 = 'after',
    ident.2 = 'before',
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
ATC12_markers %>% dim()
write.csv(ATC12_markers, file.path(out_dir, 'ATC12_af_be_DEGs.csv'))


# enrichment -----------
ATC12_de_genes_for_enrichmenet <- 
  ATC12_markers%>% 
  rownames_to_column('gene') %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1-pct.2 > 0.01) %>%
  dplyr::filter(abs(avg_log2FC) > 0.05) %>% 
  pull("gene")
ATC12_de_genes_for_enrichmenet <- list(ATC12_de_genes_for_enrichmenet)
names(ATC12_de_genes_for_enrichmenet) <- 'ATC12_af_vs_be'


enrich_ls <- 
  ATC12_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output <- file.path(out_dir, 'enrichment') 
dir.create(er_output)


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

ATC12_de_genes_id_for_enrichmenet <-
  ATC12_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  ATC12_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = T
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(ATC12, file.path(ATC12_dir, 'ATC12.rds'))