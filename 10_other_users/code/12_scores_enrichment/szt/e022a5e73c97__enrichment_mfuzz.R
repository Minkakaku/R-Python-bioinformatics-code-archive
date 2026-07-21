#enrichment---------------
type <- unique(all_data$annotation)
type<-c("M2")
library(EnhancedVolcano)
library("org.Hs.eg.db")
library(clusterProfiler)
type %>%
  setNames(.,.) %>%
  imap(function(celltype,name){ #browser()
    Idents(Macrophage)<-'subtype'
    diffmarkers <-
      FindMarkers(
        Macrophage,
        ident.1 = celltype,
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.01,
        group.by = 'subtype'
      )
    dir.create(file.path(str_glue('diffgene_enrichment/{celltype}')),recursive = TRUE)
    write.csv(diffmarkers,file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} group diffmarker.csv')))
    enrich_ls <- enrichGO(
      rownames(diffmarkers),
      OrgDb = org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL')
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} Go enriment table.csv')),row.names = FALSE)
    enrich_dot <- enrichplot::dotplot(
      enrich_ls,
      showCategory = 10,
      label_format = 50)
    enrich_bar <- barplot(
      enrich_ls,
      showCategory = 10,
      label_format = 50)
    ggsave(
      enrich_dot,
      filename = file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} Go dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} Go barplot.pdf')),
      height = 6,
      width = 8
    )
    de_genes_id_for_enrichmenet <- 
      clusterProfiler::bitr(
        rownames(diffmarkers),
        OrgDb = "org.Hs.eg.db",
        fromType = "SYMBOL",
        toType = "ENTREZID"
      ) %>%
      .$ENTREZID
    enrich_ls <- 
      enrichKEGG(
        de_genes_id_for_enrichmenet,
        organism = "hsa", 
        keyType = "ncbi-geneid", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        use_internal_data = TRUE
      )
    enrich_ls <- setReadable(enrich_ls, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} KEGG enriment table.csv')),row.names = FALSE)
    enrich_dot <- enrichplot::dotplot(
      enrich_ls,
      showCategory = 10,
      label_format = 50)
    enrich_bar <- barplot(
      enrich_ls,
      showCategory = 10,
      label_format = 50)
    ggsave(
      enrich_dot,
      filename = file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} KEGG dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/{celltype}/{celltype} KEGG barplot.pdf')),
      height = 6,
      width = 8
    )
    
    p_volcano <- EnhancedVolcano(diffmarkers,
                                 lab = rownames(diffmarkers),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 title = str_glue('{celltype} volcano plot'),
                                 pCutoff = 10e-2,
                                 FCcutoff = 0,
                                 pointSize = 3.0,
                                 labSize = 6.0)
    ggsave(p_volcano,
           filename = file.path(str_glue("diffgene_enrichment/{celltype}/{celltype} diffgene valcano.pdf")),
           height = 6,
           width = 8)
  })



Idents(Epithelial) <- "RNA_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    Epithelial,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )


#聚类图
library(Mfuzz)
group_average <- AverageExpression(all_data,group.by = c('annotation','group.y'))
type<-c("Epithelial","Macrophage","Fibroblasts","MergeDC","pDCs")
type<-c("T cells")
type %>%
  setNames(.,.) %>%
  imap(function(celltype,name){ #browser()
    marker <- read.csv(file.path(str_glue("diffgene_enrichment/{celltype}/{celltype} group diffmarker.csv"))) %>% pull(X)
    CA_temp <- group_average %>%
      as.data.frame() %>% 
      select(contains('CA')) %>% 
      # select(contains({celltype})) %>%
      select(contains('T.cells')) %>%
      # rename_with(~ gsub(str_glue("RNA\\.{celltype}_"), "", .), everything()) %>%  # 移除前缀
      rename_with(~ gsub(str_glue("RNA\\.T.cells_"), "", .), everything()) %>%  # 移除前缀
      rename_with(~ gsub("\\.", "_", .), everything()) %>% 
      filter(rownames(.) %in% marker)%>%
      select(N0_CA, Ni_CA, Nmi_CA, N1_CA)
    
    protein <- as.matrix(CA_temp)
    
    mfuzz_class <- new('ExpressionSet',exprs = protein)
    mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
    mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
    mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
    
    mfuzz_class <- standardise(mfuzz_class)
    
    set.seed(123)
    cluster_num <- 10
    mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

    pdf(
      file.path('mfuzz',str_glue("CA_{celltype}.pdf")),
      height = 6,
      width  = 15
    )
    mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = colnames(protein), x11 = FALSE)
    dev.off()
    
    saveRDS(mfuzz_class,str_glue("mfuzz/CA_{celltype}_mfuzz_class.rds"))
    saveRDS(mfuzz_cluster,str_glue("mfuzz/CA_{celltype}_mfuzz_cluster.rds"))
  })



T_average <- AverageExpression(All_T,group.by = c('subtype','group.y'))
type<-c("Treg")
type %>%
  setNames(.,.) %>%
  imap(function(celltype,name){ #browser()
    marker <- read.csv(file.path(str_glue("diffgene_enrichment/{celltype}/{celltype} group diffmarker.csv"))) %>% pull(X)
    CA_temp <- T_average %>%
      as.data.frame() %>% 
      select(contains('CA')) %>% 
      select(contains({celltype})) %>%
      rename_with(~ gsub(str_glue("RNA\\.{celltype}_"), "", .), everything()) %>%  # 移除前缀
      rename_with(~ gsub("\\.", "_", .), everything()) %>% 
      filter(rownames(.) %in% marker)%>%
      select(N0_CA, Ni_CA, Nmi_CA, N1_CA)
    
    protein <- as.matrix(LN_temp)
    
    mfuzz_class <- new('ExpressionSet',exprs = protein)
    mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
    mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
    mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
    
    mfuzz_class <- standardise(mfuzz_class)
    
    set.seed(123)
    cluster_num <- 10
    mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
    
    pdf(
      file.path('mfuzz',str_glue("CA_{celltype}.pdf")),
      height = 6,
      width  = 15
    )
    mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = colnames(protein), x11 = FALSE)
    dev.off()
    
    saveRDS(mfuzz_class,str_glue("mfuzz/CA_{celltype}_mfuzz_class.rds"))
    saveRDS(mfuzz_cluster,str_glue("mfuzz/CA_{celltype}_mfuzz_cluster.rds"))
  })


M_average <- AverageExpression(Macrophage,group.by = c('subtype','group.y'))
type<-c("M2")
type %>%
  setNames(.,.) %>%
  imap(function(celltype,name){ #browser()
    marker <- read.csv(file.path(str_glue("diffgene_enrichment/{celltype}/{celltype} group diffmarker.csv"))) %>% pull(X)
    CA_temp <- M_average %>%
      as.data.frame() %>% 
      select(contains('CA')) %>% 
      select(contains({celltype})) %>%
      rename_with(~ gsub(str_glue("RNA\\.{celltype}_"), "", .), everything()) %>%  # 移除前缀
      rename_with(~ gsub("\\.", "_", .), everything()) %>% 
      filter(rownames(.) %in% marker)%>%
      select(N0_CA, Ni_CA, Nmi_CA, N1_CA)
    
    protein <- as.matrix(CA_temp)
    
    mfuzz_class <- new('ExpressionSet',exprs = protein)
    mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
    mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
    mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
    
    mfuzz_class <- standardise(mfuzz_class)
    
    set.seed(123)
    cluster_num <- 10
    mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
    
    pdf(
      file.path('mfuzz',str_glue("CA_{celltype}.pdf")),
      height = 6,
      width  = 15
    )
    mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(2, 5), time.labels = colnames(protein), x11 = FALSE)
    dev.off()
    
    saveRDS(mfuzz_class,str_glue("mfuzz/CA_{celltype}_mfuzz_class.rds"))
    saveRDS(mfuzz_cluster,str_glue("mfuzz/CA_{celltype}_mfuzz_cluster.rds"))
  })

