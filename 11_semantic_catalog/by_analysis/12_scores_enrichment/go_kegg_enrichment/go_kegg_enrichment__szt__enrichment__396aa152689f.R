type <- unique(all_data$annotation2)

library("org.Hs.eg.db")
library(clusterProfiler)
type %>%
  setNames(.,.) %>%
  imap(function(celltype,name){ #browser()
    Idents(all_data)<-'annotation2'
    diffmarkers <-
      FindMarkers(
        JoinLayers(all_data),
        ident.1 = 'Post_treatment',
        ident.2 = 'During_treatment',
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        group.by = 'group_analysis',
        subset.ident = celltype
      )
    dir.create(file.path(str_glue('diffgene_enrichment/post_during/{celltype}')))
    write.csv(diffmarkers,file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} group diffmarker.csv')))
    enrich_ls <- enrichGO(
      rownames(diffmarkers),
      OrgDb = org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL')
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} Go enriment table.csv')),row.names = FALSE)
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
      filename = file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} Go dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} Go barplot.pdf')),
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
          use_internal_data = FALSE
        )
    enrich_ls <- setReadable(enrich_ls, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} KEGG enriment table.csv')),row.names = FALSE)
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
      filename = file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} KEGG dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/post_during/{celltype}/{celltype} KEGG barplot.pdf')),
      height = 6,
      width = 8
    )
  })


type %>%
  setNames(.,.) %>%
  imap(function(celltype,name){ #browser()
    Idents(all_data)<-'annotation2'
    diffmarkers <-
      FindMarkers(
        JoinLayers(all_data),
        ident.1 = 'During_treatment',
        ident.2 = 'Pre_treatment',
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        group.by = 'group_analysis',
        subset.ident = celltype
      )
    dir.create(file.path(str_glue('diffgene_enrichment/during_pre/{celltype}')))
    write.csv(diffmarkers,file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} group diffmarker.csv')))
    enrich_ls <- enrichGO(
      rownames(diffmarkers),
      OrgDb = org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL')
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} Go enriment table.csv')),row.names = FALSE)
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
      filename = file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} Go dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} Go barplot.pdf')),
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
        use_internal_data = FALSE
      )
    enrich_ls <- setReadable(enrich_ls, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} KEGG enriment table.csv')),row.names = FALSE)
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
      filename = file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} KEGG dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/during_pre/{celltype}/{celltype} KEGG barplot.pdf')),
      height = 6,
      width = 8
    )
  })







type <- c(unique(All_T$subtype),
          unique(Macrophage$subtype),
          unique(Fibroblast$subtype)) 


type %>%
  setNames(.,c(rep('All_T',11),rep('Macrophage',7),rep('Fibroblast',7))) %>%
  imap(function(celltype,name){ #browser()
    temp <- get(name)
    Idents(temp)<-'subtype'
    diffmarkers <-
      FindMarkers(
        JoinLayers(temp),
        ident.1 = 'Post_treatment',
        ident.2 = 'During_treatment',
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        group.by = 'group_analysis',
        subset.ident = celltype
      )
    dir.create(file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/')), recursive = TRUE)
    write.csv(diffmarkers,file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} group diffmarker.csv')))
    enrich_ls <- enrichGO(
      rownames(diffmarkers),
      OrgDb = org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL')
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} Go enriment table.csv')),row.names = FALSE)
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
      filename = file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} Go dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} Go barplot.pdf')),
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
        use_internal_data = FALSE
      )
    enrich_ls <- setReadable(enrich_ls, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.csv(enrich_ls@result,file = file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} KEGG enriment table.csv')),row.names = FALSE)
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
      filename = file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} KEGG dotplot.pdf')),
      height = 6,
      width = 8
    )
    ggsave(
      enrich_bar,
      filename = file.path(str_glue('diffgene_enrichment/{name}/post_during/{celltype}/{celltype} KEGG barplot.pdf')),
      height = 6,
      width = 8
    )
  })



type<-c('LC_BM_cancer_cells','M_BM_cancer_cells','OM_BM_cancer_cells','UTC_BM_cancer_cells')
type%>%
  walk(function(celltype){
    DefaultAssay(harmony) <- "peaks"
    diffpeak <- 
      FindMarkers(harmony,
                  ident.1 = celltype,
                  ident.2 = "LC_cancer_cells",
                  test.use = "LR",
                  min.pct = 0.3, 
                  logfc.threshold = 1,
                  latent.vars = 'nCount_peaks') 
    write.csv(diffpeak,str_glue('/home/linrui/ATACnew/diff_gene_peaks/Primary_each_metastatic/{celltype}_LC_diffpeaks.csv'))
    top <- rbind(diffpeak %>% top_n(n = 25, wt = avg_log2FC),diffpeak %>% top_n(n = -25, wt = avg_log2FC))
    a<-harmony[,harmony$sample_type=='LC_cancer_cells'|harmony$sample_type==celltype]
    a[['RNA']]<-a[['RNA']] %>% as('Assay5')
    a.genes <- rownames(a) 
    a <- ScaleData(a, features =a.genes) 
    p<-DoHeatmap(a, label = F,
                 features = rownames(top),
                 group.by = "group1",
                 disp.max= 1,
                 disp.min= -1,
                 assay = 'peaks')+ #设置组别颜色
      scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
    pdf(str_glue("/home/linrui/ATACnew/diff_gene_peaks/Primary_each_metastatic/{celltype}_LC_diffpeaks_heatmap.pdf"), height = 8, width =10)
    print(p)
    dev.off()
  })