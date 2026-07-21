library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

T_out_dir <- 'T_subtype/annotation/'
Myeloid_out_dir <- 'Myeloid_subtype/'
Fibro_out_dir <- 'Fibroblast_subtype/'

All_T_anno_remove_Tumor3 <- readRDS('T_subtype/All_T_anno_remove_Tumor3.rds')
Myeloid <- readRDS(file.path(Myeloid_out_dir, 'Myeloid.rds'))
Fibroblast <- readRDS(file.path(Fibro_out_dir, 'Fibroblast.rds'))

all_data <- readRDS(file.path("R_OUTPUT", "all_data_remove_Tumor3.rds"))

# 映射---------------------------------------------------------------------------
All_T_anno_remove_Tumor3$subtype %>% table() 
# DimPlot(All_T_anno_remove_Tumor3, group.by = 'subtype', reduction = 'harmonyumap')
Myeloid$subtype %>% table()
# DimPlot(Myeloid, group.by = 'subtype', reduction = 'harmonyumap', label = T) + NoLegend()
Fibroblast$subtype %>% table()
# DimPlot(Fibroblast, group.by = 'subtype', reduction = 'harmonyumap', label = T) + NoLegend()


library(tibble)
library(data.table)
subtype <- rbind(All_T_anno_remove_Tumor3[['subtype']] %>% rownames_to_column('cells'), 
                 Myeloid[['subtype']] %>% rownames_to_column('cells'), 
                 Fibroblast[['subtype']] %>% rownames_to_column('cells'))
meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_full = ifelse((as.character(annotation) %in% c('T cell', 'Myeloid cell','Fibroblast')),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation_full) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data, meta.data.t['annotation_full'])
all_data$annotation %>% table() %>% sum()
all_data$annotation_full %>% table() %>% sum()
all_data2 <- subset(all_data, annotation_full != 'NA')
all_data2$annotation <- str_replace_all(all_data2$annotation, c('Pancreatic β-cell' = 'Pancreatic-beta cell'))
all_data2$annotation_full <- str_replace_all(all_data2$annotation_full, c('Pancreatic β-cell' = 'Pancreatic-beta cell'))

all_data$group <- factor(all_data$group, levels = c('N', 'CP', 'T'))
all_data2$group <- factor(all_data2$group, levels = c('N', 'CP', 'T'))
all_data2$sample <- factor(all_data2$sample, levels = c(
  'Normal1', 'Normal2', 'Normal3', 'Normal4', 
  'CP1', 'CP2', 'CP3', 'CP4', 
  'Tumor1', 'Tumor2', 'Tumor4'
  ))
saveRDS(all_data, file.path("R_OUTPUT", "all_data_remove_Tumor3.rds"))
saveRDS(all_data2, 'all_data_remove_Tumor3_annotation_full.rds')




# enrichment analysis -----------------------------------------------------------
all_data2$annotation_full %>% table()
## T cell -----------------------------------------------------------------------
All_T <- subset(all_data2, annotation_full %in% c(
  'CD4 Treg', 
  'CD4 Trm', 
  'CD8 T', 
  'CD8 Tem', 
  'CD8 Tnaive', 
  'Cycling T', 
  'gdT'
))




### cell props barplot  --------
All_T_meta <- All_T[[]]

p_col_proj <-
  All_T_meta %>%
  group_by(sample, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation_full)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(T_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_group <-
  All_T_meta %>%
  group_by(group, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_group,
  filename = file.path(T_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


### volcano-------
Idents(All_T) <- 'group'
groups <- c('T', 'CP', 'N')
all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
all_combns_ls <- list(all_combns_ls[[1]], all_combns_ls[[3]])
names(all_combns_ls) <- 
  all_combns_ls %>% 
  map(function(x){ # browser()
    grp <- groups[x]
    str_glue('{grp[1]} vs {grp[2]}')
  }) %>% 
  unlist()
  
T_subtype_ls <- All_T$annotation_full %>% unique()
names(T_subtype_ls) <- T_subtype_ls
  
T_markers_ls <- 
  T_subtype_ls %>% 
  map(function(T_subtype){ # browser()
    marker_ls <- 
      all_combns_ls %>% 
      map(function(x){  
        grp <- groups[x]
        FindMarkers(
          subset(JoinLayers(All_T), annotation_full == T_subtype),
          # All_T, 
          group.by = 'group',
          test.use = "wilcox", 
          ident.1 = grp[1],
          ident.2 = grp[2],
          logfc.threshold = 0.25,
          min.pct = 0.1,
          only.pos = TRUE
        ) %>% 
          rownames_to_column('gene')
      })
      
    marker_ls
  })

T_markers_ls %>% 
  iwalk(function(T_marker, T_subtype){
    openxlsx::write.xlsx(T_marker, file.path(T_out_dir, str_glue('{T_subtype}_diff_gene.xlsx')))
  })



library(EnhancedVolcano)

volcano_ls <- 
  T_markers_ls %>%
  imap(function(T_markers, subtype){
    T_markers %>% 
      imap(function(T_marker, group){
        p <- EnhancedVolcano(T_marker,
                             lab = T_marker$gene,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = str_glue('{subtype}: {group} volcano plot'),
                             pCutoff = 10e-2,
                             FCcutoff = 0,
                             pointSize = 3.0,
                             labSize = 6.0)
        p
      })
  })

volcano_ls %>% 
  iwalk(function(x, subtype){
    x %>% 
      imap(function(valcano, group){
        ggsave(
          valcano,
          filename = file.path(T_out_dir, str_glue("volcano/{subtype}_{group} diffgene valcano.pdf")),
          height = 6,
          width = 8
        )
      })
  })




### enrichment ----------------------------
T_de_genes_for_enrichmenet_ls <- 
  T_markers_ls %>% 
  map(function(T_markers){
    T_markers %>% 
      map(function(T_marker){
        T_marker %>% 
          dplyr::filter(p_val_adj < 0.05) %>%
          dplyr::filter(pct.1 - pct.2 > 0.1) %>%
          dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
          nest() %>%
          deframe() %>%
          map("gene")
      })
  })


enrich_ls <- 
  T_de_genes_for_enrichmenet_ls %>% 
  map(function(T_de_genes_for_enrichmenet_subtype){ 
    T_de_genes_for_enrichmenet_subtype %>% 
      map(function(T_de_genes_for_enrichmenet_group){ # browser()
        T_de_genes_for_enrichmenet_group %>% 
          map(function(T_de_genes_for_enrichmenet){
            clusterProfiler::enrichGO(
              T_de_genes_for_enrichmenet,
              OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
              ont = 'ALL',
              pAdjustMethod = 'BH',
              pvalueCutoff = 0.05, 
              qvalueCutoff = 0.05,
              keyType = 'SYMBOL'
            )
          })
      })
  })
  

er_output <- file.path(T_out_dir, 'enrichment') # dir.create(er_output)

# expot enrichment table
enrich_ls %>% 
  iwalk(function(enrich, subtype){ # browser()
    dir.create(file.path(er_output, subtype))
    
    enrich %>% 
      iwalk(function(enrich_result, name){ # browser()
        enrich_result %>% 
          discard(function(x) length(x) < 1) %>% 
          walk(function(x){
            write.csv(
              x@result,
              file = file.path(er_output, subtype, str_glue("{name} Go enrichment table.csv")),
              row.names = FALSE
            )
          })
      })
})

# dotplot
enrich_dot_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })

enrich_bar_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })


# export dotplot
enrich_dot_ls %>% 
  iwalk(function(enrich_dot, subtype){
    enrich_dot %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} go dotplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })
  

enrich_bar_ls %>% 
  iwalk(function(enrich_bar, subtype){
    enrich_bar %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} go barplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


T_de_genes_id_for_enrichmenet_ls <-
  T_de_genes_for_enrichmenet_ls %>% 
  map(function(T_de_genes_for_enrichmenet_subtype){
    T_de_genes_for_enrichmenet_subtype %>% 
      map(function(T_de_genes_for_enrichmenet_group){
        T_de_genes_for_enrichmenet_group %>% 
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
      })
  })
 


enrich_ls <- 
  T_de_genes_id_for_enrichmenet_ls %>% 
  map(function(T_de_genes_id_for_enrichmenet_subtype){
    T_de_genes_id_for_enrichmenet_subtype %>% 
      map(function(T_de_genes_id_for_enrichmenet_group){
        T_de_genes_id_for_enrichmenet_group %>% 
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
              use_internal_data = FALSE
            )
          })
      })
    
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>% 
  iwalk(function(enrich, subtype){ # browser()
    dir.create(file.path(er_output, subtype))
    
    enrich %>% 
      iwalk(function(enrich_result, name){ # browser()
        enrich_result %>% 
          discard(function(x) length(x) < 1) %>% 
          walk(function(x){
            write.csv(
              x@result,
              file = file.path(er_output, subtype, str_glue("{name} KEGG enrichment table.csv")),
              row.names = FALSE
            )
          })
      })
  })


# dotplot
enrich_dot_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })

enrich_bar_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })


# export dotplot
enrich_dot_ls %>% 
  iwalk(function(enrich_dot, subtype){
    enrich_dot %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          purrr::discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} kegg dotplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


enrich_bar_ls %>% 
  iwalk(function(enrich_bar, subtype){
    enrich_bar %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          purrr::discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} kegg barplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })




### subtype group freq diff -----------------------------------------------------
Tsub_boxplot_out_dir <- file.path(T_out_dir, 'Percentage_diff_boxplot')
dir.create(Tsub_boxplot_out_dir)
All_T$annotation_full %>% table()
All_T_count <- 
  All_T@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(sample, group) %>% 
  summarise(count_T = n())
Tsub_meta_ls <- 
  All_T@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(annotation_full) %>% 
  nest() %>% 
  deframe() %>% 
  map(function(T_meta){
    T_meta %>% 
      group_by(sample, group) %>% 
      summarise(count_Tsub = n()) %>% 
      left_join(All_T_count, by = c('sample', 'group')) %>% 
      mutate(Percentage = count_Tsub/count_T*100) %>% 
      group_by(group) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat){
        qValue <- 
          dat[['Percentage']] %>% 
          quantile()
        
        dat %>% 
          mutate(Q2 = qValue[2], 
                 Q4 = qValue[4])
      }) %>% 
      enframe(name = 'group') %>% 
      unnest() %>% 
      mutate(group = factor(group, levels = c('N', 'CP', 'T')))
  })

Tsub_meta_ls %>% 
  imap(function(Tsub_meta, Tsub){
    print(Tsub)
    print(shapiro.test(Tsub_meta$Percentage))
    print(bartlett.test(Percentage~group, data = Tsub_meta))
  })
saveRDS(Tsub_meta_ls, file.path(T_out_dir, 'Tsub_meta_ls.rds'))

Tsub_meta_tTest_ls <- Tsub_meta_ls[c('gdT', 'CD8 T', 'CD4 Treg')]
Tsub_meta_wilcox_ls <- Tsub_meta_ls[c('CD8 Tem', 'CD8 Tnaive', 'CD4 Trm', 'Cycling T')]
  
box_p_tTest_ls <- 
  Tsub_meta_tTest_ls %>% 
  imap(function(Tsub_meta, Tsub){
    Tsub_meta %>% 
      ggplot(aes(x = group, y = Percentage, fill = group)) +
      geom_errorbar(aes(ymin = Q2, ymax = Q4), color = 'darkgrey') + 
      geom_jitter(
        width = 0.2, 
        height = 0.1, 
        alpha = 1, 
        size = 2, 
        aes(color = group)
      ) +
      ggpubr::geom_pwc(
        tip.length = 0, 
        hide.ns = F, 
        label = 'p', 
        # method.args = list(paired = F, exact = T), 
        # p.adjust.method = "BH", 
        method = 't_test',
        # method = 'wilcox_test',
        # group.by = "x.var"
        group.by = "legend.var"
      ) +  
      scale_color_manual(values = c(
        'N' = 'lightgreen', 
        'CP' = 'lightblue', 
        'T' = 'pink')
      ) +
      ggtitle(Tsub) + 
      theme_bw()  
  })

box_p_tTest_ls %>% 
  iwalk(function(box_p, Tsub){ # browser()
    pdf(file.path(Tsub_boxplot_out_dir, str_glue('{Tsub} by group.pdf')), width = 6, height = 6)
    print(box_p)
    dev.off()
  })


box_p_wilcox_ls <- 
  Tsub_meta_wilcox_ls %>% 
  imap(function(Tsub_meta, Tsub){
    Tsub_meta %>% 
      ggplot(aes(x = group, y = Percentage, fill = group)) +
      geom_errorbar(aes(ymin = Q2, ymax = Q4), color = 'darkgrey') + 
      geom_jitter(
        width = 0.2, 
        height = 0.1, 
        alpha = 1, 
        size = 2, 
        aes(color = group)
      ) +
      ggpubr::geom_pwc(
        tip.length = 0, 
        hide.ns = F, 
        label = 'p', 
        # method.args = list(paired = F, exact = T), 
        # p.adjust.method = "BH", 
        # method = 't_test',
        method = 'wilcox_test',
        # group.by = "x.var"
        group.by = "legend.var"
      ) +  
      scale_color_manual(values = c(
        'N' = 'lightgreen', 
        'CP' = 'lightblue', 
        'T' = 'pink')
      ) +
      ggtitle(Tsub) + 
      theme_bw()  
  })

box_p_wilcox_ls %>% 
  iwalk(function(box_p, Tsub){ # browser()
    pdf(file.path(Tsub_boxplot_out_dir, str_glue('{Tsub} by group.pdf')), width = 6, height = 6)
    print(box_p)
    dev.off()
  })




## Myeloid -----------------------------------------------------------------------
Myeloid <- subset(all_data2, annotation_full %in% c(
  'cDC2', 
  'IFN-responsive macrophage', 
  'LAMP3 DC', 
  'M1 macrophage', 
  'M2 macrophage', 
  'Neutropil', 
  'Proliferating macrophage'
))




### cell props barplot  --------
Myeloid_meta <- Myeloid[[]]

p_col_proj <-
  Myeloid_meta %>%
  group_by(sample, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation_full)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Myeloid_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_group <-
  Myeloid_meta %>%
  group_by(group, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_group,
  filename = file.path(Myeloid_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)




### volcano-------
Idents(Myeloid) <- 'group'
groups <- c('T', 'CP', 'N')
all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
all_combns_ls <- list(all_combns_ls[[1]], all_combns_ls[[3]])
names(all_combns_ls) <- 
  all_combns_ls %>% 
  map(function(x){ # browser()
    grp <- groups[x]
    str_glue('{grp[1]} vs {grp[2]}')
  }) %>% 
  unlist()

Myeloid_subtype_ls <- Myeloid$annotation_full %>% unique()
names(Myeloid_subtype_ls) <- Myeloid_subtype_ls

Myeloid_markers_ls <- 
  Myeloid_subtype_ls %>% 
  map(function(Myeloid_subtype){ # browser()
    marker_ls <- 
      all_combns_ls %>% 
      map(function(x){  
        grp <- groups[x]
        FindMarkers(
          subset(JoinLayers(Myeloid), annotation_full == Myeloid_subtype),
          # Myeloid, 
          group.by = 'group',
          test.use = "wilcox", 
          ident.1 = grp[1],
          ident.2 = grp[2],
          logfc.threshold = 0.25,
          min.pct = 0.1,
          only.pos = TRUE
        ) %>% 
          rownames_to_column('gene')
      })
    
    marker_ls
  })

Myeloid_markers_ls %>% 
  iwalk(function(Myeloid_marker, Myeloid_subtype){
    openxlsx::write.xlsx(Myeloid_marker, file.path(Myeloid_out_dir, str_glue('{Myeloid_subtype}_diff_gene.xlsx')))
  })



library(EnhancedVolcano)

volcano_ls <- 
  Myeloid_markers_ls %>%
  imap(function(Myeloid_markers, subtype){
    Myeloid_markers %>% 
      imap(function(Myeloid_marker, group){
        p <- EnhancedVolcano(Myeloid_marker,
                             lab = Myeloid_marker$gene,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = str_glue('{subtype}: {group} volcano plot'),
                             pCutoff = 10e-2,
                             FCcutoff = 0,
                             pointSize = 3.0,
                             labSize = 6.0)
        p
      })
  })

volcano_ls %>% 
  iwalk(function(x, subtype){
    x %>% 
      imap(function(valcano, group){
        ggsave(
          valcano,
          filename = file.path(Myeloid_out_dir, str_glue("volcano/{subtype}_{group} diffgene valcano.pdf")),
          height = 6,
          width = 8
        )
      })
  })




### enrichment ----------------------------
Myeloid_de_genes_for_enrichmenet_ls <- 
  Myeloid_markers_ls %>% 
  map(function(Myeloid_markers){
    Myeloid_markers %>% 
      map(function(Myeloid_marker){
        Myeloid_marker %>% 
          dplyr::filter(p_val_adj < 0.05) %>%
          dplyr::filter(pct.1 - pct.2 > 0.1) %>%
          dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
          nest() %>%
          deframe() %>%
          map("gene")
      })
  })


enrich_ls <- 
  Myeloid_de_genes_for_enrichmenet_ls %>% 
  map(function(Myeloid_de_genes_for_enrichmenet_subtype){ 
    Myeloid_de_genes_for_enrichmenet_subtype %>% 
      map(function(Myeloid_de_genes_for_enrichmenet_group){ # browser()
        Myeloid_de_genes_for_enrichmenet_group %>% 
          map(function(Myeloid_de_genes_for_enrichmenet){
            clusterProfiler::enrichGO(
              Myeloid_de_genes_for_enrichmenet,
              OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
              ont = 'ALL',
              pAdjustMethod = 'BH',
              pvalueCutoff = 0.05, 
              qvalueCutoff = 0.05,
              keyType = 'SYMBOL'
            )
          })
      })
  })


er_output <- file.path(Myeloid_out_dir, 'enrichment') # dir.create(er_output)

# expot enrichment table
enrich_ls %>% 
  iwalk(function(enrich, subtype){ # browser()
    dir.create(file.path(er_output, subtype))
    
    enrich %>% 
      iwalk(function(enrich_result, name){ # browser()
        enrich_result %>% 
          discard(function(x) length(x) < 1) %>% 
          walk(function(x){
            write.csv(
              x@result,
              file = file.path(er_output, subtype, str_glue("{name} Go enrichment table.csv")),
              row.names = FALSE
            )
          })
      })
  })

# dotplot
enrich_dot_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })

enrich_bar_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })


# export dotplot
enrich_dot_ls %>% 
  iwalk(function(enrich_dot, subtype){
    enrich_dot %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} go dotplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


enrich_bar_ls %>% 
  iwalk(function(enrich_bar, subtype){
    enrich_bar %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} go barplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


Myeloid_de_genes_id_for_enrichmenet_ls <-
  Myeloid_de_genes_for_enrichmenet_ls %>% 
  map(function(Myeloid_de_genes_for_enrichmenet_subtype){
    Myeloid_de_genes_for_enrichmenet_subtype %>% 
      map(function(Myeloid_de_genes_for_enrichmenet_group){
        Myeloid_de_genes_for_enrichmenet_group %>% 
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
      })
  })



enrich_ls <- 
  Myeloid_de_genes_id_for_enrichmenet_ls %>% 
  map(function(Myeloid_de_genes_id_for_enrichmenet_subtype){
    Myeloid_de_genes_id_for_enrichmenet_subtype %>% 
      map(function(Myeloid_de_genes_id_for_enrichmenet_group){
        Myeloid_de_genes_id_for_enrichmenet_group %>% 
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
              use_internal_data = FALSE
            )
          })
      })
    
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>% 
  iwalk(function(enrich, subtype){ # browser()
    dir.create(file.path(er_output, subtype))
    
    enrich %>% 
      iwalk(function(enrich_result, name){ # browser()
        enrich_result %>% 
          discard(function(x) length(x) < 1) %>% 
          walk(function(x){
            write.csv(
              x@result,
              file = file.path(er_output, subtype, str_glue("{name} KEGG enrichment table.csv")),
              row.names = FALSE
            )
          })
      })
  })


# dotplot
enrich_dot_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })

enrich_bar_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })


# export dotplot
enrich_dot_ls %>% 
  iwalk(function(enrich_dot, subtype){
    enrich_dot %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          purrr::discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} kegg dotplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


enrich_bar_ls %>% 
  iwalk(function(enrich_bar, subtype){
    enrich_bar %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          purrr::discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} kegg barplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })




### subtype group freq diff -----------------------------------------------------
Myeloid_sub_boxplot_out_dir <- file.path(Myeloid_out_dir, 'Percentage_diff_boxplot')
dir.create(Myeloid_sub_boxplot_out_dir)
Myeloid$annotation_full %>% table()
Myeloid_count <- 
  Myeloid@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(sample, group) %>% 
  summarise(count_Myeloid = n())
Myeloid_sub_meta_ls <- 
  Myeloid@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(annotation_full) %>% 
  nest() %>% 
  deframe() %>% 
  map(function(Myeloid_meta){
    Myeloid_meta %>% 
      group_by(sample, group) %>% 
      summarise(count_Myeloid_sub = n()) %>% 
      left_join(Myeloid_count, by = c('sample', 'group')) %>% 
      mutate(Percentage = count_Myeloid_sub/count_Myeloid*100) %>% 
      group_by(group) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat){
        qValue <- 
          dat[['Percentage']] %>% 
          quantile()
        
        dat %>% 
          mutate(Q2 = qValue[2], 
                 Q4 = qValue[4])
      }) %>% 
      enframe(name = 'group') %>% 
      unnest() %>% 
      mutate(group = factor(group, levels = c('N', 'CP', 'T')))
  })

Myeloid_sub_meta_ls %>% 
  imap(function(Myeloid_sub_meta, Myeloid_sub){
    print(Myeloid_sub)
    print(shapiro.test(Myeloid_sub_meta$Percentage))
    print(bartlett.test(Percentage~group, data = Myeloid_sub_meta))
  })
saveRDS(Myeloid_sub_meta_ls, file.path(Myeloid_out_dir, 'Myeloid_sub_meta_ls.rds'))

Myeloid_sub_meta_tTest_ls <- Myeloid_sub_meta_ls[c("M1 macrophage", "cDC2", "Proliferating macrophage", "LAMP3 DC", "M2 macrophage", "IFN-responsive macrophage")]
Myeloid_sub_meta_wilcox_ls <- Myeloid_sub_meta_ls[c("Neutropil")]

box_p_tTest_ls <- 
  Myeloid_sub_meta_tTest_ls %>% 
  imap(function(Myeloid_sub_meta, Myeloid_sub){
    Myeloid_sub_meta %>% 
      ggplot(aes(x = group, y = Percentage, fill = group)) +
      geom_errorbar(aes(ymin = Q2, ymax = Q4), color = 'darkgrey') + 
      geom_jitter(
        width = 0.2, 
        height = 0.1, 
        alpha = 1, 
        size = 2, 
        aes(color = group)
      ) +
      ggpubr::geom_pwc(
        tip.length = 0, 
        hide.ns = F, 
        label = 'p', 
        # method.args = list(paired = F, exact = T), 
        # p.adjust.method = "BH", 
        method = 't_test',
        # method = 'wilcox_test',
        # group.by = "x.var"
        group.by = "legend.var"
      ) +  
      scale_color_manual(values = c(
        'N' = 'lightgreen', 
        'CP' = 'lightblue', 
        'T' = 'pink')
      ) +
      ggtitle(Myeloid_sub) + 
      theme_bw()  
  })

box_p_tTest_ls %>% 
  iwalk(function(box_p, Myeloid_sub){ # browser()
    pdf(file.path(Myeloid_sub_boxplot_out_dir, str_glue('{Myeloid_sub} by group.pdf')), width = 6, height = 6)
    print(box_p)
    dev.off()
  })


box_p_wilcox_ls <- 
  Myeloid_sub_meta_wilcox_ls %>% 
  imap(function(Myeloid_sub_meta, Myeloid_sub){
    Myeloid_sub_meta %>% 
      ggplot(aes(x = group, y = Percentage, fill = group)) +
      geom_errorbar(aes(ymin = Q2, ymax = Q4), color = 'darkgrey') + 
      geom_jitter(
        width = 0.2, 
        height = 0.1, 
        alpha = 1, 
        size = 2, 
        aes(color = group)
      ) +
      ggpubr::geom_pwc(
        tip.length = 0, 
        hide.ns = F, 
        label = 'p', 
        # method.args = list(paired = F, exact = T), 
        # p.adjust.method = "BH", 
        # method = 't_test',
        method = 'wilcox_test',
        # group.by = "x.var"
        group.by = "legend.var"
      ) +  
      scale_color_manual(values = c(
        'N' = 'lightgreen', 
        'CP' = 'lightblue', 
        'T' = 'pink')
      ) +
      ggtitle(Myeloid_sub) + 
      theme_bw()  
  })

box_p_wilcox_ls %>% 
  iwalk(function(box_p, Myeloid_sub){ # browser()
    pdf(file.path(Myeloid_sub_boxplot_out_dir, str_glue('{Myeloid_sub} by group.pdf')), width = 6, height = 6)
    print(box_p)
    dev.off()
  })




## Fibroblast -----------------------------------------------------------------------
Fibroblast <- subset(all_data2, annotation_full %in% c(
  'dCAF', 
  'iCAF', 
  'mCAF', 
  'tCAF', 
  'Normal fibroblast'
))




### cell props barplot  --------
Fibroblast_meta <- Fibroblast[[]]

p_col_proj <-
  Fibroblast_meta %>%
  group_by(sample, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation_full)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Fibro_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_group <-
  Fibroblast_meta %>%
  group_by(group, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_group,
  filename = file.path(Fibro_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)




### volcano-------
Idents(Fibroblast) <- 'group'
groups <- c('T', 'CP', 'N')
all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
all_combns_ls <- list(all_combns_ls[[1]], all_combns_ls[[3]])
names(all_combns_ls) <- 
  all_combns_ls %>% 
  map(function(x){ # browser()
    grp <- groups[x]
    str_glue('{grp[1]} vs {grp[2]}')
  }) %>% 
  unlist()

Fibroblast_subtype_ls <- Fibroblast$annotation_full %>% unique()
names(Fibroblast_subtype_ls) <- Fibroblast_subtype_ls

Fibroblast_markers_ls <- 
  Fibroblast_subtype_ls %>% 
  map(function(Fibroblast_subtype){ # browser()
    tryCatch({
      marker_ls <- 
        all_combns_ls %>% 
        map(function(x){  
          grp <- groups[x]
          FindMarkers(
            subset(JoinLayers(Fibroblast), annotation_full == Fibroblast_subtype),
            # Fibroblast, 
            group.by = 'group',
            test.use = "wilcox", 
            ident.1 = grp[1],
            ident.2 = grp[2],
            logfc.threshold = 0.25,
            min.pct = 0.1,
            only.pos = TRUE
          ) %>% 
            rownames_to_column('gene')
        })
      
      marker_ls
    }, 
    error = function(x){
      print(x$message)
      NULL
    }
    )  
  }) 

Fibroblast_markers_ls %>% 
  iwalk(function(Fibroblast_marker, Fibroblast_subtype){
    openxlsx::write.xlsx(Fibroblast_marker, file.path(Fibro_out_dir, str_glue('{Fibroblast_subtype}_diff_gene.xlsx')))
  })



library(EnhancedVolcano)

volcano_ls <- 
  Fibroblast_markers_ls %>%
  imap(function(Fibroblast_markers, subtype){
    Fibroblast_markers %>% 
      imap(function(Fibroblast_marker, group){
        p <- EnhancedVolcano(Fibroblast_marker,
                             lab = Fibroblast_marker$gene,
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = str_glue('{subtype}: {group} volcano plot'),
                             pCutoff = 10e-2,
                             FCcutoff = 0,
                             pointSize = 3.0,
                             labSize = 6.0)
        p
      })
  })

volcano_ls %>% 
  iwalk(function(x, subtype){
    x %>% 
      imap(function(valcano, group){
        ggsave(
          valcano,
          filename = file.path(Fibro_out_dir, str_glue("volcano/{subtype}_{group} diffgene valcano.pdf")),
          height = 6,
          width = 8
        )
      })
  })




### enrichment ----------------------------
Fibroblast_de_genes_for_enrichmenet_ls <- 
  Fibroblast_markers_ls %>% 
  map(function(Fibroblast_markers){
    Fibroblast_markers %>% 
      map(function(Fibroblast_marker){
        Fibroblast_marker %>% 
          dplyr::filter(p_val_adj < 0.05) %>%
          dplyr::filter(pct.1 - pct.2 > 0.1) %>%
          dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
          nest() %>%
          deframe() %>%
          map("gene")
      })
  })


enrich_ls <- 
  Fibroblast_de_genes_for_enrichmenet_ls %>% 
  map(function(Fibroblast_de_genes_for_enrichmenet_subtype){ 
    Fibroblast_de_genes_for_enrichmenet_subtype %>% 
      map(function(Fibroblast_de_genes_for_enrichmenet_group){ # browser()
        Fibroblast_de_genes_for_enrichmenet_group %>% 
          map(function(Fibroblast_de_genes_for_enrichmenet){
            clusterProfiler::enrichGO(
              Fibroblast_de_genes_for_enrichmenet,
              OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
              ont = 'ALL',
              pAdjustMethod = 'BH',
              pvalueCutoff = 0.05, 
              qvalueCutoff = 0.05,
              keyType = 'SYMBOL'
            )
          })
      })
  })


er_output <- file.path(Fibro_out_dir, 'enrichment') # dir.create(er_output)

# expot enrichment table
enrich_ls %>% 
  iwalk(function(enrich, subtype){ # browser()
    dir.create(file.path(er_output, subtype))
    
    enrich %>% 
      iwalk(function(enrich_result, name){ # browser()
        enrich_result %>% 
          discard(function(x) length(x) < 1) %>% 
          walk(function(x){
            write.csv(
              x@result,
              file = file.path(er_output, subtype, str_glue("{name} Go enrichment table.csv")),
              row.names = FALSE
            )
          })
      })
  })

# dotplot
enrich_dot_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })

enrich_bar_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })


# export dotplot
enrich_dot_ls %>% 
  iwalk(function(enrich_dot, subtype){
    enrich_dot %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} go dotplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


enrich_bar_ls %>% 
  iwalk(function(enrich_bar, subtype){
    enrich_bar %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} go barplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


Fibroblast_de_genes_id_for_enrichmenet_ls <-
  Fibroblast_de_genes_for_enrichmenet_ls %>% 
  map(function(Fibroblast_de_genes_for_enrichmenet_subtype){
    Fibroblast_de_genes_for_enrichmenet_subtype %>% 
      map(function(Fibroblast_de_genes_for_enrichmenet_group){
        Fibroblast_de_genes_for_enrichmenet_group %>% 
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
      })
  })



enrich_ls <- 
  Fibroblast_de_genes_id_for_enrichmenet_ls %>% 
  map(function(Fibroblast_de_genes_id_for_enrichmenet_subtype){
    Fibroblast_de_genes_id_for_enrichmenet_subtype %>% 
      map(function(Fibroblast_de_genes_id_for_enrichmenet_group){
        Fibroblast_de_genes_id_for_enrichmenet_group %>% 
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
              use_internal_data = FALSE
            )
          })
      })
    
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>% 
  iwalk(function(enrich, subtype){ # browser()
    dir.create(file.path(er_output, subtype))
    
    enrich %>% 
      iwalk(function(enrich_result, name){ # browser()
        enrich_result %>% 
          discard(function(x) length(x) < 1) %>% 
          walk(function(x){
            write.csv(
              x@result,
              file = file.path(er_output, subtype, str_glue("{name} KEGG enrichment table.csv")),
              row.names = FALSE
            )
          })
      })
  })


# dotplot
enrich_dot_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })

enrich_bar_ls <- 
  enrich_ls %>% 
  map(function(enrich){
    enrich %>% 
      map(function(enrich_result){
        enrich_result %>% 
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
      })
  })


# export dotplot
enrich_dot_ls %>% 
  iwalk(function(enrich_dot, subtype){
    enrich_dot %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          purrr::discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} kegg dotplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })


enrich_bar_ls %>% 
  iwalk(function(enrich_bar, subtype){
    enrich_bar %>% 
      iwalk(function(enrich, group){
        enrich %>% 
          purrr::discard(function(x) length(x) < 1) %>%
          iwalk(function(x, name){
            # cluster_name <- paste0("cluster_", name)
            # dir.create(file.path(er_output, cluster_name))
            
            ggsave(
              x,
              filename = file.path(er_output, subtype, str_glue("{group} kegg barplot.pdf")),
              height = 6,
              width = 8
            )
          })
      })
  })




### subtype group freq diff -----------------------------------------------------
Fibroblast_sub_boxplot_out_dir <- file.path(Fibro_out_dir, 'Percentage_diff_boxplot')
dir.create(Fibroblast_sub_boxplot_out_dir)
Fibroblast$annotation_full %>% table()
Fibroblast_count <- 
  Fibroblast@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(sample, group) %>% 
  summarise(count_Fibroblast = n())
Fibroblast_sub_meta_ls <- 
  Fibroblast@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  group_by(annotation_full) %>% 
  nest() %>% 
  deframe() %>% 
  map(function(Fibroblast_meta){
    Fibroblast_meta %>% 
      group_by(sample, group) %>% 
      summarise(count_Fibroblast_sub = n()) %>% 
      left_join(Fibroblast_count, by = c('sample', 'group')) %>% 
      mutate(Percentage = count_Fibroblast_sub/count_Fibroblast*100) %>% 
      group_by(group) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat){
        qValue <- 
          dat[['Percentage']] %>% 
          quantile()
        
        dat %>% 
          mutate(Q2 = qValue[2], 
                 Q4 = qValue[4])
      }) %>% 
      enframe(name = 'group') %>% 
      unnest() %>% 
      mutate(group = factor(group, levels = c('N', 'CP', 'T')))
  })

Fibroblast_sub_meta_ls %>% 
  imap(function(Fibroblast_sub_meta, Fibroblast_sub){ 
    tryCatch({
      print(Fibroblast_sub)
      print(shapiro.test(Fibroblast_sub_meta$Percentage))
      print(bartlett.test(Percentage~group, data = Fibroblast_sub_meta))
    }, 
    error = function(x){
      print(x$message)
      NULL
    })
  })
saveRDS(Fibroblast_sub_meta_ls, file.path(Fibro_out_dir, 'Fibroblast_sub_meta_ls.rds'))

Fibroblast_sub_meta_tTest_ls <- Fibroblast_sub_meta_ls[c("iCAF")]
Fibroblast_sub_meta_wilcox_ls <- Fibroblast_sub_meta_ls[c("mCAF", "tCAF", "Normal fibroblast", "dCAF")]

box_p_tTest_ls <- 
  Fibroblast_sub_meta_tTest_ls %>% 
  imap(function(Fibroblast_sub_meta, Fibroblast_sub){
    Fibroblast_sub_meta %>% 
      ggplot(aes(x = group, y = Percentage, fill = group)) +
      geom_errorbar(aes(ymin = Q2, ymax = Q4), color = 'darkgrey') + 
      geom_jitter(
        width = 0.2, 
        height = 0.1, 
        alpha = 1, 
        size = 2, 
        aes(color = group)
      ) +
      ggpubr::geom_pwc(
        tip.length = 0, 
        hide.ns = F, 
        label = 'p', 
        # method.args = list(paired = F, exact = T), 
        # p.adjust.method = "BH", 
        method = 't_test',
        # method = 'wilcox_test',
        # group.by = "x.var"
        group.by = "legend.var"
      ) +  
      scale_color_manual(values = c(
        'N' = 'lightgreen', 
        'CP' = 'lightblue', 
        'T' = 'pink')
      ) +
      ggtitle(Fibroblast_sub) + 
      theme_bw()  
  })

box_p_tTest_ls %>% 
  iwalk(function(box_p, Fibroblast_sub){ # browser()
    pdf(file.path(Fibroblast_sub_boxplot_out_dir, str_glue('{Fibroblast_sub} by group.pdf')), width = 6, height = 6)
    print(box_p)
    dev.off()
  })


box_p_wilcox_ls <- 
  Fibroblast_sub_meta_wilcox_ls %>% 
  imap(function(Fibroblast_sub_meta, Fibroblast_sub){
    Fibroblast_sub_meta %>% 
      ggplot(aes(x = group, y = Percentage, fill = group)) +
      geom_errorbar(aes(ymin = Q2, ymax = Q4), color = 'darkgrey') + 
      geom_jitter(
        width = 0.2, 
        height = 0.1, 
        alpha = 1, 
        size = 2, 
        aes(color = group)
      ) +
      ggpubr::geom_pwc(
        tip.length = 0, 
        hide.ns = F, 
        label = 'p', 
        # method.args = list(paired = F, exact = T), 
        # p.adjust.method = "BH", 
        # method = 't_test',
        method = 'wilcox_test',
        # group.by = "x.var"
        group.by = "legend.var"
      ) +  
      scale_color_manual(values = c(
        'N' = 'lightgreen', 
        'CP' = 'lightblue', 
        'T' = 'pink')
      ) +
      ggtitle(Fibroblast_sub) + 
      theme_bw()  
  })

box_p_wilcox_ls %>% 
  iwalk(function(box_p, Fibroblast_sub){ # browser()
    pdf(file.path(Fibroblast_sub_boxplot_out_dir, str_glue('{Fibroblast_sub} by group.pdf')), width = 6, height = 6)
    print(box_p)
    dev.off()
  })






# all data cell prop ------------------------------------------------------------
all_data_out_dir <- '2. dimension reduction and cluster'
all_data2_meta <- all_data2[[]]

all_data2_meta$sample <- factor(all_data2_meta$sample, 
                                levels=c(all_data2_meta$sample %>% unique())
)
p_col_proj <-
  all_data2_meta %>%
  group_by(sample, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(all_data_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_group <-
  all_data2_meta %>%
  group_by(group, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_group,
  filename = file.path(all_data_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)




# Annotation info output --------------------------------------------------------
dir.create('anno_info')
all_data2$annotation %>% table()
cellType_anno_ls <- 
  all_data2@meta.data[c('project', 'annotation')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

cellType_anno_ls %>% 
  iwalk(function(cellType, project){
    write.csv(cellType, str_glue('anno_info/cell type {project}.csv'), row.names = F)
  })


all_data2$annotation_full %>% table()
subType_anno_ls <- 
  all_data2@meta.data[c('project', 'annotation_full')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*_', ''), 
    group = annotation_full
  ) %>% 
  dplyr::select(cell_id, project, group) %>% 
  group_by(project) %>%
  nest() %>% 
  deframe()

subType_anno_ls %>% 
  iwalk(function(subType, project){
    write.csv(subType, str_glue('anno_info/sub cell type {project}.csv'), row.names = F)
  })

