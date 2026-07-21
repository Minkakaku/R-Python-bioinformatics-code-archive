library(tidyverse)
library(Seurat)

# define default-------------
out_dir <- "CN1-REST_tumor_vs_CN7-REST_tumor_DEGs_enrichment" 
dir.create(out_dir)

all_data <- readRDS("all_data_annotation_filtered.rds")
all_data$annotation_full %>% table()

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') %>% 
  left_join(all_data_cn_meta[c('cell_id', 'CN')], by = 'cell_id') %>% 
  column_to_rownames('cell_id')
all_data_meta %>% head()
all_data <- AddMetaData(all_data, all_data_meta['CN'])
all_data@meta.data %>% head()




# Tumor DEGs ---------------
Tumor <- subset(JoinLayers(all_data), annotation_full == 'REST_tumor cell' & CN %in% c('CN1', 'CN7'))
Tumor$annotation_full %>% table()
Tumor$CN %>% table()
Idents(Tumor) <- 'CN'

## vocanol -------------------
Tumor_markers <- 
  FindAllMarkers(
    JoinLayers(Tumor),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Tumor_markers, file.path(out_dir, 'Tumor_diff_gene.csv'))

library(EnhancedVolcano)
type<-unique(Tumor_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Tumor_markers[Tumor_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })

volcano_ls %>% 
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


Tumor_markers_major <- list(Tumor = Tumor_markers)
# scatter
p_volcanoplot_major_ls <- 
  Tumor_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>% 
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 6,
      width = 10
    )
  })


## enrrichment ---------------------------------------
Tumor_de_genes_for_enrichmenet <-  
  Tumor_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 - pct.2 > 0.05) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

Tumor_er_output <- file.path(out_dir, 'enrichment') 
dir.create(Tumor_er_output)

enrich_ls <-
  Tumor_de_genes_for_enrichmenet %>%
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


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Tumor_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Tumor_er_output, name, "Go enriment table.csv"),
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
      filename = file.path(Tumor_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(Tumor_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Tumor_de_genes_id_for_enrichmenet <-
  Tumor_de_genes_for_enrichmenet %>%
  map(function(x){
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
  Tumor_de_genes_id_for_enrichmenet %>%
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

# saveRDS(Tumor_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Tumor_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Tumor_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(Tumor_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(Tumor_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(Tumor, file.path(out_dir, 'all_Tumor.rds'))
