# /data/rawdata/rawdata-sp/FW2024-389/FW2024-389-01/客户反馈数据
# sc <- readRDS('BP_AD_Pso-v2.rds')
sc <- readRDS('BP_AD_Pso_20240930.rds')

pdf('dotplot singlecell.pdf', height = 6, width = 6)
DotPlot(
  sc,
  features = c("CD68", "CD86", "CD80", "MARCO")
)
dev.off()


sc$orig.ident %>% unique


## find markers -------
# Find all markers using FindConservedMarkers
FindAllConservedMarkers <-
  function(obj, 
           ident.1 = NULL, 
           grouping.var,
           ...) {
    
    stopifnot(is(obj, "Seurat"))
    
    all_ident <- unique(Idents(obj))
    stopifnot(length(all_ident) > 1)
    
    
    if(is.null(ident.1)){
      message("ident.1 is null, will find all markers of each ident in all idents.")
    }else{
      all_ident <-  ident.1
    }
    
    all_ident %>%
      set_names(., .) %>%
      map(function(ident){
        FindConservedMarkers(
          object = obj,
          ident.1 = ident,
          grouping.var = grouping.var,
          ...
        )
      })
  }

# Idents(sc_fibro) <- "RNA_snn_res.0.2"
# DefaultAssay(sc_fibro) <- "RNA"

# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_fibro),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
all_markers <-
  FindAllConservedMarkers(
    sc,
    grouping.var = "ident",
    test.use = "wilcox",
    # logfc.threshold = 0.25,
    # min.pct = 0.1,
    only.pos = TRUE
  ) %>%
  map(rownames_to_column, "gene") %>%
  enframe("cluster") %>%
  unnest(value) %>%
  rowwise() %>%
  mutate(
    avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
    pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
  ) %>%
  mutate(p_val_adj  = minimump_p_val)



all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # arrange(p_val_adj) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

avg_expr <- AverageExpression(sc)[[1]]

avg_expr[all_de_genes %>% unique, -(7:9)] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )