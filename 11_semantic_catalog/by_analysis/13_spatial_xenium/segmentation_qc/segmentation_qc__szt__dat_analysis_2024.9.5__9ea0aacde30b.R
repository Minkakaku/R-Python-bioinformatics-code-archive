library(Seurat)
library(tidyverse)
 
# define function -------------
#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}

##### nestSingleList
nestSingleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  name <- names(x)
  id_ <- seq_along(name)
  stopifnot(all(str_detect(name, sep)))
  name_1 <- str_split(name, pattern = sep) %>% map_chr(~ .x[1])
  name_2 <- str_split(name, pattern = sep) %>% map(~ .x[2:length(.x)]) %>% map_chr(str_c, collapse = sep)
  
  res_ls <- list()
  for(i in id_){    
    n1 <- name_1[i]
    n2 <- name_2[i]
    res <- x[[i]]
    res_ls[[n1]][[n2]] <- res
  }
  res_ls
}


out_dir <- 'dat_analysis_cycle3'
dir.create(out_dir)


# project to drop out -----
project_to_remove <- c(
  '09975N5',
  '09925E1',
  '77054'
)


# read data --------------
all_data <- readRDS("dat_analysis_cycle3/all_data_cycle3.rds") 
all_data_cn_meta <- 
  readRDS(file.path(out_dir, "all_data_cn_meta.rds")) %>%
  dplyr::filter(! project %in% project_to_remove)


all_data_cn <- 
  all_data %>%
  AddMetaData(all_data_cn_meta[c('CN', 'CN2', 'patient', 'group3', 'group4')])

all_data_cn <- all_data_cn %>% JoinLayers()

all_data_cn <- all_data_cn %>% subset(cells = rownames(all_data_cn_meta))

# export rds
# 已经导出数据，为免误操作，注释saveRDS代码
# saveRDS(all_data_cn, file.path('data_for_analysis', "all_data_cn_remove_samples.rds"))
# saveRDS(all_data_cn[[]], file.path('data_for_analysis', "all_data_meta_cn_remove_samples.rds"))


# !import
# read data2 -- replace read data --------------
all_data_cn <- readRDS(file.path('data_for_analysis', "all_data_cn_remove_samples.rds"))
all_data_cn_meta <- readRDS(file.path('data_for_analysis', "all_data_meta_cn_remove_samples.rds"))


# add sub group ---------------
sub_group <- c(
  A_00037G2 = 'A1',
  A_00037G3 = 'A1',
  A_06685N2 = 'A1',
  A_06685N3 = 'A1',
  A_09925H4 = 'A1',
  A_09975J1 = 'A1',
  A_10858   = 'A2',
  A_13638   = 'A2',
  A_47983A3 = 'A2',
  A_73479   = 'A1',
  B_06938F1 = 'B4',
  B_06938I6 = 'B1',
  B_12876F7 = 'B4',
  B_12876F9 = 'B1',
  B_15201H1 = 'B2',
  B_28519G1 = 'B2',
  B_28519M4 = 'B1',
  B_43158G1 = 'B4',
  B_44994A1 = 'B3',
  B_45827   = 'B3',
  B_57986I1 = 'B5',
  B_57986L1 = 'B5',
  B_58198   = 'B2',
  B_67718F4 = 'B1',
  B_67718G1 = 'B2',
  B_71402A6 = 'B2',
  B_72696   = 'B2',
  C_09791H1 = 'C',
  C_47634I1 = 'C'
)
names(sub_group) <- sub_group %>% names() %>% str_remove("^[ABC]_")
sub_group_df <- 
  sub_group %>%
  enframe('project','sub_group')

add_sub_group_df <-
  all_data_cn[['project']] %>%
  rownames_to_column('cell_id') %>%
  left_join(sub_group_df, by = "project") %>%
  column_to_rownames("cell_id")

all_data_cn <-
  all_data_cn %>%
  AddMetaData(add_sub_group_df) 

# cell props ----------------
cell_props_df <-
  all_data_cn[[]] %>%
  mutate(project = str_c(group, project, sep = "_")) %>%
  dplyr::select(project, CN2, annotation2) %>%
  group_by(project, CN2, annotation2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq = count/sum(count) * 100)

p_cell_props_cn_cell_project <-
  cell_props_df %>%
  separate(project, into = c("group", "project"), sep = "_") %>%
  group_by(group, CN2, annotation2) %>%
  summarise(freq = mean(freq)) %>%
  group_by(group, CN2)  %>%
  mutate(freq = freq/sum(freq)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation2)) +
  facet_wrap(~ CN2) +
  geom_col() +
  ggsci::scale_fill_igv()

ggsave(
  p_cell_props_cn_cell_project,
  filename = file.path(out_dir, "cell props barplot_group, cell and CN.pdf"),
  height = 18,
  width  = 15
)


p_cell_props_cn_project <-
  cell_props_df %>%
  ggplot(aes(x = CN2, y = freq, fill = annotation2)) +
  geom_col() +
  facet_wrap(~ project, scales = 'free_y') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  ggsci::scale_fill_igv()

ggsave(
  p_cell_props_cn_project,
  filename = file.path(out_dir, "cell props barplot_project and CN.pdf"),
  height = 12,
  width  = 15
)

# cn per sample
cell_props_df_only_cn <-
  all_data_cn[[]] %>%
  mutate(project = str_c(group, project, sep = "_")) %>%
  dplyr::select(project, CN2) %>%
  group_by(project, CN2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq = count/sum(count) * 100)

p_cell_props_df_only_cn <-
  cell_props_df_only_cn %>%
  ggplot(aes(x = project, y = freq, fill = CN2)) +
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  ggsci::scale_fill_igv()

ggsave(
  p_cell_props_df_only_cn,
  filename = file.path(out_dir, "cell props barplot_only CN.pdf"),
  height = 8,
  width  = 12
)


# All cells in all CNs, de analysis between A and B --------

all_de_ls2 <-
  all_data_cn[[]] %>%
  rownames_to_column('cell_id') %>%
  dplyr::filter(CN2 %in% c(
    "Tumor boundary",
    "Tumor core",
    "CD8T + pan-immune",
    "CD4T + pan-immune"
  )) %>%
  group_by(CN2) %>%
  nest() %>%
  deframe() %>%
  map(function(meta){
    meta %>%
      group_by(annotation2) %>%
      nest() %>%
      deframe() %>%
      map(function(sub_meta){
        cell.1 <- sub_meta$cell_id[sub_meta$group == "A"]
        cell.2 <- sub_meta$cell_id[sub_meta$group == "B"]
        
        # FindMarkers
        de <- tryCatch(
          
            FindMarkers(
              all_data_cn[['RNA']],
              cells.1 = cell.1,
              cells.2 = cell.2,
              fc.slot = "data",
              pseudocount.use = 1,
              base = 2,
              logfc.threshold = 0.1,
              min.pct = 0.1, # 0.01
              only.pos = FALSE
            ),
          error = function(e) return(NULL)
        )
        
        de
      })
  })
  
# saveRDS(all_de_ls, file = file.path(out_dir, "all_de_ls.rds"))

# export de table
all_de_ls2 %>%
  iwalk(function(dat, name){
    dat %>%
      map(as.data.frame) %>%
      map(rownames_to_column, "Gene") %>%
      enframe("Cell", "value") %>%
      unnest(value) %>%
      openxlsx::write.xlsx(
        file = file.path(out_dir, str_glue("de table - CN - {name}.xlsx")),
        asTable = TRUE
      )
  })

# enrichmn analysis -------------
cells_to_analysis <- c(
  "Epithelial",
  "CD4T",
  "CD8T",
  "Epi0_CD4T",
  "Epi0_CD83_T",
  "Epi0_T"
)

all_enrich2 <-
  all_de_ls2 %>%
  map(~ .x[names(.x) %in% cells_to_analysis]) %>%
  flatDoubleList() %>%
  map(function(de){ # browser()
    de <-
      de %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(abs(avg_log2FC) > 1) %>%
      dplyr::filter(pct.1 > 0.1)
    gene_up   <- rownames(de)[de$avg_log2FC > 0]
    gene_down <- rownames(de)[de$avg_log2FC < 0]
    
    # id transfer
    gene_up_id <- gene_up %>% clusterProfiler::bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    gene_down_id <- gene_down %>% clusterProfiler::bitr(fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
    
    # enrichemnt 
    ernich_up <-
      clusterProfiler::enrichGO(
        gene_up_id,
        ont = "BP",
        OrgDb = "org.Hs.eg.db",
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        qvalueCutoff = 0.2, 
        minGSSize = 10,       
        maxGSSize = 500
      )
    ernich_down <-
      clusterProfiler::enrichGO(
        gene_up_id,
        ont = "BP",
        OrgDb = "org.Hs.eg.db",
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        qvalueCutoff = 0.2, 
        minGSSize = 10,       
        maxGSSize = 500
      )
    list(
      UP = ernich_up,
      Down = ernich_down
    )
  })


# saveRDS(all_enrich2, file = file.path(out_dir, "all_enrich2.rds"))


all_enrich_dot2 <-
  all_enrich2 %>% 
  map(function(enrich){
    up_dot <- NULL
    down_dot <- NULL
    
    if(!is.null(enrich$UP)){
      up_dot <- 
        enrich$UP %>%
        enrichplot::dotplot(showCategory = 15)
    }
    
    if(!is.null(enrich$Down)){
      down_dot <- 
        enrich$Down %>%
        enrichplot::dotplot(showCategory = 15)
    }
    
    # return
    list(
      UP = up_dot,
      Down = down_dot
    )
  })

# export
pdf(
  file.path(out_dir, "enrichment analysis CN.pdf"),
  height = 6,
  width  = 15
)
all_enrich_dot2 %>%
  iwalk(function(enrich, name){
    null_idx <- map_lgl(enrich, is.null)
    if(all(null_idx)) return()
    if(any(null_idx)) {
      print(enrich[!null_idx])
      return()
    }
    # browser()
    p <- 
      enrich %>%
      imap(~ .x + ggtitle(.y)) %>%
      patchwork::wrap_plots(nrow = 1) +
      patchwork::plot_annotation(title = name)
    print(p)
  })
dev.off()


# de gene of CN 
all_de_gene_ls2 <-
  all_de_ls2 %>%
  map(~ .x[names(.x) %in% cells_to_analysis]) %>%
  flatDoubleList() %>%
  map(function(de){ # browser()
    de <-
      de %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(abs(avg_log2FC) > 1) %>%
      dplyr::filter(pct.1 > 0.1)
    gene_up   <- rownames(de)[de$avg_log2FC > 0]
    gene_down <- rownames(de)[de$avg_log2FC < 0]
    
    list(
      Up = gene_up,
      Down = gene_down
    )
  })

# geneset ------
'Chemokine Receptors' <- c(
  "ACKR3",
  'CX3CR1',
  'CCR5',
  'CCRL2',
  'CCR1',
  'CXCR6',
  'ACKR1',
  'ACKR4',
  'CXCR4',
  'CCR2',
  'CCR10',
  'CCR7',
  'CCR9',
  'CXCR5',
  'ACKR2',
  'PITPNM3',
  'CCR4',
  'CCR6',
  'CXCR6',
  'CXCR3',
  'CCR3',
  'CCR8'
)

Chemokines <- c(
  'CCL8',
  'CXCL10',
  'CCL22',
  'CXCL13',
  'CCL20',
  'CCL27',
  'CCL25',
  'CXCL9',
  'CXCL8',
  'CXCL17',
  'CCL21',
  'CCL2',
  'CXCL1',
  'CXCL2',
  'CXCL12',
  'CCL11',
  'CXCL11',
  'CXCL14',
  'XCL2',
  'CCL18',
  'XCL1',
  'CCL4',
  'CCL4L2',
  'CCL5',
  'CXCL16',
  'CCL28',
  'CCL3',
  'CCL23',
  'CXCL3',
  'CX3CL1',
  'CCL16'
)


'Interleukins & interferons' <- c(
  'IL34',
  'IL24',
  'IL16',
  'IL1A',
  'IL18',
  'IL21',
  'IL10',
  'IL22',
  'IL17F',
  'IL1B',
  'IL15',
  'IL2',
  'IL23A',
  'IL12A',
  'IL6',
  'IL11',
  'IFNG',
  'IFNL1',
  'IL4',
  'IL17C',
  'IL1RN',
  'IL5',
  'IL3',
  'IL32',
  'IL19',
  'IL17',
  'IL26',
  'IL13'
)


'Interleukin Receptors' <- c(
  'IL9R',
  'IL2RB',
  'IL18RAP',
  'IL12RB1',
  'IL15RA1',
  'IL21R',
  'IL2RG',
  'IL5RA',
  'CXCR2',
  'CXCR1',
  'IL10RA',
  'IL11RA',
  'IL17RB',
  'IL4R',
  'IL7R',
  'IL13RA1',
  'IL1RAP',
  'IL17RA',
  'IL3RA',
  'IL17RC',
  'IL17RE',
  'IL6R',
  'IL27RA',
  'IL1RL1',
  'IL1R1',
  'IL1R2',
  'IL2RA',
  'IL10RB',
  'IL12RB2',
  'IL20RB',
  'IL18R1',
  'IL6ST'
)


'TNF Superfamily' <- c(
  'TNFSF15',
  'TNFSF13',
  'CD40LG',
  'TNFSF11',
  'LTA',
  'LTB',
  'TNFSF12',
  'TNFSF8',
  'TNFSF13B',
  'TNF',
  'EDA',
  'TNFSF14',
  'FASLG',
  'TNFSF9',
  'TNFSF4',
  'CD70',
  'TNFSF10'
)

'TNF Receptors' <- c(
  'TNFRSF1A',
  'NGFR',
  'TNFRSF9',
  'TNFRSF19',
  'TNFRSF1B',
  'TNFRSF10B',
  'TNFRSF17',
  'TNFRSF12A',
  'RELT',
  'EDAR',
  'TNFRSF10A',
  'TNFRSF10D',
  'LTBR',
  'CD40',
  'EDA2R',
  'TNFRSF11A',
  'TNFRSF21',
  'TNFRSF13B',
  'CD27',
  'TNFRSF14',
  'TNFRSF8',
  'TNFRSF4',
  'FAS',
  'TNFRSF18',
  'TNFRSF25',
  'TNFRSF13C'
)

pathway_ls <- list(
  'Chemokine Receptors'        = `Chemokine Receptors`,
  Chemokines                   = `Chemokines`,
  'Interleukins & interferons' = `Interleukins & interferons`,
  'Interleukin Receptors'      = `Interleukin Receptors`,
  'TNF Superfamily'            = `TNF Superfamily`,
  'TNF Receptors'              = `TNF Receptors`
)

# gene expr
avg_expr <-
  all_data_cn %>%
  AverageExpression(group.by = c('group', "CN2", "annotation2"))

avg_expr_df <- 
  avg_expr$RNA %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "group", values_to = "expr") %>%
  separate(group, into = c("group", "CN", "Cell"), sep = "_")




anno_col <- 
  avg_expr$RNA %>%
  colnames() %>%
  data.frame(group = ., group2 = .) %>%
  separate(group, into = c("group", "CN", "Cell"), sep = "_") %>%
  column_to_rownames("group2")

anno_colors <- list()
anno_colors$CN <- structure(
  ggsci::pal_aaas()(n_distinct(anno_col$CN)),
  names = anno_col$CN %>% unique()
)
anno_colors$Cell <- structure(
  ggsci::pal_igv()(n_distinct(anno_col$Cell)),
  names = anno_col$Cell %>% unique()
)

# avg_expr$RNA %>%
#   as.matrix() %>%
#   as.data.frame() %>%
#   .[rownames(.) %in% pathway_ls$`Chemokine Receptors`, ] %>%
#   
#   pheatmap::pheatmap(
#     scale = 'row',
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     annotation_col = anno_col,
#     show_colnames = FALSE,
#     cluster_cols = FALSE
#   )

pathway_heatmap_ls <-
  pathway_ls %>%
  imap(function(x, path_name){

    dat <- avg_expr$RNA %>%
      as.matrix() %>%
      as.data.frame() %>%
      .[rownames(.) %in% x, ] %>%
      log1p() 
    
    anno_col%>%
      rownames_to_column("id") %>%
      group_by(group) %>%
      nest() %>%
      deframe() %>%
      imap(function(anno_col_sub, group){
        anno_col_sub <- anno_col_sub %>% column_to_rownames("id")
        p <- dat[, rownames(anno_col_sub)] %>%
          pheatmap::pheatmap(
            scale = 'row',
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
            annotation_col = anno_col,
            annotation_colors = anno_colors,
            show_colnames = FALSE,
            main = str_c(path_name, group, sep = "_"),
            cluster_cols = FALSE,
            silent = TRUE
          )
        cowplot::plot_grid(p$gtable)
      })
  })


# export
pdf(
  file = file.path(out_dir, "pathway heatmap.pdf"),
  height = 8,
  width  = 12
)
pathway_heatmap_ls  %>% walk(print)
dev.off()


pathway_cols <-
  pathway_ls %>%
  imap(function(x, path_name){# browser()
    avg_expr_df %>%
      group_by(Cell) %>% nest %>% deframe() %>%
      imap(function(avg_expr_sub, cell){
        avg_expr_sub %>%
          mutate(expr = log1p(expr)) %>%
          dplyr::filter(Gene %in% x) %>%
          ggplot(aes(x = group, y = expr, fill = group)) +
          geom_col() +
          facet_grid(Gene ~ CN) +
          ggsci::scale_fill_igv() +
          ggtitle(str_c(path_name, cell, sep = "_")) +
          theme(axis.text.x = element_text(angle = 90, size = 3, hjust = 1, vjust = 0.5))
      })
      
  })

# export
pdf(
  file = file.path(out_dir, "pathway barplot.pdf"),
  height = 12,
  width  = 12
)
pathway_cols  %>% walk(print)
dev.off()

# geneset2 -------
# group CN
# gene expr
avg_expr2 <-
  all_data_cn %>%
  AverageExpression(group.by = c('group','project', "CN2"))

avg_expr_df2 <- 
  avg_expr2$RNA %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "group", values_to = "expr") %>%
  separate(group, into = c("group", "project", "CN"), sep = "_")


pathway_cols2 <-
  pathway_ls %>%
  imap(function(x, path_name){# browser()
    avg_expr_df2 %>%
      mutate(expr = log1p(expr)) %>%
      dplyr::filter(Gene %in% x) %>%
      ggplot(aes(x = group, y = expr, fill = group)) +
      # geom_col() +
      geom_boxplot() +
      facet_grid(Gene ~ CN) +
      ggsci::scale_fill_igv() +
      ggtitle(str_c(path_name, sep = "_")) +
      theme(axis.text.x = element_text(angle = 90, size = 3, hjust = 1, vjust = 0.5))
    
  })

# export
pdf(
  file = file.path(out_dir, "pathway boxplot - Only CN.pdf"),
  height = 12,
  width  = 12
)
pathway_cols2  %>% walk(print)
dev.off()

# geneset3 ------
# gene expr
avg_expr3 <-
  all_data_cn %>%
  AverageExpression(group.by = c('sub_group', "CN2", "annotation2"))

avg_expr_df3 <- 
  avg_expr3$RNA %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "group", values_to = "expr") %>%
  separate(group, into = c("group", "CN", "Cell"), sep = "_")


pathway_cols3 <-
  pathway_ls %>%
  imap(function(x, path_name){# browser()
    avg_expr_df3 %>%
      group_by(Cell) %>% nest %>% deframe() %>%
      imap(function(avg_expr_sub, cell){
        avg_expr_sub %>%
          mutate(expr = log1p(expr)) %>%
          dplyr::filter(Gene %in% x) %>%
          ggplot(aes(x = group, y = expr, fill = group)) +
          geom_col() +
          facet_grid(Gene ~ CN) +
          ggsci::scale_fill_igv() +
          ggtitle(str_c(path_name, cell, sep = "_")) +
          theme(axis.text.x = element_text(angle = 90, size = 3, hjust = 1, vjust = 0.5))
      })
    
  })

# export
pdf(
  file = file.path(out_dir, "pathway barplot - sub_group.pdf"),
  height = 12,
  width  = 12
)
pathway_cols3  %>% walk(print)
dev.off()


# pathway score ----
p_pathway_score_ls <-
  pathway_ls %>%
  imap(function(geneset, path_name){ # browser()
    geneset <- geneset[geneset %in% rownames(all_data_cn)]
    dat <-
      all_data_cn %>%
      AddModuleScore(
        features = list(geneset),
        nbin = 24,
        ctrl = 10) %>%
      .[[]]  
      # dplyr::rename({{path_name}} := Cluster1)
    
    p <- dat %>%
      group_by(project, group, CN2) %>% 
      summarise(Cluster1 = mean(Cluster1)) %>%
      
      ggplot(aes(x = group, y = Cluster1, fill = group)) +
      geom_boxplot() +
      facet_wrap(~ CN2) +
      labs(title = path_name, y = path_name) +
      ggpubr::geom_pwc(tip.length = 0) +
      scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]))
    p
  })

# export
pdf(
  file = file.path(out_dir, "pathway score boxplot.pdf"),
  height = 9,
  width  = 9
)
p_pathway_score_ls  %>% walk(print)
dev.off()


# pathway score - T cell and tumor --------------
T_cell <- 
  readxl::read_excel("tumor and T related genes.xlsx", sheet = "T cell") %>%
  dplyr::select(Gene, Function) %>%
  with(split(Gene, Function))
  
  
  
tumor <- 
  readxl::read_excel("tumor and T related genes.xlsx", sheet = "tumor") %>%
  dplyr::select(Gene, Function) %>%
  with(split(Gene, Function))


p_pathway_score_t_cell_ls <-
  T_cell %>%
  imap(function(geneset, path_name){ # browser()
    geneset <- geneset[geneset %in% rownames(all_data_cn)]
    dat <-
      all_data_cn %>%
      AddModuleScore(
        features = list(geneset),
        nbin = 24,
        ctrl = 10) %>%
      .[[]]  
    # dplyr::rename({{path_name}} := Cluster1)
    
    p <- dat %>%
      group_by(project, group, CN2) %>% 
      summarise(Cluster1 = mean(Cluster1)) %>%
      
      ggplot(aes(x = group, y = Cluster1, fill = group)) +
      geom_boxplot() +
      facet_wrap(~ CN2) +
      labs(title = path_name, y = path_name) +
      ggpubr::geom_pwc(tip.length = 0) +
      scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]))
    p
  })

p_pathway_score_tumor_ls <-
  tumor %>%
  imap(function(geneset, path_name){ # browser()
    geneset <- geneset[geneset %in% rownames(all_data_cn)]
    dat <-
      all_data_cn %>%
      AddModuleScore(
        features = list(geneset),
        nbin = 24,
        ctrl = 10) %>%
      .[[]]  
    # dplyr::rename({{path_name}} := Cluster1)
    
    p <- dat %>%
      group_by(project, group, CN2) %>% 
      summarise(Cluster1 = mean(Cluster1)) %>%
      
      ggplot(aes(x = group, y = Cluster1, fill = group)) +
      geom_boxplot() +
      facet_wrap(~ CN2) +
      labs(title = path_name, y = path_name) +
      ggpubr::geom_pwc(tip.length = 0) +
      scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]))
    p
  })

# export
pdf(
  file = file.path(out_dir, "pathway score boxplot - T cell.pdf"),
  height = 9,
  width  = 9
)
p_pathway_score_t_cell_ls  %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir, "pathway score boxplot - tumor.pdf"),
  height = 9,
  width  = 9
)
p_pathway_score_tumor_ls  %>% walk(print)
dev.off()


# pathway score - -----
immune_checkpoint <- c(
  'PDCD1',
  'CD274',
  'LAG3',
  'TIGIT',
  'CTLA4',
  'TNFRSF9',
  'TNFRSF18',
  'HAVCR2',
  'PDCD1LG2',
  'CD80',
  'CD86',
  'CD276',
  'ICOS',
  'ICOSLG',
  'TMIGD2',
  'TNFRSF14',
  'TNFRSF25',
  'TNFRSF4',
  'TNFRSF8',
  'LGALS9',
  'CD244',
  'CD27',
  'CD40',
  'CD40LG',
  'CD44',
  'CD48',
  'CD70',
  'CD160',
  'CD200',
  'CD200R1',
  'NRP1',
  'TNFSF9',
  'TNFSF14',
  'LAIR1'
)

p_pathway_score_immune_checkpoint_ls <-
  list(immune_checkpoint = immune_checkpoint) %>%
  imap(function(geneset, path_name){ # browser()
    geneset <- geneset[geneset %in% rownames(all_data_cn)]
    dat <-
      all_data_cn %>%
      AddModuleScore(
        features = list(geneset),
        nbin = 24,
        ctrl = 10) %>%
      .[[]]  
    # dplyr::rename({{path_name}} := Cluster1)
    
    p <- dat %>%
      group_by(project, group, CN2) %>% 
      summarise(Cluster1 = mean(Cluster1)) %>%
      
      ggplot(aes(x = group, y = Cluster1, fill = group)) +
      geom_boxplot() +
      facet_wrap(~ CN2) +
      labs(title = path_name, y = path_name) + 
      scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]))
    p
  })

pdf(
  file = file.path(out_dir, "pathway score boxplot - immune checkpoint.pdf"),
  height = 9,
  width  = 9
)
p_pathway_score_immune_checkpoint_ls  %>% walk(print)
dev.off()


# 2024.9.25 -------------------
# T and Tumor gene, A vs B in same CN --------------

T_cell <- 
  readxl::read_excel("tumor and T related genes.xlsx", sheet = "T cell") %>%
  dplyr::select(Gene, Function) %>%
  with(split(Gene, Function))



tumor <- 
  readxl::read_excel("tumor and T related genes.xlsx", sheet = "tumor") %>%
  dplyr::select(Gene, Function) %>%
  with(split(Gene, Function))


bulk_all_data <-
  all_data_cn %>%
  AggregateExpression(group.by = c("project", "CN2"), return.seurat = TRUE) %>%
  NormalizeData()


bulk_all_data$project <-
  bulk_all_data$project %>% 
  str_remove("^g")

  

project_group <- 
  all_data_cn_meta[c('project', 'group')] %>%
  distinct()


plot_data <-
  LayerData(bulk_all_data, "data") %>%
  t %>%
  cbind(bulk_all_data[[]]) %>%
  as.data.frame() %>%
  mutate(project = str_remove(project, "^g")) %>%
  left_join(project_group, by = "project")


# plot
p_AvsB_cn_t_ls <-
  T_cell %>%
  map(function(gene_ls){
    gene_ls %>%
      set_names(., .) %>%
      map(function(gene){
        var_cols <- c(gene, "project", "CN2", "group")
        
        plot_data[var_cols] %>%
          dplyr::filter(group != "C") %>%
          ggplot(aes(x = CN2, y = .data[[gene]], fill = group)) +
          geom_boxplot() +
          ggpubr::stat_pwc(tip.length = FALSE) + #, hide.ns = TRUE) +
          ggpubr::theme_pubclean() +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))
      })
  })

p_AvsB_cn_tumor_ls <-
  tumor %>%
  map(function(gene_ls){
    gene_ls %>%
      set_names(., .) %>%
      map(function(gene){
        var_cols <- c(gene, "project", "CN2", "group")
        
        plot_data[var_cols] %>%
          dplyr::filter(group != "C") %>%
          ggplot(aes(x = CN2, y = .data[[gene]], fill = group)) +
          geom_boxplot() +
          ggpubr::stat_pwc(tip.length = FALSE) + #, hide.ns = TRUE) +
          ggpubr::theme_pubclean() +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))
      })
  })

# ggsave
pdf(
  file.path(out_dir, "T gene boxplot A vs B in one CN.pdf"),
  height = 6, 
  width  = 8
)
p_AvsB_cn_t_ls %>%
  flatDoubleList() %>%
  imap(~ .x + ggtitle(.y)) %>%
  walk(print)

dev.off()

pdf(
  file.path(out_dir, "Tumor gene boxplot A vs B in one CN.pdf"),
  height = 6, 
  width  = 8
)
p_AvsB_cn_tumor_ls %>%
  flatDoubleList() %>%
  imap(~ .x + ggtitle(.y)) %>%
  walk(print)

dev.off()


# v2 T and Tumor gene, A vs B in same CN - v2 --------------

bulk_all_data <-
  all_data_cn %>%
  AggregateExpression(group.by = c("project", "CN2"), return.seurat = TRUE) %>%
  NormalizeData()

bulk_meta_1 <-
  bulk_all_data[['project']] %>%
  rownames_to_column('id') %>%
  mutate(project = str_remove(project, "^g")) %>%
  left_join(distinct(all_data_cn_meta[c("project", "group")]), by = "project") %>%
  column_to_rownames('id')
bulk_all_data <- bulk_all_data %>% AddMetaData(bulk_meta_1) 
  
Idents(bulk_all_data) <- "group"
de_AvsB_cn_ls <-
  bulk_all_data$CN2 %>%
  unique() %>%
  set_names(., .) %>%
  map(function(cn){
    sc <- bulk_all_data[, bulk_all_data$CN2 == cn]
    FindMarkers(
      sc,
      ident.1 = 'A',
      ident.2 = 'B',
      test.use = "DESeq2"
    )
  })

p_de_AvsB_cn_ls <-
  de_AvsB_cn_ls %>%
  imap(function(dat, name){# browser()
    dat <- 
      dat %>%
      rownames_to_column("Gene") %>%
      mutate(
        `-log10pval` = -log10(p_val_adj),
        Group = case_when(
          avg_log2FC >  1 & p_val_adj < 0.05 ~ "UP",
          avg_log2FC < -1 & p_val_adj < 0.05 ~ "DOWN",
          TRUE ~ "non-sig"
        )
      ) 
    
    sig_dat <- 
      dat %>%
      dplyr::filter(Group != "non-sig") %>%
      group_by(Group) %>%
      slice_min(p_val, n = 10)
    
    dat %>%
      
      ggplot(aes(x = avg_log2FC, y = `-log10pval`)) +
      geom_point(aes(color = Group)) +
      geom_hline(yintercept = 1.30103, linetype = 2) +
      geom_vline(xintercept = c(-1, 1), linetype = 2) +
      
      ggrepel::geom_text_repel(data = sig_dat, aes(label = Gene)) +
      
      scale_color_manual(values = c(
        UP       = "red",
        DOWN     = "blue",
        "non-sig" = "gray"
      )) +
      theme_bw() +
      theme(
        legend.position = "none"
      ) +
      labs(
        title = name,
        subtitle = "A vs B"
      )
    
  })

# export
pdf(
  file.path(out_dir, "All gene volcanoplot A vs B in one CN.pdf"),
  height = 6, 
  width  = 8
)
p_de_AvsB_cn_ls %>%
  walk(print)

dev.off()




# T and Tumor gene, pairwise compare between CN in one project --------------

# bulk_all_data <-
#   all_data_cn %>%
#   AggregateExpression(group.by = c("project", "CN2"), return.seurat = TRUE) %>%
#   NormalizeData()

p_pwc_cn_t_ls <-
  T_cell %>%
  map(function(gene_ls){
    gene_ls %>%
      set_names(., .) %>%
      map(function(gene){
        var_cols <- c(gene, "project", "CN2", "group")
        
        plot_data[var_cols] %>%
          mutate(project = str_c(project, group, sep = "_")) %>%
          ggplot(aes(x = CN2, y = .data[[gene]], fill = CN2)) +
          geom_col() +
          # ggpubr::stat_pwc(tip.length = FALSE) + #, hide.ns = TRUE) +
          ggsci::scale_fill_igv() +
          ggpubr::theme_pubclean() +
          facet_wrap(~ project) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))
      })
  })

p_pwc_cn_tumor_ls <-
  tumor %>%
  map(function(gene_ls){
    gene_ls %>%
      set_names(., .) %>%
      map(function(gene){
        var_cols <- c(gene, "project", "CN2", "group")
        
        plot_data[var_cols] %>%
          mutate(project = str_c(project, group, sep = "_")) %>%
          ggplot(aes(x = CN2, y = .data[[gene]], fill = CN2)) +
          geom_col() +
          # ggpubr::stat_pwc(tip.length = FALSE) + #, hide.ns = TRUE) +
          ggsci::scale_fill_igv() +
          ggpubr::theme_pubclean() +
          facet_wrap(~ project) +
          theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))
      })
  })

# ggsave
pdf(
  file.path(out_dir, "T gene colplot of CN in one sample.pdf"),
  height = 6 * 2.5, 
  width  = 8 * 3
)
p_pwc_cn_t_ls %>%
  flatDoubleList() %>%
  imap(~ .x + ggtitle(.y)) %>%
  walk(print)

dev.off()

pdf(
  file.path(out_dir, "Tumor gene colplot of CN in one sample.pdf"),
  height = 6 * 2.5, 
  width  = 8 * 3
)
p_pwc_cn_tumor_ls %>%
  flatDoubleList() %>%
  imap(~ .x + ggtitle(.y)) %>%
  walk(print)

dev.off()

# v2 T and Tumor gene, pairwise compare between CN in one project --------------

Idents(bulk_all_data) <- "CN2"
de_cn_pwc_group_ls <-
  bulk_all_data$group %>%
  unique() %>%
  set_names(., .) %>%
  map(function(g){
    sc <- bulk_all_data[, bulk_all_data$group == g]
    FindAllMarkers(
      sc,
      
      logfc.threshold = 0,
      min.pct = 0,
      min.cells.feature = 0,
      min.cells.group = 0,
      return.thresh = 1,
      
      group.by = "CN2",
      test.use = "DESeq2"
    )
  })


p_de_cn_pwc_group_ls <-
  de_cn_pwc_group_ls %>%
  .[! names(.) %in% 'C'] %>%
  imap(function(dat, name){# browser()
    set.seed(1234)
    
    dat <- dat%>% 
      slice_sample(n = nrow(.)) %>%
      mutate(group = ifelse(p_val_adj < 0.05, 'p < 0.05', 'p > 0.05')) %>%
      mutate(cluster_x = as.numeric(factor(cluster)) - 1) %>%
      group_by(cluster) %>%
      mutate(order = 1:n()/n()) %>%
      mutate(order = order * 0.8 + 0.1) %>%
      mutate(x = order + cluster_x) %>%
      mutate(avg_log2FC = ifelse(avg_log2FC > 0, avg_log2FC + 0.5, avg_log2FC - 0.5))
    
    bar_dat <- dat %>% group_by(cluster) %>% summarise(x = unique(cluster_x))
    sig_dat <- 
      dat %>%
      mutate(fc_group = case_when(
        avg_log2FC > 1  ~ "UP",
        avg_log2FC < -1 ~ "DOWN",
        TRUE ~ "non-sig"
      )) %>%
      dplyr::filter(fc_group != "non-sig", p_val_adj < 0.05) %>%
      group_by(cluster, fc_group) %>%
      # slice_min(p_val, n = 10)
      slice_max(abs(avg_log2FC), n = 10)
    
    
    dat %>%
      ggplot(aes(x = x, y = avg_log2FC)) +
      geom_point(aes(color = group), size = 0.5) +
      geom_tile(
        data = bar_dat,
        aes(x = x + 0.5, y = 0, fill = cluster),
        height = 0.5,
        width  = 1,
        inherit.aes = FALSE
      ) +
      geom_text(
        data =  bar_dat,
        aes(x = x + 0.5, y = 0, label = cluster),
        inherit.aes = FALSE
      ) +
      ggrepel::geom_text_repel(data = sig_dat, aes(label = gene)) +
      scale_color_manual(values = c("p < 0.05" = "red", "p > 0.05" = "gray")) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        legend.position = "none"
      ) +
      labs(
        title = name,
        # subtitle = "A vs B"
      )
    
  })

# export
pdf(
  file.path(out_dir, "All gene multi-de-plot pwc CN in one group A and B.pdf"),
  height = 6, 
  width  = 18
)
p_de_cn_pwc_group_ls %>%
  walk(print)

dev.off()



# 2024.10.8 -------------------
gene_grps <- c(
  "Fibroblast"           = "other",
  "Plasma"               = "immune",
  "Mono-like Macrophage" = "immune",
  "CD83-dim B"           = "immune",
  "Epi0_B"               = "immune",
  "CD8T"                 = "immune",
  "Endothelial"          = "other",
  "CD4T"                 = "immune",
  "Smooth muscle cells"  = "other",
  "Epi0_Fibroblast"      = "other",
  "CD4T-Macrophage"      = "immune",
  "Epi0_T"               = "immune",
  "B"                    = "immune",
  "pDCs"                 = "immune",
  "CD4 Endothelial"      = "immune",
  "Epi0_Mast cell"       = "immune",
  "Epi0_CD83_T"          = "immune",
  "Epithelial"           = "tumor",
  "Fat cells"            = "other",
  "CD83 B"               = "immune",
  "Epi0_CD4T"            = "immune"
)



bulk_all_data_for_AvsB <-
  all_data_cn %>%
  AggregateExpression(
    group.by = c("project", "CN2", "annotation2"), 
    return.seurat = TRUE
  ) 

bulk_all_data_for_AvsB_group <-
  all_data_cn[[c(c("project", "CN2", "annotation2", "group"))]] %>%
  distinct() %>%
  mutate(
    project = str_c("g", project, sep = ""), 
    id = str_c(project, CN2, annotation2, sep = "_")
  ) %>%
  dplyr::select(id, group) %>%
  remove_rownames() %>%
  column_to_rownames("id") 

bulk_all_data_for_AvsB <-
  bulk_all_data_for_AvsB %>%
  AddMetaData(bulk_all_data_for_AvsB_group) %>%
  NormalizeData()

# T and Tumor gene, A vs B in same CN and same cell --------------


# T_cell <- 
#   readxl::read_excel("tumor and T related genes.xlsx", sheet = "T cell") %>%
#   dplyr::select(Gene, Function) %>%
#   with(split(Gene, Function))
# 
# 
# 
# tumor <- 
#   readxl::read_excel("tumor and T related genes.xlsx", sheet = "tumor") %>%
#   dplyr::select(Gene, Function) %>%
#   with(split(Gene, Function))


# all_data_cn$annotation_full %>% table
# all_data_cn$annotation2 %>% table


all_celltypes <- all_data_cn$annotation2 %>% unique() %>% as.character()
all_cns <- all_data_cn$CN2 %>% unique() %>% as.character()

# all_de_df_ls <-
#   all_celltypes %>%
#   set_names(., .) %>%
#   map(function(cell){
#     all_cns %>%
#       set_names(., .) %>%
#       map(function(cn){ # browser()
#         cell.1 <- 
#           all_data_cn$annotation2 == cell &
#           all_data_cn$CN2 == cn &
#           all_data_cn$group == "A"
#         cell.1 <- colnames(all_data_cn)[cell.1]
#         
#         cell.2 <- 
#           all_data_cn$annotation2 == cell &
#           all_data_cn$CN2 == cn &
#           all_data_cn$group == "B"
#         cell.2 <- colnames(all_data_cn)[cell.2]
#         
#         # FindMarkers
#         de <- tryCatch(
#           
#           FindMarkers(
#             all_data_cn[['RNA']],
#             cells.1 = cell.1,
#             cells.2 = cell.2,
#             fc.slot = "data",
#             pseudocount.use = 1,
#             base = 2,
#             logfc.threshold = 0.1,
#             min.pct = 0.1, # 0.01
#             only.pos = FALSE
#           ),
#           error = function(e) return(NULL)
#         )
#         
#         de
#       })
#   })

# saveRDS(all_de_df_ls, 'all_de_df_ls.rds')

all_de_df_ls <-
  all_celltypes %>%
  set_names(., .) %>%
  map(function(cell){
    all_cns %>%
      set_names(., .) %>%
      map(function(cn){ # browser()
        cell.1 <- 
          bulk_all_data_for_AvsB$annotation2 == cell &
          bulk_all_data_for_AvsB$CN2 == cn &
          bulk_all_data_for_AvsB$group == "A"
        cell.1 <- colnames(bulk_all_data_for_AvsB)[cell.1]
        
        cell.2 <- 
          bulk_all_data_for_AvsB$annotation2 == cell &
          bulk_all_data_for_AvsB$CN2 == cn &
          bulk_all_data_for_AvsB$group == "B"
        cell.2 <- colnames(bulk_all_data_for_AvsB)[cell.2]
        
        # FindMarkers
        de <- tryCatch(
          
          FindMarkers(
            bulk_all_data_for_AvsB[['RNA']],
            cells.1 = cell.1,
            cells.2 = cell.2,
            fc.slot = "data",
            test.use = "DESeq2",
            pseudocount.use = 1,
            base = 2,
            logfc.threshold = 0, # 0.1,
            min.pct = 0.01,
            only.pos = FALSE
          ),
          error = function(e) NULL # return(e$message)
        )
        
        de
      })
  })

p_de_AvsB_same_cell_cn_ls <-
  all_de_df_ls %>%
  flatDoubleList(sep = " - (Cell - CN) - ") %>%
  imap(function(dat, name){ # browser()
    
    if(is.null(dat)) return(NULL)
    
    dat <- 
      dat %>%
      as.data.frame() %>%
      rownames_to_column("Gene") %>%
      mutate(
        `-log10pval` = -log10(p_val_adj) %>% pmin(12),
        Group = case_when(
          avg_log2FC >  0.585 & p_val_adj < 0.05 ~ "UP",
          avg_log2FC < -0.585 & p_val_adj < 0.05 ~ "DOWN",
          TRUE ~ "non-sig"
        )
      ) 
    
    sig_dat <- 
      dat %>%
      dplyr::filter(Group != "non-sig") %>%
      group_by(Group) %>%
      slice_min(p_val, n = 10)
    
    dat %>%
      
      ggplot(aes(x = avg_log2FC, y = `-log10pval`)) +
      geom_point(aes(color = Group)) +
      geom_hline(yintercept = 1.30103, linetype = 2) +
      geom_vline(xintercept = c(-0.585, 0.585), linetype = 2) +
      
      ggrepel::geom_text_repel(data = sig_dat, aes(label = Gene)) +
      
      scale_color_manual(values = c(
        UP       = "red",
        DOWN     = "blue",
        "non-sig" = "gray"
      )) +
      theme_bw() +
      theme(
        legend.position = "none"
      ) +
      labs(
        title = name,
        subtitle = "A vs B"
      )
  })

# export
pdf(
  file.path(out_dir, "10.9 All gene volcanoplot A vs B in one cell and CN - FC1.5.pdf"),
  height = 6,
  width  = 8
)
p_de_AvsB_same_cell_cn_ls %>%
  walk(print)

dev.off()


saveRDS(all_de_df_ls, "all_de_df_ls-bulk.rds")

# T and Tumor gene, pwc CN in same group and same cell --------------

all_groups <- all_data_cn$group %>% unique() %>% as.character()
all_celltypes <- all_data_cn$annotation2 %>% unique() %>% as.character()

all_de_df_ls2 <-
  all_groups %>%
  setdiff("C") %>%
  set_names(., .) %>%
  map(function(group){
    all_celltypes %>%
      set_names(., .) %>%
      map(function(cell){ # browser()
        
        sc <- all_data_cn[,
          all_data_cn$group == group & all_data_cn$annotation2 == cell
        ]
        
        # FindMarkers
        Idents(sc) <- 'CN2'
        de <- tryCatch(
          
          FindAllMarkers(
            sc,
            
            logfc.threshold = 0,
            min.pct = 0,
            min.cells.feature = 0,
            min.cells.group = 0,
            return.thresh = 1
          ),
          error = function(e) return(NULL)
        )
        
        de
      })
  })

# saveRDS(all_de_df_ls2, 'all_de_df_ls2.rds')

plot_multi_de_scatter <- 
  function(dat, name){# browser()
    set.seed(1234)
    
    dat <- dat%>% 
      slice_sample(n = nrow(.)) %>%
      mutate(group = ifelse(p_val_adj < 0.05, 'p < 0.05', 'p > 0.05')) %>%
      mutate(cluster_x = as.numeric(factor(cluster)) - 1) %>%
      group_by(cluster) %>%
      mutate(order = 1:n()/n()) %>%
      mutate(order = order * 0.8 + 0.1) %>%
      mutate(x = order + cluster_x) %>%
      mutate(avg_log2FC = ifelse(avg_log2FC > 0, avg_log2FC + 0.5, avg_log2FC - 0.5))
    
    bar_dat <- dat %>% group_by(cluster) %>% summarise(x = unique(cluster_x))
    sig_dat <- 
      dat %>%
      mutate(fc_group = case_when(
        avg_log2FC > 1  ~ "UP",
        avg_log2FC < -1 ~ "DOWN",
        TRUE ~ "non-sig"
      )) %>%
      dplyr::filter(fc_group != "non-sig", p_val_adj < 0.05) %>%
      group_by(cluster, fc_group) %>%
      # slice_min(p_val, n = 10)
      slice_max(abs(avg_log2FC), n = 10)
    
    
    dat %>%
      ggplot(aes(x = x, y = avg_log2FC)) +
      geom_point(aes(color = group), size = 0.5) +
      geom_tile(
        data = bar_dat,
        aes(x = x + 0.5, y = 0, fill = cluster),
        height = 1,
        width  = 1,
        inherit.aes = FALSE
      ) +
      geom_text(
        data =  bar_dat,
        aes(x = x + 0.5, y = 0, label = cluster),
        inherit.aes = FALSE
      ) +
      ggrepel::geom_text_repel(data = sig_dat, aes(label = gene)) +
      scale_color_manual(values = c("p < 0.05" = "red", "p > 0.05" = "gray")) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        legend.position = "none"
      ) +
      labs(
        title = name,
        # subtitle = "A vs B"
      )
  }


# mult de plot

p_all_de_df_cn_pwc_ls <-
  all_de_df_ls2  %>%
  imap(function(dat_ls, grp){
    dat_ls %>%
      imap(function(dat, name){
        plot_multi_de_scatter(dat, paste(grp, name, sep = "_"))
      })
  })

# export
pdf(
  file.path(out_dir, "10.9 All gene multi-de-plot pwc CN in one group and one cell - v2.pdf"),
  height = 6, 
  width  = 25
)
p_all_de_df_cn_pwc_ls %>%
  flatDoubleList() %>%
  walk(print)

dev.off()


## gene expr plot - T gene ----

bulk_all_data_avg <-
  all_data_cn %>%
  AverageExpression(
    group.by = c("project", "CN2", "annotation2"), 
    return.seurat = FALSE
  ) 
bulk_all_data_avg <- bulk_all_data_avg$RNA

# saveRDS(bulk_all_data_avg, 'bulk_all_data_avg.rds')


bulk_all_data_avg <-
  bulk_all_data_avg %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  pivot_longer(
    -gene,
    names_to  = "group",
    values_to = "value"
  ) %>%
  separate(
    group,
    into = c("project", "CN", "cell"),
    sep  = "_"
  )


bulk_all_data_avg$project <-
  bulk_all_data_avg$project %>% 
  str_remove("^g")

plot_data_avg_expr <- 
  bulk_all_data_avg %>%
  left_join(
    distinct(all_data_cn[[c('project', 'group')]] ),
    by = "project"
  )


# t_gene
p_t_gene_avg_expr_ls <-
  plot_data_avg_expr %>% 
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  .[! names(.) %in% "C"] %>%
  imap(function(dat, grp){ # <-- group
    T_cell %>%
      imap(function(gene, cell_name){
        dat %>%
          group_by(cell) %>%
          nest %>%
          deframe() %>%
          imap(function(sub_dat, annotation){ # <-- cell
            # browser()
            gene %>%
              set_names(., .) %>%
              map(function(gene){ # <-- gene
                sub_dat2 <- sub_dat[sub_dat$gene == gene, ]
                
                avg_expr <- 
                  sub_dat2 %>%
                  group_by(CN) %>%
                  summarise(value = mean(value)) %>%
                  ungroup() %>%
                  arrange(desc(value))
                
                sub_dat2 %>%
                  mutate(CN = factor(CN, levels = avg_expr$CN)) %>%
                  ggplot(aes(x = CN, y = value, fill = CN)) +
                  geom_boxplot(outliers = FALSE) +
                  geom_point(position = position_jitter(width = 0.15)) +
                  ggpubr::theme_pubr() +
                  ggsci::scale_fill_igv() +
                  theme(
                    legend.position = "right",
                    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5)
                  ) +
                  labs(
                    x = NULL,
                    y = 'gene expression',
                    title = str_glue("{cell_name} marker {gene} in {annotation} cell"),
                    subtitle = "marker gene in cell"
                  )
              })

          })
          
      })
  })



# export
dir.create(file.path(out_dir, "gene expression in CN", "t gene"), recursive = TRUE)

p_t_gene_avg_expr_ls %>%
  flatDoubleList(sep = " group -- ") %>%
  flatDoubleList(sep = " marker gene - in - cells ") %>% 
  iwalk(~{
    pdf(
      file.path(
        out_dir, "gene expression in CN", "t gene", 
        str_glue("{.y}.pdf")
      ),
      height = 6, 
      width  = 8
    )
    print(.x)
    dev.off()
  })



## gene expr plot - tumor gene ----

p_tumor_gene_avg_expr_ls <-
  plot_data_avg_expr %>% 
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  .[! names(.) %in% "C"] %>%
  imap(function(dat, grp){ # <-- group
    tumor %>%
      imap(function(gene, cell_name){
        dat %>%
          group_by(cell) %>%
          nest %>%
          deframe() %>%
          imap(function(sub_dat, annotation){ # <-- cell
            # browser()
            gene %>%
              set_names(., .) %>%
              map(function(gene){ # <-- gene
                sub_dat2 <- sub_dat[sub_dat$gene == gene, ]
                
                avg_expr <- 
                  sub_dat2 %>%
                  group_by(CN) %>%
                  summarise(value = mean(value)) %>%
                  ungroup() %>%
                  arrange(desc(value))
                
                sub_dat2 %>%
                  mutate(CN = factor(CN, levels = avg_expr$CN)) %>%
                  ggplot(aes(x = CN, y = value, fill = CN)) +
                  geom_boxplot(outliers = FALSE) +
                  geom_point(position = position_jitter(width = 0.15)) +
                  ggpubr::theme_pubr() +
                  ggsci::scale_fill_igv() +
                  theme(
                    legend.position = "right",
                    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5)
                  ) +
                  labs(
                    x = NULL,
                    y = 'gene expression',
                    title = str_glue("{cell_name} marker {gene} in {annotation} cell"),
                    subtitle = "marker gene in cell"
                  )
              })
            
          })
        
      })
  })



# export
dir.create(file.path(out_dir, "gene expression in CN", "tumor gene"), recursive = TRUE)

p_tumor_gene_avg_expr_ls %>%
  flatDoubleList(sep = " group -- ") %>%
  flatDoubleList(sep = " marker gene - in - cells ") %>% 
  iwalk(~{
    pdf(
      file.path(
        out_dir, "gene expression in CN", "tumor gene", 
        str_glue("{.y}.pdf")
      ),
      height = 6, 
      width  = 8
    )
    print(.x)
    dev.off()
  })



# ACTG2  --------------
all_data_cn %>%
  VlnPlot(features = "ACTG2", group.by = "annotation2", pt.size = 0)



