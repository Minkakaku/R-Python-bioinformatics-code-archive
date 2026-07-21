library(tidyverse)
library(Seurat)

out_dir_0107 <- "dat_analysis_250107"
dir.create(out_dir_0107)


all_data <- readRDS("all_data.rds")

# parse roi -------------
all_rois_path <- list.files("cell coords_20250107/", full.names = TRUE) 
names(all_rois_path) <- 
  all_rois_path %>%
  basename() %>%
  # str_remove("^.*_") %>%
  str_remove("\\..*$")

all_rois_ls <-
  all_rois_path %>%
  map(read.csv) %>%
  imap(function(dat, name){
    name <- ifelse(str_detect(name, "DZZ_Saline"), "DZZ_Saline", "SYZ_NH4")
    dat$Barcode <- str_c(name, dat$Barcode, sep = "_")
    dat
  })

all_rois_df <-
  all_rois_ls %>%
  map(function(x){
    x %>%
      dplyr::rename(group = 2) %>%
      dplyr::filter( group != "")
  }) %>%
  bind_rows()

# some cells of tow ROIs are duplicated ---------
dup_cells <- all_rois_df$Barcode %>% .[duplicated(.)]
dup_roi_dat <- all_rois_df[all_rois_df$Barcode %in% dup_cells, ]
dup_roi_dat$group %>% unique() # "NH4-100um-sp5-R" "NH4-100um-VRG-R"

keep_dup_group <- "NH4-100um-VRG-R"

all_rois_df_remove_dup <- all_rois_df[
  ! (all_rois_df$Barcode %in% dup_cells & all_rois_df$group != keep_dup_group),
]

all_rois_df_remove_dup$group[all_rois_df_remove_dup$Barcode %in% dup_cells] %>% unique
# "NH4-100um-VRG-R"



# all_data_roi ----------------
table(all_rois_df_remove_dup$Barcode %in% colnames(all_data) )

all_data_roi <- all_data %>% subset(cells = intersect(all_rois_df_remove_dup$Barcode, colnames(all_data)))

all_data_roi <- 
  all_data_roi %>%
  AddMetaData(column_to_rownames(remove_rownames(all_rois_df_remove_dup), "Barcode"))

all_data_roi$group2 <- 
  all_data_roi$group %>%
  str_remove("^.*um-") %>%
  str_remove("-.*$")

# cell annotation according to gene expression ---------------
choose_genes_ls <-
  list(
    sp5 = c("Calb1", "Pdyn", "Pvalb", "Sst", "Tac1", "Th"),
    VRG = c("Chat", "Dbh", "Phox2b", "Pvalb", "Sst", "Th"),
    NTS = c("Calb1", "Chat", "Dbh", "Pdyn", "Phox2b", "Tac1", "Th")
  )

choose_genes_ls %>%
  map(~ . %in% rownames(all_data_roi) )



all_data_roi_ls <-
  all_data_roi %>% 
  SplitObject("group2")


cells_df_ls <-
  all_data_roi_ls %>% 
  imap(function(x, name){ # browser()
    expr_df <- x %>% FetchData(vars = choose_genes_ls[[name]]) 
    T_F_df  <- expr_df > 0 
    
    all_zeros <- rowSums(T_F_df) == 0
    
    cells <- T_F_df[!all_zeros, ] %>% apply( 1, function(y) str_c(colnames(expr_df)[y], collapse = "_") ) 
    
    cells_full <- rep("", nrow(T_F_df))
    cells_full[!all_zeros] <- cells
      
    data.frame(
      row.names = rownames(T_F_df),
      cells = cells_full,
      number = rowSums(T_F_df)
    ) 
  }) 

cells_df <-
  cells_df_ls %>%
  bind_rows()

all_data_roi <- all_data_roi %>% AddMetaData(cells_df)



# 同一区域 de analysis ---------
all_data_roi_only_one <- all_data_roi[, all_data_roi$number == 1]

Idents(all_data_roi_only_one) <- "group2"

all_de_df_ls <-
  all_data_roi_only_one$cells %>% 
  unique() %>% 
  set_names(., .) %>%
  map(function(x){
    all_data_roi_only_one$group2 %>%
      unique() %>%
      set_names(., .) %>%
      map(function(y){
        tryCatch(
          FindMarkers(
            JoinLayers(all_data_roi_only_one),
            group.by        = "cells",
            ident.1         = x,
            subset.ident    = y,
            logfc.threshold = 0,
            min.pct         = 0
          ),
          error = function(e) return(NULL)
        )
        
      })
  })



all_de_df <-
  all_de_df_ls %>%
  map(function(de_ls){
    de_ls <- de_ls[!map_lgl(de_ls, is.null)]
    de_ls %>%
      map(~ rownames_to_column(., "gene")) %>%
      enframe(name = "group") %>%
      unnest(value)
  }) %>%
  enframe(name = "cells") %>%
  unnest(value)


## enrichment analysis ----------------
all_de_gene_ls <-
  all_de_df %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  mutate(id = str_c(group, cells, sep = "_")) %>%
  dplyr::select(id, gene) %>%
  group_by(id) %>%
  nest() %>%
  deframe() %>%
  map(unlist)

all_de_gene_id_ls <-
  all_de_gene_ls %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
    .$ENTREZID
  })

all_enrich_ls <-
  all_de_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      qvalueCutoff = 0.2, 
      minGSSize = 10,       
      maxGSSize = 500
    )
  })


all_enrich_dot_ls <-
  all_enrich_ls %>% 
  imap(function(enrich, name){
    if(nrow(as.data.frame(enrich)) < 1) return(NULL)
    
    enrichplot::dotplot(enrich, showCategory = 15, label_format = 50) +
      ggtitle(name)
    
  })



# export
pdf(
  file.path(out_dir_0107, 'enrichment analysis - highly experssed cells - dotplot.pdf'),
  height = 6,
  width  = 10
)
all_enrich_dot_ls %>% walk(print)
dev.off()


## dotplot------------
marker_gene_ls <-
  all_de_df %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  mutate(id = str_c(group, cells, sep = "_")) %>%
  group_by(id) %>%
  slice_max(avg_log2FC, n = 5) %>%
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  map('gene') 
  
# export
pdf(
  file.path(out_dir_0107, 'marker gene - highly experssed cells - dotplot.pdf'),
  height = 6,
  width  = 10
)
marker_gene_ls %>%
  imap(~{
    DotPlot(
      all_data_roi_only_one[, all_data_roi_only_one$group2 == .y],
      features = .x,
      group.by = 'cells'
    ) +
    ggtitle(.y) +
    RotatedAxis()
  })
dev.off()

# 区域间 de analysis -----------------

all_de_area_pwc_ls <-
  choose_genes_ls %>%
  names() %>%
  combn(2) %>%
  apply(2, c, simplify = FALSE) %>%
  setNames(., map(., str_c, collapse = "_")) %>%
  
  map(function(x){
    genes <- intersect(choose_genes_ls[[ x[1] ]], choose_genes_ls[[ x[2] ]])
    
    genes %>%
      set_names(., .) %>%
      map(function(g){
        cells.1 <- colnames(all_data_roi_only_one)[all_data_roi_only_one$cells == g & all_data_roi_only_one$group2 %in% x[1]]
        cells.2 <- colnames(all_data_roi_only_one)[all_data_roi_only_one$cells == g & all_data_roi_only_one$group2 %in% x[2]]
        
        
        FindMarkers(
          JoinLayers(all_data_roi_only_one)[['Spatial.008um']],
          cells.1 = cells.1,
          cells.2 = cells.2,
          logfc.threshold = 0,
          min.pct         = 0
        )
      })
  })


all_de_area_pwc_df <-
  all_de_area_pwc_ls %>%
  map(function(de_ls){
    de_ls <- de_ls[!map_lgl(de_ls, is.null)]
    de_ls %>%
      map(~ rownames_to_column(., "gene")) %>%
      enframe(name = "cells") %>%
      unnest(value)
  }) %>%
  enframe(name = "group") %>%
  unnest(value)
 
## enrichment analysis ----------------
all_de_area_pwc_df2 <-
  all_de_area_pwc_df %>%
  mutate(group = group %>% str_split("_") %>% map_chr(~ str_c(rev(.x), collapse = "_")) ) %>%
  mutate(avg_log2FC  = - avg_log2FC )

all_de_gene_area_pwc_ls <-
  bind_rows(
    all_de_area_pwc_df,
    all_de_area_pwc_df2
  ) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  mutate(id = str_c(group, cells, sep = "-")) %>%
  dplyr::select(id, gene) %>%
  group_by(id) %>%
  nest() %>%
  deframe() %>%
  map(unlist)

all_de_gene_id_area_pwc_ls <-
  all_de_gene_area_pwc_ls %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
      .$ENTREZID
  })

all_enrich_area_pwc_ls <-
  all_de_gene_id_area_pwc_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      qvalueCutoff = 0.2, 
      minGSSize = 10,       
      maxGSSize = 500
    )
  })


all_enrich_area_pwc_dot_ls <-
  all_enrich_area_pwc_ls %>% 
  imap(function(enrich, name){
    if(nrow(as.data.frame(enrich)) < 1) return(NULL)
    
    enrichplot::dotplot(enrich, showCategory = 15, label_format = 50) +
      ggtitle(name)
    
  })



# export
pdf(
  file.path(out_dir_0107, 'enrichment analysis - brain area comparison - dotplot.pdf'),
  height = 6,
  width  = 10
)
all_enrich_area_pwc_dot_ls %>% walk(print)
dev.off()

  
## dotplot------------
marker_gene_ls2 <-
  all_de_area_pwc_df %>%
  dplyr::filter(abs(avg_log2FC) > 1, p_val_adj < 0.05) %>%
  mutate(id = str_c(group, cells, sep = "-")) %>%
  group_by(id) %>%
  {bind_rows(
    slice_max(., avg_log2FC, n = 5),
    slice_min(., avg_log2FC, n = 5) 
  )} %>%
  group_by(id) %>%
  nest() %>%
  deframe() %>%
  map('gene') 
  
  
# export
pdf(
  file.path(out_dir_0107, 'marker gene - brain area comparison - dotplot.pdf'),
  height = 6,
  width  = 10
)
marker_gene_ls2 %>% # .[1] %>%
  imap(function(marker_gene, name){
    
    # name: VRG_NTS-Chat
    name2 <- name %>% str_split_1("-")
    name3 <- name2[1] %>% str_split_1("_")
    
    DotPlot(
      all_data_roi_only_one[,
        all_data_roi_only_one$group2 %in% name3 &
          all_data_roi_only_one$cells %in% name2[2]
        ], 
        features = marker_gene %>% unique,
        group.by = 'group2'
      ) +
      ggtitle(name) +
      RotatedAxis()
  })

dev.off()


# 各区域 组别间 de analysis ------------

all_data_roi_only_one$group3 <-
  all_data_roi_only_one$project %>% 
  fct_recode(
     Control = 'DZZ_Saline',
     Test    = 'SYZ_NH4'
  )


all_de_grp_ls <-
  choose_genes_ls %>%
  
  imap(function(cells, grp){
    
    cells %>%
      set_names(., .) %>%
      map(function(cell){
        cells.1 <- all_data_roi_only_one %>% {colnames(.)[.$cells == cell & .$group2 == grp & .$group3 == "Test"]}
        cells.2 <- all_data_roi_only_one %>% {colnames(.)[.$cells == cell & .$group2 == grp & .$group3 == "Control"]}
        
        
        FindMarkers(
          JoinLayers(all_data_roi_only_one)[['Spatial.008um']],
          cells.1 = cells.1,
          cells.2 = cells.2,
          logfc.threshold = 0,
          min.pct         = 0
        )
      })
  })


all_de_grg_df <-
  all_de_grp_ls %>% 
  map(function(de_ls){
    de_ls %>%
      map(~ rownames_to_column(., "gene")) %>%
      enframe(name = "cells") %>%
      unnest(value)
  }) %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group = group %>% str_c(cells, sep = "_")) %>%
  mutate(de = "Test vs Control")




## enrichment analysis ---------------- 
all_de_grg_df2 <-
  all_de_grg_df %>%
  mutate(de = "Control vs Test") %>%
  mutate(avg_log2FC  = - avg_log2FC ) %>%
  mutate(temp = pct.1, pct.1 = pct.2, pct.2 = temp) %>% 
  dplyr::select(-temp)

all_de_gene_grg_ls <- 
  bind_rows(
    all_de_grg_df,
    all_de_grg_df2
  ) %>%
  dplyr::filter(avg_log2FC > 0.25, p_val_adj < 1) %>%
  mutate(id = str_c(group, de, sep = " - ")) %>%
  dplyr::select(id, gene) %>%
  group_by(id) %>%
  nest() %>%
  deframe() %>%
  map(unlist)

all_de_gene_id_grg_ls <-
  all_de_gene_grg_ls %>%
  map(function(x){
    tryCatch({
      clusterProfiler::bitr(
        x,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = "org.Mm.eg.db"
      ) %>%
      .$ENTREZID
    },
    error = function(e) return(NULL)
  )
})
 

all_enrich_grp_ls <-
  all_de_gene_id_grg_ls %>%
  discard(is.null) %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      qvalueCutoff = 0.2, 
      minGSSize = 10,       
      maxGSSize = 500
    )
  })


all_enrich_grg_dot_ls <-
  all_enrich_grp_ls %>% 
  imap(function(enrich, name){
    if(nrow(as.data.frame(enrich)) < 1) return(NULL)
    
    enrichplot::dotplot(enrich, showCategory = 15, label_format = 50) +
      ggtitle(name)
    
  })



# export
pdf(
  file.path(out_dir_0107, 'enrichment analysis - brain area comparison - Test vs Control - dotplot.pdf'),
  height = 6,
  width  = 10
)
all_enrich_grg_dot_ls %>% walk(print)
dev.off()


# export table 
all_enrich_grp_ls %>%
  map('result') %>%
  enframe("group") %>%
  unnest("value") %>%
  openxlsx::write.xlsx(
    file.path(out_dir_0107, "enrichment table - brain area comparison - Test vs Control.xlsx"),
    asTable = TRUE
  )




## dotplot------------
marker_gene_ls3 <-
  all_de_grg_df %>%
  dplyr::filter(abs(avg_log2FC) > 1, p_val_adj < 0.05) %>%
  mutate(id = str_c(group, "Test vs Control", sep = "-")) %>%
  group_by(id) %>%
  {bind_rows(
    slice_max(., avg_log2FC, n = 5),
    slice_min(., avg_log2FC, n = 5) 
  )} %>%
  group_by(id) %>%
  nest() %>%
  deframe() %>%
  map('gene') 


# export
pdf(
  file.path(out_dir_0107, 'marker gene - brain area comparison - Test vs Control - dotplot.pdf'),
  height = 6,
  width  = 10
)
marker_gene_ls3 %>% # .[1] %>%
  imap(function(marker_gene, name){ # browser()
    
    # name: sp5_Pvalb-Test vs Control
    name2 <- name %>% str_split_1("-")
    name3 <- name2[1] %>% str_split_1("_")
    name4 <- name2[2] %>% str_split_1("vs") %>% map_chr(str_trim)
    
    DotPlot(
      all_data_roi_only_one[,
                            all_data_roi_only_one$group2 %in% name3[1] &
                            all_data_roi_only_one$group3 %in% name4 &
                              all_data_roi_only_one$cells %in% name3[2]
      ], 
      features = marker_gene %>% unique,
      group.by = 'group3'
    ) +
      ggtitle(name) +
      RotatedAxis()
  })

dev.off()









