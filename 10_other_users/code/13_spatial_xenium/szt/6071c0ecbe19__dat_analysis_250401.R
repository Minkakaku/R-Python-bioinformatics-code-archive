library(tidyverse)
library(Seurat)

# 需求： 20250305肿瘤细胞分析任务.docx
# 0. ----------------
out_dir <- "dat_analysis_cycle3/result_20250401" 
dir.create(out_dir, recursive = TRUE)

cn_out_dir <- "dat_analysis_cycle3/cn_analysis_250219"
# dir.create(cn_out_dir)


all_data <- readRDS(file.path("dat_analysis_cycle3", "subtype", "all_data_cycle3_250219.rds"))
all_data_meta_cn <- readRDS(file.path(cn_out_dir, "meta_data_cn.rds"))


# subset samples
group_A <- c(
  '13638',
  '10858',
  '47983A3'
)

group_B <- c(
  '45827',
  '72696',
  '58198',
  '71402A6',
  '06938F1',
  '15201H1',
  '28519G1',
  '44994A1',
  '43158G1',
  '67718G1'
)

stopifnot( all(group_A %in% all_data$project) )
stopifnot( all(group_B %in% all_data$project) )

all_data_subset <- all_data[, all_data$project %in% c(group_A, group_B)]
all_data_meta_cn_subset <- all_data_meta_cn[all_data$project %in% c(group_A, group_B), ]

stopifnot(all(colnames(all_data_subset) == all_data_meta_cn_subset$cell_id))

# 1. tumor cells -----------------

## 1.1 A vs B - epithelials ---------
de_AvsB_epi <- 
  FindMarkers(
    all_data_subset,
    # method = "DESeq2",
    ident.1 = "A",
    ident.2 = "B",
    group.by = "group",
    subset.ident = "Epithelial",
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0
  )

p_val_adj_min <- de_AvsB_epi$p_val_adj %>% .[. != 0] %>% min
de_AvsB_epi$p_val_adj[de_AvsB_epi$p_val_adj == 0] <- p_val_adj_min/10
  
de_AvsB_epi_tidy <-
  de_AvsB_epi %>%
  rownames_to_column("gene") %>% 
  mutate(group = case_when(
    avg_log2FC >  1 & p_val_adj < 0.05 ~ "Pos",
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "Neg",
    TRUE ~ "Non_Sig"
  )) %>%
  mutate(
    `-log10(p.value)` = -log10(p_val_adj)
  ) 


de_AvsB_epi_tidy_label_up <- 
  de_AvsB_epi_tidy %>%
  dplyr::filter(group == "Pos") %>%
  dplyr::slice_min(p_val, n = 10) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10)


de_AvsB_epi_tidy_label_down <- 
  de_AvsB_epi_tidy %>%
  dplyr::filter(group == "Neg") %>%
  dplyr::slice_min(p_val, n = 10) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) 



p_AvsB_epi_add_text <-
  de_AvsB_epi_tidy %>%
  ggplot(aes(x = avg_log2FC, y = `-log10(p.value)`, color = group)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  
  ggrepel::geom_text_repel(
    data = de_AvsB_epi_tidy_label_up,
    aes(label = gene),
    max.overlaps = 100,
    color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = de_AvsB_epi_tidy_label_down,
    aes(label = gene),
    max.overlaps = 100,
    color = "black"
  ) +
  
  scale_color_manual(values = c(
    Neg = "blue",
    Pos = "red",
    Non_Sig = "gray"
  )) + 
  theme_bw() +
  theme(
    legend.position = "none"
  )

# export
pdf(
  file = file.path(out_dir, "1.1 volcano plot of A vs B of Epithelial cells.pdf"),
  height = 5,
  width  = 8
)
p_AvsB_epi_add_text %>% print()
dev.off()


de_AvsB_epi %>% 
  rownames_to_column("gene") %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "1.1 de analysis table of A vs B of Epithelial cells.xlsx"),
    asTable = TRUE
  )

## 1.2 A vs B - one CN ---------
all_de_cn_ls <-
  all_data_meta_cn_subset$CN2 %>% 
  unique() %>% 
  as.character() %>% 
  set_names(., .) %>%
  map(function(cn){
    # 比较AB组间同一CN中肿瘤细胞的基因表达差异
    
    cell_id <- 
      all_data_meta_cn_subset %>%
      dplyr::filter(annotation5 == 'Epithelial' & CN2 == cn) %>%
      dplyr::pull("cell_id")
    
    de <- 
      FindMarkers(
        all_data_subset[, cell_id],
        # method = "DESeq2",
        ident.1 = "A",
        ident.2 = "B",
        only.pos = FALSE,
        group.by = "group",
        min.pct = 0,
        logfc.threshold = 0
      )
    
    p_val_adj_min <- de$p_val_adj %>% .[. != 0] %>% min
    de$p_val_adj[de$p_val_adj == 0] <- p_val_adj_min/10
    
    de
  })


p_AvsB_epi_all_cn_ls <-
  all_de_cn_ls %>%
  imap(function(dat, name){
    dat_tidy <-
      dat %>%
      rownames_to_column("gene") %>% 
      mutate(group = case_when(
        avg_log2FC >  1 & p_val_adj < 0.05 ~ "Pos",
        avg_log2FC < -1 & p_val_adj < 0.05 ~ "Neg",
        TRUE ~ "Non_Sig"
      )) %>%
      mutate(
        `-log10(p.value)` = -log10(p_val_adj)
      ) 
    
    dat_tidy_label_up <- 
      de_AvsB_epi_tidy %>%
      dplyr::filter(group == "Pos") %>%
      dplyr::slice_min(p_val, n = 10) %>%
      dplyr::slice_max(abs(avg_log2FC), n = 10)
    
    
    dat_tidy_label_down <- 
      de_AvsB_epi_tidy %>%
      dplyr::filter(group == "Neg") %>%
      dplyr::slice_min(p_val, n = 10) %>%
      dplyr::slice_max(abs(avg_log2FC), n = 10) 
    
    dat_tidy %>% 
      ggplot(aes(x = avg_log2FC, y = `-log10(p.value)`, color = group)) +
      geom_point() +
      geom_vline(xintercept = c(-1, 1), linetype = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) +
      ggrepel::geom_text_repel(
        data = dat_tidy_label_up,
        aes(label = gene),
        max.overlaps = 100,
        color = "black"
      ) +
      ggrepel::geom_text_repel(
        data = dat_tidy_label_down,
        aes(label = gene),
        max.overlaps = 100,
        color = "black"
      ) +
      scale_color_manual(values = c(
        Neg = "blue",
        Pos = "red",
        Non_Sig = "gray"
      )) + 
      theme_bw() +
      theme(
        legend.position = "none"
      ) +
      labs(
        title = name
      )
  })



# export
pdf(
  file = file.path(out_dir, "1.2 volcano plot of A vs B of Epithelial cells in each CN.pdf"),
  height = 5,
  width  = 6
)
p_AvsB_epi_all_cn_ls %>% walk(print)
dev.off()

all_de_cn_ls %>% 
  map(rownames_to_column, "gene") %>% 
  openxlsx::write.xlsx(
    file = file.path(out_dir, "1.2 de analysis table of A vs B of Epithelial cells in each CN.xlsx"),
    asTable = TRUE
  )



## 1.3 pairwise CN - A and B ---------
all_de_pairwise_cn_ls <-
  all_data_meta_cn_subset$group %>% 
  unique() %>% 
  as.character() %>% 
  set_names(., .) %>%
  map(function(grp){ # browser()
    
    # 比较AB组内不同CN中肿瘤细胞的基因表达差异
    cell_idx <- all_data_meta_cn_subset$group == grp & all_data_meta_cn_subset$annotation5 == "Epithelial"
    
    all_data_sub <- all_data_subset[, cell_idx]
    meta_data    <- all_data_meta_cn_subset[cell_idx, ]
    
    all_data_sub$CN2 <- meta_data$CN2
    
    Idents(all_data_sub) <- "CN2"
    
    de <- 
      FindAllMarkers(
        all_data_sub,
        only.pos = FALSE,
        min.pct = 0,
        return.thresh = 1.1,
        logfc.threshold = 0
      )
    
    de 
  })

p_pwc_cn_epi_all_A_B_ls <-
  all_de_pairwise_cn_ls %>%
  imap(function(dat, name){ # browser()
    # dat
    # p_val avg_log2FC pct.1 pct.2     p_val_adj cluster       gene
    # MTRNR2L12  1.259519e-114  2.0699800 0.421 0.139 3.364427e-110    CD-F  MTRNR2L12
    # FP236383.3 1.647682e-110  0.5967926 1.000 1.000 4.401288e-106    CD-F FP236383.3
    # FP236383.5  4.720774e-79  0.9386205 0.750 0.558  1.261013e-74    CD-F FP236383.5
    # FP236383.4  8.018641e-69  0.9385583 0.673 0.449  2.141939e-64    CD-F FP236383.4
    # FP671120.7  3.710999e-63  0.6419110 0.977 0.949  9.912819e-59    CD-F FP671120.7
    # SLCO5A1     1.146982e-47  0.9227888 0.460 0.258  3.063819e-43    CD-F    SLCO5A1
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    # scatter
    df_scatter <-
      dat %>%
      mutate(cluster = factor(cluster), x_base = as.numeric(cluster)) %>%
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
    
    
    # text label
    df_text_grp <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(x = x_base - 1) %>%
      mutate(x = (x_base + x)/2)
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_text_grp,
        aes(x = x, y = 0, label = cluster)
      ) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      ggsci::scale_fill_igv() +
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



# export
pdf(
  file = file.path(out_dir, "1.3 multi volcano plot of epithelial cells of pairwise comparsion of CN in each A or B.pdf"),
  height = 7,
  width  = 18
)
p_pwc_cn_epi_all_A_B_ls %>% walk(print)
dev.off()

all_de_pairwise_cn_ls %>% 
  # map(rownames_to_column, "gene") %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir, "1.3 de analysis table of pairwise comparsion of CN of Epithelial cells in A or B.xlsx"),
    asTable = TRUE
  )




# 2. 肿瘤细胞各种通路打分 --------------
pathway_ls <- 
  c(
    "Pathway_Gene_Intersection_gmt.xlsx",
    "Pathway_Gene_Intersection_modified.xlsx"
  ) %>%
  map(function(x){
    sheets <- readxl::excel_sheets(x) %>% setNames(., .)
    map(sheets, ~ readxl::read_excel(x, sheet = .x))
  })

pathway_tidy_ls <-
  pathway_ls %>%
  map(~ .x %>% map("gene")) %>%
  purrr::list_c()

names(pathway_tidy_ls) %in% colnames( all_data_subset[[]] )

names(pathway_tidy_ls) <- c(
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING',
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
  'KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_I_MOLECULES',
  'GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL',
  'HALLMARK_ANGIOGENESIS',
  'GOBP_EPITHELIAL_CELL_PROLIFERATION',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE',
  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
  
  "inflammation",
  "metastasis",
  "HLA",
  "immune_checkpoint",
  "proliferation",
  "IFN" 
)

pathway_grp <- list(
  gmt      = names(pathway_tidy_ls) %>% str_subset("^(HALLMARK|KEGG|GO)"),
  modified = names(pathway_tidy_ls) %>% str_subset("^(HALLMARK|KEGG|GO)", negate = TRUE)
)



all_data_subset_score <-
  all_data_subset %>%
  AddModuleScore(
    features = pathway_tidy_ls,
    nbin = 5,
    ctrl = 50
  )

all_data_meta_subset_score <-
  all_data_subset_score[[]] %>%
  dplyr::select(!all_of(c(
    'HLA', 'tumor_cell_LN_metastasis', 
    'IFN', 'cancer_stem_cell_and_EMT', 
    'Inflammation'
  ))) %>%
  dplyr::rename(!!! structure(
    str_c("Cluster", seq_along(pathway_tidy_ls), sep = ""),
    names = names(pathway_tidy_ls)
  ))



## 2.1 pathway score A vs B - epithelials ----------
path_de_df <- 
  names(pathway_tidy_ls) %>%
  set_names(., .) %>%
  map(function(pathway){
    dat <- 
      all_data_meta_subset_score %>%
      dplyr::filter(annotation5 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group')] 
    
    dat <- dat %>%
      group_by(project, group) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    pathway_a <- dat[[pathway]][dat$group == "A"] 
    pathway_b <- dat[[pathway]][dat$group == "B"] 
    
    delta <- mean(pathway_a) - mean(pathway_b)
    wtest <- wilcox.test(pathway_a, pathway_b)
    
    w_stat <- wtest$statistic
    pvalue <- wtest$p.value
    
    data.frame(
      statistic = w_stat,
      pvalue = pvalue,
      delta = delta
    )
  }) %>%
  enframe("pathway") %>%
  unnest(value)


# plot
p_path_de_df_ls <-
  pathway_grp %>%
  map(function(path_grp){
    path_de_df %>%
      dplyr::filter(pathway %in% path_grp) %>%
      # mutate(pvalue = ifelse(pvalue < 1e-30, 0, pvalue)) %>%
      mutate(pathway = reorder(pathway, delta)) %>%
      ggplot(aes(y = pathway, x = delta, fill = statistic, size = pvalue)) +
      geom_point(shape = 21, color = "black") +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds")) +
      ggpubr::theme_pubclean() +
      guides(
        color = guide_colorbar(theme = theme(
          legend.key.width  = unit(10, "lines") 
        )),
        size = guide_legend(reverse = TRUE)
      ) +
      scale_size(range = c(2, 6), transform = "reverse") + 
      theme(
        # legend.spacing = unit(2, "lines")
        # legend.position = "right"
      ) +
      labs(
        fill = "wilcox statistic",
        x = "delta (pathway score A - pathway score  B)",
        y = NULL
      )
    
  })
  
ggsave(
  p_path_de_df_ls$gmt,
  filename = file.path(out_dir, '2.1 dotplot plot of pathway score - A vs B - epithelial - gmt.pdf'),
  height = 5,
  width  = 15
)
ggsave(
  p_path_de_df_ls$modified,
  filename = file.path(out_dir, '2.1 dotplot plot of pathway score - A vs B - epithelial - modified.pdf'),
  height = 5,
  width  = 7
)

## 2.1 pathway score A vs B - epithelial - boxplot -----

p_path_de_df_box_ls  <-
  names(pathway_tidy_ls) %>%
  set_names(., .) %>%
  # .[1] %>%
  map(function(pathway){
    dat <- 
      all_data_meta_subset_score %>%
      dplyr::filter(annotation5 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group')] 
    
    dat <- dat %>%
      group_by(project, group) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    set.seed(1234)
    dat %>%
      ggplot(aes(x = group, y = .data[[pathway]])) +
      geom_boxplot(width = 0.35, outlier.colour = NA) +
      geom_point(position = position_jitter(width = 0.15), color = "blue") +
      ggpubr::stat_compare_means(
        label = 'p.format',
        comparisons = list(c("A", "B")),
        tip.length = 0) +
      theme_classic()
  }) 


pdf(
  file = file.path(out_dir, '2.1 boxplot plot of pathway score - A vs B - epithelial.pdf'),
  height = 5,
  width  = 7
)
p_path_de_df_box_ls %>% walk(print)
dev.off()

## 2.2 A vs B - one CN ---------


path_de_df_cn <- 
  names(pathway_tidy_ls) %>%
  set_names(., .) %>%
  map(function(pathway){ # browser()
    stopifnot(all(rownames(all_data_meta_subset_score) == all_data_meta_cn_subset$cell_id))
    dat <- 
      all_data_meta_subset_score %>%
      cbind(all_data_meta_cn_subset['CN2']) %>% 
      dplyr::filter(annotation5 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group', 'CN2')]
    
    dat2 <-
      dat %>%
      group_by(project, group, CN2) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    # print(head(dat2))
    
    dat2 %>%
      group_by(CN2) %>%
      nest() %>%
      deframe() %>%
      map(function(x){
        pathway_a <- x[[pathway]][x$group == "A"] 
        pathway_b <- x[[pathway]][x$group == "B"] 
        
        delta <- mean(pathway_a) - mean(pathway_b)
        ttest <- wilcox.test(pathway_a, pathway_b)
        
        t_stat <- ttest$statistic
        pvalue <- ttest$p.value
        
        data.frame(
          statistic = t_stat,
          pvalue = pvalue,
          delta = delta
        )
      }) %>%
      enframe("CN2") %>%
      unnest(value)
    
  }) %>%
  enframe("pathway") %>%
  unnest(value)


# plot
p_path_de_df_cn_ls <-
  pathway_grp %>%
  map(function(path_grp){
    p_path_de_df_cn <-
      path_de_df_cn %>%
      
      dplyr::filter(pathway %in% path_grp) %>%
      
      mutate(pathway = factor(pathway, levels = path_de_df_cn %>% group_by(pathway) %>% summarise(s = sum(delta)) %>% arrange(s) %>% .[[1]])) %>%
      ggplot(aes(y = pathway, x = delta, fill = statistic, size = pvalue)) +
      geom_point(shape = 21, color = "black") +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds")) + 
      facet_wrap(~ CN2) +
      ggpubr::theme_pubclean() +
      guides(
        fill = guide_colorbar(theme = theme(
          legend.key.width  = unit(10, "lines") 
        )),
        size = guide_legend(reverse = TRUE)
      ) +
      scale_size(range = c(2, 6), transform = "reverse") + 
      theme(
        # legend.spacing = unit(2, "lines")
      ) +
      labs(
        fill = "wilcox statistic",
        x = "delta (pathway score A - pathway score  B)",
        y = NULL
      )
    
    p_path_de_df_cn
  })


ggsave(
  p_path_de_df_cn_ls$gmt,
  filename = file.path(out_dir, '2.2 dotplot plot of pathway score - A vs B - split by CN - epithelial - gmt.pdf'),
  height = 10,
  width  = 18
)
ggsave(
  p_path_de_df_cn_ls$modified,
  filename = file.path(out_dir, '2.2 dotplot plot of pathway score - A vs B - split by CN - epithelial - modified.pdf'),
  height = 8,
  width  = 10
)

## 2.2 A vs B - one CN - boxplot -----
p_path_de_df_cn_box_ls <-
  names(pathway_tidy_ls) %>%
  set_names(., .) %>% 
  # .[1] %>%
  map(function(pathway){ # browser()
    stopifnot(all(rownames(all_data_meta_subset_score) == all_data_meta_cn_subset$cell_id))
    dat <- 
      all_data_meta_subset_score %>%
      cbind(all_data_meta_cn_subset['CN2']) %>% 
      dplyr::filter(annotation5 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group', 'CN2')]
    
    dat2 <-
      dat %>%
      group_by(project, group, CN2) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    set.seed(1234)
    dat2 %>%
      ggplot(aes(x = CN2, y = .data[[pathway]], fill = group)) +
      geom_boxplot(width = 0.35, outlier.colour = NA) +
      geom_point(position = position_dodge( width = 0.35), size = 1) +
      ggpubr::geom_pwc(
        label = 'p.format',
        # comparisons = list(c("A", "B")),
        tip.length = 0) +
      theme_classic()
  })



pdf(
  file = file.path(out_dir, '2.2 boxplot plot of pathway score - A vs B - split by CN - epithelial.pdf'),
  height = 5,
  width  = 13
)
p_path_de_df_cn_box_ls %>% walk(print)
dev.off()

## 2.3 pairwise CN - A and B ---------

path_de_df_group_ls <- 
  names(pathway_tidy_ls) %>%
  set_names(., .) %>%
  map(function(pathway){
    stopifnot(all(rownames(all_data_meta_subset_score) == all_data_meta_cn_subset$cell_id))
    dat <- 
      all_data_meta_subset_score %>%
      cbind(all_data_meta_cn_subset['CN2']) %>% 
      dplyr::filter(annotation5 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group', 'CN2')]
    
    dat2 <-
      dat %>%
      group_by(project, group, CN2) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    dat2 %>%
      group_by(group) %>%
      nest() %>%
      deframe() %>%
      map(function(x){ # browser()
        x$CN2 %>%
          unique() %>%
          set_names(., .) %>% 
          map(function(cn){
            pathway_a <- x[[pathway]][x$CN2 == cn] 
            pathway_b <- x[[pathway]][x$CN2 != cn] 
            
            delta <- mean(pathway_a) - mean(pathway_b)
            ttest <- wilcox.test(pathway_a, pathway_b)
            
            t_stat <- ttest$statistic
            pvalue <- ttest$p.value
            
            data.frame(
              statistic = t_stat,
              pvalue = pvalue,
              delta = delta
            )
          }) %>%
          enframe("CN2") %>%
          unnest(value)
        
      }) %>%
      enframe("group") %>%
      unnest(value)
  })  


# plot
p_path_de_df_group_ls <-
  path_de_df_group_ls %>%
  imap(function(dat, name){
    dat %>%
      ggplot(aes(y = CN2, x = delta, fill = statistic, size = pvalue)) + 
      geom_point(shape = 21, color = "black") +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds")) + 
      facet_wrap(~ group) +
      ggpubr::theme_pubclean() +
      guides(
        fill = guide_colorbar(theme = theme(
          legend.key.width  = unit(10, "lines") 
        )),
        size = guide_legend(reverse = TRUE)
      ) +
      scale_size(range = c(2, 6), transform = "reverse") + 
      theme(
        # legend.spacing = unit(2, "lines")
      ) +
      labs(
        title = name, 
        fill  = "wilcox statistic",
        x = "delta (pathway score of CN - average pathway score of Others CNs)",
        y = NULL
      )
  })

# 
pdf(
  file = file.path(out_dir, '2.3 dotplot plot of pathway score - CN - split by group - epithelial.pdf'),
  height = 6,
  width  = 15
)
p_path_de_df_group_ls %>% walk(print)
dev.off()

## 2.3 pairwise CN - A and B - boxplot -----

p_path_de_df_group_box_ls <-
  names(pathway_tidy_ls) %>%
  set_names(., .) %>% 
  # .[1] %>%
  map(function(pathway){ # browser()
    stopifnot(all(rownames(all_data_meta_subset_score) == all_data_meta_cn_subset$cell_id))
    dat <- 
      all_data_meta_subset_score %>%
      cbind(all_data_meta_cn_subset['CN2']) %>% 
      dplyr::filter(annotation5 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group', 'CN2')]
    
    dat2 <-
      dat %>%
      group_by(project, group, CN2) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    set.seed(1234)
    dat2 %>%
      ggplot(aes(x = group, y = .data[[pathway]], fill = CN2)) +
      geom_boxplot(width = 0.35, outlier.colour = NA) +
      geom_point(position = position_dodge( width = 0.35), size = 1) +
      # ggpubr::geom_pwc(
      #   label = 'p.format',
      #   hide.ns = TRUE,
      #   # remove.bracket = TRUE,
      #   tip.length = 0) +
      theme_classic()
  })



pdf(
  file = file.path(out_dir, '2.3 boxplot plot of pathway score - CN - split by group - epithelial.pdf'),
  height = 5,
  width  = 10
)
p_path_de_df_group_box_ls %>% walk(print)
dev.off()

# # 3. 各基因对该通路打分的contribution ------
# pathway_tidy_ls %>%
#   imap(function(pathway, name){
#     all_data_subset
#   })
# 
# all_data_subset[unique(unlist(pathway_tidy_ls)), all_data_subset$annotation5 == "Epithelial"] %>%
#   { cbind(LayerData(., layer = "data"), .[['group']]) } %>%
#   group_by(group) %>%
#   summ





