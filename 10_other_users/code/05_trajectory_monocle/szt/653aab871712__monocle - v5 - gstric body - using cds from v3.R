detach("package:monocle3")

library(Seurat)
library(tidyverse)
library(ggpubr)
library(monocle)
source("monocle/monocle fix bugs.R")
source("monocle/color_pallete.R")

all_data_full <- readRDS('all_data_full_v6.rds')



# set default ------- 
monocle_out_dir <- "monocle - v5 - gastric body" 
dir.create(monocle_out_dir) 

# cds --------
cds_body <- readRDS("monocle - v3 - gastric body/cds.rds")

meta_cds_body <- 
  data.frame(cell_id = colnames(cds_body)) %>%
  left_join(
    rownames_to_column(all_data_full[[c('sub_celltype', 'sub_celltype3')]], "cell_id"),
    by = 'cell_id'
  )
table(is.na(meta_cds_body$sub_celltype))
stopifnot(all(meta_cds_body$cell_id == colnames(cds_body)))

pData(cds_body)$sub_celltype  <- meta_cds_body$sub_celltype %>% fct_drop()
pData(cds_body)$sub_celltype3 <- meta_cds_body$sub_celltype3 %>% fct_drop()
pData(cds_body)$annotation_full <- pData(cds_body)$sub_celltype3

remove_idx_body <- (!is.na(cds_body$sub_celltype3)) & (cds_body$sub_celltype3 != "endocrine_cells")


# select origin -- antrm - PMC

plot_cell_trajectory(cds_body, color_by="sub_celltype3", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds_body, color_by="State", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds_body, color_by="Pseudotime", size=1,show_backbone=TRUE) 

cds_body <- orderCells(cds_body, root_state = 3)
saveRDS(cds_body, file = file.path(monocle_out_dir, "cds2.rds"))




# plot ---------- 
color_used_temp <- c(
  enterocytes   = "red",
  goblet_cells  = "red",
  
  # 'TUMOR+MKI67'   = "#0000FF33",
  'tumor+S100A6'  = "#0000FF88",
  'tumor+OLFM4'  = "#0000FFDD",
  
  enteroendocrine_cells  = 'yellow',
  endocrine_cell         = 'yellow',
  
  gmc                    = "#00FF0033",
  pmc                    = "#00FF0044",
  parietal_cells         = "#00FF00aa",
  chief_cells            = "#00FF00dd"
)



cds_body %>% 
  Biobase::pData() %>%
  # mutate(annotation_full = Location) %>%
  dplyr::select(annotation_full, State) %>%
  group_by(annotation_full, State) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_full, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )

cds_body %>% 
  Biobase::pData() %>%
  # mutate(annotation_full = Location) %>%
  dplyr::select(annotation_full, State) %>%
  group_by(State, annotation_full) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_full)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )


library(ggpubr)

## density plot of annotation  ##
pdf(
  file.path(monocle_out_dir, "cell density.pdf"),
  height = 6,
  width = 18
)
pData(cds_body) %>%
  .[remove_idx_body, ] %>% 
  # dplyr::filter(! State %in% c("5")) %>%
  ggplot(aes(x = Pseudotime, y = annotation_full, fill = annotation_full, color = annotation_full)) +
  ggridges::geom_density_ridges(panel_scaling = FALSE) +
  ggtitle('IM') +
  ggsci::scale_color_igv() +
  ggsci::scale_fill_igv() +
  ggridges::theme_ridges() +
  
  
  pData(cds_body)[remove_idx_body, ] %>%
  # dplyr::filter(! State %in% c("5")) %>%
  ggplot(aes(x = Pseudotime, y = annotation_full, color = annotation_full)) +
  geom_jitter(width = 0, height = 0.2)+
  ggridges::theme_ridges() +
  ggsci::scale_color_igv() 
dev.off()


 
## density plot of group ##
pdf(
  file.path(monocle_out_dir, "IM cell density - group.pdf"),
  height = 6,
  width = 18
)
pData(cds_body) %>%
  .[remove_idx_body, ] %>% 
  # dplyr::filter(! State %in% c("5")) %>%
  ggplot(aes(x = Pseudotime, y = group, fill = group, color = group)) +
  ggridges::geom_density_ridges(panel_scaling = FALSE) +
  ggtitle('IM') +
  ggsci::scale_color_igv() +
  ggsci::scale_fill_igv() +
  ggridges::theme_ridges() +
  
  
  pData(cds_body)[remove_idx_body, ] %>%
  # dplyr::filter(! State %in% c("5")) %>%
  ggplot(aes(x = Pseudotime, y = group, color = group)) +
  geom_jitter(width = 0, height = 0.2)+
  ggridges::theme_ridges() +
  ggsci::scale_color_igv() 
dev.off()


 

### 
pdf(
  file.path(monocle_out_dir, "barplot of State group by annotation.pdf"),
  height = 6,
  width = 6
)
cds_body %>% 
  Biobase::pData() %>%
  .[remove_idx_body, ] %>% 
  # mutate(annotation_full = Location) %>%
  dplyr::select(annotation_full, State) %>%
  group_by(annotation_full, State) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_full, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
dev.off()


pdf(
  file.path(monocle_out_dir, "barplot of annotation group by State.pdf"),
  height = 6,
  width  = 6
)
cds_body %>% 
  Biobase::pData() %>%
  .[remove_idx_body, ] %>% 
  # mutate(annotation_full = Location) %>%
  dplyr::select(annotation_full, State) %>%
  group_by(State, annotation_full) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_full)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory.pdf"), width = 8, height = 4)

plot_cell_trajectory(cds_body,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds_body,color_by="group", cell_size =0.1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

p <- plot_cell_trajectory(cds_body,color_by="annotation_full", cell_size =0.1, show_backbone=TRUE) + 
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))
p$data <- p$data[remove_idx_body, ]
p
plot_cell_trajectory(cds_body, color_by = "State",cell_size =1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 4)
local({
  p <- plot_cell_trajectory(cds_body, color_by = "annotation_full",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_1) +
    facet_wrap(~ group) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
  p$data <- p$data[remove_idx_body, ]
  p
})
dev.off()



pdf(file.path(monocle_out_dir, "cell trajectory - tree.pdf"), width = 12, height = 8)
local({
  p <- 
    plot_complex_cell_trajectory(cds_body, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_full") +
    scale_color_manual(values = col_plan_1) +
    facet_wrap(~ annotation_full)
  p$data <- p$data[remove_idx_body, ]
  p
})
dev.off()


pdf(file.path(monocle_out_dir, "barplot tof cell props group by State - group.pdf"), width = 8, height = 5)
cds_body %>%
  Biobase::pData() %>%
  .[remove_idx_body, ] %>% 
  group_by(State, project) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = project))
dev.off()

pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation.pdf"), width = 8, height = 5)
cds_body %>%
  Biobase::pData() %>%
  .[remove_idx_body, ] %>% 
  group_by(State, annotation_full ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_full ))
dev.off()


# cell barplot 
p_col_cell_state <-
  cds_body %>% 
  Biobase::pData() %>%
  .[remove_idx_body, ] %>%
  dplyr::select(annotation_full, State) %>%
  group_by(annotation_full, State) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_full, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_out_dir, "cell stat barplot.pdf"), height = 4, width = 8)

p_col_cell_state2 <-
  cds_body %>% 
  Biobase::pData() %>%
  .[remove_idx_body, ] %>% 
  dplyr::select(annotation_full, State) %>%
  group_by(State, annotation_full) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_full)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state2, filename = file.path(monocle_out_dir, "cell stat barplot2.pdf"), height = 4, width = 8)


# plot pseudotime by geneset signature -----------------
meta_geneset_score  <- readRDS('data/adata_gensets_score.rds')
stopifnot(all(colnames(cds_body) %in% meta_geneset_score$cell_id))
table(colnames(cds_body) %in% meta_geneset_score$cell_id)

meta_geneset_score2 <- 
  data.frame(cell_id = colnames(cds_body)) %>%
  left_join(meta_geneset_score, by = "cell_id")

colnames(meta_geneset_score2) <- meta_geneset_score2 %>%
  colnames() %>% 
  str_trim() %>%
  make.names()

cds_body2 <- cds_body
pData(cds_body2) <- cbind(
  pData(cds_body2),
  meta_geneset_score2[-1]
)

p_pseudo_geneset_ls <- 
  meta_geneset_score2 %>%
  colnames() %>%
  .[-1] %>%
  setNames(., .) %>% # .[1] %>%
  map(function(x){
    p <- plot_cell_trajectory(
      cds_body2,
      color_by=x, 
      size=1,
      show_backbone=TRUE
    ) 
    p$data <- p$data %>% .[!is.na(.[[x]]), ]
    p$layers <- p$layers[c(2, 1, 3, 4)]
    
    
    p2 <- p
    p2$data <- arrange(p2$data, .data[[x]])
    
    p <- p + ggtitle('original')
    p2 <- p2 + ggtitle('small value located in background')
    
    p + p2
  })

# export 
pdf(file.path(monocle_out_dir, "cell trajectory - geneset score - v2 .pdf"), width = 10, height = 4)
p_pseudo_geneset_ls %>%  walk(print)
dev.off()


# plot smooth pseudotime by geneset signature -----------------
meta_geneset_score  <- readRDS('data/adata_gensets_score.rds')
stopifnot(all(colnames(cds_body) %in% meta_geneset_score$cell_id))
table(colnames(cds_body) %in% meta_geneset_score$cell_id)

meta_geneset_score2 <- 
  data.frame(cell_id = colnames(cds_body)) %>%
  left_join(meta_geneset_score, by = "cell_id")

colnames(meta_geneset_score2) <- meta_geneset_score2 %>%
  colnames() %>% 
  str_trim() %>%
  make.names()

cds_body2 <- cds_body
pData(cds_body2) <- cbind(
  pData(cds_body2),
  meta_geneset_score2[-1]
)

p_pseudo_smooth_geneset_ls <-
  meta_geneset_score2 %>%
  colnames() %>%
  .[-1] %>%
  setNames(., .) %>% # .['Oncogenic.TFs'] %>%
  map(function(x){
    dat <- pData(cds_body2) %>% slice_sample(n = 5000) %>% dplyr::filter(!is.na(.data[[x]]))
    
    cor_val <- cor(dat[[x]], dat$Pseudotime)
    rng_path <- dat[[x]] %>% pmin(quantile(., 0.99)) %>% pmax(quantile(., 0.01)) %>% range(na.rm = T)
    # rng_path <- loess(as.formula(str_glue("Pseudotime ~ {x}")), data = dat) %>% .$fitted %>% range(na.rm = T)
    
    print(cor_val)
    print(rng_path)
    
    set.seed(1234)
    dat %>%
      ggplot(aes(x = Pseudotime, y = .data[[x]])) +
      geom_smooth(method = "loess", size = 0.8, color = "black") +
      ggpubr::stat_cor(
        label.sep = "\n",
        label.y = ifelse(
          cor_val > 0, 
          rng_path[2] - 0.1*(diff(rng_path)), 
          rng_path[1] + 0.1*(diff(rng_path))
        ),
        size = 4,
        vjust = ifelse(cor_val > 0, 1, 0)
      ) +
      theme_classic() 
  })

# export 
pdf(file.path(monocle_out_dir, "cell trajectory - geneset score - smooth line.pdf"), width = 6, height = 4)
p_pseudo_smooth_geneset_ls %>%  walk(print)
dev.off()


# # de analysis -----------
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性

Time_diff <- differentialGeneTest(
  cds_body[rownames(fData(cds_body))[fData(cds_body)$use_for_ordering], ],
  cores = 1,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
saveRDS(Time_diff, file.path(monocle_out_dir, "Time_diff.rds"))
# Time_diff <- readRDS(file.path(monocle_out_dir, "Time_diff.rds"))

# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

# Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
Time_genes <- dplyr::filter(Time_diff, qval < 0.05) %>% pull(gene_short_name) %>% as.character()


p = plot_pseudotime_heatmap(cds_body[Time_genes,], num_clusters=5, show_rownames=F, return_heatmap=T)
ggsave(file.path(monocle_out_dir, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)


library(ClusterGVis)

p$data %>% .[1:5, 1:5]


newdata <- data.frame(
  Pseudotime = seq(
    min(pData(cds_body)$Pseudotime), 
    max(pData(cds_body)$Pseudotime), 
    length.out = 100
  )
)
m <- genSmoothCurves(
  cds_body[Time_genes,], 
  cores = 1, 
  trend_formula = "~sm.ns(Pseudotime, df=3)", 
  relative_expr = T, 
  new_data = newdata
)
pseu_time <- cds_body[Time_genes,]$Pseudotime %>% range

set.seed(1234)
ck <- clusterData(
  obj = m,
  clusterMethod = "kmeans",
  clusterNum = 4
)

enrich <- enrichCluster(object = ck,
                        OrgDb = 'org.Hs.eg.db',
                        type = "BP",
                        idTrans = T,
                        toType = "ENTREZID",
                        readable = F,
                        pvalueCutoff = 0.5,
                        topn = 5)
write.csv(enrich, file = file.path(monocle_out_dir, "enrichment table of heatmap of pseudotime.csv"))


# export cluster
ck$cluster.list %>%
  enframe('cluster', 'gene') %>%
  unnest(gene) %>%
  write.csv(
    file = file.path(monocle_out_dir, "gene in each cluster in pseudotime heatmap.csv"),
    row.names = FALSE
  )


# export heatmap
pdf(
  file.path(monocle_out_dir, "heatmap of pseudotime.pdf"),
  height = 8,
  width  = 12
)
visCluster(
  object = ck,
  plotType = "both",
  lineSide = "left",
  addSampleAnno = FALSE,
  annnoblockText = T,
  markGenes = ck$cluster.list %>% map(head) %>% unlist(),
  annoTermData = enrich %>% dplyr::select(id = group, term = Description, pval = pvalue),
  ctAnnoCol = ggsci::pal_npg()(4),
  # byGo = "anno_block",
  # heatmapAnnotation = ComplexHeatmap::HeatmapAnnotation(
  #   Pseudotime = structure(
  #     seq(pseu_time[1], pseu_time[2], length.out = 100),
  #     names = colnames(m)
  #   )
  #   # col = list(
  #   #   Pseudotime =
  #   # )
  # )
)
dev.off()


plot_genes_jitter(
  cds['CEACAM5', ],
  color_by = "group"
)

plot_genes_in_pseudotime(
  cds['FABP1', ],
  color_by = "group"
) +
  ggsci::scale_color_igv()


plot_genes_in_pseudotime(
  cds[c("CDC73"), ],
  color_by = "group"
) +
  ggsci::scale_color_igv()


plot_genes_in_pseudotime(
  cds[c("CDC73", "CDK13", "CCNK", "LEO1", "MED14", "VHL"), ],
  color_by = "group"
) +
  ggsci::scale_color_igv()

## C2
top_path_gene_c2 <-
  all_enrich_de_ls[['2']]@result %>% 
  slice_min(p.adjust, n = 10) %>%
  .$geneSymbol %>% 
  str_split(",") %>%
  unlist %>% 
  unique() %>% 
  sort %>% 
  walk(cat, " ")

gene_mean <- cds %>% exprs %>% rowMeans()
gene_mean[top_path_gene_c2] %>% sort(decreasing = TRUE) %>% names %>% dput

p_pseudotime_c2 <- 
  plot_genes_in_pseudotime(
    cds[c("ANXA2", "FABP1", "PTTG1IP", "DGAT1", "TSPO", "ROMO1"), ],
    color_by = "group"
  ) +
  ggsci::scale_color_igv()

p_pseudotime_c2_full <- 
  plot_genes_in_pseudotime(
    cds[top_path_gene_c2, ],
    color_by = "group"
  ) +
  ggsci::scale_color_igv()

ggsave(
  p_pseudotime_c2,
  filename = file.path(monocle_out_dir, "C2 scatter plot of pseudotime.pdf"),
  width = 8,
  height = 25
)

ggsave(
  p_pseudotime_c2_full,
  filename = file.path(monocle_out_dir, "C2 scatter plot of pseudotime - full.pdf"),
  width = 5,
  height = 50,
  limitsize = FALSE
)

## C1
top_path_gene_c1 <-
  all_enrich_de_ls[['1']]@result %>% 
  dplyr::filter(p.adjust < 0.05) %>%
  .$geneSymbol %>% 
  str_split(",") %>%
  unlist %>% 
  unique() %>% 
  sort %>% 
  walk(cat, " ")

gene_mean <- cds %>% exprs %>% rowMeans()
gene_mean[top_path_gene_c1] %>% sort(decreasing = TRUE) %>% names %>% dput


p_pseudotime_c1_full <- 
  plot_genes_in_pseudotime(
    cds[top_path_gene_c1, ],
    color_by = "group"
  ) +
  ggsci::scale_color_igv()

ggsave(
  p_pseudotime_c1_full,
  filename = file.path(monocle_out_dir, "C1 scatter plot of pseudotime.pdf"),
  width = 8,
  height = 25
)



# branch gene enrichment 
de_genes_ls <- 
  p$tree_row %>% cutree(k = 5) %>%
  split(names(.), .) 

# write gene
p$tree_row %>% 
  cutree(k = 5) %>%
  enframe("gene", "cluster") %>% 
  write.csv(file.path(monocle_out_dir, "de cluster gene.csv"))




de_gene_id_ls <-
  de_genes_ls %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Hs.eg.db"
    ) %>%
      .$ENTREZID
  })


all_enrich_de_ls <- 
  de_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2
    )
  })

# enrichment 
all_enrich_de_ls <-
  all_enrich_de_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })



# enrichment plot
p_all_enrich_de_ls <-
  all_enrich_de_ls %>%
  imap(function(enrich, name){
    if(nrow(as.data.frame(enrich)) < 1) return(NULL)
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 20,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export
all_enrich_de_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file.path(monocle_out_dir, "monocle de gene enrichment analysis table - full.xlsx"),
    asTable = TRUE
  )

pdf(
  file.path(monocle_out_dir, "monocle de gene enrichment analysis dotplot - full.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_de_ls %>%  walk(print)
dev.off()



# BEAM analysis ------------------
undebug(BEAM)
# gene_for_beam <- VariableFeatures(sc_epi2)
gene_for_beam <- rownames(fData(cds))[fData(cds)$use_for_ordering]
BEAM_res <- BEAM(cds[gene_for_beam,], branch_point = 2, cores = 10) 

saveRDS(BEAM_res, file.path(monocle_out_dir, "BEAM_res.rds"))

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res, file.path(monocle_out_dir, "BEAM_res.csv"), row.names = F)


#选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

undebug(plot_genes_branched_heatmap)
p <- 
  cds[BEAM_genes,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 4, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )
ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full.pdf"), p$ph_res, width = 6.5, height = 10)


# branch gene enrichment 
beam_genes_ls <- 
  p$ph_res$tree_row %>% cutree(k = 4) %>%
  split(names(.), .) 


beam_gene_id_ls <-
  beam_genes_ls %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Hs.eg.db"
    ) %>%
      .$ENTREZID
  })


all_enrich_beam_ls <- 
  beam_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2
    )
  })

# enrichment 
all_enrich_beam_ls <-
  all_enrich_beam_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    xx2 <-
      xx %>%
      map_chr(function(g){
        clusterProfiler::bitr(
          g,
          fromType = "ENTREZID",
          toType   = "SYMBOL",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$SYMBOL %>%
          str_c(collapse = "/")
      })
    x@result$geneSymbol <- xx2
    x
  })



# enrichment plot
p_all_enrich_beam_ls <-
  all_enrich_beam_ls %>%
  imap(function(enrich, name){
    if(nrow(as.data.frame(enrich)) < 1) return(NULL)
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export
all_enrich_beam_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file.path(monocle_out_dir, "monocle BEAM gene enrichment analysis table - full.xlsx"),
    asTable = TRUE
  )

pdf(
  file.path(monocle_out_dir, "monocle BEAM gene enrichment analysis dotplot - full.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_beam_ls %>%  walk(print)
dev.off()


