library(Seurat)
library(tidyverse)
library(monocle)
source("monocle/monocle fix bugs.R")
source("monocle/color_pallete.R")

all_data_full <- readRDS('all_data_full_v4.rds')

sc_epi <- all_data_full[, all_data_full$sub_celltype2 %in% c(
  'chief_cells',
  # 'endocrine_cells',
  # 'enteroendocrine_cells',
  'gmc',
  'parietal_cells',
  'pmc',
  'enterocytes',
  'goblet_cells',
  'tumor_cell'
)]
DefaultAssay(sc_epi) <- "Spatial.008um" 

# sc_epi$site %>% unique
# [1] distal_stomach Gastric_antrum body

# CN 
meta_data_cn <- readRDS("CN_analysis - v4 - remove late_cancer/all_data_meta_cn2.rds")
cn7_cell <- meta_data_cn$cell_id[meta_data_cn$CN2 == "CN7_Intestinal metaplasia"]

# QC
sc_epi2 <- sc_epi[, 
  sc_epi$nCount_Spatial.008um > 100 & 
    sc_epi$site %in% c('body') &
    (! sc_epi$group %in% c("TLS", "other")) &
    (colnames(sc_epi) %in% cn7_cell)
]

# factor 
sc_epi2$annotation_full <- sc_epi2$sub_celltype2 %>% fct_drop()

# set default ------- 
monocle_out_dir <- "monocle - v4 - gastric body - CN7" 
dir.create(monocle_out_dir) 


sc_epi2 <- FindVariableFeatures(sc_epi2, nfeatures = 3000)

saveRDS(sc_epi2, file = file.path(monocle_out_dir, "sc_epi2.rds"))

# cds --------

# subset
set.seed(1234)
meta_sc_epi2_subset <-
  sc_epi2[[]] %>%
  rownames_to_column("cell_id") %>%
  # group_by( annotation_full, project ) %>%
  # slice_sample(n = 200)
  slice_sample(n = 1e4)
# meta_sc_epi2_subset <- sc_epi2[['annotation_full']] %>% rownames_to_column("cell_id") 

table(meta_sc_epi2_subset$project, meta_sc_epi2_subset$annotation_full )


expr_dat_monocle <- as(LayerData(JoinLayers(sc_epi2), "counts"), "sparseMatrix")
p_data_monocle <- sc_epi2[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sc_epi2),
  row.names       = rownames(sc_epi2)
)

cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_sc_epi2_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_sc_epi2_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# saveRDS(cds, "cds.rds")


# pData(cds)$group <- all_data_meta_cn2$group[match(colnames(cds), all_data_meta_cn2$cell_id)]



##使用seurat选择的高变基因⚠️ ------
# express_genes_monocle <- VariableFeatures(sc_epi2)
# drop_gene <- dispersionTable(cds) %>% {.$gene_id[.$dispersion_empirical == 0]}
# express_genes_monocle <- express_genes_monocle %>% setdiff(drop_gene)
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)

# 
# disp_gene_eq_zero <- dispersionTable(cds) %>% {.$gene_id[.$dispersion_empirical == 0]}
# express_genes_monocle2 <- fData(cds) %>% {rownames(.)[.$use_for_ordering]} %>%  setdiff(disp_gene_eq_zero)
# cds <- setOrderingFilter(cds, express_genes_monocle2)
# plot_ordering_genes(cds)

# ##使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(JoinLayers(sc_epi2))
# express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)
# 
##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

# #这一步输入的expressed_genes来自于步骤4。
# fData(cds)$num_cells_expressed <- rowSums(exprs(cds) > 0)
# expressed_genes2 <- row.names(subset(fData(cds), num_cells_expressed >= 10))
# diff <- differentialGeneTest(cds[expressed_genes2, ],fullModelFormulaStr="~annotation",cores=1) 
# 
# head(diff)
# 
# ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
# deg <- subset(diff, qval < 0.01) #选出2829个基因
# deg <- deg[order(deg$qval,decreasing=F),]
# head(deg)
# 
# ##差异基因的结果文件保存
# write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
# 
# ## 轨迹构建基因可视化
# ordergene <- rownames(deg) 
# cds <- setOrderingFilter(cds, ordergene)  
# plot_ordering_genes(cds)


cds <- reduceDimension(cds)
cds <- orderCells(cds)

saveRDS(cds, file = file.path(monocle_out_dir, "cds.rds"))


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

plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds,color_by="annotation_full", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="group", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds,color_by="Location", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)

plot_cell_trajectory(cds,color_by="sub_celltype", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="sub_celltype", size=1,show_backbone=TRUE) + scale_color_manual(values = color_used_temp)
plot_cell_trajectory(cds,color_by="sub_celltype", size=1,show_backbone=TRUE) + facet_wrap(~ sub_celltype)
plot_cell_trajectory(cds,color_by="sub_celltype", size=0.1,show_backbone=TRUE) +
  facet_wrap(~ sub_celltype) +
  geom_bin2d(bins = 100) +
  scale_color_manual(values = col_plan_1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "OrRd"), transform = "log10")


plot_cell_trajectory(cds, color_by="Pseudotime", size=1,show_backbone=TRUE) 

# cds <- orderCells(cds, root_state = 3)

saveRDS(cds, file = file.path(monocle_out_dir, "cds2.rds"))

# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by="annotation_full", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by="group", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_full")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "group", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)



cds %>% 
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

cds %>% 
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
pData(cds) %>%
  # dplyr::filter(! State %in% c("5")) %>%
  ggplot(aes(x = Pseudotime, y = annotation_full, fill = annotation_full, color = annotation_full)) +
  ggridges::geom_density_ridges(panel_scaling = FALSE) +
  ggtitle('IM') +
  ggsci::scale_color_igv() +
  ggsci::scale_fill_igv() +
  ggridges::theme_ridges() +
  
  
  pData(cds) %>%
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
pData(cds) %>%
  # dplyr::filter(! State %in% c("5")) %>%
  ggplot(aes(x = Pseudotime, y = group, fill = group, color = group)) +
  ggridges::geom_density_ridges(panel_scaling = FALSE) +
  ggtitle('IM') +
  ggsci::scale_color_igv() +
  ggsci::scale_fill_igv() +
  ggridges::theme_ridges() +
  
  
  pData(cds) %>%
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
cds %>% 
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
dev.off()


pdf(
  file.path(monocle_out_dir, "barplot of annotation group by State.pdf"),
  height = 6,
  width  = 6
)
cds %>% 
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
dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory.pdf"), width = 8, height = 4)

plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="group", cell_size =0.1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="annotation_full", cell_size =0.1, show_backbone=TRUE) + 
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 4)
local({
  plot_cell_trajectory(cds, color_by = "annotation_full",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_1) +
    facet_wrap(~ group) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()



pdf(file.path(monocle_out_dir, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_full") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_full)
dev.off()


pdf(file.path(monocle_out_dir, "barplot tof cell props group by State - group.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, project) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = project))
dev.off()

pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_full ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_full ))
dev.off()


# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
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
  cds %>% 
  Biobase::pData() %>%
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




# # de analysis -----------
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性

Time_diff <- differentialGeneTest(
  cds[rownames(fData(cds))[fData(cds)$use_for_ordering], ],
  cores = 10,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)
saveRDS(Time_diff, file.path(monocle_out_dir, "Time_diff.rds"))

# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

# Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
Time_genes <- dplyr::filter(Time_diff, qval < 0.05) %>% pull(gene_short_name) %>% as.character()


p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=5, show_rownames=F, return_heatmap=T)
ggsave(file.path(monocle_out_dir, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)




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


