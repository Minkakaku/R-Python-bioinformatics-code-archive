library(monocle)
source("monocle/monocle fix bugs.R")


all_data_full2 <- readRDS('all_data_full2')

sc_monocle <- all_data_full2[, all_data_full2$celltype %in% c(
  'Chief_cell',
  'Endocrine',
  'Enterocytes',
  'Enteroendocrine',
  'GMC',
  'Goblet_cell',
  'PMC',
  'Tumor',
  'parietal_cells',
  'tuft_cells'
)]



# monocle out
monocle_out_dir <- "monocle"
dir.create(monocle_out_dir)

DimPlot(sc_monocle, reduction = 'harmonyumap') + scale_color_manual(values = col_plan_1)

VariableFeatures(sc_monocle)

expr_dat_monocle <- as(LayerData(JoinLayers(sc_monocle), "counts"), "sparseMatrix")
p_data_monocle <- sc_monocle[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sc_monocle),
  row.names = rownames(sc_monocle)
)
cds <- 
  newCellDataSet(
    expr_dat_monocle,
    phenoData   = new("AnnotatedDataFrame", p_data_monocle),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


# cds <- clusterCells(cds)


##使用seurat选择的高变基因⚠️
express_genes_monocle <- VariableFeatures(sc_monocle)
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

# ##使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(sc_monocle)
# express_genes_monocle <- subset(deg.cluster,p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)
# 
# ##使用monocle选择的高变基因⚠️
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# cds <- setOrderingFilter(cds, disp.genes)
# plot_ordering_genes(cds)

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

saveRDS(cds, file.path(monocle_out_dir, "cds.rds"))



# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds,color_by="annotation", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation")

# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

pdf(file.path(monocle_out_dir, "cell trajectory.pdf"), width = 8, height = 4)
  plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) + theme(legend.position = "right")
  plot_cell_trajectory(cds,color_by="annotation", size=1,show_backbone=TRUE) + 
    scale_color_manual(values = col_plan_1) +
    theme(legend.position = "right")
  plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory - annotation - split by group.pdf"), width = 21, height = 4)
  local({
    pData(cds)$group <- pData(cds)$group %>% factor(levels = c("Ctrl", "CD", "CD-F"))
    plot_cell_trajectory(cds, color_by = "annotation",size=1,show_backbone=TRUE) +
      scale_color_manual(values = col_plan_1) +
      facet_wrap(~ group) +
      theme(legend.position = "right")
  })
dev.off()



pdf(file.path(monocle_out_dir, "cell trajectory - tree.pdf"), width = 12, height = 10)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation)
dev.off()

# cell barplot 
p_col_cell_state <-
  cds %>% 
  pData() %>%
  dplyr::select(annotation, State) %>%
  group_by(annotation, State) %>%
  mutate(annotation = as.character(annotation)) %>%
  mutate(annotation = factor(annotation, levels = c(
    'S1 fibroblasts', 
    'Peri-cryptal S2 fibroblasts',
    'Immediate early-response (IER) fibroblasts',
    'S3 fibroblasts',
    'Myofibroblasts'
  ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq)) +
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

# de analysis
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 1, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(monocle_out_dir, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)  
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)


# BEAM analysis
BEAM_res <- BEAM(cds[express_genes_monocle,], branch_point = 1, cores = 2) 

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
write.csv(BEAM_res, "BEAM_res.csv", row.names = F)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T) # 有632个gene，太多了
#选前100个基因可视化
BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(
    cds[BEAM_genes,],  
    branch_labels = c("Myofibroblasts", "S3 fibroblasts"),
    branch_point = 1, 
    num_clusters = 3, 
    show_rownames = T, 
    return_heatmap = T
  )
ggsave(file.path(monocle_out_dir, "BEAM_heatmap.pdf"), p$ph_res, width = 6.5, height = 10)



plot_genes_branched_pseudotime(
  cds[c("CCDC80", "DCN", "ACTG2"),],
  color_by = "annotation"
)

plot_genes_jitter(
  cds[c("CCDC80", "DCN", "ACTG2"),],
  color_by = c("State")
)

# branch gene enrichment 
beam_genes_ls <- 
  p$ph_res$tree_row %>% cutree(k = 3) %>%
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

# enrichment plot
p_all_enrich_beam_ls <-
  all_enrich_beam_ls %>%
  imap(function(enrich, name){
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })



pdf(
  file.path(monocle_out_dir, "monocle BEAM gene enrichment analysis dotplot.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_beam_ls %>%  walk(print)
dev.off()



