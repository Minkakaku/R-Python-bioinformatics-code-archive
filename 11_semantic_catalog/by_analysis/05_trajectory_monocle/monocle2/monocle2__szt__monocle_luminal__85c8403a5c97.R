library(Seurat)
library(tidyverse)
library(monocle)
source("monocle fix bugs.R")
source('color_pallete.R')

# set default ------- 
monocle_out_dir <- "monocle_result" 
# dir.create(monocle_out_dir) 


# ADPC luminal ------------------------------------------------------------------
monocle_out_dir <- "monocle_result/ADPC" 
# dir.create(monocle_out_dir)
Luminal4 <- readRDS('result_8um/Luminal_reannotation3/Luminal4_final.rds')
Luminal4$clusters %>% table()

Luminal4$annotation_monocle <- 
  Luminal4$clusters %>% 
  str_replace_all(c('0' = 'LE_0', 
                    '1' = 'LE_1', 
                    '2' = 'LE_2', 
                    '3' = 'LE_3', 
                    '4' = 'LE_4', 
                    '5' = 'LE_5'))

Luminal4$annotation_monocle %>% table()


# subset ---------------------
set.seed(1234)
meta_Luminal4_subset <-
  Luminal4[['annotation_monocle']] %>%
  rownames_to_column("cell_id") %>%
  group_by(annotation_monocle) %>%
  slice_sample(n = 3000)


## cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(Luminal4), "counts"), "sparseMatrix")
p_data_monocle <- Luminal4[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(Luminal4),
  row.names       = rownames(Luminal4)
)


cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_Luminal4_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_Luminal4_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)




##使用seurat选择的高变基因⚠️ ------
express_genes_monocle <- VariableFeatures(Luminal4)
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'seuratHVGs')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_seuratHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_seuratHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_seuratHVGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Luminal4[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()
saveRDS(cds, file = file.path(monocle_out_dir, "cds_seuratHVGs.rds"))




##使用clusters差异表达基因 --------------
Idents(Luminal4) <- "annotation_monocle"
deg.cluster <- FindAllMarkers(JoinLayers(Luminal4))
express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'clusterDEGs')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_clusterDEGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_clusterDEGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_clusterDEGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_clusterDEGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_clusterDEGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Luminal4[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()
saveRDS(cds, file = file.path(monocle_out_dir, "cds_clusterDEGs.rds"))




##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.001 & dispersion_empirical >= 0.6 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'monocleHVGs')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_monocleHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_monocleHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_monocleHVGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Luminal4[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()
saveRDS(cds, file = file.path(monocle_out_dir, "cds_monocleHVGs.rds"))




## dpFeature ----------------
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.1)) ## 不要也写0.1 ，而是要写0.01。
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'dpFeature')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_dpFeature.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_dpFeature.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_dpFeature.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_dpFeature.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_dpFeature.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Luminal4[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()
saveRDS(cds, file = file.path(monocle_out_dir, "cds_dpFeature.rds"))





















#cds <- orderCells(cds, root_state = 4)









