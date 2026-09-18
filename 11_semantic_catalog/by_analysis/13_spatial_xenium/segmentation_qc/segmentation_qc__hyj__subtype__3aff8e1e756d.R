library(Seurat)
library(tidyverse)
library(patchwork)
library(sf)
library(ggplot2)
library(harmony)

all_data <- readRDS("data_for_analysis/all_data_projected.rds")

# PT subtype-----------------------
PT_out_dir <- "PT_sub"
dir.create(PT_out_dir)

# PT <- all_data %>% JoinLayers() %>% .[, .$first_type %in% "PT"]
PT <- all_data[, all_data$annotation_full %in% c('PT')]

pdf(file.path(PT_out_dir, 'nCount_Spatial_Polygons plot.pdf'), height = 6, width = 8)
VlnPlot(PT, features = "nCount_Spatial.Polygons", group.by = 'project', pt.size = 0)
dev.off()

PT_temp <- PT

# PT[['Spatial.Polygons']] <- split(PT[['Spatial.Polygons']], PT$sample, layers = "counts")

DefaultAssay(PT) <- "Spatial.Polygons"
Idents(PT) <- colnames(PT)

# cellcycle
PT <- 
  PT %>% 
  CellCycleScoring(
    s.features = cc.genes.updated.2019$s.genes %>% str_to_title,
    g2m.features = cc.genes.updated.2019$g2m.genes %>% str_to_title
  )


# preprocess
PT <- 
  PT %>% 
  NormalizeData(scale.factor = floor(mean(PT$nCount_Spatial.Polygons)))%>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData(
    vars.to.regress = c(
      'nCount_Spatial.Polygons'
    )
  )

saveRDS(PT, file.path(PT_out_dir, "PT_scale.rds"))
#PT <- readRDS(file.path(PT_out_dir, "PT_scale.rds"))

# PT <- PT %>% ScaleData(vars.to.regress = c('nCount_RNA', 'S.Score', 'G2M.Score'))
# VariableFeaturePlot(PT)

# VariableFeatures(PT) <- VariableFeatures(PT) %>% setdiff(c(ig_genes, tcr_genes, mt_genes, ribo_genes))



## pca -----
PT <- PT %>% RunPCA()

pdf(file.path(PT_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(PT, ndims = 30)
dev.off()

# set.seed(1234)
# suset_cells_t <- 
#   PT[[]] %>%
#   rownames_to_column("cell_id") %>%
#   group_by(project) %>%
#   slice_sample(prop = 50000/ncol(PT) ) %>%
#   dplyr::pull("cell_id")


pdf(file.path(PT_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(PT, dims = 1:5)
dev.off()
pdf(file.path(PT_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(PT, dims = 6:10)
dev.off()
pdf(file.path(PT_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(PT, dims = 11:15)
dev.off()
pdf(file.path(PT_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(PT, dims = 16:20)
dev.off()


## dimension reduction --------
PT_dims <- 15


library(future)
plan(multisession, workers = 24)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


PT <- RunUMAP(PT, dims = 1:PT_dims)


pdf(file.path(PT_out_dir, str_glue("Dimplot dim{PT_dims} project.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "project", label = TRUE) + NoLegend()
dev.off()


PT <- PT %>% 
  harmony::RunHarmony(
    dims.use = 1:PT_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

PT <- RunUMAP(
  PT, 
  dims = 1:PT_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(PT_out_dir, str_glue("Dimplot dim{PT_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(PT, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# plot QC figure
pdf(file.path(PT_out_dir, str_glue("qc - Dimplot dim{PT_dims}.pdf")), height = 6, width = 18)
FeaturePlot(PT, features = c('nCount_Spatial.Polygons', 'nFeature_Spatial.Polygons', 'log10GenesPerUMI'), ncol = 3) + NoLegend()
FeaturePlot(PT, features = c('nCount_Spatial.Polygons', 'nFeature_Spatial.Polygons', 'log10GenesPerUMI'), reduction = "harmonyumap", ncol = 3) + NoLegend()
dev.off()

pdf(file.path(PT_out_dir, str_glue("qc - Vlnplot dim{PT_dims}.pdf")), height = 6, width = 12)
VlnPlot(PT, features = c('nCount_Spatial.Polygons', 'nFeature_Spatial.Polygons', 'log10GenesPerUMI'), group.by = 'project', pt.size = 0) + NoLegend()
dev.off()

pdf(file.path(PT_out_dir, str_glue("FeaturePlot - dim{PT_dims} - CellCycle.pdf")), height = 6, width = 12)
FeaturePlot(PT, features = c('S.Score', 'G2M.Score'), reduction = "harmonyumap")
dev.off()

pdf(file.path(PT_out_dir, str_glue("DimPlot - dim{PT_dims} - Phase.pdf")), height = 5, width = 6)
DimPlot(PT, group.by = c('Phase'), reduction = "harmonyumap")
dev.off()

pdf(file.path(PT_out_dir, str_glue("DimPlot - dim{PT_dims} - Phase split.pdf")), height = 5, width = 15)
DimPlot(PT, group.by = c('Phase'), reduction = "harmonyumap", split.by = "Phase")
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
PT <-
  FindNeighbors(
    PT,
    reduction = reduction.use,
    dims = 1:PT_dims,
    k.param = 10
  )

PT <-
  FindClusters(
    PT,
    resolution = c(0.1, 0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{PT_dims} cluster - {reduction.use}.pdf")), height = 5, width = 6)
# DimPlot(PT, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(PT, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(PT, group.by = "Spatial.Polygons_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(PT, group.by = "Spatial.Polygons_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(PT, group.by = "Spatial.Polygons_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(PT, group.by = "Spatial.Polygons_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(PT, group.by = "Spatial.Polygons_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# Dotplot -------
# epi cluster
major_marker_ls <- c(
  'epithelials'  = 'Epcam',
  'epithelials'  = 'Krt5',
  'epithelials'  = 'Krt14',
  
  "Foveolar cells" = "Gkn1", #ACT
  "Foveolar cells" = "Muc5ac",
  "Foveolar cells" = "Tff1",
  
  'endothelials' = 'Pecam1', 
  'endothelials' = 'Vwf', 
  
  "Fibroblasts" = "Pdgfra", 
  "Fibroblasts" = "Ccdc80",
  "Fibroblasts" = "Sfrp4",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  "Fibroblasts" = "Dcn",
  "Fibroblasts" = "Acta2",
  
  "B cells" = "Cd19",
  "B cells" = "Cd79a",
  "B cells" = "Ms4a1",
  "Plasma" = "Tnfrsf17",
  "Plasma" = "Jchain",
  "Plasma" = "Mzb1",
  "Plasma" = "Xbp1",
  
  'T cells'           = 'Cd3d', 
  'T cells'           = 'Cd3e', 
  'T cells'           = 'Cd3g', 
  'T cells'           = 'Trac', 
  
 # "NK" = "Gnly",
  "NK" = "Nkg7",
  "NK" = "Ncam1",
  
  'Proliferating' = 'Mki67',
  'Proliferating' = 'Top2a',
  
  "Pericytes" = "Pdgfrb",
  "Pericytes" = "Rgs5",
  "Pericytes" = "Mcam",
  
  "Monocytes" = "Cd14", 
  #"Monocytes" = "Fcgr3a",
  
  "Macrophage" = "Cd68",
  "Macrophage" = "Cd163",
  "Macrophage" = "C1qa",
  "Macrophage" = "C1qc", 
  'Macrophage'  = 'Csf1r',
  
  "Dendritic cell" = "Cd83",
  #'Dendritic cell' = "Lilra4", 
  'pDCs' = "Irf4", 
  'pDCs' = "Irf7", 
  'pDCs' = "Irf8", 
  'Dendritic cell' = "Clec9a",
  'Dendritic cell' = "Ido1",
  #'Dendritic cell' = "Cd1c",
  'Dendritic cell' = "Fcer1a",
  'Lamp3_DCs' = "Lamp3",
  
  'netrophil' = "Csf3r",
  'netrophil' = "G0s2",
  'netrophil' = "Retnlg",
  
  "Myocytes" = "Tnnc1", 
  "Myocytes" = "Myh7", 
  
  "Mast cell" = "Cpa3",
  "Mast cell" = "Ctsg",
  "Mast cell" = "Tpsab1",
  "Mast cell" = "Kit"
)


PT_major_marker_ls1 <- c(
  'Healthy'='Slc5a12',
  'Healthy'='Slc2a8',
  'Healthy'='Cyp2e1',
  'Healthy'='Slc7a13',
  'Inj'='Havcr1',
  'Inj'='Krt20',
  'Failed_repair'='C3',
  #'Failed_repair'='Ccl74',
  'Failed_repair'='Vcam1'
)


PT_major_marker_ls2 <- c(
  'Healthy'='Slc5a2',
  'Healthy'='Slc34a1',
  'Healthy'='Slco4c1',
  #'Healthy'='Fgg',
  #'Healthy'='Ml2',
  'Stressed'='Fgg',
  #'Stressed'='Ml2',
  #'Injured'='Ml2',
  'Injured'='Krt20',
  #'Injured'='Gm20481',
  #'Severely injured PT'='Gm20481',
  'Failed repair PT'='Vcam1'
)


table(PT$Spatial.Polygons_snn_res.0.6)
# 指定cluster顺序
PT@meta.data$Spatial.Polygons_snn_res.0.6 <- factor(
  PT@meta.data$Spatial.Polygons_snn_res.0.6, levels = sort(as.numeric(levels(PT$Spatial.Polygons_snn_res.0.6)))
)

pdf(
  file.path(PT_out_dir, str_glue("Dotplot dim{PT_dims} cluster - {reduction.use} - major marker - Spatial.Polygons_snn_res.0.6.pdf")), 
  height = 5,
  width = 15
)
# major marker
DotPlot(
  PT,
  features = major_marker_ls,
  group.by = "Spatial.Polygons_snn_res.0.6"
) + 
  RotatedAxis()
  DotPlot(
  PT,
  features = PT_major_marker_ls1,
  group.by = "Spatial.Polygons_snn_res.0.6"
) + 
  RotatedAxis()
DotPlot(
  PT,
  features = PT_major_marker_ls2,
  group.by = "Spatial.Polygons_snn_res.0.6"
) + 
  RotatedAxis()
dev.off()


### FindAllMarkers --------
Idents(PT) <- "Spatial.Polygons_snn_res.0.6"
DefaultAssay(PT) <- "Spatial.Polygons"

enrich_outdir <- file.path(PT_out_dir, "enrichment analysis")
dir.create(enrich_outdir)


library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Mm.eg.db")
# library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 


diffmarkers <- FindAllMarkers(
  JoinLayers(PT), 
  test.use = "wilcox", 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  only.pos = TRUE
)

write.csv(diffmarkers, file.path(enrich_outdir, 'cluster diffmarker.csv'))

all_top5_gene <- 
  diffmarkers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene, 
    names = as.character(cluster)
  ))
all_top5_gene_unique <- all_top5_gene[!duplicated(all_top5_gene)]

p_dot_top5 <- 
  DotPlot(
    PT, 
    features = all_top5_gene_unique, 
    dot.scale = 12
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(enrich_outdir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 18
)

all_top10_gene <- 
  diffmarkers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 10)

write.csv(all_top10_gene,file = file.path(enrich_outdir, "top 10 gene.csv"))


# pdf(
#   file.path(PT_out_dir, str_glue("Dotplot dim{PT_dims} cluster - {reduction.use} - hep marker - RNA_snn_res.0.4.pdf")), 
#   height = 5,
#   width = 10
# )
# # liver hep marker
# DotPlot(
#   PT,
#   features = marker_Hep,
#   group.by = "Spatial.008um_snn_res.0.4"
# ) +
#   RotatedAxis()
# dev.off()


## annotation ----
Idents(PT) <- "Spatial.Polygons_snn_res.0.6"

cell_cluster <- c(
  '0'    = "C0", 
  '1'    = "C1", 
  '2'    = "C2", 
  '3'    = "C3", 
  '4'    = "C4", 
  '5'    = "C5"
)


PT <- RenameIdents(PT, cell_cluster)
PT$annotation_PT <- Idents(PT)


# 保存中间文件
saveRDS(PT, file.path(PT_out_dir, "PT.rds"))
PT <- readRDS(file.path(PT_out_dir, "PT.rds"))

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

# 注释后umap
pdf(file.path(PT_out_dir, str_glue("Dimplot dim{PT_dims} annotation_PT {reduction.use} .pdf")), height = 6, width = 8)
# DimPlot(sc_fb, group.by = "annotation_PT", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(PT, group.by = "annotation_PT", label = TRUE, reduction = 'harmonyumap', cols = plan) + NoLegend()
dev.off()

pdf(file.path(PT_out_dir, str_glue("dotplot dims{PT_dims} {reduction.use} annotation_PT major_marker.pdf")), width = 15, height = 6)
DotPlot(PT, group.by = "annotation_PT", features = major_marker_ls) + theme_90angle
dev.off()

# pdf(file.path(PT_out_dir, str_glue("dotplot dims{PT_dims} {reduction.use} annotation_PT fb_marker.pdf")), width = 15, height = 6)
# DotPlot(PT, group.by = "annotation_PT", features = marker_fb) + theme_90angle
# dev.off()



# annotation 空间图
dir.create("PT_subcluster/spatial_plot", recursive = T)

for (i in Images(PT)){
  p1 <- SpatialDimPlot(PT, 
  group.by = "annotation_PT", images = i, image.scale = "hires", plot_segmentations = TRUE, image.alpha = 0.2, stroke = 0, pt.size.factor = 0.5)
  p2 <- p1 + scale_fill_manual(values = plan, na.value = NA, labels = function(x) str_replace_na(x, "")) + 
  guides(fill = guide_legend(override.aes = list(size = 5, color = NA))) + 
  theme_bw() + xlab("")+ ylab("") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())
  
  ggsave(p2, filename = paste0("PT_subcluster/spatial_plot/", str_remove(i, "\\.polygons$"), "_annotation_PT_spatial.pdf"), width = 20, height = 20)
  ggsave(p2, filename = paste0("PT_subcluster/spatial_plot/", str_remove(i, "\\.polygons$"), "_annotation_PT_spatial.png"), width = 20, height = 20, dpi = 300)
}






### FindAllMarkers --------
Idents(PT) <- "annotation_PT"
DefaultAssay(PT) <- "Spatial.Polygons"

enrich_outdir <- file.path(PT_out_dir, "enrichment analysis")
dir.create(enrich_outdir)


library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library('clusterProfiler')
library("org.Mm.eg.db")
# library("org.Hs.eg.db")
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 



diffmarkers <- FindAllMarkers(
  JoinLayers(PT), 
  test.use = "wilcox", 
  logfc.threshold = 0.25, 
  min.pct = 0.1, 
  only.pos = TRUE
)

write.csv(diffmarkers, file.path(enrich_outdir, 'celltype diffmarker.csv'))



all_top5_gene <- 
  diffmarkers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 > 0.25) %>% 
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice_head(n = 5) %>% 
  with(structure(
    gene, 
    names = as.character(cluster)
  ))
all_top5_gene_unique <- all_top5_gene[!duplicated(all_top5_gene)]

p_dot_top5 <- 
  DotPlot(
    PT, 
    features = all_top5_gene_unique, 
    dot.scale = 12
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(enrich_outdir, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 18
)



PT <- 
  PT %>% 
  AddMetaData(roi_df)

PT$ROI[is.na(PT$ROI)] <- "Other"
table(PT$ROI)




# cell props barplot  --------

PT_meta <- 
  PT[[]] %>% 
  rownames_to_column("Barcode")

# PT_meta$project <- factor(PT_meta$project, levels = c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  PT_meta %>%
  group_by(ROI, annotation_PT) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = ROI, y = freq, fill = annotation_PT)) + 
  geom_col() +
  # ggsci::scale_fill_igv() + 
  scale_fill_manual(values = plan) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  PT_meta %>%
  group_by(ROI, project, annotation_PT) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_PT, y = freq, fill = ROI)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  # ggsci::scale_fill_igv() +
  scale_fill_manual(values = plan) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(PT_out_dir, "cell_freq_project.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(PT_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)


saveRDS(PT, file.path(PT_out_dir, "PT.rds"))





# rostral vs caudal diff analysis -------
epi_enrich_outdir <- file.path(PT_out_dir, "enrichment analysis")

PT <- readRDS(file.path(PT_out_dir, "PT.rds"))
Idents(PT) <- "ROI"
DefaultAssay(PT) <- "Spatial.Polygons"

PT_Control <- PT[, PT$project %in% c("Control")]
PT_P2 <- PT[, PT$project %in% c("P2")]

## PT_Control ------------
diffmarkers <- 
    FindMarkers(
        PT_Control,
        ident.1 = 'Control_rostral',
        ident.2 = 'Control_caudal',
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        group.by = 'ROI'
        )

write.csv(diffmarkers, file.path(epi_enrich_outdir, 'Control diffmarker.csv'))

enrichGO_ls <- enrichGO(
  rownames(diffmarkers), 
  OrgDb = org.Mm.eg.db, 
  ont = "ALL", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  keyType = "SYMBOL"
)
write.csv(enrichGO_ls@result, file = file.path(epi_enrich_outdir, "Control GO enrichment table.csv"), row.names = FALSE)

enrichGO_dot <- enrichplot::dotplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_dot, filename = file.path(epi_enrich_outdir, "Control GO dotplot.pdf"), height = 6, width = 8)

enrichGO_bar <- barplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_bar, filename = file.path(epi_enrich_outdir, "Control GO barplot.pdf"), height = 6, width = 8)


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(epi_enrich_outdir, 'Control KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(epi_enrich_outdir, 'Control KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(epi_enrich_outdir, 'Control KEGG barplot.pdf'),
  height = 6,
  width = 8)




## PT_P2 ------------
diffmarkers <- 
    FindMarkers(
        PT_P2,
        ident.1 = 'P2_rostral',
        ident.2 = 'P2_caudal',
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        group.by = 'ROI'
        )

write.csv(diffmarkers, file.path(epi_enrich_outdir, 'P2 diffmarker.csv'))

enrichGO_ls <- enrichGO(
  rownames(diffmarkers), 
  OrgDb = org.Mm.eg.db, 
  ont = "ALL", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  keyType = "SYMBOL"
)
write.csv(enrichGO_ls@result, file = file.path(epi_enrich_outdir, "P2 GO enrichment table.csv"), row.names = FALSE)

enrichGO_dot <- enrichplot::dotplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_dot, filename = file.path(epi_enrich_outdir, "P2 GO dotplot.pdf"), height = 6, width = 8)

enrichGO_bar <- barplot(
  enrichGO_ls, 
  showCategory = 10, 
  label_format = 50
)
ggsave(enrichGO_bar, filename = file.path(epi_enrich_outdir, "P2 GO barplot.pdf"), height = 6, width = 8)


de_genes_id_for_enrichmenet <- 
  clusterProfiler::bitr(
    rownames(diffmarkers),
    OrgDb = "org.Mm.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID

enrichKEGG_ls <- enrichKEGG(
  de_genes_id_for_enrichmenet,
  organism = "mmu", 
  keyType = "ncbi-geneid",
  # keyType = "kegg", 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
  # use_internal_data = TRUE
)
enrichKEGG_ls <- setReadable(enrichKEGG_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(enrichKEGG_ls@result,file = file.path(epi_enrich_outdir, 'P2 KEGG enriment table.csv'),row.names = FALSE)

enrichKEGG_dot <- enrichplot::dotplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
enrichKEGG_bar <- barplot(
  enrichKEGG_ls,
  showCategory = 10,
  label_format = 50)
ggsave(
  enrichKEGG_dot,
  filename = file.path(epi_enrich_outdir, 'P2 KEGG dotplot.pdf'),
  height = 6,
  width = 8)
ggsave(
  enrichKEGG_bar,
  filename = file.path(epi_enrich_outdir, 'P2 KEGG barplot.pdf'),
  height = 6,
  width = 8)











