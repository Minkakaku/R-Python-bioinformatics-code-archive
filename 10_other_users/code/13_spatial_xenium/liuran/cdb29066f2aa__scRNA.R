setwd("scRNA")

library(Seurat)
library(tidyverse)
library(dplyr)

sc_outdir <- file.path("scRNA/")

all_data_sc.data <- Read10X(data.dir = "scRNA/Matrix/")
all_data_sc <- CreateSeuratObject(counts = all_data_sc.data)

all_data_sc <- NormalizeData(all_data_sc, normalization.method = "LogNormalize", scale.factor = 10000)
all_data_sc <- FindVariableFeatures(all_data_sc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(all_data_sc)
all_data_sc <- ScaleData(all_data_sc, features = all.genes, vars.to.regress = NULL)

all_data_sc <- RunPCA(object = all_data_sc)
pdf(file.path(sc_outdir, "pca plot.pdf"), height = 6, width = 8)
ElbowPlot(all_data_sc, ndims = 30)
dev.off()

all_data_sc <- FindNeighbors(all_data_sc, reduction = "pca", dims = 1:20)
all_data_sc <- FindClusters(all_data_sc, resolution = 0.8, algorithm = 1)
all_data_sc <- RunUMAP(all_data_sc, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
all_data_sc <- RunTSNE(all_data_sc, dims = 1:20, perplexity = 30)


pdf(file.path(sc_outdir, str_glue("UMAP Dimplot.pdf")), height = 5, width = 6)
DimPlot(all_data_sc, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


# add annotation -------
annotation <- read.csv("annotations.csv", header = T)
annotation <- annotation %>% column_to_rownames("barcode")

all_data_sc <- AddMetaData(all_data_sc, metadata = annotation)


pdf(file.path(sc_outdir, str_glue("UMAP Dimplot annotation.pdf")), height = 5, width = 6)
DimPlot(all_data_sc, group.by = "annotation", label = TRUE, reduction = "umap", repel = T) + NoLegend()
dev.off()

saveRDS(all_data_sc, "all_data_sc.rds")




Idents(all_data_sc) <- "annotation"

all_markers <- FindAllMarkers(
  all_data_sc,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = TRUE
)

all_markers %>% 
  dplyr::select(gene, everything()) %>% 
  write.csv(
    file = file.path(sc_outdir, "marker.csv"),
    row.names = FALSE
  )

all_top5_gene <- 
  all_markers %>% 
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
    all_data_sc,
    features = all_top5_gene_unique,
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  p_dot_top5, 
  filename = file.path(sc_outdir, "Top 5 gene dotplot.pdf"),
  height = 6, 
  width = 20
  )




old.size <- object.size(all_data_sc)
format(old.size, units = "Gb")

all_data_sc@assays$RNA@layers$scale.data <- matrix()
new_size <- object.size(all_data_sc)
format(new_size, units = "Gb")

saveRDS(all_data_sc, "all_data_sc_plot.rds")








# cell marker ----
all_data_sc$annotation <- factor(all_data_sc$annotation, levels = c(
  "Basal keratinocytes",
  "Epithelial Stem Cell",
  "Sweat Gland Cell",
  "Melanocyte",
  "Endothelial Cell",
  "Lymphatic Endothelial Cell",
  "Vascular Endothelial Cell",
  "Fibroblast",
  "Pericyte",
  "Pericyte_or_Smooth Muscle Cell",
  "T Cell",
  "Macrophage_or_Dendritic Cell",
  "Mast Cell",
  "Neuronal cell"
))
major_marker_ls <- c(
  # "Basal_cell" = "KRT5", 
  "Basal_cell" = "KRT14",
  "Basal_cell" = "COL17A1", 
  "KC" = "KRT5", 
  "KC" = "CDH1", 
  # "Keratinocyte_progenitor_cell" = "ABCG2", 
  "Keratinocyte_progenitor_cell" = "CD34", 
  "Keratinocyte_progenitor_cell" = "TP63", 
  "Spinous_cell" = "KRT1", 
  "Spinous_cell" = "KRT10", 
  # "Spinous_cell" = "KRT5", 
  # "Spinous_cell" = "KRT14", 
  "Sebaceous_gland" = "MGST1", 
  "Sebaceous_gland" = "FASN", 
  "Granular_cell" = "SPINK5", 
  "Granular_cell" = "FLG", 
  # "Granular_cell" = "KRT1", 
  # "Granular_cell" = "KRT10", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "UBE2C", 
  "Proliferating_keratinocyte_Mitotic_cell" = "TOP2A", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "PCNA", 
  "Proliferating_keratinocyte_Mitotic_cell" = "MKI67", 
  "Hair_follicle" = "KRT6B", 
  "Hair_follicle" = "KRT17", 
  # "Hair_follicle" = "CD34", 
  # "Hair_follicle" = "RUNX1", 
  # "Sweat_Gland_Cel" = "SCGB1D2", 
  "Sweat_Gland_Cel" = "SCGB2A2", 
  # "Sweat_Gland_Cel" = "DCD", 
  "Sweat_Gland_Cel" = "AQP5",
  "Endothelial_cell" = "PECAM1", 
  "Endothelial_cell" = "SELE", 
  "Endothelial_cell" = "CLDN5", 
  "Endothelial_cell" = "CDH5", 
  # "Endothelial_cell" = "VWF", 
  "VEC" = "VWF", 
  "VEC" = "CD36", 
  "VEC" = "ACKR1", 
  # "VEC" = "FABP4", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "COL1A2", 
  # "Fibroblast" = "SFRP2", 
  # "Fibroblast" = "COL6A5", 
  # "Fibroblast" = "COL18A1", 
  "LEC" = "CCL21", 
  "LEC" = "LYVE1", 
  "Langerhans_cells" = "CD1A", 
  "Langerhans_cells" = "CD207", 
  "Dendritic_cell" = "CD1C", 
  "Dendritic_cell" = "CD80", 
  "Dendritic_cell" = "CD83", 
  "Dendritic_cell" = "CD86", 
  # "Dendritic_cell" = "HLA-A", 
  # "Dendritic_cell" = "HLA-B", 
  # "Dendritic_cell" = "HLA-C", 
  # "Dendritic_cell" = "HLA-DRA", 
  # "Dendritic_cell" = "HLA-DPA1", 
  # "Dendritic_cell" = "HLA-DQA1", 
  "Dendritic_cell" = "LAMP3", 
  "T_cell" = "CD3E", 
  "T_cell" = "CD4", 
  "T_cell" = "CD8A", 
  "NK_cell" = "FCGR3A", 
  "NK_cell" = "CD2", 
  # "NK_cell" = "TNFRSF8", 
  # "NK_cell" = "CD44", 
  "Mast_cell" = "TPSAB1", 
  "Mast_cell" = "TPSB2", 
  "Mast_cell" = "CPA3", 
  "Melanocyte_cell" = "MLANA", 
  "Melanocyte_cell" = "TYRP1", 
  "Melanocyte_cell" = "PMEL", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "FABP4", 
  "vSMC" = "MYL9", 
  "vSMC" = "TPM2", 
  "vSMC" = "RERGL", 
  # "Stem_cell" = "LGR5", 
  "Stem_cell" = "SOX9", 
  # "Stem_cell" = "CD34", 
  "Stem_cell" = "KRT15", 
  "Stem_cell" = "TCF4", 
  "Stem_cell" = "CTNNB1", 
  "Stem_cell" = "LRIG1", 
  "Stem_cell" = "CD200", 
  "Stem_cell" = "KRT19", 
  # "Osteogenic_precursor_like_cell" = "BMPR1B", 
  "Dermal_papilla_cell" = "PROM1", 
  # "Embryonic_stem_cell" = "NANOG", 
  # "Embryonic_stem_cell" = "POU5F1", 
  "Merkel_cell" = "KRT18", 
  # "Merkel_cell" = "KRT20", 
  "Merkel_cell" = "NCAM1", 
  "Pluripotent_stem_cell" = "ANPEP", 
  "Progenitor_cell" = "EPCAM", 
  "Progenitor_cell" = "NCAM2",
  # "Neuron" = "CAMK2A",
  # "Neuron" = "CBLN2",
  "Neuron" = "LDB2",
  # "Neuron" = "GAD1",
  "Neuron" = "LHFPL3",
  "Neuron" = "L1CAM",
  "Neuron" = "TUBB3"#,
  # "Neuron" = "PCDH15"
  # "Progenitor_cell" = "PROM1", 
  # "Treg" = "LAG3", 
  # "Treg" = "ITGA2", 
  # "Treg" = "FOXP3", 
  # "Treg" = "IKZF2", 
  # "Treg" = "NRP1", 
  # "Treg" = "NT5E", 
  # "Treg" = "IL10", 
  # "Treg" = "TGFB1", 
)


pdf(file.path("scRNA", str_glue("dotplot major_marker.pdf")), width = 30, height = 8)
DotPlot(all_data_sc, group.by = "seurat_clusters", features = major_marker_ls) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(angle = 45, hjust = 0, vjust = 0))
DotPlot(all_data_sc, group.by = "annotation", features = major_marker_ls) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(angle = 45, hjust = 0, vjust = 0))
dev.off()








