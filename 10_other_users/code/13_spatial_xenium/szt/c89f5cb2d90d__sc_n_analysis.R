library(tidyverse)
library(Seurat)

visim_ls_pass <- readRDS("visim_ls_pass.rds")


# N 样本 -------------------------------

sc_n <- visim_ls_pass$N
sc_n <- sc_n %>% 
  ScaleData(vars.to.regress = c(
    # str_glue('nCount_{assay}')
    "nCount_Spatial.008um"
  ))

sc_n <-  sc_n %>% FindVariableFeatures() %>% RunPCA()


sc_n %>% ElbowPlot()

DimHeatmap(sc_n)


## dimension reduction --------
dims <- 15

sc_n <- RunUMAP(sc_n, dims = 1:dims, return.model = TRUE, assay  =  "Spatial.008um") 


# pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
# DimPlot(sc_n, group.by = "ROI", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
# dev.off()


## cluster ------
sc_n <-
  FindNeighbors(
    sc_n,
    reduction = "pca",
    dims = 1:dims,
    k.param = 20
  )

sc_n <-
  FindClusters(
    sc_n,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("sc_n Dimplot dim{dims} cluster - pca.pdf")), height = 6, width = 8)

DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

dev.off()


# spatial-dimplot
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} cluster - pca.pdf")), height = 12, width = 12)

str_c('Spatial.008um_snn_res.', seq(0.2, 0.8, by = 0.2) ) %>%
  map(function(res){
    SpatialDimPlot(
      sc_n, 
      group.by = res,
      label = T, 
      repel = T,
      alpha = 0.5,
      shape = 22,
      label.size = 4,
      pt.size.factor = 6,
      image.scale = "hires"
      # images = name
    ) +
      ggsci::scale_fill_igv()
  })

dev.off()



## dotplot -------
major_marker_ls <- c(
  'Myoid cell' = 'ACTA2',
  'Myoid cell' = 'MYH11', 
  'Myoid cell' = 'PTCH1',
  'Myoid cell' = 'PTCH2',
  
  'Sertoli cell' = 'WT1', 
  'Sertoli cell' = 'SOX9', 
  'Sertoli cell' = 'BEX2',
  'Sertoli cell' = 'GATA4',

  'Macrophage cell' = 'CD68',
  'Macrophage cell' = 'MS4A6A', 
  'Macrophage cell' = 'LYZ', 
  'Macrophage cell' = 'IL3RA',
  
  'Leydig cell' = 'CYP11A1',
  'Leydig cell' = 'STAR', 
  'Leydig cell' = 'IGF1', 
  'Leydig cell' = 'IGFBP5',
  
  
  'Endothelial cell' = 'VWF',
  'Endothelial cell' = 'PECAM1', 
  'Endothelial cell' = 'CD34', 
  'Endothelial cell' = 'TGFBR2',
  
  
  'T cell' = 'CD3D',
  'T cell' = 'CD3E', 
  'T cell' = 'CD3G', 
  'T cell' = 'CD2',
  
  Spermatogonia = 'FGFR3', 
  Spermatogonia = 'STRA8', 
  Spermatogonia = 'DMRT1', 
  Spermatogonia = 'UCHL1', 
  Spermatogonia = 'TCF3',
  
  Spermatocytes = 'SYCP1',
  Spermatocytes = 'SMC1B', 
  Spermatocytes = 'SYCP3', 
  Spermatocytes = 'OVOL2',
  
  Spermatids = 'TNP1',
  Spermatids = 'PRM1', 
  Spermatids = 'TXNDC2', 
  Spermatids = 'PRM3'
  
)


theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("sc_n dotplot dims{dims} pca main_marker.pdf")), width = 25, height = 6)

DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle

dev.off()


## pathway score -------

major_marker_ls2 <-
  major_marker_ls %>%
  enframe("pathway", "gene") %>%
  group_by(pathway) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE) 


pathway_score_ls <-
  major_marker_ls2 %>%
  imap(function(x, name){# browser()
    sc_n %>%
      AddModuleScore(features = list(x) ) %>%
      FeaturePlot(
        features = "Cluster1"
      ) +
      ggtitle(name)
  })
  
pdf(
  file.path(dr_outdir, str_glue("cell annotation pathway scores dims{dims}.pdf")),
  height = 18,
  width  = 24
)
pathway_score_ls %>% patchwork::wrap_plots()
dev.off()

## find markers -------
Idents(sc_n) <- "Spatial.008um_snn_res.0.6"
DefaultAssay(sc_n) <- "Spatial.008um"
all_markers <-
  FindAllMarkers(
    JoinLayers(sc_n),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

# all_de_genes <-
#   all_markers %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   dplyr::filter(pct.1 > 0.25) %>%
#   dplyr::filter(abs(avg_log2FC) > 1) %>% ##
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 10) %>%
#   pull("gene")

## annotation ------
Idents(sc_n) <- "Spatial.008um_snn_res.0.6"
cell_cluster <- c(
  '0'   = "Spermatids",
  '1'   = "Myoid cell", ##
  '2'   = "Spermatocytes",
  '3'   = "Sertoli cell", # 
  '4'   = "Leydig cell", # 
  '5'   = "Spermatids",
  '6'   = "Myoid cell",
  '7'   = "Leydig cell",
  '8'   = "Spermatids",
  '9'   = "Leydig cell", ##
  '10'  = "Myoid cell",
  '11'  = "Myoid cell", ##
  '12'  = "Myoid cell",
  '13'  = "Myoid cell",
  '14'  = "Myoid cell"
)


sc_n <- RenameIdents(sc_n, cell_cluster)
sc_n$annotation <- Idents(sc_n)


# spatial-dimplot
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} annotation - pca.pdf")), height = 12, width = 12)
SpatialDimPlot(
  sc_n, 
  group.by = "annotation",
  label = T, 
  repel = T,
  alpha = 0.5,
  shape = 22,
  label.size = 4,
  pt.size.factor = 6,
  image.scale = "hires"
  # images = name
) +
  ggsci::scale_fill_igv()
dev.off()


# BANKSY ----
# BANKSY
sc_n_ls <- list(N = visim_ls$N) %>%
  map(function(visium){

    DefaultAssay(visium) <- "Spatial.008um"
    all_data2 <- visium %>% FindVariableFeatures()
    
    DefaultAssay(all_data2) <- "Spatial.008um"
    all_data2 <- RunBanksy2(all_data2 %>% JoinLayers(),
                            lambda = 0.6, 
                            verbose = TRUE,
                            assay = "Spatial.008um", 
                            # slot = "data",
                            layer = "data",
                            features = "variable",
                            k_geom = 10
    )
    
    
    
    DefaultAssay(all_data2) <- "BANKSY"
    all_data2 <- RunPCA(all_data2, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(all_data2), npcs = 30)
    ElbowPlot(all_data2, reduction = "pca.banksy")
    all_data2 <- FindNeighbors(all_data2, reduction = "pca.banksy", dims = 1:10)
    
    resolution <- 0.5
    all_data2  <- FindClusters(all_data2, cluster.name = "banksy_cluster", resolution = resolution)
    
    while (n_distinct(all_data2$banksy_cluster) < 10 & resolution < 1.3) {
      resolution <- resolution + 0.1
      all_data2  <- FindClusters(all_data2, cluster.name = "banksy_cluster", resolution = resolution)
    }
    
    all_data2
  })


SpatialDimPlot(
  sc_n_ls$N, 
  group.by = 'banksy_cluster',
  label = T, 
  repel = T,
  alpha = 0.5,
  shape = 22,
  label.size = 4,
  pt.size.factor = 6,
  image.scale = "hires"
  # images = name
) +
  ggsci::scale_fill_igv()



# CARD free -------
# devtools::install_github('YingMa0107/CARD')

library(CARD)

spatial_count <-
  sc_n %>%
  LayerData(layer = "counts") 
  
spatial_location <- 
  sc_n %>%
  GetTissueCoordinates() %>%
  dplyr::select(x, y)

CARDfree_obj <-
  createCARDfreeObject(
    markerList = major_marker_ls2,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    minCountGene = 10,
    minCountSpot = 2
  ) 


CARDfree_obj <- CARD_refFree(CARDfree_obj)

## STOP! The average number of unique marker genes for each cell type is less than 20 ...

# colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
#            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77",
#            "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
# 
# CARDfree_obj@Proportion_CARD = CARDfree_obj@Proportion_CARD[,c(8,10,14,2,1,6,12,18,7,13,20,19,16,17,11,15,4,9,3,5)]
# colnames(CARDfree_obj@Proportion_CARD) = paste0("CT",1:20)
# p1 <- CARD.visualize.pie(CARDfree_obj@Proportion_CARD,CARDfree_obj@spatial_location,colors = colors)
# p1

# STdeconvolve ----
library(STdeconvolve)

data(mOB)
mOB$annot %>% head

## remove pixels with too few genes
counts <- cleanCounts(counts = spatial_count,
                      min.lib.size = 100,
                      min.reads = 1,
                      min.detected = 1)

## feature select for genes
corpus <- restrictCorpus(counts,
                         removeAbove=1,
                         removeBelow = 0.05,
                         alpha = 0.1,
                         plot = TRUE,
                         verbose = TRUE)

all_zero_idx <- colSums(corpus) == 0

## Note: the input corpus needs to be an integer count matrix of pixels x genes
ldas <- fitLDA(t(as.matrix(corpus[, !all_zero_idx])), Ks = seq(2, 9, by = 1), plot=TRUE, verbose=TRUE)


## select model with minimum perplexity
optLDA <- optimalModel(models = ldas, opt = "min")
## extract pixel cell-type proportions (theta) and cell-type gene expression profiles (beta) for the given dataset
results <- getBetaTheta(optLDA, corpus = t(as.matrix(corpus[, !all_zero_idx])))
deconProp <- results$theta
deconGexp <- results$beta*1000
vizAllTopics(deconProp, 
             spatial_location[!all_zero_idx, ], 
             # groups = annot, 
             # group_cols = rainbow(length(levels(annot))),
             r=0.4)


deconProp_for_plot <-
  deconProp %>% 
  as.data.frame() %>%
  merge(spatial_location[!all_zero_idx, ], by = "row.names") 

p_deconProp_for_plot <-
  deconProp_for_plot %>% 
  # dplyr::filter(x < 15000, y < 15000) %>%
  ggplot(aes(x = x, y = y)) +
  # coord_equal(
  #   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff), 
  #   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
  # ) +
  coord_fixed(ratio = 1/1) +
  scatterpie::geom_scatterpie(
    aes(x = x, y = y, group = Row.names, r = 15,
        # color = Pixel.Groups
        ), 
    pie_scale = 1,
        # lwd = lwd, 
    data = deconProp_for_plot, # %>% dplyr::filter(x < 15000, y < 15000), 
    cols = as.character(1:4), # topicColumns, 
    legend_name = "CellTypes"
  ) +
  theme_void()


ggsave(
  p_deconProp_for_plot,
  file = file.path(dr_outdir, "STdeconvolve spatial plot.pdf"),
  height = 32,
  width  = 32
)




# p <- ggplot2::ggplot(mapping = ggplot2::aes(x = 0:dim(overlay)[2], 
#       y = 0:dim(overlay)[1])) + ggplot2::coord_equal(xlim = c(0, 
#       dim(overlay)[2]), ylim = c(0, dim(overlay)[1]), 
#       expand = FALSE) + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
#       axis.line = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), 
#       axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
#       axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), 
#       panel.background = ggplot2::element_blank(), plot.background = ggplot2::element_blank(), 
#       legend.text = ggplot2::element_text(size = 12, colour = "black"), 
#       legend.title = ggplot2::element_text(size = 12, 
#         colour = "black")) + ggplot2::annotation_raster(overlay, 
#       xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
#       scatterpie::geom_scatterpie(ggplot2::aes(x = x, 
#         y = y, group = Row.names, r = r, color = Pixel.Groups), 
#         lwd = lwd, data = theta_ordered_pos, cols = topicColumns, 
#         legend_name = "CellTypes") + ggplot2::scale_fill_manual(values = as.vector(topicCols)) + 
#       ggplot2::scale_color_manual(values = group_cols)



