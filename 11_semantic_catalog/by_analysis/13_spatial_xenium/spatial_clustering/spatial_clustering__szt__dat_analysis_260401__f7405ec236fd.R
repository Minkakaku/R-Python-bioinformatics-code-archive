# subtype code: dat_analysis_cycle3/subtype/subtype_code.R

library(tidyverse)
library(Seurat)

CD8T <- readRDS("dat_analysis_cycle3/subtype/CD8T/CD8_subtype.rds")
Idents(CD8T) <- 'subtype_full'
# # marker
# CD8_Naive/Tcm：TCF7, IL7R, CCR7, SELL, CD27, BACH2 高；PDCD1/TOX 低 
# CD8_Tpex：PDCD1 阳性，同时保留 TCF7/IL7R/CXCR5/XCL1/CCR7/SELL 程序；terminal exhaustion markers 不高 
# CD8_Tex：PDCD1, TOX, HAVCR2, LAG3, TIGIT, ENTPD1, CXCL13 高；TCF7/IL7R 低 
# CD8_Tem / cytotoxic-like：CCL3/4/5, GZMB, PRF1, IFNG, NKG7 高；耗竭程序不主导 
# CD8_Trm：ITGAE, CXCR6, CD69 高 
# Cycling CD8 T：MKI67, TOP2A, CENPF, PCLAF 高


cd8_marker <- c(
  'CD8_Naive/Tcm' = 'TCF7',
  'CD8_Naive/Tcm' = 'IL7R',
  'CD8_Naive/Tcm' = 'CCR7',
  'CD8_Naive/Tcm' = 'SELL',
  'CD8_Naive/Tcm' = 'CD27',
  'CD8_Naive/Tcm' = 'BACH2',
  
  'CD8_Tpex' = 'PDCD1', # TCF7/IL7R/CXCR5/XCL1/CCR7/SELL
  
  'CD8_Tex' = 'TOX', 
  'CD8_Tex' = 'HAVCR2', 
  'CD8_Tex' = 'LAG3', 
  'CD8_Tex' = 'TIGIT', 
  'CD8_Tex' = 'ENTPD1', 
  'CD8_Tex' = 'CXCL13',
  
  'CD8_Tem' = 'CCL3',
  'CD8_Tem' = 'CCL4',
  'CD8_Tem' = 'CCL5',
  'CD8_Tem' = 'GZMB',
  'CD8_Tem' = 'PRF1',
  'CD8_Tem' = 'IFNG',
  'CD8_Tem' = 'NKG7',
  
  'CD8_Trm' = 'ITGAE',
  'CD8_Trm' = 'CXCR6',
  'CD8_Trm' = 'CD69',
  
  'Cycling CD8 T' = 'MKI67',
  'Cycling CD8 T' = 'TOP2A', 
  'Cycling CD8 T' = 'CENPF', 
  'Cycling CD8 T' = 'PCLAF'
  
)



CD8T %>% dim


DimPlot(
  CD8T,
  group.by = "subtype_full",
  reduction = "full.harmonyumap.1"
)

DefaultAssay(CD8T) <- "sketch"

p_dot_marker <-
  DotPlot(
    CD8T,
    features = cd8_marker,
    group.by = "sub.cluster"
  ) +
  scale_color_gradientn(colours = c("white", "white", "blue")) +
  RotatedAxis()

CD8T$subtype_full %>% is.na %>% table


# ll <- CellsByIdentities(object = CD8T, cells = colnames(CD8T[['RNA']]))

# export 
p_dot_marker$theme$strip.clip <- "off"
ggsave(
  p_dot_marker,
  filename = file.path("dat_analysis_cycle3/", "subtype", "dotplot of CD8T using markers provided by customer.pdf"),
  height = 6,
  width = 12
)

