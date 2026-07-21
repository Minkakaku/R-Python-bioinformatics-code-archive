library(spacexr)

rctd_outdir <- "RCTD"
dir.create(rctd_outdir)

# FORMAT CONVERSION ---------
# config python
library(reticulate)
use_condaenv(condaenv = "/home/liuran/software/miniconda3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")

##1.单细胞数据-----
sc <- import("scanpy")
adata_sc <- sc$read_h5ad("/mnt/NexData_file/home/hyj/Rstudio/FWSP20241023_shangzhiqun/RCTD/ad_TME_Celltype_final.h5ad")

###counts-----
#counts_x <- adata_sc$layers$count %>% t %>% as.matrix()
counts_x <- t( as.matrix( adata_sc$layers[["counts"]] ) )
colnames(counts_x) <- rownames(adata_sc$obs)
rownames(counts_x) <- rownames(adata_sc$var)

counts_x[1:10,1:10]

adata_sc_meta %>% head()

###metadata-----
adata_sc_meta <- adata_sc$obs
table(adata_sc_meta$Celltype_final,adata_sc_meta$sample)
adata_sc_meta$annotation <- adata_sc_meta$Celltype_final

anno_counts <- table(adata_sc_meta$annotation)
print(anno_counts)

scRNA2 <- CreateSeuratObject(
  counts_x,
  meta.data = adata_sc_meta
)

scRNA2@meta.data <- adata_sc_meta
saveRDS(scRNA2, "RCTD/scRNA2.rds")

anno_counts <- table(scRNA2@meta.data$annotation)
print(anno_counts)

# 设为默认 Assay
DefaultAssay(scRNA2) <- "RNA"  

Idents(scRNA2) <- "annotation"
counts_matrix <- GetAssayData(scRNA2, "RNA","counts")
cluster_sc <- as.factor(scRNA2$annotation)
nUMI_sc    <- scRNA2$total_counts

#sc_cell_filter_idx <-  cluster_sc != "Unknown"

## 2.create the RCTD reference object------
reference <- Reference(
  counts_matrix, 
  droplevels(cluster_sc), 
  nUMI_sc
)


#xenium数据
all_data <- readRDS("all_data_annotation.rds")
## 3. create the RCTD query object-----
all_data_sp_rctd_ls <- all_data %>% 
  JoinLayers() %>% 
  SplitObject('project')

cells_position_ls <- readRDS("~/data_for_analysis/cells_position_ls.rds")

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords    <- 
      cells_position_ls[[name]] %>% 
      dplyr::select(
        cell_id,
        x = x_centroid, 
        y = y_centroid
      ) %>%
      column_to_rownames('cell_id')
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })
saveRDS(rctd_query_ls,"RCTD/rctd_query_ls.rds")


##4. run RCTD-----
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, 'RCTD/rctd_result_ls.rds')

#rctd_result_ls <- readRDS("RCTD/rctd_result_ls.rds")

##5. addMetaData-----
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

saveRDS(all_data_sp_rctd_anno_ls,'RCTD/rctd.rds')

#all_data_sp_rctd_anno_ls <- readRDS("rctd.rds")


#6. dimplot 
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )

meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_sp_rctd_anno <-
  all_data %>% 
  AddMetaData(meta_rctd_result)


all_data2 <- JoinLayers(all_data_sp_rctd_anno)


xenium <- readRDS("result/xenium_subtype.rds")
# all_data2_metadata <- all_data2[[]]
xenium_RCTD <- xenium
xenium_RCTD@meta.data %>% head()
xenium_RCTD@meta.data$cell_id1 %>% head()
xenium_RCTD@meta.data$cell_id <- xenium_RCTD@meta.data$cell_id1 
xenium_RCTD@meta.data$sample_cell_id <- rownames(xenium_RCTD@meta.data)

all_data2@meta.data %>% head()
all_data2@meta.data$cell_id <- gsub("1_","",rownames(all_data2@meta.data))

source_meta <- all_data2@meta.data[, c("cell_id", "first_type")]

xenium_RCTD@meta.data <- merge(
  xenium_RCTD@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  # 保留所有目标对象细胞
)
rownames(xenium_RCTD@meta.data) <- xenium_RCTD@meta.data$sample_cell_id 

xenium_RCTD@meta.data %>% head()

saveRDS(xenium_RCTD,"RCTD/xenium_RCTD.rds")

xenium_RCTD$first_type <- as.character(xenium_RCTD$first_type)
xenium_RCTD$first_type[is.na(xenium_RCTD$first_type)] <- "other"

xenium_RCTD$first_type <- factor(xenium_RCTD$first_type)

xenium_RCTD$first_type %>% table()

#反卷积结果验证--------
cell_type_ordered <- c(
  # 内皮细胞 Endothelial
  "ArtECs", "CapECs", "EndoMT", "LECs", "TipECs", "VenECs",
  
  # 成纤维细胞 / CAF
  "Fibro_COL15A1", "Fibro_HSPA1B", "Fibro_ICAM1", "Fibro_IL6",
  "Fibro_PI16_MEOX2", "Fibro_PI16_MFAP4", "Fibro_POSTN",
  "Fibro_RPL36", "Fibro_like_ACTG2",
  
  # 基质细胞 Stromal
  "Pericyte_RGS5", "SMC_MYH11",
  
  # T 细胞
  "CD4_Tem_GZMK", "CD4_TfhTh1_CXCL13", "CD4_Th17_IL17A",
  "CD4_Tm_ANXA1", "CD4_Tm_CREM", "CD4_Tm_TXNIP", "CD4_Tn_SELL", "CD4_Treg_FOXP3",
  "CD8_ISG", "CD8_MAIT_SLC4A10", "CD8_MT1X", "CD8_Tem_GZMK",
  "CD8_Temra_CX3CR1", "CD8_Tex_CXCL13", "CD8_Tm_ANXA1", "CD8_Tm_IL7R",
  "CD8_Trm_GZMB", "CD8_Trm_ITGA1", "CD8_Trm_ZNF683", "CD8_Tstr_HSPA1A", "CD8_p-Tex_TOX",
  
  # NK / ILC
  "ILC3_KIT", "NK_CD16hi_CX3CR1", "NK_CD16hi_GZMB",
  "NK_CD16hi_RUNX3", "NK_CD16lo_CCL3", "NK_CD16lo_CREM",
  
  # B / Plasma
  "Bgc_RGS13", "Bn_TCL1A", "Plasma",
  
  # 髓系 Myeloid
  "Bm_CCR7", "Bm_CD74", "Bm_GPR183", "Bm_Stress_response",
  "cDC1_CLEC9A", "cDC2_CD1C", "cDC_LAMP3", "Mast",
  "Mo_CD14", "Mo_FCGR3A", "Mph_C1QC", "Mph_CX3CR1",
  "Mph_CXCL10", "Mph_FOLR2", "Mph_GPR183", "Mph_MT1X", "Mph_SPP1",
  "Neu_CSF3R", "pDC_LILRA4",
  
  # 其他
  "Hypoxia", "other"
)

xenium_RCTD$first_type <- as.character(xenium_RCTD$first_type)
xenium_RCTD$first_type <- factor(xenium_RCTD$first_type,levels = cell_type_ordered)

cell_color <- c(
  "ArtECs" = "#1F77B4",
  "CapECs" = "#FF7F0E",
  "EndoMT" = "#2CA02C",
  "LECs" = "#D62728",
  "TipECs" = "#9467BD",
  "VenECs" = "#8C564B",
  "Fibro_COL15A1" = "#E377C2",
  "Fibro_HSPA1B" = "#7F7F7F",
  "Fibro_ICAM1" = "#BCBD22",
  "Fibro_IL6" = "#17BECF",
  "Fibro_PI16_MEOX2" = "#AEC7E8",
  "Fibro_PI16_MFAP4" = "#FFBB78",
  "Fibro_POSTN" = "#98DF8A",
  "Fibro_RPL36" = "#FF9896",
  "Fibro_like_ACTG2" = "#C5B0D5",
  "Pericyte_RGS5" = "#C49C94",
  "SMC_MYH11" = "#F7B6D2",
  "CD4_Tem_GZMK" = "#A4A3E5",
  "CD4_TfhTh1_CXCL13" = "#DBDB8D",
  "CD4_Th17_IL17A" = "#9EDAE5",
  "CD4_Tm_ANXA1" = "#393B79",
  "CD4_Tm_CREM" = "#637939",
  "CD4_Tm_TXNIP" = "#8C6D31",
  "CD4_Tn_SELL" = "#7B4173",
  "CD4_Treg_FOXP3" = "#8E2800",
  "CD8_ISG" = "#5D6970",
  "CD8_MAIT_SLC4A10" = "#191970",
  "CD8_MT1X" = "#006400",
  "CD8_Tem_GZMK" = "#8B008B",
  "CD8_Temra_CX3CR1" = "#FF1493",
  "CD8_Tex_CXCL13" = "#00CED1",
  "CD8_Tm_ANXA1" = "#B22222",
  "CD8_Tm_IL7R" = "#228B22",
  "CD8_Trm_GZMB" = "#4682B4",
  "CD8_Trm_ITGA1" = "#DAA520",
  "CD8_Trm_ZNF683" = "#32CD32",
  "CD8_Tstr_HSPA1A" = "#FF4500",
  "CD8_p-Tex_TOX" = "#9932CC",
  "ILC3_KIT" = "#00FA9A",
  "NK_CD16hi_CX3CR1" = "#CD5C5C",
  "NK_CD16hi_GZMB" = "#F08080",
  "NK_CD16hi_RUNX3" = "#3CB371",
  "NK_CD16lo_CCL3" = "#D2691E",
  "NK_CD16lo_CREM" = "#FF6347",
  "Bgc_RGS13" = "#4B0082",
  "Bn_TCL1A" = "#5F9EA0",
  "Plasma" = "#DC143C",
  "Bm_CCR7" = "#00BFFF",
  "Bm_CD74" = "#FFD700",
  "Bm_GPR183" = "#FF8C00",
  "Bm_Stress_response" = "#ADFF2F",
  "cDC1_CLEC9A" = "#FF00FF",
  "cDC2_CD1C" = "#1E90FF",
  "cDC_LAMP3" = "#00FF7F",
  "Mast" = "#FF69B4",
  "Mo_CD14" = "#8A2BE2",
  "Mo_FCGR3A" = "#2E8B57",
  "Mph_C1QC" = "#DDA0DD",
  "Mph_CX3CR1" = "#F4A460",
  "Mph_CXCL10" = "#808000",
  "Mph_FOLR2" = "#FF0000",
  "Mph_GPR183" = "#00FF00",
  "Mph_MT1X" = "#0000FF",
  "Mph_SPP1" = "#FFFF00",
  "Neu_CSF3R" = "#00FFFF",
  "pDC_LILRA4" = "#808080",
  "Hypoxia" = "#BD2870",
  "other" = "#C7C7C7"
)

p <- DimPlot(xenium_RCTD,reduction = "harmonyumap",group.by = "first_type",cols =cell_color )
ggsave(p, file="RCTD/RCTD_umap.pdf",width = 25, height = 10)

major_marker_ls <- c(
  # 内皮细胞 Endothelial
  "Endothelials" = "PECAM1",
  "Endothelials" = "ACKR1",
  "Endothelials" = "CLDN5",
  "Endothelials" = "RAMP2",
  "Endothelials" = "TM4SF1",
  "Endothelials" = "RNASE1",
  "Endothelials" = "EGFL7",
  "Endothelials" = "RAMP3",
  "Endothelials" = "PLVAP",
  "Endothelials" = "AQP1",
  "Endothelials" = "ECSCR",
  "Endothelials" = "FKBP1A",
  "Endothelials" = "AC011526.1",
  "Endothelials" = "EMP1",
  "Endothelials" = "DARC",
  "Endothelials" = "VWF",
  "Endothelials" = "EMCN",
  
  #  成纤维细胞 / CAF 
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "LUM",
  "Fibroblasts" = "VIM",
  "Fibroblasts" = "IGFBP7",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  "Fibroblasts" = "PDGFRA",
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "IGF1",
  "Fibroblasts" = "APOD",
  "Fibroblasts" = "FBLN1",
  "Fibroblasts" = "MEG3",
  "Fibroblasts" = "CXCL12",
  
  #  基质细胞 Stromal 
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  "Pericytes" = "ACTA2",
  "Pericytes" = "MT1M",
  "Pericytes" = "FRZB",
  "Pericytes" = "MT1A",
  "Pericytes" = "NDUFA4L2",
  "Pericytes" = "PPP1R14A",
  "Pericytes" = "MYLK",
  "Pericytes" = "PHLDA1",
  
  "Smooth muscle cells" = "MYH11",
  "Smooth muscle cells" = "TAGLN",
  "Smooth muscle cells" = "CNN1",
  "Smooth muscle cells" = "DES",
  "Smooth muscle cells" = "ACTG2",
  
  #  T 细胞 
  "T cells" = "CD3D",
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "TRAC",
  "T cells" = "NKG7",
  "T cells" = "KLRC1",
  "T cells" = "CCR7",
  "T cells" = "FOXP3",
  "T cells" = "CLTA4",
  "T cells" = "CD8B",
  "T cells" = "XCR6",
  "T cells" = "IL32",
  "T cells" = "CXCR6",
  
  #  NK / ILC 
  "NK cells" = "GNLY",
  "NK cells" = "NCAM1",
  "NK cells" = "GZMB",
  "NK cells" = "KLRD1",
  
  #  B / Plasma 
  "B cells" = "MS4A1",
  "B cells" = "CD79A",
  "B cells" = "CD19",
  "B cells" = "BTG1",
  "B cells" = "VPREB3",
  "B cells" = "BANK1",
  "B cells" = "IGLL5",
  "B cells" = "CD79B",
  "B cells" = "IGHM",
  "B cells" = "CD22",
  
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  "Plasma cells" = "TNFRSF17",
  "Plasma cells" = "JCHAIN",
  "Plasma cells" = "MZB1",
  "Plasma cells" = "XBP1",
  "Plasma cells" = "SEC11C",
  "Plasma cells" = "PRDX4",
  "Plasma cells" = "SPCS2",
  "Plasma cells" = "SSR3",
  "Plasma cells" = "SDF2L1",
  "Plasma cells" = "C19orf10",
  "Plasma cells" = "MANF",
  "Plasma cells" = "TMEM258",
  "Plasma cells" = "DNAJB9",
  "Plasma cells" = "CD138",
  "Plasma cells" = "CD38",
  "Plasma cells" = "BCMA",
  "Plasma cells" = "IRF4",
  "Plasma cells" = "BLIMP1",
  
  #  髓系 Myeloid 
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC",
  "Macrophages" = "CSF1R",
  "Macrophages" = "C1QB",
  "Macrophages" = "GPR34",
  "Macrophages" = "MS4A4A",
  
  "Monocytes" = "CD14",
  "Monocytes" = "FCGR3A",
  "Monocytes" = "S100A9",
  "Monocytes" = "FCN1",
  "Monocytes" = "S100A8",
  "Monocytes" = "CSTA",
  "Monocytes" = "EREG",
  "Monocytes" = "NEAT1",
  "Monocytes" = "TKT",
  "Monocytes" = "THBS1",
  "Monocytes" = "VCAN",
  "Monocytes" = "TSPO",
  "Monocytes" = "LYZ",
  
  "Dendritic cell" = "CD83",
  "Dendritic cell" = "LILRA4",
  "Dendritic cell" = "CLEC9A",
  "Dendritic cell" = "IDO1",
  "Dendritic cell" = "CD1C",
  "Dendritic cell" = "FCER1A",
  
  "Mast cells" = "KIT",
  "Mast cells" = "TPSAB1",
  "Mast cells" = "CPA3",
  "Mast cells" = "HPGDS",
  "Mast cells" = "RGS13",
  "Mast cells" = "IL1RL1",
  "Mast cells" = "GATA2",
  "Mast cells" = "FOSB",
  "Mast cells" = "RP11-354E11.2",
  "Mast cells" = "VWA5A",
  "Mast cells" = "SLC18A2",
  "Mast cells" = "HDC",
  "Mast cells" = "CAPG",
  "Mast cells" = "FCER1G",
  
  "Netrophils" = "CSF3R",
  "Netrophils" = "G0S2",
  "Netrophils" = "RETNLG",
  
  #  其他 
  "Epithelial cells" = "EPCAM",
  "Epithelial cells" = "KRT5",
  "Epithelial cells" = "KRT14",
  
  "Luminal cells" = "KRT8",
  "Luminal cells" = "KRT18",
  "Luminal cells" = "NKX3-1",
  "Luminal cells" = "KLK4",
  "Luminal cells" = "KLK3",
  "Luminal cells" = "KLK2",
  "Luminal cells" = "MSMB",
  "Luminal cells" = "ACPP"
)

# dotplot
pdf(file.path(rctd_outdir, str_glue("dotplot major_marker.pdf")), width = 60, height = 25)
DotPlot(xenium_RCTD, group.by = "first_type", features = major_marker_ls)+RotatedAxis()
dev.off()


# find markers -------
Idents(xenium_RCTD) <- "first_type"
all_markers <-
  FindAllMarkers(
    xenium_RCTD,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

top5_markers <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 30) %>% 
  arrange(desc(pct.1)) %>% 
  slice_head(n=5) %>% 
  nest() %>% 
  deframe() %>% 
  map(pull,'gene') %>% 
  imap(function(genes,name){
    re <- genes %>% 
      map(function(gene){
        names(gene) = str_glue("{name}")
        gene
      })
    re
  }) %>% 
  unlist() %>% 
  .[!duplicated(.)]

names(top5_markers) <- sub('.*\\.','',names(top5_markers))

# dotplot
pdf(file.path(rctd_outdir, str_glue("dotplot top5_markers.pdf")), width = 100, height = 25)
DotPlot(xenium_RCTD, group.by = "first_type", features = top5_markers)+RotatedAxis()
dev.off()


#关注细胞类型空间映射---------
Imagedimplot_dir <- 'result/Imagedimplot'
dir.create(Imagedimplot_dir)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928"
)

all_colors <- c(
  "Fibro_POSTN" = "#FF6B6B",
  "TipECs" = "#A6CEE3",
  "CD4_Treg_FOXP3" = "#00F5A7", 
  "CD4_TfhTh1_CXCL13"  = "#FDBB69",
  "CD8_Tex_CXCL13" = "#feee7d", 
  "Mph_SPP1" = "#00BBF9"
)

cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1"
)

sample_label <- unique(xenium_RCTD@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_RCTD@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_RCTD,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.05,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("RCTD/Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("RCTD/Imagedimplot/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
}

#邻近度分析-------
cells_boundary_ls <- readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "RCTD/proximity_out_dir"
dir.create(proximity_out_dir)

xenium_RCTD@meta.data %>% head()
xenium_RCTD@meta.data$cell_id %>% head()
rownames(xenium_RCTD@meta.data) %>% head()
xenium_RCTD@meta.data$cell_id <- paste0("1_",xenium_RCTD@meta.data$cell_id1)

table(xenium_RCTD$group1,xenium_RCTD$sample)
###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
## Proximity analysis -------
## 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- xenium_RCTD[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(first_type, sample, "group1","cell_id") %>% 
      rownames_to_column("sample_cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(first_type = factor(first_type)) %>%
      dplyr::rename(Classification1 = first_type,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group1
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,"RCTD/proximity_out_dir/all_data_meta_pos.csv")
table(all_data_meta_pos$Slide_ID,all_data_meta_pos$Group)

all_data_meta_pos <- read.csv("RCTD/proximity_out_dir/all_data_meta_pos.csv")

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


## 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 5)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path("RCTD/proximity_out_dir/", "proximity analysis.pdf"),
  height = 10,
  width  = 20
)
p_res_pa %>% print()
dev.off()

table(res_pa_df$Target,res_pa_df$Group)


all_colors <- c(
  "Fibro_POSTN" = "#FF6B6B",
  "TipECs" = "#A6CEE3",
  "CD4_Treg_FOXP3" = "#00F5A7", 
  "CD4_TfhTh1_CXCL13"  = "#FDBB69",
  "CD8_Tex_CXCL13" = "#feee7d", 
  "Mph_SPP1" = "#00BBF9"
)

cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1"
)

###关注细胞类型-------
lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    dplyr::filter(Others %in% cell_list) %>% 
    dplyr::filter(Target == cell) %>%
    
    ggplot(aes(x = Others, y = cell_props, fill = Group
               )) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path("RCTD/proximity_out_dir/", paste0(cell,"_proximity analysis.pdf")),
    height = 5,
    width  = 5
  )
  p_res_pa %>% print()
  dev.off()
})



#0409analysis--------
result_dir <- "RCTD/result_noEPI"
dir.create(result_dir)
#剔除大类注释的上皮-------
xenium_RCTD@meta.data$annotation %>% table()

xenium_RCTD_noEPI <- subset(xenium_RCTD, annotation != "Epithelial")

xenium_RCTD_noEPI$first_type %>% unique()

##1、反卷积结果验证--------
###反卷积亚类排序-------
first_type_ordered <- c(
  # 内皮细胞 Endothelial
  "ArtECs", "CapECs", "EndoMT", "LECs", "TipECs", "VenECs",
  
  # 成纤维细胞 / CAF
  "Fibro_COL15A1", "Fibro_HSPA1B", "Fibro_ICAM1", "Fibro_IL6",
  "Fibro_PI16_MEOX2", "Fibro_PI16_MFAP4", "Fibro_POSTN",
  "Fibro_RPL36", "Fibro_like_ACTG2",
  
  # 基质细胞 Stromal
  "Pericyte_RGS5", "SMC_MYH11",
  
  # T 细胞
  "CD4_Tem_GZMK", "CD4_TfhTh1_CXCL13", "CD4_Th17_IL17A",
  "CD4_Tm_ANXA1", "CD4_Tm_CREM", "CD4_Tm_TXNIP", "CD4_Tn_SELL", "CD4_Treg_FOXP3",
  "CD8_ISG", "CD8_MAIT_SLC4A10", "CD8_MT1X", "CD8_Tem_GZMK",
  "CD8_Temra_CX3CR1", "CD8_Tex_CXCL13", "CD8_Tm_ANXA1", "CD8_Tm_IL7R",
  "CD8_Trm_GZMB", "CD8_Trm_ITGA1", "CD8_Trm_ZNF683", "CD8_Tstr_HSPA1A", "CD8_p-Tex_TOX",
  
  # NK / ILC
  "ILC3_KIT", "NK_CD16hi_CX3CR1", "NK_CD16hi_GZMB",
  "NK_CD16hi_RUNX3", "NK_CD16lo_CCL3", "NK_CD16lo_CREM",
  
  # B / Plasma
  "Bgc_RGS13", "Bn_TCL1A", "Plasma",
  
  # 髓系 Myeloid
  "Bm_CCR7", "Bm_CD74", "Bm_GPR183", "Bm_Stress_response",
  "cDC1_CLEC9A", "cDC2_CD1C", "cDC_LAMP3", "Mast",
  "Mo_CD14", "Mo_FCGR3A", "Mph_C1QC", "Mph_CX3CR1",
  "Mph_CXCL10", "Mph_FOLR2", "Mph_GPR183", "Mph_MT1X", "Mph_SPP1",
  "Neu_CSF3R", "pDC_LILRA4",
  
  # 其他
  "Hypoxia", "other"
)

xenium_RCTD_noEPI$first_type <- as.character(xenium_RCTD_noEPI$first_type)
xenium_RCTD_noEPI$first_type <- factor(xenium_RCTD_noEPI$first_type,levels = first_type_ordered)
xenium_RCTD_noEPI$first_type %>% unique()

###合并大类------
xenium_RCTD_noEPI$celltype <- case_when(
  # 1. 内皮
  xenium_RCTD_noEPI$first_type %in% c("ArtECs","CapECs","EndoMT","LECs","TipECs","VenECs") ~ "Endothelial",
  
  # 2. 成纤维
  xenium_RCTD_noEPI$first_type %in% c("Fibro_COL15A1","Fibro_HSPA1B","Fibro_ICAM1","Fibro_IL6",
                    "Fibro_PI16_MEOX2","Fibro_PI16_MFAP4","Fibro_POSTN",
                    "Fibro_RPL36","Fibro_like_ACTG2") ~ "Fibroblast",
  
  # 3. 周细胞
  xenium_RCTD_noEPI$first_type %in% c("Pericyte_RGS5") ~ "Pericyte",
  
  # 4. 平滑肌
  xenium_RCTD_noEPI$first_type %in% c("SMC_MYH11") ~ "SMC",
  
  # 5. T 细胞（全部 CD4/CD8/NK/ILC 都归为 Tcell）
  xenium_RCTD_noEPI$first_type %in% c(
    "CD4_Tem_GZMK","CD4_TfhTh1_CXCL13","CD4_Th17_IL17A",
    "CD4_Tm_ANXA1","CD4_Tm_CREM","CD4_Tm_TXNIP","CD4_Tn_SELL","CD4_Treg_FOXP3",
    "CD8_ISG","CD8_MAIT_SLC4A10","CD8_MT1X","CD8_Tem_GZMK",
    "CD8_Temra_CX3CR1","CD8_Tex_CXCL13","CD8_Tm_ANXA1","CD8_Tm_IL7R",
    "CD8_Trm_GZMB","CD8_Trm_ITGA1","CD8_Trm_ZNF683","CD8_Tstr_HSPA1A","CD8_p-Tex_TOX",
    "ILC3_KIT","NK_CD16hi_CX3CR1","NK_CD16hi_GZMB",
    "NK_CD16hi_RUNX3","NK_CD16lo_CCL3","NK_CD16lo_CREM"
  ) ~ "Tcell",
  
  # 6. B 细胞 / 浆细胞
  xenium_RCTD_noEPI$first_type %in% c("Bgc_RGS13","Bn_TCL1A","Plasma",
                    "Bm_CCR7","Bm_CD74","Bm_GPR183","Bm_Stress_response") ~ "Bcell",
  
  # 7. 髓系
  xenium_RCTD_noEPI$first_type %in% c(
    "cDC1_CLEC9A","cDC2_CD1C","cDC_LAMP3","Mast",
    "Mo_CD14","Mo_FCGR3A","Mph_C1QC","Mph_CX3CR1",
    "Mph_CXCL10","Mph_FOLR2","Mph_GPR183","Mph_MT1X","Mph_SPP1",
    "Neu_CSF3R","pDC_LILRA4"
  ) ~ "Myeloid",
  
  # 8. 上皮 + 其他未归类
  TRUE ~ "other"
)

# 设置因子顺序（画图顺序固定）
xenium_RCTD_noEPI$celltype <- factor(xenium_RCTD_noEPI$celltype, 
                               levels = c("Endothelial","Fibroblast",
                                         "Pericyte","SMC","Tcell","Bcell", "Myeloid","other"))

####大类亚类umap-------
p <- DimPlot(xenium_RCTD_noEPI,reduction = "harmonyumap",group.by = "celltype",#cols =cell_color
)
ggsave(
  plot = p,
  filename = file.path(result_dir, "RCTD_celltype_umap.pdf"), 
  width = 15,
  height = 10
)

p <- DimPlot(xenium_RCTD_noEPI,reduction = "harmonyumap",group.by = "first_type",#cols =cell_color
             )
ggsave(
  plot = p,
  filename = file.path(result_dir, "RCTD_first_type_umap.pdf"), 
  width = 25,
  height = 10
)

####marker大类亚类dotplot-------
major_marker_ls <- c(
  # 内皮细胞 Endothelial
  "Endothelials" = "PECAM1",
  "Endothelials" = "ACKR1",
  "Endothelials" = "CLDN5",
  "Endothelials" = "RAMP2",
  "Endothelials" = "TM4SF1",
  "Endothelials" = "RNASE1",
  "Endothelials" = "EGFL7",
  "Endothelials" = "RAMP3",
  "Endothelials" = "PLVAP",
  "Endothelials" = "AQP1",
  "Endothelials" = "ECSCR",
  "Endothelials" = "FKBP1A",
  "Endothelials" = "AC011526.1",
  "Endothelials" = "EMP1",
  "Endothelials" = "DARC",
  "Endothelials" = "VWF",
  "Endothelials" = "EMCN",
  
  #  成纤维细胞 / CAF 
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "LUM",
  "Fibroblasts" = "VIM",
  "Fibroblasts" = "IGFBP7",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  "Fibroblasts" = "PDGFRA",
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "IGF1",
  "Fibroblasts" = "APOD",
  "Fibroblasts" = "FBLN1",
  "Fibroblasts" = "MEG3",
  "Fibroblasts" = "CXCL12",
  
  #  基质细胞 Stromal 
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  "Pericytes" = "ACTA2",
  "Pericytes" = "MT1M",
  "Pericytes" = "FRZB",
  "Pericytes" = "MT1A",
  "Pericytes" = "NDUFA4L2",
  "Pericytes" = "PPP1R14A",
  "Pericytes" = "MYLK",
  "Pericytes" = "PHLDA1",
  
  "Smooth muscle cells" = "MYH11",
  "Smooth muscle cells" = "TAGLN",
  "Smooth muscle cells" = "CNN1",
  "Smooth muscle cells" = "DES",
  "Smooth muscle cells" = "ACTG2",
  
  #  T 细胞 
  "T cells" = "CD3D",
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "TRAC",
  "T cells" = "NKG7",
  "T cells" = "KLRC1",
  "T cells" = "CCR7",
  "T cells" = "FOXP3",
  "T cells" = "CLTA4",
  "T cells" = "CD8B",
  "T cells" = "XCR6",
  "T cells" = "IL32",
  "T cells" = "CXCR6",
  
  #  NK / ILC 
  "NK cells" = "GNLY",
  "NK cells" = "NCAM1",
  "NK cells" = "GZMB",
  "NK cells" = "KLRD1",
  
  #  B / Plasma 
  "B cells" = "MS4A1",
  "B cells" = "CD79A",
  "B cells" = "CD19",
  "B cells" = "BTG1",
  "B cells" = "VPREB3",
  "B cells" = "BANK1",
  "B cells" = "IGLL5",
  "B cells" = "CD79B",
  "B cells" = "IGHM",
  "B cells" = "CD22",
  
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  "Plasma cells" = "TNFRSF17",
  "Plasma cells" = "JCHAIN",
  "Plasma cells" = "MZB1",
  "Plasma cells" = "XBP1",
  "Plasma cells" = "SEC11C",
  "Plasma cells" = "PRDX4",
  "Plasma cells" = "SPCS2",
  "Plasma cells" = "SSR3",
  "Plasma cells" = "SDF2L1",
  "Plasma cells" = "C19orf10",
  "Plasma cells" = "MANF",
  "Plasma cells" = "TMEM258",
  "Plasma cells" = "DNAJB9",
  "Plasma cells" = "CD138",
  "Plasma cells" = "CD38",
  "Plasma cells" = "BCMA",
  "Plasma cells" = "IRF4",
  "Plasma cells" = "BLIMP1",
  
  #  髓系 Myeloid 
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC",
  "Macrophages" = "CSF1R",
  "Macrophages" = "C1QB",
  "Macrophages" = "GPR34",
  "Macrophages" = "MS4A4A",
  
  "Monocytes" = "CD14",
  "Monocytes" = "FCGR3A",
  "Monocytes" = "S100A9",
  "Monocytes" = "FCN1",
  "Monocytes" = "S100A8",
  "Monocytes" = "CSTA",
  "Monocytes" = "EREG",
  "Monocytes" = "NEAT1",
  "Monocytes" = "TKT",
  "Monocytes" = "THBS1",
  "Monocytes" = "VCAN",
  "Monocytes" = "TSPO",
  "Monocytes" = "LYZ",
  
  "Dendritic cell" = "CD83",
  "Dendritic cell" = "LILRA4",
  "Dendritic cell" = "CLEC9A",
  "Dendritic cell" = "IDO1",
  "Dendritic cell" = "CD1C",
  "Dendritic cell" = "FCER1A",
  
  "Mast cells" = "KIT",
  "Mast cells" = "TPSAB1",
  "Mast cells" = "CPA3",
  "Mast cells" = "HPGDS",
  "Mast cells" = "RGS13",
  "Mast cells" = "IL1RL1",
  "Mast cells" = "GATA2",
  "Mast cells" = "FOSB",
  "Mast cells" = "RP11-354E11.2",
  "Mast cells" = "VWA5A",
  "Mast cells" = "SLC18A2",
  "Mast cells" = "HDC",
  "Mast cells" = "CAPG",
  "Mast cells" = "FCER1G",
  
  "Netrophils" = "CSF3R",
  "Netrophils" = "G0S2",
  "Netrophils" = "RETNLG"
)

# dotplot
pdf(file.path(result_dir, str_glue("dotplot major_marker celltype.pdf")), width = 20, height = 10)
DotPlot(xenium_RCTD_noEPI, group.by = "celltype", features = major_marker_ls)+RotatedAxis()
dev.off()

# dotplot
pdf(file.path(result_dir, str_glue("dotplot major_marker first_type.pdf")), width = 60, height = 25)
DotPlot(xenium_RCTD_noEPI, group.by = "first_type", features = major_marker_ls)+RotatedAxis()
dev.off()

####topgene 大类亚类dotplot-------
#find markers
Idents(xenium_RCTD_noEPI) <- "celltype"
all_markers <-
  FindAllMarkers(
    xenium_RCTD_noEPI,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

top5_markers <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 30) %>% 
  arrange(desc(pct.1)) %>% 
  slice_head(n=5) %>% 
  nest() %>% 
  deframe() %>% 
  map(pull,'gene') %>% 
  imap(function(genes,name){
    re <- genes %>% 
      map(function(gene){
        names(gene) = str_glue("{name}")
        gene
      })
    re
  }) %>% 
  unlist() %>% 
  .[!duplicated(.)]

names(top5_markers) <- sub('.*\\.','',names(top5_markers))

# dotplot
pdf(file.path(result_dir, str_glue("dotplot top5_markers celltype.pdf")), width = 15, height = 10)
DotPlot(xenium_RCTD_noEPI, group.by = "celltype", features = top5_markers)+RotatedAxis()
dev.off()

#find markers
Idents(xenium_RCTD_noEPI) <- "first_type"
all_markers <-
  FindAllMarkers(
    xenium_RCTD_noEPI,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

top5_markers <- 
  all_markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 30) %>% 
  arrange(desc(pct.1)) %>% 
  slice_head(n=5) %>% 
  nest() %>% 
  deframe() %>% 
  map(pull,'gene') %>% 
  imap(function(genes,name){
    re <- genes %>% 
      map(function(gene){
        names(gene) = str_glue("{name}")
        gene
      })
    re
  }) %>% 
  unlist() %>% 
  .[!duplicated(.)]

names(top5_markers) <- sub('.*\\.','',names(top5_markers))

# dotplot
pdf(file.path(result_dir, str_glue("dotplot top5_markers first_type.pdf")), width = 100, height = 25)
DotPlot(xenium_RCTD_noEPI, group.by = "first_type", features = top5_markers)+RotatedAxis()
dev.off()


##2、关注细胞类型空间映射---------
Imagedimplot_dir <- 'RCTD/result_noEPI/Imagedimplot/'
dir.create(Imagedimplot_dir)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928"
)

all_colors <- c(
  "Fibro_POSTN" = "#FF6B6B",
  "TipECs" = "#A6CEE3",
  "CD4_Treg_FOXP3" = "#00F5A7", 
  "CD4_TfhTh1_CXCL13"  = "#FDBB69",
  "CD8_Tex_CXCL13" = "#feee7d", 
  "Mph_SPP1" = "#00BBF9"
)

cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1"
)

sample_label <- unique(xenium_RCTD_noEPI@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_RCTD_noEPI@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_RCTD_noEPI,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.05,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0(Imagedimplot_dir,i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0(Imagedimplot_dir,i,"_anno_spatial.pdf"),width = 12 ,height = 12)
}

##3、细胞比例统计--------
library(ggplot2)
library(ggpubr)   
library(dplyr)    

plot_data <- xenium_RCTD_noEPI[[]] %>% 
  dplyr::select("sample", "group1", "first_type", "celltype")

###大类箱线图--------
boxplot_data <- plot_data %>%
  filter(celltype != "other") %>% 
  group_by(sample, group1, celltype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>% 
  mutate(group = factor(group1, levels = c("PERI","BPH","ADPC","CRPC")))

group_colors <- c(
  "ADPC" = "#1f77b4",    
  "BPH" = "#ff7f0e",   
  "CRPC" = "#2ca02c",
  "PERI" = "#B294C7"
)

p <- boxplot_data %>%
  ggplot(aes(x = celltype, y = freq, fill = group)) +
  geom_boxplot(width = 0.6) +
  #facet_wrap(~celltype, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    title = "",
    x = "celltype",
    y = "freq"
  ) +
  scale_fill_manual(values = group_colors)  

print(p)

ggsave(file.path(result_dir,paste0("boxplot_celltype/","celltype_boxplot.pdf")), p, width = 10, height = 5)
ggsave(file.path(result_dir,paste0("boxplot_celltype/","celltype_boxplot.png")), p, width = 10, height = 5)


# 循环绘制箱线图
for (i in unique(boxplot_data$celltype)) {
  
  dat <- boxplot_data %>% 
    filter(celltype == i)

  p <- dat %>%
    ggplot(aes(x = group, y = freq, fill = group)) +
    geom_boxplot(width = 0.6) +
    #facet_wrap(~celltype, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "",
      x = "group",
      y = "freq"
    ) +
    scale_fill_manual(values = group_colors)  
  
  print(p)
  
  ggsave(file.path(result_dir,paste0("boxplot_celltype/",i,"_boxplot.pdf")), p, width = 6, height = 5)
  ggsave(file.path(result_dir,paste0("boxplot_celltype/",i,"_boxplot.png")), p, width = 6, height = 5)

}

###亚箱线图--------
boxplot_data <- plot_data %>%
  filter(celltype != "other") %>% 
  group_by(sample, group1, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>% 
  mutate(group = factor(group1, levels = c("PERI","BPH","ADPC","CRPC")))

group_colors <- c(
  "ADPC" = "#1f77b4",    
  "BPH" = "#ff7f0e",   
  "CRPC" = "#2ca02c",
  "PERI" = "#B294C7"
)

p <- boxplot_data %>%
  ggplot(aes(x = first_type, y = freq, fill = group)) +
  geom_boxplot(width = 0.6) +
  #facet_wrap(~celltype, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    title = "",
    x = "subtype",
    y = "freq"
  ) +
  scale_fill_manual(values = group_colors)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1#, vjust = 0.5
                               )
  )

print(p)

ggsave(file.path(result_dir,paste0("boxplot_subtype/","subtype_boxplot.pdf")), p, width = 40, height = 8)
ggsave(file.path(result_dir,paste0("boxplot_subtype/","subtype_boxplot.png")), p, width = 40, height = 8)


# 循环绘制箱线图
for (i in unique(boxplot_data$first_type)) {
  
  dat <- boxplot_data %>% 
    filter(first_type == i)
  
  p <- dat %>%
    ggplot(aes(x = group, y = freq, fill = group)) +
    geom_boxplot(width = 0.6) +
    #facet_wrap(~celltype, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "",
      x = "group",
      y = "freq"
    ) +
    scale_fill_manual(values = group_colors)  
  
  print(p)
  
  ggsave(file.path(result_dir,paste0("boxplot_subtype/",i,"_boxplot.pdf")), p, width = 6, height = 5)
  ggsave(file.path(result_dir,paste0("boxplot_subtype/",i,"_boxplot.png")), p, width = 6, height = 5)
  
}


##4、邻近度分析------
cells_boundary_ls <- readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "RCTD/result_noEPI/proximity_out_dir"
dir.create(proximity_out_dir)

xenium_RCTD_noEPI@meta.data %>% head()
xenium_RCTD_noEPI@meta.data$cell_id %>% head()
rownames(xenium_RCTD_noEPI@meta.data) %>% head()
#xenium_RCTD_noEPI@meta.data$cell_id <- paste0("1_",xenium_RCTD_noEPI@meta.data$cell_id1)

table(xenium_RCTD_noEPI$group1,xenium_RCTD_noEPI$sample)

saveRDS(xenium_RCTD_noEPI,"RCTD/result_noEPI/xenium_RCTD_noEPI.rds")

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
### Proximity analysis -------
#### 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- xenium_RCTD_noEPI[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(first_type, sample, "group1","cell_id") %>% 
      rownames_to_column("sample_cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(first_type = factor(first_type)) %>%
      dplyr::rename(Classification1 = first_type,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group1
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path(proximity_out_dir,"all_data_meta_pos.csv"))
table(all_data_meta_pos$Slide_ID,all_data_meta_pos$Group)

#all_data_meta_pos <- read.csv(file.path(proximity_out_dir,"all_data_meta_pos.csv"))

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


#### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞


res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 5)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis.pdf"),
  height = 10,
  width  = 20
)
p_res_pa %>% print()
dev.off()

table(res_pa_df$Target,res_pa_df$Group)


all_colors <- c(
  "Fibro_POSTN" = "#FF6B6B",
  "TipECs" = "#A6CEE3",
  "CD4_Treg_FOXP3" = "#00F5A7", 
  "CD4_TfhTh1_CXCL13"  = "#FDBB69",
  "CD8_Tex_CXCL13" = "#feee7d", 
  "Mph_SPP1" = "#00BBF9"
)



####以Fibro_POSTN为核心-------
cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs"#,
               #"Mph_SPP1"
)

lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    dplyr::filter(Others %in% cell_list) %>% 
    dplyr::filter(Others != "Fibro_POSTN") %>% 
    dplyr::filter(Target == cell) %>%
    
    ggplot(aes(x = Others, y = cell_props, fill = Group
    )) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis.pdf")),
    height = 5,
    width  = 5
  )
  p_res_pa %>% print()
  dev.off()
})


####以Mph_SPP1为核心-------
cell_list <- c(#"Fibro_POSTN",
  "CD4_Treg_FOXP3",
  "CD4_TfhTh1_CXCL13",
  "CD8_Tex_CXCL13",
  "TipECs",
  "Mph_SPP1"
)

lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    dplyr::filter(Others %in% cell_list) %>%
    dplyr::filter(Others != "Mph_SPP1") %>% 
    dplyr::filter(Target == cell) %>%
    
    ggplot(aes(x = Others, y = cell_props, fill = Group
    )) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis.pdf")),
    height = 5,
    width  = 5
  )
  p_res_pa %>% print()
  dev.off()
})


#0410analysis------
xenium_RCTD_noEPI <- readRDS("RCTD/result_noEPI/xenium_RCTD_noEPI.rds")

source("spatial analysis.R")
##1、邻域分析--------
CN_out_dir <- "RCTD/result_noEPI/CN analysis"
dir.create(CN_out_dir)

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

## CN analysis ---------
# options(warn = 0)
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  xenium_RCTD_noEPI[[c("sample", "group1", "first_type","cell_id")]] %>% 
  rownames_to_column("sample_cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    first_type = factor(first_type,level = first_type_ordered), # import!
    obj_id = cell_id, # import!
    group = group1
  )

all_data_submeta_for_cn_analysis %>% head()

write.csv(all_data_submeta_for_cn_analysis,file.path(CN_out_dir,"all_data_submeta_for_cn_analysis.csv"))
#all_data_submeta_for_cn_analysis<- read.csv(file.path(CN_out_dir,"all_data_submeta_for_cn_analysis.csv")

### CN knn10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "first_type",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_knn10.rds"))

#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


#### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  xenium_RCTD_noEPI[[]] %>%
  #rownames_to_column('no_cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(CN_out_dir,'all_data_meta_cn10.rds'))

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


#### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(sample, group1, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count),
         group = factor(group1, levels = c("PERI","BPH","ADPC","CRPC"))
  )

p_cn10_props <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### CN knn50 ---------
res_cn_50 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "first_type",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_50, file.path(CN_out_dir, "res_cn_50.rds"))

#### heatmap ----------------------
# FoldChange cluster/allCell
cn50_cluster_mean <- res_cn_50$cluster_centers
cn50_baseMean     <- res_cn_50$base_mean

# log2(FC)
cn50_cluster_FC <- log2(1 + sweep(cn50_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn50_baseMean))

# heatmap
p_cn50_heat <-
  pheatmap::pheatmap(
    cn50_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn50_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_50_heatmap.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn50_heat$gtable)
dev.off()


#### CN50 spatial plot -----------------
cn50_sp_out_dir <- file.path(CN_out_dir, 'CN50 spatial plot') 
dir.create(cn50_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn50 <- 
  xenium_RCTD_noEPI[[]] %>%
  #rownames_to_column('no_cell_id') %>%
  # left_join(rownames_to_column(res_cn_50$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_50$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn50, file.path(CN_out_dir,'all_data_meta_cn50.rds'))

cn50_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn50[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn50, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn50_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn50_sp_out_dir, str_glue("cn50 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


#### cell props ------
All_type_cn50_count_df <-
  All_type_meta_cn50 %>% 
  group_by(sample, group1, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count),
         group = factor(group1, levels = c("PERI","BPH","ADPC","CRPC"))
  )

p_cn50_props <-
  All_type_cn50_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn50_props_group <-
  All_type_cn50_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "cn50_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn50_props %>% print
p_cn50_props_group %>% print
dev.off()

##2、邻近度分析-----------
cells_boundary_ls <- readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "RCTD/result_noEPI/proximity_out_dir"
dir.create(proximity_out_dir)

xenium_RCTD_noEPI@meta.data %>% head()
xenium_RCTD_noEPI@meta.data$cell_id %>% head()
rownames(xenium_RCTD_noEPI@meta.data) %>% head()
#xenium_RCTD_noEPI@meta.data$cell_id <- paste0("1_",xenium_RCTD_noEPI@meta.data$cell_id1)

table(xenium_RCTD_noEPI$group1,xenium_RCTD_noEPI$sample)

cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1",
               "Mph_FOLR2"
)

xenium_RCTD_noEPI$subtype7 <- case_when(
  xenium_RCTD_noEPI$first_type %in% cell_list ~ xenium_RCTD_noEPI$first_type,
  TRUE ~ xenium_RCTD_noEPI$celltype
) 

subtype7_level <- c("Fibro_POSTN","CD4_Treg_FOXP3","CD4_TfhTh1_CXCL13","CD8_Tex_CXCL13","TipECs","Mph_SPP1","Mph_FOLR2","Endothelial","Fibroblast","Pericyte","SMC","Tcell","Bcell", "Myeloid","other")
xenium_RCTD_noEPI$subtype7 <- factor(xenium_RCTD_noEPI$subtype7,
                                     levels = subtype7_level)
xenium_RCTD_noEPI@meta.data$subtype7 %>% unique()

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
### Proximity analysis -------
#### 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- xenium_RCTD_noEPI[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(subtype7, sample, "group1","cell_id") %>% 
      rownames_to_column("sample_cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(subtype7 = factor(subtype7, levels = subtype7_level)) %>%
      dplyr::rename(Classification1 = subtype7,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group1
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path(proximity_out_dir,"all_data_meta_pos_subtype7.csv"))
table(all_data_meta_pos$Slide_ID,all_data_meta_pos$Group)

#all_data_meta_pos <- read.csv(file.path(proximity_out_dir,"all_data_meta_pos_subtype7.csv"))

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


#### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞
#### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 5)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis.pdf"),
  height = 10,
  width  = 20
)
p_res_pa %>% print()
dev.off()

table(res_pa_df$Target,res_pa_df$Group)


all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928"
)

all_colors <- c(
  "Fibro_POSTN" = "#FF6B6B",
  "TipECs" = "#A6CEE3",
  "CD4_Treg_FOXP3" = "#00F5A7", 
  "CD4_TfhTh1_CXCL13"  = "#FDBB69",
  "CD8_Tex_CXCL13" = "#feee7d", 
  "Mph_SPP1" = "#00BBF9",
  "Mph_FOLR2" = "#2D82AF"
)



#####以七个特定亚群为核心-------
cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1",
               "Mph_FOLR2"
)

res_pa_df$Group <- factor(res_pa_df$Group,levels = c("PERI","BPH","ADPC","CRPC"))

lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    #dplyr::filter(Others %in% cell_list) %>% 
    dplyr::filter(Others != "other") %>% 
    dplyr::filter(Target == cell) %>%
    
    ggplot(aes(x = Others, y = cell_props, fill = Group
    )) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis.pdf")),
    height = 8,
    width  = 8
  )
  p_res_pa %>% print()
  dev.off()
})


#### 2.2 proximity analysis - series analysis -----
init_cells <- c(
  '10'  = 750,
  '50'  = 800,
  '100' = 1000,
  '150' = 1200,
  '200' = 1500
)

res_pa_ls <-
  c(10, 50, 100, 150, 200) %>%
  set_names(., .) %>%
  map(function(x){
    k = init_cells[as.character(x)]
    res_pa <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos,
        query = all_data_meta_pos,
        x = "X",
        y = "Y",
        group.by = "Slide_ID",
        classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
        k = 'auto', #k, # 'auto',
        searchtype = "radius",
        radius = x, # 邻域半径
        mpp = 1, # !important, 1 for xenium
        returntype = "count",
        cell_enough_ratio_cutoff = 90
      )
    res_pa
  })

# res_pa_ls <- readRDS("data/res_pa_ls.rds")

# 将同心圆统计 转换成 环状同心的统计
# 注意一定要倒序操作
# 200 = 200 - 150
# 150 = 150 - 100 
# ...
# 
res_pa_ls2 <- res_pa_ls
diff_idx_ls <-
  c(10, 50, 100, 150, 200) %>% 
  as.character() %>%
  rev() %>%
  {matrix(c(.[-length(.)], .[-1]), ncol = 2)} %>%
  apply(1, c, simplify = FALSE)

for(i in diff_idx_ls){
  first  <- i[1]
  second <- i[2]
  res_pa_ls2[[ first ]] <- res_pa_ls2[[ first ]] - res_pa_ls2[[ second ]] 
}

# count to freq
res_pa_ls2 <- 
  res_pa_ls2 %>%
  map(function(x){# browser()
    x_colsum <- colSums(x)
    xx <- sweep(x, 2, x_colsum, FUN = "/")
    if(any(x_colsum == 0)){
      xx[, x_colsum == 0] <- 0
    }
    xx
  })

# 每个样本中的每种Target的其他细胞的计数取中位数
res_pa_count_tidy <- 
  res_pa_ls2 %>%
  imap(function(x, distance){# browser()
    # 分组
    pa_groups <-
      with(
        all_data_meta_pos,
        str_c(Classification1, Group, Slide_ID, sep = '..')
      )
    # 分组内取中位数
    res <- apply(x, 1, function(y){tapply(y, INDEX = pa_groups, FUN = "median")})
    
    # data tidy
    res_df  <-
      res %>%
      as.data.frame() %>%
      rownames_to_column('group') %>%
      pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
      separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
    res_df$distance <- distance
    
    # return
    res_df
  }) %>%
  bind_rows()



# plot 
# only plot "Epithelial cells" ....
res_pa_count_tidy_epi <-
  res_pa_count_tidy %>%
  mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
  dplyr::filter(Target == "Mph_SPP1") 
  
res_pa_count_tidy_epi_median <-
  res_pa_count_tidy_epi %>%
  group_by(Others, Group, distance) %>% 
  summarise(cell_props = median(cell_props))

p_proximity_series_epi <-
  res_pa_count_tidy_epi %>%
  ggplot(aes(x = distance, y = cell_props, color = Group)) +
  geom_boxplot( ) +
  ggsci::scale_color_aaas() +
  # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
  geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
  geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
  facet_wrap(~ Others, scales = 'free_y') +
  ggpubr::theme_pubclean()


# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis_series distance.pdf"),
  height = 14,
  width  = 18
)
p_proximity_series_epi %>% print()
dev.off()

#####循环七个亚群-------
cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1",
               "Mph_FOLR2"
)
res_pa_count_tidy$Group <- factor(res_pa_count_tidy$Group,levels = c("PERI","BPH","ADPC","CRPC"))

lapply(cell_list,function(cell){
  
  res_pa_count_tidy_epi <-
    res_pa_count_tidy %>%
    mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
    dplyr::filter(Target == cell) %>% 
    dplyr::filter(Others != "other") 
  
  res_pa_count_tidy_epi_median <-
    res_pa_count_tidy_epi %>%
    group_by(Others, Group, distance) %>% 
    summarise(cell_props = median(cell_props))
  
  p_proximity_series_epi <-
    res_pa_count_tidy_epi %>%
    ggplot(aes(x = distance, y = cell_props, color = Group)) +
    geom_boxplot( ) +
    ggsci::scale_color_aaas() +
    # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
    geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
    geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
    facet_wrap(~ Others, scales = 'free_y') +
    ggpubr::theme_pubclean()
  
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis_series distance.pdf")),
    height = 14,
    width  = 18
  )
  p_proximity_series_epi %>% print()
  dev.off()
})

#####循环三个亚群-------
cell_list3 <- c("Fibro_POSTN",
               "TipECs",
               "Mph_SPP1"
)
res_pa_count_tidy$Group <- factor(res_pa_count_tidy$Group,levels = c("PERI","BPH","ADPC","CRPC"))

res_pa_count_tidy1 <- res_pa_count_tidy %>% 
                      dplyr::filter(Others %in% cell_list3)

#cell <- c("Fibro_POSTN")
lapply(cell_list3,function(cell){
  
  res_pa_count_tidy_epi <-
    res_pa_count_tidy1 %>%
    mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
    dplyr::filter(Target == cell) 
     
  res_pa_count_tidy_epi_median <-
    res_pa_count_tidy_epi %>%
    group_by(Others, Group, distance) %>% 
    summarise(cell_props = median(cell_props))
  
  p_proximity_series_epi <-
    res_pa_count_tidy_epi %>%
    ggplot(aes(x = distance, y = cell_props, color = Group)) +
    geom_boxplot( ) +
    ggsci::scale_color_aaas() +
    # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
    geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
    geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
    facet_wrap(~ Others, scales = 'free_y') +
    ggpubr::theme_pubclean()
  
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis_series distance_3types.pdf")),
    height = 6,
    width  = 18
  )
  p_proximity_series_epi %>% print()
  dev.off()
})

##3、最近邻分析-----
Nearest_neighbor_out_dir <- "RCTD/result_noEPI/Nearest_neighbor_out_dir"
dir.create(Nearest_neighbor_out_dir)

all_data_meta_pos 
res_nn_ls <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos, # 临近度分析中计算出来的
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = NULL, # 若为NULL: 提前返回window的细胞近邻的list
    k = 10, # 'auto',
    searchtype = "radius",
    radius = 100,
    mpp = 1, # !important, 1 for xenium
    cell_enough_ratio_cutoff = 20
  )

saveRDS(res_nn_ls, file = file.path(Nearest_neighbor_out_dir, "res_nn_ls.rds"))
# res_nn_ls <- readRDS(file.path(All_type_Nearest_neighbor_out_dir, "res_nn_ls.rds"))


res_nn_df <-
  res_nn_ls %>%
  map(function(x){  #browser()
    cell_idx  <- x$nn.idx[, 2]
    cell_idx[cell_idx == 0] <- NA # some cells have no nn cells
    cell_id   <- x$total_cell_id[cell_idx]
    cell_dists <- x$nn.dists[, 2]
    
    all_data_used <- all_data_meta_pos %>% column_to_rownames("cell_id")
    cell_anno <- all_data_used[cell_id, 'Classification1']
    
    res <-
      data.frame(
        Others = cell_anno,
        distance = cell_dists
      ) %>%
      cbind(all_data_used[x$query_cell_id, c('Classification1', 'Group', 'Slide_ID')]) %>%
      rownames_to_column("cell_id") %>%
      dplyr::rename(Target = Classification1) 
    
    res
  }) %>%
  bind_rows()

if(any(is.na(res_nn_df$Others) )){
  message(str_glue("{sum(is.na(res_nn_df$Others))} cells have no nearest neighbor cell"))
  res_nn_df <- drop_na(res_nn_df)
}


##Target细胞：Mph_SPP1、Fibro_POSTN、TipECs-----
cell_list3 <- c("Fibro_POSTN",
                "TipECs",
                "Mph_SPP1"
)
res_nn_df$Group <- factor(res_nn_df$Group,levels = c("PERI","BPH","ADPC","CRPC"))

### plot 密度图------
p_nn_ls <-
  res_nn_df %>%
  dplyr::filter(Target %in% cell_list3) %>% 
  dplyr::filter(Others != "other") %>% 
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){ # browser()
    x %>%
      ggplot(aes(x = distance, color = Group, fill = Group)) + 
      geom_density(alpha = 0.5) +
      facet_wrap(~ Others) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        x = "distance(um)"
      )
    
  })

# export plot
pdf(
  file.path(Nearest_neighbor_out_dir, "nearest neighbor analysis_density.pdf"),
  height = 14,
  width  = 18
)
p_nn_ls %>% print()
dev.off()


### plot 箱线图-----
p_nn_ls <-
  res_nn_df %>%
  dplyr::filter(Target %in% cell_list3) %>% 
  dplyr::filter(Others != "other") %>% 
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){  
    x %>%
      group_by(Slide_ID, Group, Others) %>%
      summarise(count = n(), .groups = "drop") %>%  
      mutate(freq = count / sum(count)) %>%
      
      ggplot(aes(x = Group, y = freq, fill = Group)) + 
      geom_boxplot(width = 0.6, outlier.size = 1) +  # 箱线图
      facet_wrap(~ Others, 
                 scales = "free_y"
                 ) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        y = "freq"
      )
  })

# export plot
pdf(
  file.path(Nearest_neighbor_out_dir, "nearest_neighbor_analysis_boxplot_free_y.pdf"),
  height = 14,
  width  = 18
)
print(p_nn_ls)
dev.off()


##4、banksy------
cells_boundary_ls <- readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))





xenium_RCTD_noEPI@meta.data %>% head()
xenium_RCTD_noEPI@meta.data$cell_id %>% head()
rownames(xenium_RCTD_noEPI@meta.data) %>% head()


library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)


# HD示例数据：/mnt/NexData_file/home/ztt/test_data/Visium_HD/Colorectal_Cancer/output/data

out_dir <- "RCTD/result_noEPI/banksy"
dir.create(out_dir)


meta_data_anno <- xenium_RCTD_noEPI@meta.data  #read.csv('adata_celltype_obs.csv')
class(xenium_RCTD_noEPI[['Xenium']])

all_data <- xenium_RCTD_noEPI  #readRDS("adata_filter.rds")
#all_data[['RNA']] <- as(all_data[['RNA']], 'Assay5')



# all_data$Y <- Embeddings(all_data[['spatial']])[, 'SPATIAL_1'] 
# all_data$Y <- Embeddings(all_data[['spatial']])[, 'SPATIAL_2'] 

all_data_meta_pos <- readRDS("RCTD/proximity_out_dir/all_data_meta_pos.csv")
all_data@meta.data <- 
  all_data@meta.data %>%
  left_join(
    all_data_meta_pos %>% 
      dplyr::select(sample_cell_id, X, Y), 
    by = "sample_cell_id"
  )
rownames(all_data@meta.data) <- all_data$sample_cell_id



# # coords 平移
# rng_coord <-
#   all_data[[]] %>%
#   group_by(sample) %>%
#   mutate(
#     xmin = min(x),
#     xmax = max(x),
#     rng = xmax - xmin
#   )
# 
# all_data$xx <- all_data$x - rng_coord$xmin + as.numeric(factor(all_data$sample)) * 1.2 * max(rng_coord$rng)
# 
# data.frame(
#   all_data$xx,
#   all_data$y
# ) %>%
#   slice_sample(n = 10000) %>%
#   plot

all_data@meta.data %>% head()
all_data %>% colnames() %>% head()
all_data %>% rownames() %>% head()

### banksy ----------
all_data <- 
  all_data %>% 
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 1000) 


all_data <- 
  RunBanksy(all_data,
            lambda = 0.8, 
            verbose = TRUE,
            assay = "RNA", # "Spatial.008um", 
            slot = "data",
            features = "variable",
            k_geom = 50,
            
            dimx = "x", # use coord 
            dimy = "y",
            
            group = "sample"
  )

DefaultAssay(all_data) <- "BANKSY"

# data_own <- SeuratWrappers:::get_data(all_data, assay = "RNA", slot = 'data', features = 'variable', verbose = F)
# locs <- SeuratWrappers:::get_locs(all_data, dimx = "x", dimy = "y", dimz = NULL,
#                                   ndim = 2, data_own = data_own, group = 'sample', verbose = F)
# 
# 
# locs %>%
#   slice_sample(n = 10000) %>%
#   plot





### Run PCA and UMAP -----------
all_data <- RunPCA(
  all_data, 
  assay = 'BANKSY', 
  features = rownames(all_data[['BANKSY']]), 
  npcs = 30
)

# umap
ElbowPlot(all_data, ndims = 30)
all_data <- RunUMAP(all_data, dims = 1:20)


# harmony
all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:20,
    group.by.vars = "sample",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:20,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)




# Clustering
all_data <- FindNeighbors(all_data, dims = 1:20, reduction = "harmony")
all_data <- FindClusters(all_data, resolution = c(0.2, 0.4, 0.6, 0.8))



DimPlot(
  all_data,
  group.by = "BANKSY_snn_res.0.2"
)


cluster_lt_1000 <- all_data$BANKSY_snn_res.0.2 %>% table %>% {names(.)[. < 1000]}

all_data2 <- all_data[, !all_data$BANKSY_snn_res.0.2 %in% cluster_lt_1000]



DimPlot(
  all_data2,
  group.by = "BANKSY_snn_res.0.2"
)

DimPlot(
  all_data2,
  group.by = "sample"
)


saveRDS(all_data, 'all_data_banksy.rds')

### plot  ---------------
all_data_banksy <- readRDS('all_data_banksy.rds')
# export pdf
pdf(
  file.path(out_dir, "umap plot of banksy.pdf"),
  height = 6,
  width = 8
)
DimPlot(
  all_data_banksy,
  group.by = "BANKSY_snn_res.0.2"
)
DimPlot(
  all_data_banksy,
  group.by = "sample"
)
dev.off()


# spatial plot
annotation_col <- c("BANKSY_snn_res.0.2")

cells_position_ls <-
  all_data_banksy[[c("x", "y", "sample")]] %>%
  rownames_to_column('cell_id') %>%
  dplyr::rename(
    SPATIAL_1 = x,
    SPATIAL_2 = y
  ) %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

plot_dat <-
  FetchData(
    all_data_banksy, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')


p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 15
    width <- 15
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[annotation_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "spatial plot of celltype in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()






##5、cNMF-----

#0421analysis-------
xenium_RCTD <- readRDS("~/RCTD/xenium_RCTD.rds")
Epi_anno <- readRDS("RCTD/Epi_anno.rds")
rownames(xenium_RCTD@meta.data) %>% head()
rownames(Epi_anno@meta.data) %>% head()

Epi_anno$RCTD_Subtype <- as.character(Epi_anno$RCTD_Subtype)

xenium_RCTD@meta.data$RCTD_subtype <- as.character(xenium_RCTD@meta.data$first_type)
common_cells <- intersect(rownames(xenium_RCTD@meta.data), rownames(Epi_anno@meta.data))

xenium_RCTD@meta.data[common_cells, "RCTD_subtype"] <- Epi_anno@meta.data[common_cells, "RCTD_Subtype"]
table(xenium_RCTD$annotation,xenium_RCTD$RCTD_subtype)

ttt <- xenium_RCTD$RCTD_subtype %>% is.na()
table(ttt)

Epi_anno$RCTD_Subtype %>% table()

###反卷积亚类排序-------
RCTD_subtype_ordered <- c(
  #上皮
  "Basal1",
  "Basal2",
  "Club",
  "Hillock",
  "Luminal1",
  "Luminal2",
  "Luminal3",
  "Luminal4",
  "NE",
  
  # 内皮细胞 Endothelial
  "ArtECs", "CapECs", "EndoMT", "LECs", "TipECs", "VenECs",
  
  # 成纤维细胞 / CAF
  "Fibro_COL15A1", "Fibro_HSPA1B", "Fibro_ICAM1", "Fibro_IL6",
  "Fibro_PI16_MEOX2", "Fibro_PI16_MFAP4", "Fibro_POSTN",
  "Fibro_RPL36", "Fibro_like_ACTG2",
  
  # 基质细胞 Stromal
  "Pericyte_RGS5", "SMC_MYH11",
  
  # T 细胞
  "CD4_Tem_GZMK", "CD4_TfhTh1_CXCL13", "CD4_Th17_IL17A",
  "CD4_Tm_ANXA1", "CD4_Tm_CREM", "CD4_Tm_TXNIP", "CD4_Tn_SELL", "CD4_Treg_FOXP3",
  "CD8_ISG", "CD8_MAIT_SLC4A10", "CD8_MT1X", "CD8_Tem_GZMK",
  "CD8_Temra_CX3CR1", "CD8_Tex_CXCL13", "CD8_Tm_ANXA1", "CD8_Tm_IL7R",
  "CD8_Trm_GZMB", "CD8_Trm_ITGA1", "CD8_Trm_ZNF683", "CD8_Tstr_HSPA1A", "CD8_p-Tex_TOX",
  
  # NK / ILC
  "ILC3_KIT", "NK_CD16hi_CX3CR1", "NK_CD16hi_GZMB",
  "NK_CD16hi_RUNX3", "NK_CD16lo_CCL3", "NK_CD16lo_CREM",
  
  # B / Plasma
  "Bgc_RGS13", "Bn_TCL1A", "Plasma",
  
  # 髓系 Myeloid
  "Bm_CCR7", "Bm_CD74", "Bm_GPR183", "Bm_Stress_response",
  "cDC1_CLEC9A", "cDC2_CD1C", "cDC_LAMP3", "Mast",
  "Mo_CD14", "Mo_FCGR3A", "Mph_C1QC", "Mph_CX3CR1",
  "Mph_CXCL10", "Mph_FOLR2", "Mph_GPR183", "Mph_MT1X", "Mph_SPP1",
  "Neu_CSF3R", "pDC_LILRA4",
  
  # 其他
  "Hypoxia", "other"
)

xenium_RCTD$RCTD_subtype <- as.character(xenium_RCTD$RCTD_subtype)
xenium_RCTD$RCTD_subtype <- factor(xenium_RCTD$RCTD_subtype,levels = RCTD_subtype_ordered)
xenium_RCTD$RCTD_subtype %>% unique()


###合并大类------
xenium_RCTD$celltype <- case_when(
  # 1. 内皮
  xenium_RCTD$RCTD_subtype %in% c("Basal1",
                                "Basal2",
                                "Club",
                                "Hillock",
                                "Luminal1",
                                "Luminal2",
                                "Luminal3",
                                "Luminal4",
                                "NE") ~ "Epithelial",
  
  # 1. 内皮
  xenium_RCTD$RCTD_subtype %in% c("ArtECs","CapECs","EndoMT","LECs","TipECs","VenECs") ~ "Endothelial",
  
  # 2. 成纤维
  xenium_RCTD$RCTD_subtype %in% c("Fibro_COL15A1","Fibro_HSPA1B","Fibro_ICAM1","Fibro_IL6",
                                      "Fibro_PI16_MEOX2","Fibro_PI16_MFAP4","Fibro_POSTN",
                                      "Fibro_RPL36","Fibro_like_ACTG2") ~ "Fibroblast",
  
  # 3. 周细胞
  xenium_RCTD$RCTD_subtype %in% c("Pericyte_RGS5") ~ "Pericyte",
  
  # 4. 平滑肌
  xenium_RCTD$RCTD_subtype %in% c("SMC_MYH11") ~ "SMC",
  
  # 5. T 细胞（全部 CD4/CD8/NK/ILC 都归为 Tcell）
  xenium_RCTD$RCTD_subtype %in% c(
    "CD4_Tem_GZMK","CD4_TfhTh1_CXCL13","CD4_Th17_IL17A",
    "CD4_Tm_ANXA1","CD4_Tm_CREM","CD4_Tm_TXNIP","CD4_Tn_SELL","CD4_Treg_FOXP3",
    "CD8_ISG","CD8_MAIT_SLC4A10","CD8_MT1X","CD8_Tem_GZMK",
    "CD8_Temra_CX3CR1","CD8_Tex_CXCL13","CD8_Tm_ANXA1","CD8_Tm_IL7R",
    "CD8_Trm_GZMB","CD8_Trm_ITGA1","CD8_Trm_ZNF683","CD8_Tstr_HSPA1A","CD8_p-Tex_TOX",
    "ILC3_KIT","NK_CD16hi_CX3CR1","NK_CD16hi_GZMB",
    "NK_CD16hi_RUNX3","NK_CD16lo_CCL3","NK_CD16lo_CREM"
  ) ~ "Tcell",
  
  # 6. B 细胞 / 浆细胞
  xenium_RCTD$RCTD_subtype %in% c("Bgc_RGS13","Bn_TCL1A","Plasma",
                                      "Bm_CCR7","Bm_CD74","Bm_GPR183","Bm_Stress_response") ~ "Bcell",
  
  # 7. 髓系
  xenium_RCTD$RCTD_subtype %in% c(
    "cDC1_CLEC9A","cDC2_CD1C","cDC_LAMP3","Mast",
    "Mo_CD14","Mo_FCGR3A","Mph_C1QC","Mph_CX3CR1",
    "Mph_CXCL10","Mph_FOLR2","Mph_GPR183","Mph_MT1X","Mph_SPP1",
    "Neu_CSF3R","pDC_LILRA4"
  ) ~ "Myeloid",
  
  # 8. 上皮 + 其他未归类
  TRUE ~ "other"
)

xenium_RCTD$celltype %>% unique()

# 设置因子顺序（画图顺序固定）
xenium_RCTD$celltype <- factor(xenium_RCTD$celltype, 
                                     levels = c("Epithelial","Endothelial","Fibroblast",
                                                "Pericyte","SMC","Tcell","Bcell", "Myeloid"))

saveRDS(xenium_RCTD,"RCTD/result_addEPI/xenium_RCTD_addEPI.rds")


##1、邻近度分析-----------
cells_boundary_ls <- readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "RCTD/result_addEPI/proximity_out_dir"
dir.create(proximity_out_dir)

xenium_RCTD@meta.data %>% head()
xenium_RCTD@meta.data$cell_id %>% head()
rownames(xenium_RCTD@meta.data) %>% head()
#xenium_RCTD@meta.data$cell_id <- paste0("1_",xenium_RCTD@meta.data$cell_id1)

table(xenium_RCTD$group1,xenium_RCTD$RCTD_subtype)

cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1",
               "Mph_FOLR2",
               "Luminal1",
               "Luminal2",
               "Luminal3",
               "Luminal4"
)


###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
### Proximity analysis -------
#### 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- xenium_RCTD[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(RCTD_subtype, sample, "group1","cell_id") %>% 
      rownames_to_column("sample_cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(RCTD_subtype = factor(RCTD_subtype, levels = RCTD_subtype_ordered)) %>%
      dplyr::rename(Classification1 = RCTD_subtype,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group1
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path(proximity_out_dir,"all_data_meta_pos.csv"))
table(all_data_meta_pos$Slide_ID,all_data_meta_pos$Group)

#all_data_meta_pos <- read.csv(file.path(proximity_out_dir,"all_data_meta_pos_subtype7.csv"))

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


#### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞
#### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 5)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis.pdf"),
  height = 10,
  width  = 20
)
p_res_pa %>% print()
dev.off()

table(res_pa_df$Target,res_pa_df$Group)


all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928"
)

all_colors <- c(
  "Fibro_POSTN" = "#FF6B6B",
  "TipECs" = "#A6CEE3",
  "CD4_Treg_FOXP3" = "#00F5A7", 
  "CD4_TfhTh1_CXCL13"  = "#FDBB69",
  "CD8_Tex_CXCL13" = "#feee7d", 
  "Mph_SPP1" = "#00BBF9",
  "Mph_FOLR2" = "#2D82AF",
  "Luminal1" = "#B294C7",
  "Luminal2" = "#7D54A5",
  "Luminal3" = "#9E8099",
  "Luminal4" = "#600473"
)



#####以11个特定亚群为核心-------
cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1",
               "Mph_FOLR2",
               "Luminal1",
               "Luminal2",
               "Luminal3",
               "Luminal4"
)

res_pa_df$Group <- factor(res_pa_df$Group,levels = c("PERI","BPH","ADPC","CRPC"))

lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    dplyr::filter(Others %in% cell_list) %>% 
    #dplyr::filter(Others != "other") %>% 
    dplyr::filter(Target == cell) %>%
    
    ggplot(aes(x = Others, y = cell_props, fill = Group
    )) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis.pdf")),
    height = 8,
    width  = 8
  )
  p_res_pa %>% print()
  dev.off()
})


#### 2.2 proximity analysis - series analysis -----
init_cells <- c(
  '10'  = 750,
  '50'  = 800#,
  # '100' = 1000,
  # '150' = 1200,
  # '200' = 1500
)

res_pa_ls <-
  c(10, 50#, 100, 150, 200
    ) %>%
  set_names(., .) %>%
  map(function(x){
    k = init_cells[as.character(x)]
    res_pa <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos,
        query = all_data_meta_pos,
        x = "X",
        y = "Y",
        group.by = "Slide_ID",
        classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
        k = 'auto', #k, # 'auto',
        searchtype = "radius",
        radius = x, # 邻域半径
        mpp = 1, # !important, 1 for xenium
        returntype = "count",
        cell_enough_ratio_cutoff = 90
      )
    res_pa
  })

# res_pa_ls <- readRDS("data/res_pa_ls.rds")

# 将同心圆统计 转换成 环状同心的统计
# 注意一定要倒序操作
# 200 = 200 - 150
# 150 = 150 - 100 
# ...
# 
res_pa_ls2 <- res_pa_ls
diff_idx_ls <-
  c(10, 50#, 100, 150, 200
    ) %>% 
  as.character() %>%
  rev() %>%
  {matrix(c(.[-length(.)], .[-1]), ncol = 2)} %>%
  apply(1, c, simplify = FALSE)

for(i in diff_idx_ls){
  first  <- i[1]
  second <- i[2]
  res_pa_ls2[[ first ]] <- res_pa_ls2[[ first ]] - res_pa_ls2[[ second ]] 
}

# count to freq
res_pa_ls2 <- 
  res_pa_ls2 %>%
  map(function(x){# browser()
    x_colsum <- colSums(x)
    xx <- sweep(x, 2, x_colsum, FUN = "/")
    if(any(x_colsum == 0)){
      xx[, x_colsum == 0] <- 0
    }
    xx
  })

# 每个样本中的每种Target的其他细胞的计数取中位数
res_pa_count_tidy <- 
  res_pa_ls2 %>%
  imap(function(x, distance){# browser()
    # 分组
    pa_groups <-
      with(
        all_data_meta_pos,
        str_c(Classification1, Group, Slide_ID, sep = '..')
      )
    # 分组内取中位数
    res <- apply(x, 1, function(y){tapply(y, INDEX = pa_groups, FUN = "median")})
    
    # data tidy
    res_df  <-
      res %>%
      as.data.frame() %>%
      rownames_to_column('group') %>%
      pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
      separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
    res_df$distance <- distance
    
    # return
    res_df
  }) %>%
  bind_rows()



# plot 
# only plot "Epithelial cells" ....
res_pa_count_tidy_epi <-
  res_pa_count_tidy %>%
  mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
  dplyr::filter(Target == "Mph_SPP1") 

res_pa_count_tidy_epi_median <-
  res_pa_count_tidy_epi %>%
  group_by(Others, Group, distance) %>% 
  summarise(cell_props = median(cell_props))

p_proximity_series_epi <-
  res_pa_count_tidy_epi %>%
  ggplot(aes(x = distance, y = cell_props, color = Group)) +
  geom_boxplot( ) +
  ggsci::scale_color_aaas() +
  # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
  geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
  geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
  facet_wrap(~ Others, scales = 'free_y') +
  ggpubr::theme_pubclean()


# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis_series distance.pdf"),
  height = 14,
  width  = 18
)
p_proximity_series_epi %>% print()
dev.off()

#####循环11个亚群-------
cell_list <- c("Fibro_POSTN",
               "CD4_Treg_FOXP3",
               "CD4_TfhTh1_CXCL13",
               "CD8_Tex_CXCL13",
               "TipECs",
               "Mph_SPP1",
               "Mph_FOLR2",
               "Luminal1",
               "Luminal2",
               "Luminal3",
               "Luminal4"
)

res_pa_count_tidy$Group <- factor(res_pa_count_tidy$Group,levels = c("PERI","BPH","ADPC","CRPC"))

lapply(cell_list,function(cell){
  
  res_pa_count_tidy_epi <-
    res_pa_count_tidy %>%
    mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
    dplyr::filter(Target == cell) %>% 
    #dplyr::filter(Others != "other") %>% 
    dplyr::filter(Others %in% cell_list)
    
  
  res_pa_count_tidy_epi_median <-
    res_pa_count_tidy_epi %>%
    group_by(Others, Group, distance) %>% 
    summarise(cell_props = median(cell_props))
  
  p_proximity_series_epi <-
    res_pa_count_tidy_epi %>%
    ggplot(aes(x = distance, y = cell_props, color = Group)) +
    geom_boxplot( ) +
    ggsci::scale_color_aaas() +
    # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
    geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
    geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
    facet_wrap(~ Others, scales = 'free_y') +
    ggpubr::theme_pubclean()
  
  
  # export plot
  pdf(
    file.path(proximity_out_dir, paste0(cell,"_proximity analysis_series distance.pdf")),
    height = 14,
    width  = 18
  )
  p_proximity_series_epi %>% print()
  dev.off()
})



