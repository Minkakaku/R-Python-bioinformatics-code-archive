library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

##inferCNV----
#####set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


# 1.所有上皮细胞和down sample部分免疫细胞，跑inferCNV分析，
# 识别出正常肝细胞（肝母细胞/胆管上皮细胞？），恶性上皮细胞。
# 恶性上皮细胞跑cNMF分析，识别出上皮细胞程序。

library(infercnv)

#####epi subcluster --------------
Epi_out_dir <- "Epi_subcluster"
dir.create(Epi_out_dir)

Epi <- subset(all_data, annotation_full %in% c('Epithelial',"Luminal","NE"))
#Epi <- all_data %>% JoinLayers() %>% .[, .$annotation %in% "Epithelial"]


pdf(file.path(Epi_out_dir, 'nCount_RNA plot.pdf'), height = 6, width = 8)
VlnPlot(Epi, features = "nCount_RNA", group.by = 'sample', pt.size = 0)
dev.off()

#Epi[['RNA']] <- split(Epi[['RNA']], Epi$project, layers = "counts")


# cellcycle
Epi <- 
  Epi %>%
  CellCycleScoring(
    s.features = cc.genes.updated.2019$s.genes %>% str_to_title,
    g2m.features = cc.genes.updated.2019$g2m.genes %>% str_to_title
  )

# preprocess
Epi <- 
  Epi %>% 
  NormalizeData(scale.factor = floor(mean(all_data$nCount_RNA)))%>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(
    vars.to.regress = c(
      'nCount_RNA'
    )
  )



# sc_hep <- sc_hep %>% ScaleData(vars.to.regress = c('nCount_RNA', 'S.Score', 'G2M.Score'))

VariableFeaturePlot(Epi)

# VariableFeatures(sc_hep) <- VariableFeatures(sc_hep) %>% setdiff(c(ig_genes, tcr_genes, mt_genes, ribo_genes))



#####pca -----
Epi <- Epi %>% RunPCA()

pdf(file.path(Epi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epi, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_Epi <- 
  Epi[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(Epi) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(Epi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 1:5)
dev.off()
pdf(file.path(Epi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 6:10)
dev.off()
pdf(file.path(Epi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 11:15)
dev.off()
pdf(file.path(Epi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 16:20)
dev.off()


#####dimension reduction --------
epi_dims <-  15


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Epi <- RunUMAP(Epi, dims = 1:epi_dims)


pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{epi_dims} project.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "project", label = TRUE) + NoLegend()
dev.off()


Epi <- Epi %>% 
  harmony::RunHarmony(
    dims.use = 1:epi_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Epi <- RunUMAP(
  Epi, 
  dims = 1:epi_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{epi_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# plot QC figure
pdf(file.path(Epi_out_dir, str_glue("qc - Dimplot dim{epi_dims}.pdf")), height = 6, width = 8)
FeaturePlot(Epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI')) + NoLegend()
FeaturePlot(Epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("qc - Vlnplot dim{epi_dims}.pdf")), height = 6, width = 10)
VlnPlot(Epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), group.by = 'project', pt.size = 0) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("FeaturePlot - dim{epi_dims} - CellCycle.pdf")), height = 6, width = 12)
FeaturePlot(Epi, features = c('S.Score', 'G2M.Score'), reduction = "harmonyumap")
dev.off()

pdf(file.path(Epi_out_dir, str_glue("DimPlot - dim{epi_dims} - Phase.pdf")), height = 5, width = 6)
DimPlot(Epi, group.by = c('Phase'), reduction = "harmonyumap")
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

#####cluster ------
Epi <-
  FindNeighbors(
    Epi,
    reduction = reduction.use,
    dims = 1:epi_dims,
    k.param = 10
  )

Epi <-
  FindClusters(
    Epi,
    resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 5, width = 6)
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Epi, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

#####Dotplot -------
# epi cluster
marker_epi <- list(
  "AT1"     = c("AGER","EMP2","CAV1","RTKN2","MYL9","HOPX","PDPN","MSLN","CLIC5"),
  # "AT1-AT2" = c("EGFR","SHH","MUC1","PDPN"),
  "AT2"     = c("MUC1","SFTPC","NPC2","SFTPA1","NAPSA","SFTPA2","CTSH","PGC","SFTPB","SFTPD","ETV5","ABCA3"),
  basal     = c("KRT5","KRT17","MIR205HG","DLK2","BCL11A","KRT14","TP63","DAPL1"),
  'ciliated cells' = c("CCDC113","FOXJ1","DNAH3","DNAIC2","DNAAF1","LRRC6"),
  'club cells'     = c("SCGB1A1","SCGB3A2","CYP2F2","CCKAR")
)

# Hep cluster
marker_Hep <- list(
  "PC_hep"    = c("Oat","Glul","Cyp2c37"),
  "MZ_hep"    = "Hamp",
  "PP_hep"    = c("Hal","Sds","Gls2"),
  "Cyc_hep"   = c("Mki67","Top2a","Cenpa")
)


table(Epi$Xenium_snn_res.0.4)
pdf(
  file.path(Epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - major marker - RNA_snn_res.0.1.pdf")), 
  height = 5,
  width = 20
)
# major marker
DotPlot(
  Epi,
  features = c(
    "Immune cell" = "PTPRC",  # 小鼠Ptprc（CD45）→ 人类PTPRC
    "T cell" = "IL2RB",       # 小鼠Il2rb → 人类IL2RB
    "T cell" = "MS4A4A",      # 小鼠Ms4a4b（MS4A家族）→ 人类同源物MS4A4A
    "T cell" = "IL7R",        # 小鼠Il7r → 人类IL7R
    "T cell" = "LCK",         # 小鼠Lck → 人类LCK
    "T cell" = "CD3D",        # 小鼠Cd3d → 人类CD3D
    "T cell" = "CD3E",        # 小鼠Cd3e → 人类CD3E
    "B cell" = "CD19",        # 小鼠Cd19 → 人类CD19
    "B cell" = "CD79A",       # 小鼠Cd79a → 人类CD79A
    "B cell" = "MS4A1",       # 小鼠Ms4a1（CD20）→ 人类MS4A1
    "PC_hep" = "OAT",         # 小鼠Oat → 人类OAT
    "PC_hep" = "GLUL",        # 小鼠Glul → 人类GLUL
    "PC_hep" = "CYP2C9",      # 小鼠Cyp2c37（CYP2C家族）→ 人类同源物CYP2C9
    "MZ_hep" = "HAMP",        # 小鼠Hamp → 人类HAMP（ Hepcidin）
    "PP_hep" = "HAL",         # 小鼠Hal → 人类HAL
    "PP_hep" = "SDS",         # 小鼠Sds → 人类SDS
    "PP_hep" = "GLS2",        # 小鼠Gls2 → 人类GLS2
    "Cyc_hep" = "MKI67",      # 小鼠Mki67 → 人类MKI67
    "Cyc_hep" = "TOP2A",      # 小鼠Top2a → 人类TOP2A
    "Cyc_hep" = "CENPA",      # 小鼠Cenpa → 人类CENPA
    "HSCs" = "RELN",          # 小鼠Reln → 人类RELN
    "HSCs" = "DCN",           # 小鼠Dcn → 人类DCN
    "HSCs" = "NR1H5",         # 小鼠Nr1h5 → 人类NR1H5
    "Portal_fibroblasts" = "AEBP1",  # 小鼠Aebp1 → 人类AEBP1
    "Portal_fibroblasts" = "MGP",    # 小鼠Mgp → 人类MGP
    "Portal_fibroblasts" = "IGFBP5", # 小鼠Igfbp5 → 人类IGFBP5
    "Fibroblasts" = "PDGFRA",        # 小鼠Pdgfra → 人类PDGFRA
    "Fibroblasts" = "COL1A1",        # 小鼠Col1a1 → 人类COL1A1
    "Fibroblasts" = "COL1A2",        # 小鼠Col1a2 → 人类COL1A2
    "LSECs" = "KDR",          # 小鼠Kdr（VEGFR2）→ 人类KDR
    "LSECs" = "FLT1",         # 小鼠Flt1（VEGFR1）→ 人类FLT1
    "LVECs" = "VWF",          # 小鼠Vwf → 人类VWF
    "LVECs" = "EDN1",         # 小鼠Edn1 → 人类EDN1
    "LVECs" = "JAM2",         # 小鼠Jam2 → 人类JAM2
    "Cholangiocytes" = "SPP1",       # 小鼠Spp1 → 人类SPP1（骨桥蛋白）
    "Cholangiocytes" = "EPCAM",      # 小鼠Epcam → 人类EPCAM
    "Cholangiocytes" = "GLIS3",      # 小鼠Glis3 → 人类GLIS3
    "Kupffer" = "CLEC4G",     # 小鼠Clec4f（肝窦巨噬细胞标记）→ 人类同源物CLEC4G
    "Kupffer" = "CD5L",       # 小鼠Cd5l → 人类CD5L
    "Kupffer" = "TRPM2",      # 小鼠Trpm2 → 人类TRPM2
    "Dendritic" = "SIGLEC10", # 小鼠Siglech（sialic acid受体）→ 人类同源物SIGLEC10
    "Dendritic" = "IRF8",     # 小鼠Irf8 → 人类IRF8
    "Dendritic" = "RUNX2",    # 小鼠Runx2 → 人类RUNX2
    "Monocytes/non-Kupffer_cell_macrophages" = "CCR2",  # 小鼠Ccr2 → 人类CCR2
    "Monocytes/non-Kupffer_cell_macrophages" = "MYO1F", # 小鼠Myo1f → 人类MYO1F
    "Monocytes/non-Kupffer_cell_macrophages" = "GPNMB", # 小鼠Gpnmb → 人类GPNMB
    "Pericytes" = "RGS5",     # 小鼠Rgs5 → 人类RGS5
    "Pericytes" = "KCNJ8",    # 小鼠Kcnj8 → 人类KCNJ8
    "Pericytes" = "HIGD1B",   # 小鼠Higd1b → 人类HIGD1B
    "Pericytes" = "ABCC9",    # 小鼠Abcc9 → 人类ABCC9
    "MDSCs" = "CD33",         # 小鼠Cd33 → 人类CD33
    "MDSCs" = "ITGAM"         # 小鼠Itgam（CD11b）→ 人类ITGAM
  ),
  group.by = "Xenium_snn_res.0.4"
) + 
  RotatedAxis()
dev.off()

pdf(
  file.path(Epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - marker_epi - RNA_snn_res.0.1.pdf")), 
  height = 5,
  width = 10
)
# liver hep marker
DotPlot(
  Epi,
  features = marker_epi,
  group.by = "Xenium_snn_res.0.4"
) +
  RotatedAxis()
dev.off()

#####dotplot ------
major_marker_ls<-c(
  "epithelial state" = "CDH1",
  "epithelial state" = "EPCAM",
  "epithelial state" = "CLDN3",
  "epithelial state" = "OCLN",
  "mesenchymal state" = "VIM",
  "mesenchymal state" = "CDH2",
  "mesenchymal state" = "FN1",
  "EMT" = "SNAI1",
  "EMT" = "SNAI2",
  "EMT" = "TWIST1",
  "EMT" = "TWIST2",
  "EMT" = "ZEB1",
  "EMT" = "ZEB2",
  "EMT" = "TGFB1",
  "EMT" = "MMPS"
)

celltype_marker<-c(
  "B cell"= "CD19", 
  "B cell"= "CD79A", 
  "B cell"= "CD79B", 
  "B cell"= "MS4A1",
  "Plasma" = "TNFRSF17", 
  "Plasma" = "MZB1",
  # "Plasma" = "XBP1",
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "CD4",
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  "Endothelial" = "PECAM1",
  "Endothelial" = "CLDN5", 
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT18",
  "Epithelial" = "CDH1",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Monocytes" = "ITGAM",
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  "Monocytes" = "CD68",
  "Macrophage" = "CSF1R",
  "Macrophage" = "MARCO", 
  "DC"= "FCER1A", 
  "DC"= "CD1C", 
  "DC"= "CD83", 
  "DC"= "LILRA4", 
  "DC"= "CLEC9A", 
  "DC"= "IDO1",
  # 'pDC' = 'CLEC4C',
  # 'pDC' = "LILRA4",
  # 'pDC' = 'IL3RA',
  # 'pDC' = 'IRF8',
  # 'mregDC' = 'CCR7',
  # 'mregDC' ='CCL19',
  # 'mregDC' ='LAMP3',
  # 'mregDC' ='CD83',
  "NK" = "NCAM1",
  "Mast cell" = "CTSG", 
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = 'HDC', 
  'Neutrophil'= "CSF3R", 
  'Neutrophil'= "FCGR3B", # "S100A8", "S100A9" 
  "Melanocy" = "MITF",
  "Melanocy" = "MLANA"
)

major_marker_ls1<-c(
  "basal marker" = "ACTA2",
  "basal marker" = "DST",
  "basal marker" = "SAA1",
  "basal marker" = "NTN4",
  "basal marker" = "XBP1,",
  "basal marker" = "PROM2",
  
  "LumHR marker" = "ANKRD30A",
  "LumHR marker" = "AGR2",
  "LumHR marker" = "AREG",
  
  "LumSec marker" = "KRT15",
  "LumSec marker" = "LTF",
  "LumSec marker" = "MMP7",
  
  "luminal"= "FOLH1", 
  "luminal"= "FOXA1",  
  "luminal"= "AMACR",
  
  "tumor" = "NKX3-1", 
  "tumor" = "TMPRSS2", 
  "tumor" = "ACPP",
  "tumor" = "EEF1G", 
  "tumor" = "HSPD1,", 
  "tumor" = "TOP2A", 
  "tumor" = "MCM4,", 
  "tumor" = "BIRC5",
  "tumor" = "PDIA3", 
  "tumor" = "PKM", 
  
  "emt" = "GDF15",
  "emt" = "SPOCK1",
  "emt" = "EPCAM"
)


major_marker_ls2 <- c(
  # 基底样亚型标记 (Basal-like)
  "basal marker" = "KRT5",
  "basal marker" = "KRT14",
  "basal marker" = "KRT17",
  "basal marker" = "EGFR",
  "basal marker" = "VIM",
  
  # 腔面雄激素受体亚型标记 (Luminal Androgen Receptor)
  "LumHR marker" = "AR",
  "LumHR marker" = "FOXA1",
  "LumHR marker" = "PIP",
  
  # 腔面分泌型标记 (Luminal Secretory)
  "LumSec marker" = "KRT18",
  "LumSec marker" = "LTF",
  "LumSec marker" = "MMP7",
  
  # 管腔上皮通用标记
  "luminal" = "EPCAM", 
  "luminal" = "CDH1",  
  "luminal" = "KRT19",
  
  # 肿瘤恶性特征标记
  "tumor" = "MKI67", 
  "tumor" = "TP53", 
  "tumor" = "CCND1",
  "tumor" = "TOP2A", 
  "tumor" = "BIRC5",
  "tumor" = "PKM", 
  "tumor" = "PHGDH",
  
  # 上皮间质转化标记 (EMT)
  "emt" = "ZEB1",
  "emt" = "SNAI1",
  "emt" = "FN1"
)


major_marker_ls3 <- c(
  "ESR1",
  "PGR",
  "ERBB2",
  
  #基底样亚型标记 (Basal-like, 包含 BL-1 和 BL-2)
  "basal marker" = "TP63",
  "basal marker" = "BRCA1", # BL-1: DNA 修复缺陷相关
  "basal marker" = "EGFR", # BL-2: 生长因子通路激活
  "basal marker" = "KRT5", # BL-2: CK5 编码基因，基底上皮特征
  "basal marker" = "MET", # BL-2: 生长因子受体，促进增殖
  #腔面雄激素受体亚型标记 (Luminal Androgen Receptor, LAR)
  "LumAR marker" = "AR", # LAR: 雄激素受体，信号依赖核心
  "LumAR marker" = "FOXA1", # LAR: 雄激素信号共调节因子
  "LumAR marker" = "PIP", # LAR: 管腔上皮分泌蛋白，雄激素响应标记
  #间质型标记 (Mesenchymal, M)
  "mesenchymal marker" = "VIM", # M: 波形蛋白，EMT 特征
  "mesenchymal marker" = "CD44", # M: 细胞表面黏附分子，EMT 相关
  "mesenchymal marker" = "FN1", # M: 纤连蛋白，上皮间质转化关键分子
  #间质干细胞型标记 (Mesenchymal Stem-like, MSL)
  "MSL marker" = "PROM1", # MSL: CD133 编码基因，干细胞标记
  "MSL marker" = "ALDH1A1", # MSL: ALDH1 主要功能亚型，干细胞特性与耐药相关
  #肿瘤恶性特征标记 (涵盖高增殖、DNA 异常等)
  "tumor" = "TOP2A", # 重复标注以强调 BL-1 高增殖特征
  "tumor" = "MKI67" # 重复标注以强调 BL-1 增殖指数
)

# find markers
Idents(Epi) <- "RNA_snn_res.0.4"
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epi),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

library(purrr)
Epithelial_top5_markers <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Epithelial_top5_markers) <- sub('.*\\.','',names(Epithelial_top5_markers))


#提取deepseek
Epi_top10_markers <-
  Epithelial_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)



library(stringr)
# dotplot
pdf(file.path(Epi_out_dir, str_glue("dotplot dims{epi_dims} _0.8_major_marker.pdf")), width = 25, height = 6)

DotPlot(Epi, group.by = "Xenium_snn_res.0.8", features = celltype_marker)+RotatedAxis()
DotPlot(Epi, group.by = "Xenium_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(Epi, group.by = "Xenium_snn_res.0.8", features = major_marker_ls1)+RotatedAxis()
DotPlot(Epi, group.by = "Xenium_snn_res.0.8", features = major_marker_ls2)+RotatedAxis()
DotPlot(Epi, group.by = "Xenium_snn_res.0.8", features = major_marker_ls3)+RotatedAxis()
DotPlot(Epi, group.by = "Xenium_snn_res.0.8", features = Epithelial_top5_markers)+RotatedAxis()

DotPlot(xenium, group.by = "res_0.8", features = major_marker_ls1)+RotatedAxis()
DotPlot(xenium, group.by = "res_0.8", features = major_marker_ls2)+RotatedAxis()
DotPlot(xenium, group.by = "res_0.8", features = major_marker_ls3)+RotatedAxis()

dev.off()

#####feature-----
valid_genes <- major_marker_ls2[major_marker_ls2 %in% rownames(xenium)]
cat("有效基因:", length(valid_genes), "/", length(major_marker_ls1))

for (gene in valid_genes) {
  p <- FeaturePlot(Epi, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(result_dir, paste0("Epi_subcluster/featureplot/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


###改cluster名
Epi@meta.data$subcluster1 <- paste0("EpiC",Epi@meta.data$Xenium_snn_res.0.8)

# 保存中间文件
saveRDS(Epi, file.path(Epi_out_dir, "Epi.rds"))

Epi <- readRDS(file.path(Epi_out_dir, "Epi.rds"))

####映射——add cluster类到xenium----
library(tibble)
library(data.table)
sub <- Epi[['subcluster1']]%>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    xenium[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtypeEpi1 = ifelse((as.character(annotation) %in% c("Epithelial")),
                         as.character(subcluster1), 
                         as.character(annotation)))%>%
  dplyr::select(cells, subtypeEpi1) %>%
  column_to_rownames('cells')

xenium<-AddMetaData(xenium,meta.data.t['subtypeEpi1'])
xenium$annotation %>% table()
xenium$subtypeEpi1 %>% table()
# xenium$annotation2 <- NULL
# 确认是否加入
table(xenium$subtypeEpi1, xenium$sample, useNA = "ifany")

new_labels <- ifelse(
  xenium$subtypeEpi1 %in% paste0("EpiC", 0:13),  # 保留EpiC0-EpiC13
  as.character(xenium$subtypeEpi1),               # 保留原名
  ifelse(
    xenium$subtypeEpi1 %in% c("B", "Endothelial", "Fibroblast", "Myeloid", "Plasma", "SMC", "T"),
    "other",                                     # 指定类型标记为other
    as.character(xenium$subtypeEpi1)              # 其他类型保留（可选）
  )
)

xenium$subtypeEpi1 <- factor(new_labels)

#####Imagedimplot----
# library(RColorBrewer)
# epi_colors <- brewer.pal(12, "Set3")  # 生成12色，补充2色
# epi_colors <- c(epi_colors, "#8DD3C7", "#FFFFB3")  # 补足14色
# names(epi_colors) <- paste0("EpiC", 0:13)  # 命名
# # 设置other为灰色
# all_colors <- c(epi_colors, "other" = "grey")

all_colors <- c(
  "red", "#A6CEE3", "#5B9EC9", "#FF6B3A", "#2D82AF", #0-4
  "#98D277","#52AF43", "#6F9E4C", "#FDBB69", "#F78620",#5-9
  "#00F5A7", "#F15BB5", "#7D54A5", "#02c9c9", "#600473",#10-14
  "#B294C7", "pink", "#9E8099", "blue", "#DBB466",#15-19
  "#B15928", "#4ECDC4", "#F0EB99"          #20-22,other
)
names(all_colors) <- paste0("EpiC", 0:22)
all_colors <- c(all_colors, "other" = "grey90")

per <- c(
  "#7209B7", "#4361EE",  "#00BBF9", "#00F5D4","#7EBA98","#feee7d", "#FE982C","#FF9F1C",
  "#9B5DE5", "#FEE440", "#9C89B8", "#06D6A0" )

sample_label <- unique(xenium@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors,
                    group.by = "subtypeEpi1"
  )+
    ggtitle(paste(i,"_Epithelial_cluster")) +
    theme(plot.title = element_text(
      size = 20, 
      face = "bold", #加粗
      hjust = 0.5, #居中
      vjust = -3 #下移
    ))+
    coord_flip()+ # 交换X/Y轴
    scale_x_reverse()  # 反转x轴方向（匹配顺时针效果）
  
  ggsave(p,filename = paste0("result/Epi_subcluster/Imagedimplot_res0.8/",i,"_Epithelial_cluster.png"),width = 12 ,height = 12,dpi = 300)
}


#####加入上皮再分群cluster到all_data----
Epi_anno_meta <-
  Epi[['RNA_snn_res.0.4']] %>% 
  rownames_to_column("cell_id") 

#B cells  pDCs  T cells Endothelial  Macrophage Fibroblasts        
#Plasma  MergeDC Epithelial Smooth muscle cells

all_data %>% head()
Epi %>% head()

all_data_anno_meta <- 
  #xenium[[]] %>% dplyr::filter(nCount_RNA > 100) %>% .['annotation'] %>%
  all_data[['annotation_full']] %>%
  rownames_to_column('cell_id') %>%
  # dplyr::filter(annotation %in% c("Club cells","TRB-SC","Basal","AT2","Plasma","Macrophage","T cells")) %>%
  dplyr::filter(annotation_full %in% c('Epithelial',"Luminal","NE","T cells","Macrophage")) %>%
  left_join(Epi_anno_meta, by = "cell_id") %>% 
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(
    annotation_full %in% c('Epithelial',"Luminal","NE"),
    RNA_snn_res.0.4,
    annotation_full
  )) %>% 
  dplyr::select(cell_id, annotation = annotation2)


# all_data <- all_data %>% AddMetaData(all_data_anno_meta)  
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) #%>%
#slice_sample(prop = 0.5)

###inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹

infercnv_out_dir = "Epi_subcluster/infercnv_hepsub_variable_all500_res0.4"
dir.create(infercnv_out_dir) 


# meta_data_infer <-
#   all_data[[]] %>%
#   rownames_to_column("cell_id") %>%
#   dplyr::filter(
#     annotation %in% c("Club cells","TRB-SC","Basal","AT2",'Endothelial',"Plasma","Macrophage","T cells")
#   ) %>%
#   dplyr::select(cell_id, annotation, Spatial.008um_snn_res.0.8) %>%
#   map_dfc(as.character) 

# unique(meta_data_infer$annotation)
# "Acinar cells"     "Ductal cells"     "Endocrine cells"  "Epithelial cells" "NKT cells"


# cell is too much
# less than 50000
# 取1/20细胞： 1179497 --> 47179.88
# 1662182 --> 83000
# meta_data_infer_subset <-
#   meta_data_infer %>%
#   group_by(annotation) %>%
#   slice_sample(prop = 0.5)

# 合并annotation
# meta_data_infer_subset <-
#   meta_data_infer_subset %>%
#   mutate(Spatial.008um_snn_res.0.8 = str_c('Cluster', Spatial.008um_snn_res.0.8, sep = "_")) %>%
#   mutate(annotation2 = ifelse(
#     annotation %in% c("Plasma", "Macrophage","T cells","Endothelial"),
#     annotation,
#     Spatial.008um_snn_res.0.8
#   )) %>%
#   dplyr::select(
#     cell_id,
#     annotation = annotation2
#   )

# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )


xenium <- FindVariableFeatures(xenium)

#undebug(CreateInfercnvObject)
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(all_data)[1:500], meta_data_infer_subset$cell_id] ,
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path(infercnv_out_dir, "gene_pos47.txt"),
  
  min_max_counts_per_cell = c(30, Inf),
  
  delim             = "\t",
  ref_group_names   = c("T cells","Macrophage")
)




# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.001) } # 304 <--- 
infercnv_obj@count.data %>% rowMeans %>% { sum(. < 0.005) } # 218


infercnv_obj2  = infercnv::run(
  infercnv_obj,
  cutoff              =    0.001,
  out_dir    =    infercnv_out_dir,
  cluster_by_groups   =    TRUE,
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  png_res             =    300,
  # useRaster   = FALSE,
  output_format = 'pdf',
  no_plot = T #  figure margin too large ---
)
infercnv::plot_cnv(
  infercnv_obj2,
  out_dir           =  infercnv_out_dir, 
  output_filename   =  'infercnv',
  output_format     =  "pdf",
  png_res = 300, 
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)

saveRDS(infercnv_obj2, file.path(infercnv_out_dir, 'infercnv_obj2.rds'))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, 'infercnv_obj.rds'))

# # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
# infercnv_obj_default = infercnv::run(
#   infercnv_obj,
#   cutoff              =    0.01, 
#   out_dir    =    infercnv_out_dir,
#   cluster_by_groups   =    TRUE, 
#   plot_steps          =    FALSE,
#   denoise             =    TRUE,
#   HMM                 =    FALSE,
#   no_prelim_plot      =    TRUE,
#   png_res             =    300,
#   # useRaster   = FALSE,
#   output_format = 'pdf',
#   no_plot = T #  figure margin too large ---
# )
# 
# # C stack error
# infercnv::plot_cnv(
#   infercnv_obj_default,
#   out_dir           =  infercnv_out_dir,
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 60,
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )


# download inferCNV object from linux server
# plot inferCNV locally
# saveRDS(infercnv_obj_default, file.path(infercnv_out_dir, 'infercnv_obj_default.rds'))
# 
# infercnv_obj_default <- readRDS("inferCNV/infercnv_obj_default.rds")
# 
# infercnv::plot_cnv(
#   infercnv_obj_default,
#   out_dir           =  'inferCNV', 
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 300, 
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )
# 
infercnv_obj2_default_median_filtering <-
  infercnv_obj2 %>%
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           =  infercnv_out_dir,
    output_filename   =  'infercnv_median_filtering',
    # output_format     =  "pdf",
    output_format="pdf",
    cluster_by_groups =  TRUE,
    x.range           =  "auto",
    x.center          =  1,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering .rds'))

# expFile='expFile.txt'
# groupFiles='groupinfo.txt'
# geneFile='gene_pos.txt'
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
#                                     annotations_file=groupFiles,
#                                     delim="\t",
#                                     gene_order_file= geneFile,
#                                     ref_group_names=c('B_cells','plasma_cells','T_cells','microglial')) 
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics



# infercnv::plot_cnv(
#   infercnv_obj2,
#   out_dir           =  'infercnv_output_2',
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 300,
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.01, 
#                               out_dir= "infercnv_output_2",  
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               hclust_method="ward.D2", plot_steps=F,no_plot = T)
# infercnv::plot_cnv(infercnv_obj2, 
#                    out_dir='infercnv_output_2',
#                    output_filename='infercnv', 
#                    x.range="auto",
#                    x.center=1,
#                    title = "infercnv", 
#                    output_format="pdf",
#                    color_safe_pal = FALSE)
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.1, 
#                               out_dir= "infercnv_output_sd",  
#                               write_expr_matrix=T,
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               sd_amplifier=3,
#                               denoise=T,
#                               noise_logistic=TRUE,
#                               hclust_method="ward.D2", plot_steps=F)
# saveRDS(infercnv_obj2, file.path(infercnv_output, 'infercnv_obj2.rds'))
# saveRDS(infercnv_obj, file.path(infercnv_output, 'infercnv_obj.rds'))


#####cnv score ------------
library(scales)
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}
# cnvScore <- function(data){
#   sd_num <- apply(data, 1, sd)
#   data <- data[sd_num > median(sd_num), ]
#   
#   data <- data %>% as.matrix() %>%
#     t() %>% 
#     scale() %>% 
#     t()
#   
#   cnv_score <- as.data.frame(colSums(data * data))
#   return(cnv_score)
# }
cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <-'CNV_score'
temp <-
  cnv_score %>% 
  rownames_to_column('cell_id') %>% 
  left_join(.,as.data.frame(meta_data_infer_subset),'cell_id')  

p_infer_cnv_score <-
  temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+geom_violin(aes(fill= annotation),color="NA",drop = FALSE)+
  ggsci::scale_fill_igv()+
  # scale_fill_manual(values = color_v)+
  theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height = 6,
  width  = 12
)

cols_brewer = c("#f94144","#f25c27","#eb770a","#f89f12","#f9c74f","#90be6d",
                "#43aa8b","#4d908e","#577590","#277da1","#967aa1",
                "#F2698D","#F391AE","#F3B9CF","#AF8EB8","#7B84B7","#3F6FAF",
                "#035AA6","#0473B4","#8491B7","#737AAB","#5E5E97","gray90")

p_infer_cnv_score <-
  temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+
  geom_boxplot(aes(fill= annotation),color="NA",outlier.colour = "black")+
  # ggsci::scale_fill_igv()+
  scale_fill_manual(values = cols_brewer)+
  theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_boxplot.pdf"),
  height = 6,
  width  = 12
)



#####spatial plot --------------
# Idents(all_data) <- 'Spatial.008um'
# Idents(all_data) <- 'Spatial.qc.pass'

sc_hep_infer_anno <-
  sc_hep[['Spatial.008um_snn_res.0.4']] %>% 
  rownames_to_column("cell_id") %>%
  mutate(annotation_epi = ifelse(
    Spatial.008um_snn_res.0.4 %in% c('3', '6', '12' ,'13'),
    "Tumor cells",
    "Epithelial cells"
    # Spatial.008um_snn_res.0.4
  )) 

all_data_infer_anno <-
  all_data[['annotation']] %>%
  rownames_to_column("cell_id") %>%
  left_join(sc_hep_infer_anno, by = "cell_id") %>%
  mutate(annotation2 = ifelse(
    is.na(annotation_epi),
    as.character(annotation),
    as.character(annotation_epi)
  )) %>%
  column_to_rownames("cell_id")


sp_p_ls <-
  Images(all_data) %>%
  set_names(., .) %>%
  map(function(img){
    all_data %>%
      AddMetaData(all_data_infer_anno['annotation2']) %>%
      SpatialDimPlot(
        images = img, 
        group.by = "annotation2",
        pt.size.factor	= 2,
        image.scale = "lowres",
        crop = FALSE
      ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img)
  })


# export
pdf(
  file.path(hep_out_dir, "spatial plot of infercnv analysis.pdf"),
  height = 8,
  width  = 8
)
sp_p_ls %>% walk(print)
dev.off()



####annotation -----------
cell_cluster <- c(
  '0' = "Tumor",
  '1' = "Tumor",
  '2' = "Tumor",
  '3' = "Tumor",
  '4' = "Tumor",
  '5' = "Tumor",
  '6' = "Tumor",
  '7' = "Tumor",
  '8' = "Tumor",
  '9' = "Tumor"
)

Idents(Epi) <- "RNA_snn_res.0.4"

Epi <- RenameIdents(Epi, cell_cluster)
Epi$subtype <- Idents(Epi)

dims <- 15
pdf(file.path(Epi_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epi, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

####筛选掉group3中"verification"组的细胞----
Epi_subset <- subset(Epi, subset = group1 != "verification")  # 排除verification组
# 输出PDF并绘图（仅包含剩余的group3分组）
pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation split3_no_verification.pdf")), 
    height = 7, width = 14)
DimPlot(
  Epi_subset,  # 使用筛选后的对象
  group.by = "subtype", 
  split.by = "group3",  # 此时split仅包含非verification的组
  label = TRUE, 
  reduction = "harmonyumap", 
  repel = TRUE
) + NoLegend()
dev.off()

####cell props barplot  --------
Epithelial_meta <- Epi_subset[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample, 
                                 levels=c(Epithelial_meta$sample %>% unique())
)
p_col_proj <-
  Epithelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/Epi_subcluster/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Epithelial_meta %>%
  group_by(group3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path("result/Epi_subcluster/cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)


##映射—add subtype到xenium----
library(tibble)
library(data.table)
sub <- Epi[['subtype']]%>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    xenium[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtypeEpi = ifelse((as.character(annotation) %in% c("Epithelial")),
                        as.character(subtype), 
                        as.character(annotation)))%>%
  dplyr::select(cells, subtypeEpi) %>%
  column_to_rownames('cells')

xenium<-AddMetaData(xenium,meta.data.t['subtypeEpi'])
xenium$annotation %>% table()
xenium$subtypeEpi %>% table()
# xenium$annotation2 <- NULL
# 确认是否加入
table(xenium$subtypeEpi, xenium$annotation, useNA = "ifany")

# new_labels <- ifelse(
#   xenium$subtypeEpi1 %in% paste0("EpiC", 0:13),  # 保留EpiC0-EpiC13
#   as.character(xenium$subtypeEpi1),               # 保留原名
#   ifelse(
#     xenium$subtypeEpi1 %in% c("B", "Endothelial", "Fibroblast", "Myeloid", "Plasma", "SMC", "T"),
#     "other",                                     # 指定类型标记为other
#     as.character(xenium$subtypeEpi1)              # 其他类型保留（可选）
#   )
# )
# 
# xenium$subtypeEpi1 <- factor(new_labels)

xenium$subtypeEpi <- factor(xenium$subtypeEpi)

####Imagedimplot----
# library(RColorBrewer)
# epi_colors <- brewer.pal(12, "Set3")  # 生成12色，补充2色
# epi_colors <- c(epi_colors, "#8DD3C7", "#FFFFB3")  # 补足14色
# names(epi_colors) <- paste0("EpiC", 0:13)  # 命名
# # 设置other为灰色
# all_colors <- c(epi_colors, "other" = "grey")

all_colors <- c(
  #"red", 
  "#A6CEE3", "#5B9EC9", "#FF6B3A", "#2D82AF", #0-4
  "#98D277","#52AF43", "#6F9E4C", "#FDBB69", "#F78620",#5-9
  "#00F5A7", "#F15BB5", "#7D54A5", "#02c9c9", "#600473",#10-14
  "#B294C7", "pink", "#9E8099", 
  #"blue", 
  "#DBB466",#15-19
  "#B15928", "#4ECDC4", "#F0EB99"          #20-22,other
)
# names(all_colors) <- paste0("EpiC", 0:22)
# all_colors <- c(all_colors, "other" = "grey90")
# 
per <- c(
  "#7209B7", "#4361EE",  "#00BBF9", "#00F5D4","#7EBA98","#feee7d", "#FE982C","#FF9F1C",
  "#9B5DE5", "#FEE440", "#9C89B8", "#06D6A0" )

all_colors10 <- c(
  "pink","#7EBA98","#A6CEE3","#00F5F2", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7","#FF6B6B"
)
names(all_colors10) <- c("B","Endothelial","Fibroblast","Myeloid","Plasma","SMC","T", "Basal", "LumHR", "Tumor")


sample_label <- unique(xenium@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors10,
                    group.by = "subtypeEpi"
  )+
    ggtitle(paste(i,"_annotation")) +
    theme(plot.title = element_text(
      size = 20, 
      face = "bold", #加粗
      hjust = 0.5, #居中
      vjust = -3 #下移
    ))#+
  #coord_flip()+ # 交换X/Y轴
  #scale_x_reverse()  # 反转x轴方向（匹配顺时针效果）
  
  ggsave(p,filename = paste0("result/Epi_subcluster/Imagedimplot/",i,"_Epithelial_annotation.png"),width = 12 ,height = 12,dpi = 300)
}


####癌细胞比例变化####----
#样本比例折线图
S <- xenium
Idents(S) <- S$subtypeEpi
table(S$sample)#查看各组细胞数
prop.table(table(Idents(S)))
table(Idents(S), S$sample)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(S), S$sample), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###添加分组信息
colnames(Cellratio)[1] <- "subtypeEpi"
colnames(Cellratio)[2] <- "sample"
colnames(Cellratio)[3] <- "proportion"

Cellratio$patient <- substr(Cellratio$sample,1,2)
Cellratio <- subset(Cellratio, patient %in% paste0("P", 1:8)) 
Cellratio$group1 <- "Che_imm"
Cellratio[Cellratio$patient%in%c("P5","P6","P7","P8"),]$group1 <- "Che"
Cellratio$group2 <- "R"
Cellratio[Cellratio$patient%in%c("P3","P4","P7","P8"),]$group2 <- "NR"
Cellratio$group3 <- sub(".*_","",Cellratio$sample)

#提取cancer
Cellratio <- Cellratio[Cellratio$subtypeEpi%in%c("Tumor"),]

#排序
Cellratio$group3 <- factor(Cellratio$group3, levels = c("TN","POST1","POST2"),ordered = T)

#折线图
p <- ggplot(data = Cellratio,aes(x=group3,y=proportion,group=patient,color=patient))+
  geom_line(aes(color=group1,lty=group2))+
  geom_point(aes(color=group1),size=5)+
  scale_color_manual(values = c("Che"="#b6e0f0",
                                "Che_imm"="#ea6743"))+
  theme_classic()+ theme(legend.title=element_blank())+ 
  labs(x="state")

pdf("result/Tumor_prop/prop.pdf", width = 7, height = 5)
p
dev.off()

png(filename = "result/Tumor_prop/prop.png", 
    width= 1500, height= 900, res= 300, units = "px" )
p
dev.off()
