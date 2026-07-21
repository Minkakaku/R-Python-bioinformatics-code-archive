library(Seurat)
library(tidyverse)

# dim8 -- nFeatures 2000
# dim6 -- nFeatures 3000 5000
# dim22 30 -- 5000

# read data -------
file_list <- 
  list.files("/home/ztt/spatial/visium_HD/FWSP20240170", pattern = "A1|D1", full.names = TRUE) %>%
  str_subset('no_alignment', negate = TRUE)

# names(file_list) <- basename(file_list)


visim_ls <-
  file_list %>% 
  file.path("outs") %>%
  map(Load10X_Spatial, bin.size = 8)
names(visim_ls) <- names(file_list) %>% str_remove("_align.*$")


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 


# 'Slc29a1', 'Slc17a6' gene expression ---------

pdf(
  'Spatial feature plot.pdf',
  height = 8, 
  width = 16
)

visim_ls %>%
  imap(function(visium, name){
    SpatialFeaturePlot(
      visium,
      images = 'slice1.008um',
      feature = c('Slc29a1', 'Slc17a6'),
      shape = 21,
      alpha = c(1, 1)
    ) +
    patchwork::plot_annotation(title = name)
  })
dev.off()


pdf(
  'Spatial Vlnplot.pdf',
  height = 5, 
  width = 10
)

visim_ls %>%
  imap(function(visium, name){
    VlnPlot(
      visium,
      feature = c('Slc29a1', 'Slc17a6')
    ) +
    patchwork::plot_annotation(title = name)
  })
dev.off()



visim %>% VlnPlot(features = c('Slc29a1', 'Slc17a6'))


# qc -----
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    subset(visium, nCount_Spatial.008um > 1)
  })

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = "^mt-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^Rp[sl]")
    x$log10GenesPerUMI   <- log10(x$nFeature_Spatial.008um) / log10(x$nCount_Spatial.008um)
    x
  })


p_qc_ls <-
  visim_ls %>%
  imap(function(x, name){ # browser()
    x_sub <- subset(x, cells = Cells(x[['Spatial.008um']]), downsample = 100000)
    VlnPlot(x_sub, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um",
                                "percent.mt", "log10GenesPerUMI"), 
            ncol = 4, log = FALSE, raster = TRUE) +
      patchwork::plot_annotation(title = name)
  })

p_qc_log_ls <-
  visim_ls %>%
  imap(function(x, name){
    x_sub <- subset(x, cells = Cells(x[['Spatial.008um']]), downsample = 100000)
    VlnPlot(x_sub, features = c("nFeature_Spatial.008um", "nCount_Spatial.008um",
                                  "percent.mt", "log10GenesPerUMI"), 
            ncol = 4, log = TRUE, raster = TRUE) +
      patchwork::plot_annotation(title = str_glue('{name} + log'))
  })


# p_qc_no_point_ls <-
#   all_data_ls %>%
#   imap(function(x, name){
#     VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA",
#                             "log10GenesPerUMI"),
#             ncol = 4, log = FALSE, raster = TRUE, pt.size = 0) +
#       patchwork::plot_annotation(title = name)
#   })

# p_qc_no_point_log_ls <-
#   all_data_ls %>%
#   imap(function(x, name){
#     VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
#                             "log10GenesPerUMI"), 
#             ncol = 4, log = TRUE, raster = TRUE, pt.size = 0) +
#       patchwork::plot_annotation(title = str_glue('{name} + log'))
#   })

pdf(
  "qc plot.pdf",
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
p_qc_log_ls %>% walk(print)
dev.off()


# pdf(
#   "qc plot - no point.pdf",
#   height = 6, 
#   width = 20
# )
# p_qc_no_point_ls %>% walk(print)
# p_qc_no_point_log_ls %>% walk(print)
# dev.off()

# qc filter ----------
qc_metric <- read.csv("qc_metric.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  "nFeature_Spatial.008um", "nCount_Spatial.008um",
  "percent.mt", "log10GenesPerUMI"
)
qc_metric_ls <-
  qc_metric %>% 
  apply(1, name = names(.), simplify = FALSE, function(x, name = name){
    qc_ls <- map2(x, name, function(num, variable){
      if(!variable %in% qc_metric_variable) return()
      
      if(str_detect(num, ":")){
        num <- str_split(num, pattern = ":", n = 2) %>% .[[1]]
        qc <- str_glue("{variable} {num[1]} & {variable} {num[2]}")
      }else{
        qc <- str_glue("{variable} {num}")
      }
      qc
    })
    
    res <- qc_ls %>% discard(function(x) length(x) < 1) %>% str_c(collapse = " & ")
    res
  })

names(qc_metric_ls) <- qc_metric$sample



# filter
visim_ls_pass <-
  names(visim_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- visim_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(visim_ls) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]]
    all_data_pass <- visim_ls_pass[[sample]]
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "nFeature_Spatial.008um") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # alll_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "nFeature_Spatial.008um")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass  + patchwork::plot_layout(ncol = 3)
    
    p1 / p2
  })


# export 
pdf(
  "qc plot - pass.pdf",
  height = 12, 
  width = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()


# cells filter
cell_counts_df <-
  names(visim_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]]
    all_data_pass <- visim_ls_pass[[sample]]
    
    structure(
      c(ncol(all_data), ncol(all_data_pass)),
      names = c("Before filter", "After filter")
    )
  }) %>%
  as.data.frame(check.names = FALSE)


# forget to add sample_id
visim_ls_pass <-
  visim_ls_pass %>%
  imap(function(visium, name){
    project <- name %>% str_remove("^.*?-") %>% str_remove("-.*?$")
    visium$project <- project
    RenameCells(visium, add.cell.id = project)
  })

# merge
all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])


# normalize ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
# all_data <- ScaleData(all_data, vars.to.regress = 'nCount_Spatial.008um')

# sketch data --------
all_data <-
  SketchData(
    all_data,
    ncells         = 50000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


### dimension reduction
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 5000)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_Spatial.008um')

all_data <- RunPCA(all_data)

pdf('pca plot.pdf', height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf('dimheatmap 1-5.pdf', height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf('dimheatmap 6-10.pdf', height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf('dimheatmap 11-15.pdf', height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf('dimheatmap 16-20.pdf', height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


# dimension reduction --------
dims <- 6


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(str_glue("Dimplot dim{dims} project.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(str_glue("Dimplot dim{dims} project - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'pca' # 'harmony' 'pca'
# cluster ------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----
# cell_marker <- c(
#   'Csf1r', 'C1qa', 'Hexb', 'Cx3cr1', 'P2ry12', 'Tmem119', 'Tnf', 'Ccl4', # microglia
#   'Sox10','Cldn11', 'Mog', 'Plp1', # oligo
#   'Aqp4', 'Aldh1l1','Gfap', 'Aldoc', # astrocyte
#   'Vtn', 'Flt1', 'Pecam1','Cldn5', 'Adgrf5', # pericyte
#   'Mgp' ,'Slc47a1', 'Dapl1', 'Igf2', 'Sema3g', 'Acta2',  # vascular
#   'Pdgfra', 'Vcan', 'Cspg4', 'Olig1', # OPC
#   'Ccdc153', 'Tmem212', 'Hdc', 'Kcnj8',# ependymal
#   'Pf4', 'Cd74', 'Cxcl2', 'Lyz2', 'Ms4a7',
# )



# gene sets
# microglia -- Hexb, Csf1r, C1qa, P2ry12
# OPCs -- Pdgfra, Vcan, Cspg4, Olig1
# Endo -- Vtn, Flt1, Cldn5
# Oligo -- Plp1, Mpb, Cldn11, Mog

# InN -- Sst, Pvalb, Lamp5, and Vip

# identify striatal neurons
# sc.pl.dotplot(adata, ['Otof', 'Cacng5', 'Th','Ppp1r1b', 'Drd1','Tac1', 'Tshz1', 'Pdyn', 'Drd2','Penk','Adora2a', 'Calb1','Pthlh','Cxcl14','Chat'], groupby='final_clusts',use_raw=True)



cell_marker <- c(
  ExN = 'Slc17a7',
  ExN = 'Cux2', # ExN and some InN
  ExN = 'Rspo1',
  ExN = 'Scube1',
  ExN = 'Fezf2',
  ExN = 'Ndst4', # ExN and InN
  ExN = 'Nxph4',
  ExN = 'Hs3st4',
  ExN = 'Tshz2',
  ExN = 'Ptpru', # ExN and InN
  
  InN = 'Gad2', # InN and MSN
  InN = 'Sst',
  InN = 'Pvalb',
  InN = 'Cpne7',# InN and ExN
  InN = 'Lhx6',
  InN = 'Vip',
  InN = 'Adarb2',
  InN = 'Calb2',
  InN = 'Lamp5', # InN and ExN
  InN = 'Chat',
  
  
  MSN = 'Syt6', # MSN and ExN
  MSN = 'Otof', # MSN and ExN
  MSN = 'Drd1', # MSN and ExN and InN
  MSN = 'Adora2a',
  
  Oligo = 'Plp1',
  Oligo = 'Cldn11',
  Oligo = 'Mog',
  Oligo = 'Mbp',
  Oligo = 'Sox10',
  Oligo = 'Olig2',
  Oligo = 'Olig',
  
  Astro = 'Rorb',
  Astro = 'Aqp4',
  Astro = 'Foxj1',
  Astro = 'Gfap',
  Astro = 'Aldh1l1',
  Astro = 'S100b',

  
  # Vascular = 'Cspg4', # Vascular and Oligo
  Peri = 'Vtn',
  Peri = 'Flt1',
  Endo = 'Cldn5',
  
  Macro = 'F13a1',
  Tcell = 'Cd3e',
  
  Micro = 'Ctss',
  Micro = 'Hexb',
  Micro = 'Csf1r',
  Micro = 'C1qa',
  Micro = 'P2ry12',
  
  OPC = 'Pdgfra',
  OPC = 'Cspg4', # OPC and Peri
  OPC = 'Vcan',
  OPC = 'Olig1'
  
  
)

cell_marker2 <- c(
  'immature neuron' = 'Dcx',
  'immature neuron' = 'Tubb3',
  'immature neuron' = 'Tbr1',
  'immature neuron' = 'Stmn1',
  
  'mature neron' = 'Rbfox3',
  'mature neron' = 'Map2',
  'mature neron' = 'Syp',
  'mature neron' = 'Dlg4',
  
  # excitatory neuron
  'Glutamatergic neurons' = 'Slc17a7', # 'Vglut1',
  'Glutamatergic neurons'  = 'Slc17a6', # 'Vglut2',
  'Glutamatergic neurons'  = 'Grin1', 
  'Glutamatergic neurons'  = 'Grin2b', 
  'Glutamatergic neurons'  = 'Gls', 
  
  # inhibitory neuron
  'GABAergic neurons' = 'Slc32a1',
  'GABAergic neurons' = 'Gad2',
  'GABAergic neurons' = 'Gad1',
  
  'Dopaminergic neurons' = 'Th',
  'Dopaminergic neurons' = 'Slc6a3',
  'Dopaminergic neurons' = 'Foxa2',
  'Dopaminergic neurons' = 'Kcnj6',
  'Dopaminergic neurons' = 'Nr4a2',
  'Dopaminergic neurons' = 'Lmx1b',
  
  'Serotonergic neurons' = 'Tph1',
  'Serotonergic neurons' = 'Tph2',
  'Serotonergic neurons' = 'Slc6a4',
  'Serotonergic neurons' = 'Fev',
  
  # muscle excitatory neuron
  'Cholinergic neuron' = 'Chat',
  'Cholinergic neuron' = 'Ache'

)


major_cell_marker <- c(
  neron = 'Rbfox3', # 'NEUN',
  neron = 'Tubb3',
  neron = 'Map2',
  neron = 'Gap43',
  neron = 'Pcp4',
  
  'Glutamatergic neurons' = 'Slc17a7',
  'GABAergic neurons' = 'Gad2',
  'GABAergic neurons' = 'Gad1',
  'Dopaminergic neurons' = 'Th',
  'Dopaminergic neurons' = 'Aldh1a1',
  'Serotonergic neurons' = 'Tph1',
  'Serotonergic neurons' = 'Tph2',
  'Cholinergic neuron' = 'Chat',
  
  Astro = 'Gfap',
  Astro = 'Aldh1l1',
  Astro = 'S100b',
  
  Oligo = 'Mbp',
  Oligo = 'Plp1',
  Oligo = 'Mobp',
  Oligo = 'Cnp',
  
  'Microglial cells' = 'Dock8',
  'Microglial cells' = 'Sfmbt2',
  'Microglial cells' = 'Cx3cr1',
  'Microglial cells' = 'C1qc',
  
  'Endothelial cells' =	'Flt1',
  'Endothelial cells' =	'Vwf',
  'Endothelial cells' =	'Pecam1',
  'Endothelial cells' =	'Cdh5',
  'T cells' =	'Cd96',
  'T cells' =	'Skap1',
  'T cells' =	'Cd3g',
  'T cells' =	'Cd3e',
  'T cells' =	'Cd3d',
  Macrophages	= 'Plxdc2',
  Macrophages	= 'Cybb',
  Macrophages	= 'Apbb1ip'
)

marker_major3 <- c(
  'Endothelial cells' =	'FLT1',
  'Endothelial cells' =	'VWF',
  'Endothelial cells' =	'PECAM1',
  'Endothelial cells' =	'CDH5',
  'T cells' =	'CD96',
  'T cells' =	'SKAP1',
  'T cells' =	'CD3G',
  'T cells' =	'CD3E',
  'T cells' =	'CD3D',
  Macrophages	= 'PLXDC2',
  Macrophages	= 'CYBB',
  Macrophages	= 'APBB1IP',
  Oligodendrocytes = 'MOBP',
  Oligodendrocytes = 'PLP1',
  Oligodendrocytes = 'MBP',
  'Oligodendrocyte progenitor cells' = 'SOX10',
  'Oligodendrocyte progenitor cells' = 'PCDH15',
  'Oligodendrocyte progenitor cells' = 'LHFPL3',
  'Oligodendrocyte progenitor cells' = 'VCAN',
  'Oligodendrocyte progenitor cells' = 'PDGFRA',
  'Oligodendrocyte progenitor cells' = 'CSPG4',
  'Microglial cells' = 'DOCK8',
  'Microglial cells' = 'SFMBT2',
  'Microglial cells' = 'CX3CR1',
  'Microglial cells' = 'C1QC',
  'Inhibitory neurons' = 'ERBB4',
  'Inhibitory neurons' = 'NXPH1',
  'Inhibitory neurons' = 'GAD2',
  'Excitatory neurons' = 'SYN3',
  'Excitatory neurons' = 'RBFOX3',
  'Excitatory neurons' = 'CAMK2A',
  astrocyte	= 'ALDH1L1',
  astrocyte	= 'AQP4',
  astrocyte	= 'GJA1',
  astrocyte	= 'FGFR3',
  astrocyte	= 'GFAP'
) %>% {structure(
  str_to_title(.),
  names = names(.)
)}

# dotplot
pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = cell_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = cell_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = cell_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = cell_marker)
dev.off()

pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} neuron_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = cell_marker2)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = cell_marker2)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = cell_marker2)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = cell_marker2)
dev.off()

pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_cell_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_cell_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_cell_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_cell_marker)
dev.off()

pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} major_marker3.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = marker_major3 )
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = marker_major3 )
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = marker_major3 )
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = marker_major3 )
dev.off()


# annotation ----
Idents(all_data) <- "sketch_snn_res.0.6"


cell_cluster <- c(
  '0'  = "Pcp4-dim Oligodendrocytes",
  '1'  = "Others",
  '2'  = "Glutamatergic neurons",
  '3'  = "Oligodendrocytes",
  '4'  = "Gfap-dim Oligodendrocytes",
  '5'  = "Glutamatergic neurons",
  '6'  = "Cnpase Oligodendrocytes",
  '7'  = "Pcp4 Neron",
  '8'  = "Glutamatergic neurons",
  '9'  = "Slc17a7 Astrocyte",
  '10' = "Pcp4-dim Neron",
  '11' = "Oligodendrocytes"
)


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
dev.off()


# find markers -------
Idents(all_data) <- "sketch_snn_res.0.6"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )


all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[[1]]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

avg_expr[all_de_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )


# project data ------
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    
    umap.model = "umap", 
    dims = 1:dims, 
    
    refdata = list(annotation_full = "annotation")
  )

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "Spatial.008um"

pdf(str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.umap') + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_full", features = major_cell_marker)
dev.off()


# 需求分析 ---------------
gene_to_check_in_astro <- c('Slc29a1', 'Slc17a6')

DefaultAssay(all_data) == "Spatial.008um"
Idents(all_data) %>% unique

pdf('Spatial vlnplot - split by ident.pdf', height = 5, width = 8)
VlnPlot(
  all_data,
  features = gene_to_check_in_astro,
  # idents = levels(all_data) %>% head(1)
)
dev.off()


gene_to_check <-
  c(
    'Gfap',
    'Aif1', # 'Iba1',
    'Slc29a1', # 'Ent1',
    'Slc17a7', # 'Vglut1',
    'Slc17a6', # 'Vglut2',
    'Scd1',
    'Scd2',
    'Apoe'
  )

avg_data_by_group_and_cell <-
  all_data %>%
  JoinLayers() %>%
  AverageExpression(
    assays = 'Spatial.008um',
    group.by = c('project', 'annotation_full')
  )


p_gene_col <-
  avg_data_by_group_and_cell$Spatial.008um[gene_to_check, ] %>% 
  as.matrix() %>%
  log1p %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  pivot_longer(-gene, names_to = "cols", values_to = 'expr') %>%
  separate(col = cols, into = c("project", "annotation"), sep = "_") %>%
  
  ggplot(aes(x = gene, y = expr, fill = project)) +
  geom_col(position = "dodge") +
  facet_wrap(~ annotation, scales = 'free_y') +
  ggpubr::theme_pubclean() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = 'Normalized Log Expression'
  )
  
pdf('expression barplot of selected genes.pdf', height = 6, width = 10)
p_gene_col %>% print()
dev.off()

# spatial feature 
# pdf(
#   'Spatial feature plot - 8 genes.pdf',
#   height = 8, 
#   width = 16
# )
# 
# visim_ls %>%
#   imap(function(visium, name){
#     SpatialFeaturePlot(
#       visium,
#       images = 'slice1.008um',
#       feature = c('Slc29a1', 'Slc17a6'),
#       shape = 21,
#       alpha = c(1, 1)
#     ) +
#       patchwork::plot_annotation(title = name)
#   })
# dev.off()

images_names <- 
  Images(all_data) %>%
  map(~ all_data@images[[.x]]) %>%
  map(Cells) %>%
  map(1) %>%
  map(str_split_1, pattern = "_") %>%
  map(1) %>%
  unlist()

pdf('spatial plot - annotation.pdf', height = 10, width = 13)
# ImageDimPlot(all_data, fov = "slice1.008um")
SpatialDimPlot(all_data) &
  ggsci::scale_fill_aaas() &
  guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()


pdf('spatial plot - feature.pdf', height = 10, width = 13)
  gene_to_check %>%
  map(function(gene){
    SpatialFeaturePlot(
      all_data,
      features = gene
    )
  })
dev.off()



# slice 1
# coords_1 <- GetTissueCoordinates(all_data, image = )


# ROI分析 ---------------------
roi_out_dir <- "roi_analysis" # dir.create(roi_out_dir)

all_data <- readRDS("all_data_projected.rds")

roi_df <-
  list.files('ROI', full.names = TRUE, pattern = "csv") %>%
  setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
  map(read.csv) %>%
  enframe("group") %>%
  unnest(value) %>%
  mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
  mutate(roi_group = str_c(group, ROI, sep = "_")) %>%
  dplyr::select(roi_group, group, Barcode) %>%
  column_to_rownames("Barcode") 

roi_ls <- 
  split(
    rownames(roi_df),
    roi_df$roi_group
  )

all_data <- all_data %>% AddMetaData(roi_df)

roi_assay <- all_data[['Spatial.008um']] %>% subset(cells = rownames(roi_df))
all_data[['ROI']] <- roi_assay


DefaultAssay(all_data) <- "ROI"



# Images(all_data)
# 
# SpatialDimPlot(all_data, images = "slice1.008um.2", cells.highlight = roi_ls$KA_ROI_2)
# 
# 
# 
# GetTissueCoordinates(all_data[['slice1.008um.2']]) %>%
#   dplyr::filter(cell %in% roi_ls$KA_ROI_2) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_point()
# 

## norm annd pca ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_Spatial.008um')


all_data <- RunPCA(all_data)

pdf(file.path(roi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(roi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(roi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(roi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(roi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


## dr ------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

dims <- 15
all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "ROI_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

## cell marker -----------
major_cell_marker <- c(
  neron = 'Rbfox3', # 'NEUN',
  neron = 'Tubb3',
  neron = 'Map2',
  neron = 'Gap43',
  neron = 'Pcp4',
  
  'Glutamatergic neurons' = 'Slc17a7',
  'GABAergic neurons' = 'Gad2',
  'GABAergic neurons' = 'Gad1',
  'Dopaminergic neurons' = 'Th',
  'Dopaminergic neurons' = 'Aldh1a1',
  'Serotonergic neurons' = 'Tph1',
  'Serotonergic neurons' = 'Tph2',
  'Cholinergic neuron' = 'Chat',
  
  Astro = 'Gfap',
  Astro = 'Aldh1l1',
  Astro = 'S100b',
  
  Oligo = 'Mbp',
  Oligo = 'Plp1',
  Oligo = 'Mobp',
  Oligo = 'Cnp',
  
  'Microglial cells' = 'Dock8',
  'Microglial cells' = 'Sfmbt2',
  'Microglial cells' = 'Cx3cr1',
  'Microglial cells' = 'C1qc',
  
  'Endothelial cells' =	'Flt1',
  'Endothelial cells' =	'Vwf',
  'Endothelial cells' =	'Pecam1',
  'Endothelial cells' =	'Cdh5',
  'T cells' =	'Cd96',
  'T cells' =	'Skap1',
  'T cells' =	'Cd3g',
  'T cells' =	'Cd3e',
  'T cells' =	'Cd3d',
  Macrophages	= 'Plxdc2',
  Macrophages	= 'Cybb',
  Macrophages	= 'Apbb1ip'
)


pdf(file.path(roi_out_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "ROI_snn_res.0.2", features = major_cell_marker)
DotPlot(all_data, group.by = "ROI_snn_res.0.4", features = major_cell_marker)
DotPlot(all_data, group.by = "ROI_snn_res.0.6", features = major_cell_marker)
DotPlot(all_data, group.by = "ROI_snn_res.0.8", features = major_cell_marker)
dev.off()

## annotation ----
Idents(all_data) <- "ROI_snn_res.0.6"


cell_cluster <- c(
  '0'   = "Oligo",
  '1'   = "Slc17a7 Astrocyte",
  '2'   = "Pcp4 Glutamatergic neurons",
  '3'   = "Endothelial cells",
  '4'   = "C1qc-dim Astrocyte",
  '5'   = "Cnp Oligo",
  '6'   = "Pcp4-dim Glutamatergic neurons",
  '7'   = "Astrocyte",
  '8'   = "Glutamatergic neurons",
  '9'   = "Pcp4 Oligo",
  '10'  = "GABAergic neurons"
)


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


saveRDS(all_data, file.path(roi_out_dir, "all_data_roi.rds"))

pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(roi_out_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
dev.off()





## 需求分析 ---------------
gene_to_check_in_astro <- c('Slc29a1', 'Slc17a6')

DefaultAssay(all_data) == "ROI"
Idents(all_data) %>% unique

pdf(file.path(roi_out_dir, 'Spatial vlnplot - split by ident.pdf'), height = 5, width = 8)
gene_to_check_in_astro %>%
  map(function(gene){
    VlnPlot(
      all_data,
      features = gene,
      # idents = levels(all_data) %>% head(1)
      # pt.size = 0, 
      split.by = "group"
    )
  })
dev.off()


gene_to_check <-
  c(
    'Gfap',
    'Aif1', # 'Iba1',
    'Slc29a1', # 'Ent1',
    'Slc17a7', # 'Vglut1',
    'Slc17a6', # 'Vglut2',
    'Scd1',
    'Scd2',
    'Apoe'
  )


# TODO: Gfap表达在KA中有表达，在WT中无表达
# 需要考虑使用IntegratedData进行含量校正
pdf(file.path(roi_out_dir, 'Spatial vlnplot - split by ident - 8genes.pdf'), height = 5, width = 8)
gene_to_check %>%
  map(function(gene){
    VlnPlot(
      all_data,
      features = gene,
      # idents = levels(all_data) %>% head(1)
      pt.size = 0,
      split.by = "group"
    )
  })
dev.off()

avg_data_by_group_and_cell <-
  all_data %>%
  JoinLayers() %>%
  AverageExpression(
    assays = 'ROI',
    group.by = c('project', 'annotation')
  )

avg_data_by_group_and_cell2 <-
  all_data %>%
  JoinLayers() %>%
  AverageExpression(
    assays = 'ROI',
    group.by = c('roi_group', 'annotation')
  )

p_gene_col <-
  avg_data_by_group_and_cell$ROI[gene_to_check, ] %>%
  # avg_data_by_group_and_cell$Spatial.008um[gene_to_check, ] %>% 
  as.matrix() %>%
  log1p %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  pivot_longer(-gene, names_to = "cols", values_to = 'expr') %>%
  separate(col = cols, into = c("project", "annotation"), sep = "_") %>%
  
  ggplot(aes(x = gene, y = expr, fill = project)) +
  geom_col(position = "dodge") +
  facet_wrap(~ annotation, scales = 'free_y') +
  ggpubr::theme_pubclean() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = 'Normalized Log Expression'
  )


p_gene_col2 <-
  avg_data_by_group_and_cell2$ROI[gene_to_check, ] %>%
  # avg_data_by_group_and_cell$Spatial.008um[gene_to_check, ] %>% 
  as.matrix() %>%
  log1p %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  pivot_longer(-gene, names_to = "cols", values_to = 'expr') %>%
  separate(col = cols, into = c("project", "annotation"), sep = "_") %>%
  
  ggplot(aes(x = gene, y = expr, fill = project)) +
  geom_col(position = "dodge") +
  facet_wrap(~ annotation, scales = 'free_y') +
  ggpubr::theme_pubclean() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = 'Normalized Log Expression'
  )

pdf(file.path(roi_out_dir, 'expression barplot of selected genes.pdf'), height = 6, width = 10)
p_gene_col %>% print()
dev.off()

pdf(file.path(roi_out_dir, 'expression barplot of selected genes 2.pdf'), height = 6, width = 14)
p_gene_col2 %>% print()
dev.off()

# spatial feature 
# pdf(
#   'Spatial feature plot - 8 genes.pdf',
#   height = 8, 
#   width = 16
# )
# 
# visim_ls %>%
#   imap(function(visium, name){
#     SpatialFeaturePlot(
#       visium,
#       images = 'slice1.008um',
#       feature = c('Slc29a1', 'Slc17a6'),
#       shape = 21,
#       alpha = c(1, 1)
#     ) +
#       patchwork::plot_annotation(title = name)
#   })
# dev.off()

images_names <- 
  Images(all_data) %>%
  map(~ all_data@images[[.x]]) %>%
  map(Cells) %>%
  map(1) %>%
  map(str_split_1, pattern = "_") %>%
  map(1) %>%
  unlist()


p_spatial_dimplot <-
  map2(Images(all_data), images_names, function(image, image_name){
    SpatialDimPlot(all_data, images = image) &
      ggtitle(image_name) &
      ggsci::scale_fill_igv() &
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })

pdf(file.path(roi_out_dir, 'spatial plot - annotation.pdf'), height = 10, width = 13)
# ImageDimPlot(all_data, fov = "slice1.008um")
p_spatial_dimplot %>% walk(print)
dev.off()


p_spatial_featureplot <-
  gene_to_check %>%
  map(function(gene){
    map2(Images(all_data), images_names, function(image, image_name){
      SpatialFeaturePlot(
        all_data,
        images = image,
        features = gene
      ) &
        ggtitle(image_name)
    })
  })

pdf(file.path(roi_out_dir, 'spatial plot - feature.pdf'), height = 10, width = 13)
p_spatial_featureplot %>% walk(print)
dev.off()




# all_data_roi <- subset(all_data, cells = colnames(all_data)[!is.na(all_data$roi_group)])

# gene_to_check_in_astro <- c('Slc29a1', 'Slc17a6')
# 
# Idents(all_data) %>% unique

# VlnPlot(
#   all_data,
#   assay = "ROI",
#   features = gene_to_check_in_astro,
#   # idents = levels(all_data) %>% head(1)
# )


