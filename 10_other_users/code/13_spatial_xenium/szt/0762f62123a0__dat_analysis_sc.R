library(Seurat)
library(tidyverse)

sc_out_dir <- "sc_output" # dir.create(sc_out_dir)

file_dir <- "."

all_files_dir <- "GSE182786" # list.files(file_dir, pattern = "^[A-Za-z0-9_.]*$")
all_files <- list.files(all_files_dir, full.names = TRUE)
names(all_files) <- all_files_dir

# sample_info <-
#   data.frame(
#     sample = all_files_dir,
#     id = seq_along(all_files_dir)
#   )

# 
# # 6 samples
# stopifnot(length(all_files) == 6)


all_data_ls <- 
  all_files %>%
  map(read.table, sep = "\t") %>%
  imap(~CreateSeuratObject(
    .x,
    project = .y,
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    min.cells = 10,
    min.features = 10
  ))

all_data_ls$GSE182786$project <-
  all_data_ls$GSE182786 %>% 
  colnames() %>%
  str_split("\\.") %>%
  map(1) %>%
  unlist(use.names = FALSE)

all_data_ls <- all_data_ls$GSE182786 %>% SplitObject(split.by = "project")

# all_data_ls <- 
#   all_data_ls %>%
#   map(function(x){
#     x$project <- x@project.name
#     x
#   })

### QC
all_data_ls <- 
  all_data_ls %>%
  map(function(x){
    x[["percent.mt"]]  <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^RP[SL]")
    x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
    x
  })


p_qc_ls <-
  all_data_ls %>%
  imap(function(x, name){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
                            "percent.mt", "percent.ribo", "log10GenesPerUMI"), 
            ncol = 5, log = FALSE) +
      patchwork::plot_annotation(title = name)
  })

p_qc_log_ls <-
  all_data_ls %>%
  imap(function(x, name){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
                            "percent.mt", "percent.ribo", "log10GenesPerUMI"), 
            ncol = 5, log = TRUE) +
      patchwork::plot_annotation(title = str_glue('{name} + log'))
  })

pdf(
  file.path(sc_out_dir, "qc plot.pdf"),
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
p_qc_log_ls %>% walk(print)
dev.off()




### filter
qc_metric <- read.csv("qc_metric.csv") %>% map_dfc(str_trim)
qc_metric_variable <- c(
  'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'log10GenesPerUMI'
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
all_data_ls_pass <-
  names(all_data_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- all_data_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(all_data_ls) %>%
  map(function(sample){
    all_data <- all_data_ls[[sample]]
    all_data_pass <- all_data_ls_pass[[sample]]
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # alll_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass  + patchwork::plot_layout(ncol = 3)
    
    p1 / p2
  })

# export 
pdf(
  file.path(sc_out_dir, "qc plot - pass.pdf"),
  height = 12, 
  width = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()


# cells filter
cell_counts_df <-
  names(all_data_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- all_data_ls[[sample]]
    all_data_pass <- all_data_ls_pass[[sample]]
    
    structure(
      c(ncol(all_data), ncol(all_data_pass)),
      names = c("Before filter", "After filter")
    )
  }) %>%
  as.data.frame(check.names = FALSE)

# #### add cell_id prefix
# all_data_ls_pass <-
#   all_data_ls_pass %>%
#   names() %>%
#   set_names(., .) %>%
#   map(function(x){
#     all_data <- all_data_ls_pass[[x]]
#     ids <- sample_info$id[sample_info$sample == x]
#     
#     all_data <- RenameCells(all_data, add.cell.id = ids)
#     all_data
#   })

all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])


### normalize
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')


### dimension reduction
all_data <- RunPCA(all_data)

pdf(file.path(sc_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(sc_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()
pdf(file.path(sc_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()
pdf(file.path(sc_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()
pdf(file.path(sc_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

saveRDS(all_data, "all_data_merge.rds")

# ndims 15   ---------------------

all_data <- RunUMAP(all_data, dims = 1:15)

pdf(file.path(sc_out_dir, "Dimplot dim15 project.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:15,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:15, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_"
)


pdf(file.path(sc_out_dir, "Dimplot dim15 project - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### cluster
all_data <-
  FindNeighbors(
    all_data,
    reduction = "harmony",
    dims = 1:15,
    k.param = 15
  )
all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(sc_out_dir, "Dimplot dim15 cluster - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# dotplot
marker_ls <- c(
  Myoid = "ACTA2",
  Myoid = "MYH11",
  Myoid = "PTCH1",
  Myoid = "PTCH2",
  
  Sertoli = "WT1",
  Sertoli = "SOX9",
  Sertoli = "BEX2",
  Sertoli = "GATA4",
  
  Macrophage = "CD68",
  Macrophage = "MS4A6A",
  Macrophage = "LYZ",
  Macrophage = "IL3RA",
  
  Leydig = "CYP11A1",
  Leydig = "STAR",
  Leydig = "IGF1",
  Leydig = "IGFBP5",
  
  Endothelial = "VWF",
  Endothelial = "PECAM1",
  Endothelial = "CD34",
  Endothelial = "TGFBR2",
  
  'T cell' = "CD3D",
  'T cell' = "CD3E",
  'T cell' = "CD3G",
  'T cell' = "CD2",
  
  Spermatogonia = "FGFR3",
  Spermatogonia = "STRA8",
  Spermatogonia = "DMRT1",
  Spermatogonia = "UCHL1",
  Spermatogonia = "TCF3",
  
  Spermatocytes = "SYCP1",
  Spermatocytes = "SMC1B",
  Spermatocytes = "SYCP3",
  Spermatocytes = "OVOL2",
  
  Spermatids = "TNP1",
  Spermatids = "PRM1",
  Spermatids = "TXNDC2",
  Spermatids = "PRM3"
)


pdf(file.path(sc_out_dir, "dotplot dims15.pdf"), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls)
dev.off()

# saveRDS(all_data, 'all_data_sc_dim15.rds')

# script sc_analysis.R --dims 15 --rds all_data_dim8.rds --save --k.param 15

library(Seurat)
library(tidyverse)

all_data <- readRDS('all_data_sc_dim15.rds')


# cell_cluster <- c(
#   '0'  = "Spermatids",
#   '1'  = "PRM3-dim Spermatids",
#   '2'  = "Endothelial",
#   '3'  = "PRM3 Spermatids",
#   '4'  = "Myoid+Leydig",
#   '5'  = "Spermatocytes",
#   '6'  = "PRM3-dim TXNDC2+ Spermatids",
#   '7'  = "Spermatogonia",
#   '8'  = "Smooth muscle cells",
#   '9'  = "Sertoli",
#   '10' = "Macrophage",
#   '11' = "Unknown"
# )

cell_cluster <- c(
  '0'  = "PRM3-dim Spermatids",
  '1'  = "Spermatids",
  '2'  = "PRM3 Spermatids",
  '3'  = "Spermatids",
  '4'  = "IL3RA Endothelial",
  '5'  = "Spermatocytes",
  '6'  = "PRM3-dim TXNDC2+ Spermatids",
  '7'  = "Spermatogonia",
  '8'  = "Sertoli",
  '9'  = "Leydig",
  '10' = "Myoid",
  '11' = "Endothelial",
  '12' = "Smooth muscle cells",
  '13' = "Smooth muscle cells",
  '14' = "Macrophage",
  '15' = "Unknown",
  '16' = "Unknown"
)

# Idents(all_data) <- "RNA_snn_res.0.2"
Idents(all_data) <- "RNA_snn_res.0.4"
all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(sc_out_dir, "Dimplot dim15 annotation - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# cell props
all_data$group <- all_data$project # %>% str_remove("[0-9]$")
plot_data <- 
  all_data[[]] %>%
  dplyr::select(
    group,
    annotation
  ) %>% 
  group_by(group, annotation) %>% 
  summarise(count = n())%>%
  group_by(group)%>%
  mutate(cell_props = count/sum(count)) %>%
  mutate(cell_props_label = cumsum(rev(cell_props)))
 

p_cell_props <-
  plot_data %>%
  ggplot(aes(x = group, y = cell_props, fill = annotation)) +
  geom_col() 
  # geom_text(aes(label = cell_props, y = cell_props_label), vjust = 1)
  

pdf(file.path(sc_out_dir, "cell props - dims 15.pdf"), height = 6, width = 8)
print(p_cell_props)
dev.off()

#### 
pdf(file.path(sc_out_dir, "Dimplot dim15 group - harmony.pdf"), height = 5, width = 15)
DimPlot(all_data, group.by = "group", label = TRUE, reduction = "harmonyumap") 
DimPlot(all_data, group.by = "annotation", split.by = "group", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nCount_RNA",split.by = "group",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "nFeature_RNA",split.by = "group",  label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "percent.mt", split.by = "group", label = TRUE, reduction = "harmonyumap") 
FeaturePlot(all_data, features = "log10GenesPerUMI",split.by = "group",  label = TRUE, reduction = "harmonyumap")
dev.off()


# pdf(str_glue("dim15/Dimplot dim15 project - harmony.pdf"), height = 5, width = 15)
# DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") 
# DimPlot(all_data, group.by = "annotation", split.by = "project", label = TRUE, reduction = "harmonyumap") 
# FeaturePlot(all_data, features = "nCount_RNA",split.by = "project",  label = TRUE, reduction = "harmonyumap") 
# FeaturePlot(all_data, features = "nFeature_RNA",split.by = "project",  label = TRUE, reduction = "harmonyumap") 
# FeaturePlot(all_data, features = "percent.mt", split.by = "project", label = TRUE, reduction = "harmonyumap") 
# FeaturePlot(all_data, features = "percent.ribo", split.by = "project", label = TRUE, reduction = "harmonyumap") 
# FeaturePlot(all_data, features = "log10GenesPerUMI",split.by = "project",  label = TRUE, reduction = "harmonyumap")
# dev.off()

  
# all_data[["percent.ribo"]]  <- PercentageFeatureSet(all_data, pattern = "^Rp[sl]")

# saveRDS(all_data, 'all_data_sc_dim15_anno.rds')
