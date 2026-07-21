# installed into linux server
# mamba install conda-forge::jags
# install.packages('rjags')
# BiocManager::install("infercnv")


# BiocManager::install("infercnv")
# https://sourceforge.net/projects/mcmc-jags/files
# If you are using R 4.1.3 or earlier then install JAGS-4.3.0.exe


# 1.所有上皮细胞和down sample部分免疫细胞，跑inferCNV分析，
# 识别出正常肝细胞（肝母细胞/胆管上皮细胞？），恶性上皮细胞。
# 恶性上皮细胞跑cNMF分析，识别出上皮细胞程序。
library(Seurat)
library(tidyverse)
library(infercnv)

###  epi subcluster --------------
epi_out_dir <- "epi_subcluster"
dir.create(epi_out_dir)

all_data <- readRDS("all_data_projected.rds")
sc_epi <- all_data %>% JoinLayers() %>% .[, .$annotation_full %in% c("Hepatocyte", "Epithelial cells")]
DefaultAssay(sc_epi) <- "RNA"
sc_epi[['sketch']] <- NULL

VlnPlot(sc_epi, features = "nFeature_RNA", group.by = 'project', pt.size = 0)

sc_epi[['RNA']] <- split(sc_epi[['RNA']], sc_epi$project, layers = "counts")


# cellcycle
sc_epi <- 
  sc_epi %>%
  CellCycleScoring(
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes
  )

# preprocess
sc_epi <-
  sc_epi %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 4000) %>%
  NormalizeData(scale.factor = median(sc_epi$nCount_RNA)) %>%
  ScaleData(vars.to.regress = c('nCount_RNA')) # , 'S.Score', 'G2M.Score'))

# sc_epi <- sc_epi %>% ScaleData(vars.to.regress = c('nCount_RNA', 'S.Score', 'G2M.Score'))

VariableFeaturePlot(sc_epi)
# VariableFeatures(sc_epi) <- VariableFeatures(sc_epi) %>% setdiff(c(ig_genes, tcr_genes, mt_genes, ribo_genes))



## pca -----
sc_epi <- sc_epi %>% RunPCA()

pdf(file.path(epi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(sc_epi, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_t <- 
  sc_epi[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_epi) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(epi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_t], dims = 1:5)
dev.off()
pdf(file.path(epi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_t], dims = 6:10)
dev.off()
pdf(file.path(epi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_t], dims = 11:15)
dev.off()
pdf(file.path(epi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_t], dims = 16:20)
dev.off()


## dimension reduction --------
epi_dims <-  6 # 10


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_epi <- RunUMAP(sc_epi, dims = 1:epi_dims)


pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_epi <- sc_epi %>% 
  harmony::RunHarmony(
    dims.use = 1:epi_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_epi <- RunUMAP(
  sc_epi, 
  dims = 1:epi_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# plot QC figure
pdf(file.path(epi_out_dir, str_glue("qc - Dimplot dim{epi_dims}.pdf")), height = 6, width = 8)
FeaturePlot(sc_epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI')) + NoLegend()
FeaturePlot(sc_epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(epi_out_dir, str_glue("qc - Vlnplot dim{epi_dims}.pdf")), height = 6, width = 8)
VlnPlot(sc_epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), group.by = 'project', pt.size = 0) + NoLegend()
dev.off()

FeaturePlot(sc_epi, features = c('S.Score', 'G2M.Score'), reduction = "harmonyumap")
DimPlot(sc_epi, group.by = c('Phase'), reduction = "harmonyumap")


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
sc_epi <-
  FindNeighbors(
    sc_epi,
    reduction = reduction.use,
    dims = 1:epi_dims,
    k.param = 10
  )

sc_epi <-
  FindClusters(
    sc_epi,
    resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_epi, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# Dotplot -------
table(sc_epi$RNA_snn_res.0.05)
pdf(
  file.path(epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use}.pdf")), 
  height = 6,
  width = 8
)
epi_feature <- c(
  'epithelial cells'  = 'EPCAM',
  'epithelial cells'  = 'KRT5',
  'epithelial cells'  = 'KRT14',
  
  Hepatocyte          = 'ALB',
  Hepatocyte          = 'APOB',
  Hepatocyte          = 'CYP3A5',
  Hepatocyte          = 'CPB2'
)
DotPlot(sc_epi, features = epi_feature, group.by = "RNA_snn_res.0.05")
DotPlot(sc_epi, features = epi_feature, group.by = "RNA_snn_res.0.1")
dev.off()

# saveRDS(sc_epi, file.path(epi_out_dir, str_glue("sc_epi-{epi_dims}.rds")))
# cluster_res_used <- "RNA_snn_res.0.1"
cluster_res_used <- "RNA_snn_res.0.05"

sc_epi_anno_meta <-
  sc_epi[[cluster_res_used]] %>% 
  rownames_to_column("cell") %>%
  mutate({{cluster_res_used}} := str_c("Epi", .data[[cluster_res_used]], sep = "_"))

all_data_anno_meta <- 
  all_data[['annotation_full']] %>%
  rownames_to_column('cell') %>%
  left_join(sc_epi_anno_meta, by = "cell") %>%
  mutate(annotation2 = ifelse(
    annotation_full %in% c('Hepatocyte', 'Epithelial cells'),
    .data[[cluster_res_used]],
    annotation_full
  )) %>% 
  dplyr::select(cell, annotation2) %>%
  column_to_rownames("cell") 

all_data <- all_data %>% AddMetaData(all_data_anno_meta)  

### inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹

infercnv_out_dir = "inferCNV"
dir.create(infercnv_out_dir) 


meta_data_infer <-
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(
    annotation_full %in% c("Hepatocyte", "Epithelial cells", "B/Plasma cells", "T", "Fibroblasts","Pericytes/SMCs")
  ) %>%
  dplyr::select(cell_id, annotation_full, annotation2) %>%
  map_dfc(as.character) 

# unique(meta_data_infer$annotation)
# "Acinar cells"     "Ductal cells"     "Endocrine cells"  "Epithelial cells" "NKT cells"


# cell is too much
# less than 50000
# 取1/20细胞： 1179497 --> 47179.88
# 1662182 --> 83000
meta_data_infer_subset <-
  meta_data_infer %>%
  group_by(annotation2) %>%
  slice_sample(prop = 0.1)

# epi cells no downsampling
meta_data_infer_subset <- 
  # rbind(
  #   meta_data_infer_subset[!meta_data_infer_subset$annotation_full %in% c("Hepatocyte", "Epithelial cells"), ],
  #   meta_data_infer[meta_data_infer$annotation_full %in% c("Hepatocyte", "Epithelial cells"), ]
  # ) %>%
  meta_data_infer_subset %>%
  ungroup() %>%
  dplyr::select(cell_id, annotation = annotation2)


# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )
raw_counts_matrix_infercnv <- LayerData(JoinLayers(all_data), "counts")
infercnv_obj = infercnv::CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix_infercnv[, meta_data_infer_subset$cell_id],
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path('gtf', "GRCh38.infercnv.dropdup.gtf"),
  
  delim             = "\t",
  ref_group_names   = c("B/Plasma cells", "Fibroblasts", "T")
)


# raw_counts_matrix_infercnv %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

raw_counts_matrix_infercnv[, meta_data_infer_subset$cell_id] %>% rowMeans %>% { sum(. < 0.01) } # 87 <--- 
raw_counts_matrix_infercnv[, meta_data_infer_subset$cell_id] %>% rowMeans %>% { sum(. < 0.005) } # 34

# options()$bitmapType
# [1] "cairo"
# options(bitmapType="Xlib")
# options(bitmapType="cairo")

options(scipen = 100)
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff              =    0.01, 
  out_dir             =    infercnv_out_dir,
  cluster_by_groups   =    TRUE, 
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  write_expr_matrix   =    TRUE,
  png_res             =    300,
  no_plot = T #  figure margin too large ----
)

# C stack error
infercnv::plot_cnv(
  infercnv_obj_default,
  out_dir           =  infercnv_out_dir,
  output_filename   =  'infercnv',
  output_format     =  "pdf",
  png_res = 60,
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)


# download inferCNV object from linux server
# plot inferCNV locally
saveRDS(infercnv_obj_default, file.path(infercnv_out_dir, 'infercnv_obj_default.rds'))

infercnv_obj_default <- readRDS("inferCNV/infercnv_obj_default.rds")

infercnv::plot_cnv(
  infercnv_obj_default,
  out_dir           =  infercnv_out_dir, 
  output_filename   =  'infercnv',
  # output_format     =  "pdf",
  png_res = 300,
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)

infercnv_obj_default_median_filtering <-
  infercnv_obj_default %>%
  infercnv::apply_median_filtering()

infercnv_obj_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           =  'inferCNV', 
    output_filename   =  'infercnv_median_filtering',
    output_format     =  "pdf",
    png_res = 300,
    cluster_by_groups =  TRUE,
    x.range           =  "auto",
    x.center          =  1,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
  )

saveRDS(infercnv_obj_default_median_filtering, file.path(infercnv_out_dir, 'infercnv_obj_default_median_filtering.rds'))


# options()$bitmapType
# [1] "cairo"
# options(bitmapType="Xlib")
# options(bitmapType="cairo")

# CNV_score and CNV_correlation --------------
observation_idx <- infercnv_obj_default@observation_grouped_cell_indices
tumor_observation_idx <- observation_idx[! names(observation_idx) %in% c("Pericytes/SMCs")]

cnv_score <- 
  infercnv_obj_default@expr.data[, unlist(tumor_observation_idx)] %>%
  {. ^ 2} %>%
  colMeans()

# correlation to all non-malignant cells(gene expression)
mean_tumor_cnv <- 
  infercnv_obj_default@expr.data[, unlist(tumor_observation_idx)] %>%
  rowMeans()
  
cnv_correlation <-
  cor(
    mean_tumor_cnv,
    infercnv_obj_default@expr.data[, unlist(tumor_observation_idx)]
  )

stopifnot(length(cnv_score) == length(cnv_correlation))


data.frame(
  cnv_score = cnv_score,
  cnv_correlation = t(cnv_correlation)
) %>% 
  ggplot(aes(x = cnv_score, y = cnv_correlation)) +
  geom_point()


cell_names <- infercnv_obj_default@expr.data %>% colnames()
infercnv_obj_default@observation_grouped_cell_indices %>%
  map(~ cell_names[.x]) %>%
  enframe("cell", "cell_id") %>%
  unnest(cell_id) %>%
  left_join(
    enframe(cnv_score, name = "cell_id", value = "cnv_score"),
    by = "cell_id"
  ) %>%
  dplyr::group_by(cell) %>%
  summarise(cnv_score = mean(cnv_score))


# cNMF ----------------
library(tidyverse)
library(Seurat)

cnmf_out_dir <- "cnmf"
dir.create(cnmf_out_dir)

# all epi cells arge tumor
# sc_epi 
sc_epi$RNA_snn_res.0.05 %>% unique

sc_epi_subset_meta <-
  sc_epi[['RNA_snn_res.0.05']] %>%
  rownames_to_column("cell_id") %>%
  group_by(RNA_snn_res.0.05) %>%
  slice_sample(n = 50000)


# sc_epi_expr_counts <- sc_epi %>% FetchData(layer = "counts", vars = rownames(sc_epi))
sc_epi_expr_counts <- 
  sc_epi[, sc_epi_subset_meta$cell_id] %>%
  JoinLayers() %>%
  LayerData(layer = "counts")
sc_epi_expr_counts <- t(sc_epi_expr_counts)

# colSums
data.frame(count = rowSums(sc_epi_expr_counts)) %>% 
  ggplot(aes(x = 1, y = count)) +
  geom_violin() +
  scale_y_log10()

# hvg_gene
epi_hvg_gene <- sc_epi %>% FindVariableFeatures(assay = "RNA", nfeatures = 4000) %>% VariableFeatures()

# if error occurs: x cells have zero counts of overdispersed genes.
# remove them, rerun cnmf.
sc_epi_expr_counts_subset <- 
  sc_epi_expr_counts[
    rowSums(sc_epi_expr_counts) > 100,
    epi_hvg_gene
  ]

write.table(
  sc_epi_expr_counts_subset,
  file = file.path(cnmf_out_dir, "counts.txt"),
  quote = F,
  sep = "\t",
  row.names = T,
  col.names = T
)


##########run cnmf in linux command window ############

mkdir cnmf && cd cnmf
cnmf prepare --output-dir . --name cNMF -c counts.txt -k 5 6 7 8 9 10 11 12 13 --n-iter 100 --seed 14

cnmf factorize --output-dir . --name cNMF --worker-index 0 --total-workers 1 # 2.5h

cnmf combine --output-dir . --name cNMF

cnmf k_selection_plot --output-dir . --name cNMF

cnmf consensus --output-dir . --name cNMF --components 10 --local-density-threshold 0.01 --show-clustering
cnmf consensus --output-dir . --name cNMF --components 10 --local-density-threshold 0.1 --show-clustering
cnmf consensus --output-dir . --name cNMF --components 8 --local-density-threshold 0.1 --show-clustering
cnmf consensus --output-dir . --name cNMF --components 13 --local-density-threshold 0.1 --show-clustering











