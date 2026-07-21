# installed into linux server
# mamba install conda-forge::jags
# install.packages('rjags')
# BiocManager::install("infercnv")


# BiocManager::install("infercnv")
# https://sourceforge.net/projects/mcmc-jags/files
# If you are using R 4.1.3 or earlier then install JAGS-4.3.0.exe

library(tidyverse)
library(Seurat)



all_data <- readRDS("output/data/adata_celltype.rds")
all_data_raw <- readRDS("output/data/adata_raw.rds")


all_data_raw <-
  all_data_raw %>% 
  NormalizeData(scale.factor = median(all_data_raw$total_counts)) %>% 
  FindVariableFeatures()


# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb




# 1.所有上皮细胞和down sample部分免疫细胞，跑inferCNV分析，
# 识别出正常肝细胞（肝母细胞/胆管上皮细胞？），恶性上皮细胞。
# 恶性上皮细胞跑cNMF分析，识别出上皮细胞程序。

library(infercnv)

# ###  epi subcluster --------------
# epi_out_dir <- "epi_subcluster"
# dir.create(epi_out_dir)
# 
# 
# sc_epi <- all_data %>% JoinLayers() %>% .[, .$annotation %in% c("AT2","TRB-SC","Basal","Club cells")]
# 
# 
# VlnPlot(sc_epi, features = "nCount_Spatial.008um", group.by = 'project', pt.size = 0)
# 
# sc_epi[['RNA']] <- split(sc_epi[['RNA']], sc_epi$project, layers = "counts")
# 
# 
# # cellcycle
# sc_epi <- 
#   sc_epi %>%
#   CellCycleScoring(
#     s.features = cc.genes.updated.2019$s.genes,
#     g2m.features = cc.genes.updated.2019$g2m.genes
#   )
# 
# # preprocess
# sc_epi <- 
#   sc_epi %>% 
#   NormalizeData(scale.factor = floor(mean(all_data$nCount_Spatial.008um)))%>%
#   FindVariableFeatures(nfeatures = 3000) %>%
#   ScaleData(
#     vars.to.regress = c(
#       'nCount_Spatial.008um'
#     )
#   )
# 
# 
# 
# # sc_epi <- sc_epi %>% ScaleData(vars.to.regress = c('nCount_RNA', 'S.Score', 'G2M.Score'))
# 
# VariableFeaturePlot(sc_epi)
# # VariableFeatures(sc_epi) <- VariableFeatures(sc_epi) %>% setdiff(c(ig_genes, tcr_genes, mt_genes, ribo_genes))
# 
# 
# 
# ## pca -----
# sc_epi <- sc_epi %>% RunPCA()
# 
# pdf(file.path(epi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
# ElbowPlot(sc_epi, ndims = 30)
# dev.off()
# 
# set.seed(1234)
# suset_cells_t <- 
#   sc_epi[[]] %>%
#   rownames_to_column("cell_id") %>%
#   group_by(project) %>%
#   slice_sample(prop = 50000/ncol(sc_epi) ) %>%
#   dplyr::pull("cell_id")
# 
# 
# pdf(file.path(epi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(sc_epi[, suset_cells_t], dims = 1:5)
# dev.off()
# pdf(file.path(epi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(sc_epi[, suset_cells_t], dims = 6:10)
# dev.off()
# pdf(file.path(epi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(sc_epi[, suset_cells_t], dims = 11:15)
# dev.off()
# pdf(file.path(epi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(sc_epi[, suset_cells_t], dims = 16:20)
# dev.off()
# 
# 
# ## dimension reduction --------
# epi_dims <-  20
# 
# 
# # library(future)
# # plan(multisession, workers = 10)
# # options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
# 
# 
# sc_epi <- RunUMAP(sc_epi, dims = 1:epi_dims)
# 
# 
# pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} project.pdf")), height = 6, width = 8)
# DimPlot(sc_epi, group.by = "project", label = TRUE) + NoLegend()
# dev.off()
# 
# 
# sc_epi <- sc_epi %>% 
#   harmony::RunHarmony(
#     dims.use = 1:epi_dims,
#     group.by.vars = "project",
#     reduction.save = "harmony"
#     # theta = 3, # default 2, bigger value -- more diversity per cluster
#     # lambda = 0.5, # default, smaller value -- more integrated
#   )
# 
# sc_epi <- RunUMAP(
#   sc_epi, 
#   dims = 1:epi_dims, 
#   reduction = "harmony", 
#   reduction.name = "harmonyumap", 
#   reduction.key  = "harmonyUMAP_",
#   return.model = TRUE
# )
# 
# 
# pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} project - harmony.pdf")), height = 6, width = 8)
# DimPlot(sc_epi, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()
# 
# 
# # plot QC figure
# pdf(file.path(epi_out_dir, str_glue("qc - Dimplot dim{epi_dims}.pdf")), height = 6, width = 8)
# FeaturePlot(sc_epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI')) + NoLegend()
# FeaturePlot(sc_epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), reduction = "harmonyumap") + NoLegend()
# dev.off()
# 
# pdf(file.path(epi_out_dir, str_glue("qc - Vlnplot dim{epi_dims}.pdf")), height = 6, width = 8)
# VlnPlot(sc_epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), group.by = 'project', pt.size = 0) + NoLegend()
# dev.off()
# 
# FeaturePlot(sc_epi, features = c('S.Score', 'G2M.Score'), reduction = "harmonyumap")
# DimPlot(sc_epi, group.by = c('Phase'), reduction = "harmonyumap")
# 
# 
# reduction.use <- 'harmony' # 'harmony' 'pca'
# 
# ## cluster ------
# sc_epi <-
#   FindNeighbors(
#     sc_epi,
#     reduction = reduction.use,
#     dims = 1:epi_dims,
#     k.param = 10
#   )
# 
# sc_epi <-
#   FindClusters(
#     sc_epi,
#     resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
#   )
# 
# # dimplot
# pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
# # DimPlot(sc_epi, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
# # DimPlot(sc_epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
# # DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# # DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# # DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# # DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# 
# 
# DimPlot(sc_epi, group.by = "Spatial.008um_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()
# 
# # Dotplot -------
# # epi cluster
# marker_epi <- list(
#   "AT1"     = c("AGER","EMP2","CAV1","RTKN2","MYL9","HOPX","PDPN","MSLN","CLIC5"),
#   # "AT1-AT2" = c("EGFR","SHH","MUC1","PDPN"),
#   "AT2"     = c("MUC1","SFTPC","NPC2","SFTPA1","NAPSA","SFTPA2","CTSH","PGC","SFTPB","SFTPD","ETV5","ABCA3"),
#   basal     = c("KRT5","KRT17","MIR205HG","DLK2","BCL11A","KRT14","TP63","DAPL1"),
#   'ciliated cells' = c("CCDC113","FOXJ1","DNAH3","DNAIC2","DNAAF1","LRRC6"),
#   'club cells'     = c("SCGB1A1","SCGB3A2","CYP2F2","CCKAR")
# )
# table(sc_epi$Spatial.008um_snn_res.0.4)
# pdf(
#   file.path(epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - RNA_snn_res.0.4.pdf")), 
#   height = 6,
#   width = 12
# )
#   # epi marker
#   DotPlot(
#     sc_epi,
#     features = c(
#       'epithelial cells'  = 'EPCAM',
#       'epithelial cells'  = 'KRT5',
#       'epithelial cells'  = 'KRT14',
#       'tumor' = 'CEACAM6',
#       'tumor' = 'CEACAM5'
#     ),
#     group.by = "Spatial.008um_snn_res.0.4"
#   ) + 
#     RotatedAxis()
#   # lung epi marker
#   DotPlot(
#     sc_epi,
#     features = marker_epi,
#     group.by = "Spatial.008um_snn_res.0.4"
#   ) +
#     RotatedAxis()
# dev.off()

# saveRDS(sc_epi, file.path(epi_out_dir, "sc_epi.rds"))
sc_epi <- readRDS('sp_subcluster/sp_epi.rds')
sc_epi_anno_meta <-
  sc_epi %>%
  # .[['RNA_snn_res.0.6']] %>% # run2
  .[['RNA_snn_res.0.6']] %>% 
  rownames_to_column("cell_id") 


all_data_anno_meta <- 
  all_data[[c('celltype', 'total_counts')]] %>%
  dplyr::filter(total_counts > 100) %>%
  rownames_to_column('cell_id') %>%
  dplyr::filter(celltype %in% c('Hepatocyte', 'Tumor cell', "Kupffer cell", "B/Plasma cell")) %>%
  left_join(sc_epi_anno_meta, by = "cell_id") %>% 
  map_dfc(as.character) %>%
  
  mutate(annotation2 = ifelse(
    celltype %in% c('Hepatocyte', 'Tumor cell'),
    RNA_snn_res.0.6,
    celltype
  )) %>% 
  dplyr::select(
    cell_id,
    annotation = annotation2
  )

# all_data <- all_data %>% AddMetaData(all_data_anno_meta)  
set.seed(1234)
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.1) %>%
  
  # dplyr::filter(!annotation %in% c('4', '5', '6')) # run 2
  dplyr::filter(!annotation %in% c('13', '14', '15'))

# inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹

# infercnv_out_dir = "infercnv - run1" # cutoff 0.001
# infercnv_out_dir = "infercnv - run2" # cutoff 0.01
infercnv_out_dir = "infercnv - run3" # new cluster, group2 remove batch, res0.6
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

# undebug(CreateInfercnvObject)
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = all_data_raw@assays$RNA$counts[VariableFeatures(all_data_raw), meta_data_infer_subset$cell_id] ,
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path('infercnv', "gene_pos47.txt"),
  
  min_max_counts_per_cell = c(30, Inf),
  
  delim             = "\t",
  ref_group_names   = c("Kupffer cell", "B/Plasma cell")
)



# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.01) } # 87 <--- 
infercnv_obj@count.data %>% rowMeans %>% { sum(. < 0.005) } # 34


infercnv_obj2  = infercnv::run(
  infercnv_obj,
  cutoff              =    0.01,
  out_dir    =    infercnv_out_dir,
  cluster_by_groups   =    TRUE,
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  png_res             =    300,
  # useRaster   = FALSE,
  output_format = 'pdf',
  no_plot = T #  figure margin too large ----
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
#   no_plot = T #  figure margin too large ----
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


# cnv score ------------
library(scales)
cnvScore <- function(data){
  # sd_num <- apply(data, 1, sd)
  # data <- data[sd_num > median(sd_num), ]
  
  sum_num <- apply(data, 1, sum)
  data <- data[sum_num > median(sum_num) + 2 * mad(sum_num), ]
  
  # rng <- quantile(data, c(0.01, 0.99))
  # data <- pmax( pmin(data, rng[2]), rng[1] )
  
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



# spatial plot --------------
# Idents(all_data) <- 'Spatial.008um'
# Idents(all_data) <- 'Spatial.qc.pass'

sc_epi_infer_anno <-
  sc_epi[['Spatial.008um_snn_res.0.4']] %>% 
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
  left_join(sc_epi_infer_anno, by = "cell_id") %>%
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
  file.path(epi_out_dir, "spatial plot of infercnv analysis.pdf"),
  height = 8,
  width  = 8
)
sp_p_ls %>% walk(print)
dev.off()



