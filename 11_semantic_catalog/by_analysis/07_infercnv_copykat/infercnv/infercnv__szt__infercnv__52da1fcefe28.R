# installed into linux server
# mamba install conda-forge::jags
# install.packages('rjags')
# BiocManager::install("infercnv")


# BiocManager::install("infercnv")
# https://sourceforge.net/projects/mcmc-jags/files
# If you are using R 4.1.3 or earlier then install JAGS-4.3.0.exe


# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb




# 1.所有上皮细胞和down sample部分免疫细胞，跑inferCNV分析，
# 识别出正常肝细胞（肝母细胞/胆管上皮细胞？），恶性上皮细胞。
# 恶性上皮细胞跑cNMF分析，识别出上皮细胞程序。

library(infercnv)

###  epi subcluster --------------
epi_out_dir <- "epi_subcluster"
dir.create(epi_out_dir)


sc_epi <- all_data %>% JoinLayers() %>% .[, .$annotation %in% c("Epithelial")]

# preprocess

sc_epi <- 
  sc_epi %>% 
  NormalizeData(scale.factor = median(all_data$nCount_RNA))%>%
  FindVariableFeatures(nfeatures = 1500) %>%
  ScaleData(
    vars.to.regress = c(
      'nCount_RNA'
    )
  )


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
epi_dims <-  10


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
    resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# Dotplot -------
table(sc_epi$RNA_snn_res.0.2)
pdf(
  file.path(epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - RNA_snn_res.0.1.pdf")), 
  height = 6,
  width = 8
)
DotPlot(
  sc_epi,
  features = c(
    'epithelial cells'  = 'EPCAM',
    'epithelial cells'  = 'KRT5',
    'epithelial cells'  = 'KRT14',
    
    Hepatocyte          = 'ALB',
    Hepatocyte          = 'APOB',
    Hepatocyte          = 'CYP3A5',
    Hepatocyte          = 'CPB2'
  ),
  group.by = "RNA_snn_res.0.2"
)
dev.off()

# saveRDS(sc_epi, file.path(epi_out_dir, "sc_epi.rds"))

# Dotplot -------
table(sc_epi$Spatial.008um_snn_res.0.4)
pdf(
  file.path(epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - RNA_snn_res.0.4.pdf")), 
  height = 6,
  width = 8
)
DotPlot(
  sc_epi,
  features = c(
    'epithelial cells'  = 'EPCAM',
    'epithelial cells'  = 'KRT5',
    'epithelial cells'  = 'KRT14'
  ),
  group.by = "Spatial.008um_snn_res.0.4"
)
dev.off()
# annotation ----
Idents(sc_epi) <- "RNA_snn_res.0.2"


cell_cluster <- c(
  '0' = "tumor",
  '1' = "Epithelial",
  '2' = "tumor",
  '3' = "tumor",
  '4' = "tumor",
  '5' = "tumor",
  '6' = "Epithelial",
  '7' = "Epithelial",
  '8' = "tumor",
  '9' = "tumor",
  '10' = "tumor",
  '11' = "Epithelial"
)

sc_epi <- RenameIdents(sc_epi, cell_cluster)
sc_epi$subtype <- Idents(sc_epi)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_epi, group.by = "subtype", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# saveRDS(sc_epi, file.path(epi_out_dir, "sc_epi.rds"))
sc_epi_anno_meta <-
  sc_epi[['RNA_snn_res.0.2']] %>% 
  rownames_to_column("cell_id") 


all_data_anno_meta <- 
  all_data[['annotation']] %>%
  rownames_to_column('cell_id') %>%
  dplyr::filter(annotation %in% c("B cells","Epithelial","Macrophage","T cells")) %>%
  left_join(sc_epi_anno_meta, by = "cell_id") %>% 
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(
    annotation %in% c("Epithelial"),
    RNA_snn_res.0.2,
    annotation
  )) %>% 
  dplyr::select(cell_id, annotation = annotation2)

# all_data <- all_data %>% AddMetaData(all_data_anno_meta)  
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.5)



### inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹

infercnv_out_dir = "inferCNV_episub_1"
dir.create(infercnv_out_dir) 


# meta_data_infer <-
#   all_data[[]] %>%
#   rownames_to_column("cell_id") %>%
#   dplyr::filter(
#     annotation %in% c("B cells","Epithelial","Endothelial","Macrophage","T cells")
#   ) %>%
#   dplyr::select(cell_id, annotation, RNA_snn_res.0.4) %>%
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
meta_data_infer_subset <-
  meta_data_infer_subset %>%
  mutate(RNA_snn_res.0.4 = str_c('Cluster', RNA_snn_res.0.4, sep = "_")) %>%
  mutate(annotation2 = ifelse(
    annotation %in% c("B cells","Endothelial","Macrophage","T cells"),
    annotation,
    RNA_snn_res.0.4
  )) %>%
  dplyr::select(
    cell_id,
    annotation = annotation2
  )

# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )

undebug(CreateInfercnvObject)
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path('inferCNV', "gene_pos47.txt"),
  delim             = "\t",
  # min_max_counts_per_cell = c(50,Inf),
  ref_group_names   = c("B cells","Macrophage","T cells")
)


# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.4) } # 87 <--- 
infercnv_obj@count.data %>% rowMeans %>% { sum(. < 0.005) } # 34



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
infercnv_out_dir <- "inferCNV_episub"
infercnv_obj2  = infercnv::run(
  infercnv_obj,
  cutoff              =    0.1,
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
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.01, 
#                               out_dir= "infercnv_output_2",  
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               hclust_method="ward.D2", plot_steps=F,no_plot = T)
infercnv::plot_cnv(infercnv_obj2, 
                   out_dir='infercnv_output_2',
                   output_filename='infercnv', 
                   x.range="auto",
                   x.center=1,
                   title = "infercnv", 
                   output_format="pdf",
                   color_safe_pal = FALSE)
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, 
                              out_dir= "infercnv_output_sd",  
                              write_expr_matrix=T,
                              cluster_by_groups=T ,  
                              output_format="pdf",
                              sd_amplifier=3,
                              denoise=T,
                              noise_logistic=TRUE,
                              hclust_method="ward.D2", plot_steps=F)
saveRDS(infercnv_obj2, file.path(infercnv_out_dir, 'infercnv_obj2.rds'))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, 'infercnv_obj.rds'))


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
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score)<-'CNV_score'
cnv_score %>% 
  rownames_to_column('cell_id') %>% 
  left_join(.,as.data.frame(meta_data_infer_subset),'cell_id') ->temp

temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+geom_violin(aes(fill= annotation),color="NA")+
  ggsci::scale_fill_igv()+
  # scale_fill_manual(values = color_v)+
  theme_bw()->plot
ggsave(
  plot,
  filename = file.path(infercnv_out_dir, "CNV_score.pdf"),
  height = 6,
  width  = 8
)


infercnv::plot_cnv(infercnv_obj2, 
                   out_dir='infercnv_output',
                   output_filename='infercnv_1', 
                   x.range="auto",
                   x.center=1,
                   title = "infercnv", 
                   output_format="pdf",
                   sd_amplifier=3,
                   color_safe_pal = FALSE)

#可视化图------------
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")


expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$normal
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- test_loc$test

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("EXAMPLE_ONLY_DONT_REUSE.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

#聚类，7类，提取结果
set.seed(20241230)
kmeans.result <- kmeans(t(expr), 6)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)

#定义热图的注释，及配色
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #类别数
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

#下面是绘图
pdf("try1.pdf",width = 15,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#每一类对应的CB保存在kmeans_df_s数据框中
write.table(kmeans_df_s, file = "kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)


expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
kmeans_df_s$CB=rownames(kmeans_df_s)
CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")

CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()
ggsave("CNV level.pdf",width = 10,height = 6,units = "cm")