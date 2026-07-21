setwd("/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/")
rm(list = ls())
gc()

library(sceasy)
library(Seurat)
library(reticulate)
#############颜色#################
library(scales)
#pal= ggsci::pal_d3('category20')(16)
#names(pal) <- 
show_col(pal)

#############h5ad转变成seurat对象####################

sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/B_raw.h5ad', from="anndata", to="seurat",
                      outFile='B_raw.rds')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/T_raw.h5ad', from="anndata", to="seurat",
                      outFile='T_raw.rds')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/Myeloid_raw.h5ad', from="anndata", to="seurat",
                      outFile='Myeloid_raw.rds')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/SMC_raw.h5ad', from="anndata", to="seurat",
                      outFile='SMC_raw.rds')
for (ci in c('B','T','Myeloid','SMC')) {
  obj = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/',ci,'_raw.rds'))
  EXP = subset(obj,subset = (group == 'EXP'))
  saveRDS(file = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require8/',ci,'_EXP.rds'),EXP)
  
}
################1 构建对象#################
main_fun <-function(obj,outdir){
  setwd(outdir)
  library(dplyr)
  library(monocle)
  f_data = data.frame(gene_short_name = rownames(obj),
                      row.names = rownames(obj))
  count = as.matrix(LayerData(obj, assay="RNA", layer='counts'))
  
  pd <- new("AnnotatedDataFrame", data = obj@meta.data)
  fd <- new("AnnotatedDataFrame", data = f_data)
  HSMM <- newCellDataSet(count, phenoData = pd, featureData = fd,
                         lowerDetectionLimit = 0.5,expressionFamily=negbinomial.size())
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(HSMM),
                                      num_cells_expressed >= 10))
  
  ################2挑选celltype差异gene作为特征#####################
  ##2000左右##
  diff = differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~cell_type")
  ordering_genes <- row.names(subset(diff, qval < 0.01))
  write.csv(file = 'order_gene.csv',ordering_genes)
  
  HSMM <- setOrderingFilter(HSMM, ordering_genes[1:2000])
  plot_ordering_genes(HSMM)
  
  
  ##降维##
  HSMM <- reduceDimension(HSMM, max_components = 2,
                          reduction_method = 'DDRTree',
                          residualModelFormulaStr = "~sample")
  save(file = paste0(outdir,'_reduce2.RData'),HSMM,ordering_genes)
  ##排序##
  HSMM <- orderCells(HSMM)
  ##重新排序##
  p = plot_cell_trajectory(HSMM, color_by = "State")+
    facet_wrap(~cell_type,ncol = 4)
  ggsave(filename = "Re-order_trajectory_cell_type.png",p,w=12,h=8,dpi = 100)
  setwd("/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/")
  return(HSMM)
}


plot_fun <- function(HSMM,outdir,root_s=1){
  setwd(outdir)
  ordering_genes = read.csv(file = './order_gene.csv',row.names = 1)
  HSMM_myo <- orderCells(HSMM, root_state = root_s)
  saveRDS(file = paste0(outdir,'_order.rds'),HSMM_myo,compress = F)
  
  #################################可视化###################################
  p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")
  ggsave(filename = paste0("2.",outdir,"_trajectory_cell_type.png"),p,w=6,h=6,dpi = 100)
  
  p1 = plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
  ggsave(filename = paste0("2.",outdir,"_trajectory_Pseudotime.png"),p1,w=4,h=4,dpi = 100)
  
  p0 = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+facet_wrap(~group)
  ggsave(filename = paste0("2.",outdir,"_trajectory_cell_type_group.png"),p0,w=6,h=5,dpi = 100)
  
  p = plot_cell_trajectory(HSMM_myo, color_by = "State")
  ggsave(filename = paste0("2.",outdir,"_trajectory_state.png"),p,w=6,h=5,dpi = 100)
  ##cell_type split##
  #HC_bar = colnames(HSMM_myo)[HSMM_myo$group=='NC']
  p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+facet_grid(group~cell_type)
  ggsave(filename = paste0("2.",outdir,"_trajectory_cell_type_group2.png"),p,w=15,h=5,dpi = 100)
  
  
  ##########################4 拟时相关的差异gene############################
  ##将ordergene提出来做回归，找与时间是否显著相关##
  ##否则会计算全部gene##
  diff_test_res <- differentialGeneTest(HSMM_myo[ordering_genes[1:2000,],],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_gene_names <- diff_test_res %>% dplyr::arrange(qval) %>% head(50) %>% row.names() 
  write.csv(file = paste0(outdir,'_Pseudotime_ordergene.csv'),diff_test_res)
  p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 5,
                              cores = 1,
                              show_rownames = T,return_heatmap = TRUE)
  ggsave(filename = paste0("2.",outdir,"_top50gene_cortime.png"),p,h=8,w=6,dpi = 150)
  setwd("/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/")
  return(HSMM_myo)
}

obj <- readRDS(file = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/B_raw.rds'))
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = 'cell_type')
B_HSMM = main_fun(obj = obj,outdir = 'B')
B_HSMM1 = plot_fun(HSMM = B_HSMM,outdir = 'B',root_s = 5)

rm(obj,B_HSMM,B_HSMM1)
gc()
obj <- readRDS(file = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/T_raw.rds'))
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = 'cell_type')
T_HSMM = main_fun(obj = obj,outdir = 'T')
T_HSMM1 = plot_fun(HSMM = T_HSMM,outdir = 'T',root_s = 1) #没确定

rm(obj,T_HSMM,T_HSMM1)
gc()
obj <- readRDS(file = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/SMC_raw.rds'))
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = 'cell_type')
SMC_HSMM = main_fun(obj = obj,outdir = 'SMC')
SMC_HSMM1 = plot_fun(HSMM = SMC_HSMM,outdir = 'SMC',root_s = 1,ordering_genes = ordering_genes)

rm(obj,SMC_HSMM,SMC_HSMM1)
gc()
obj <- readRDS(file = paste0('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/monocle/Myeloid_raw.rds'))
obj <- FindVariableFeatures(obj,nfeatures = 2000)
obj$cell_type = obj$leiden
Seurat::DimPlot(obj,reduction = 'umap',group.by = 'cell_type')
Mye_HSMM = main_fun(obj = obj,outdir = 'Myeloid')
Mye_HSMM1 = plot_fun(HSMM = Mye_HSMM,outdir = 'Myeloid',root_s = 1) #没确定
