setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/subcelltype/Epithelial/monocle')
rm(list = ls())
gc()

library(cowplot)
library(Seurat)
library(RColorBrewer)

# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)
#use_condaenv('myname')

sceasy::convertFormat('../celltype/anno.h5ad', from="anndata", to="seurat",outFile='Epi_anno.rds')

obj <- readRDS(file = './Epi_anno.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type','leiden'))
# cellchat.list = SplitObject(obj,split.by = "group")
# saveRDS(file = './Myeloid_Pre.rds',cellchat.list[["pre"]])
# saveRDS(file = './Myeloid_Post.rds',cellchat.list[["post"]])

############################################reorder##########################################################
library(monocle)
library(Seurat)
library(ggplot2)
library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)

setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/subcelltype/Epithelial/monocle')


HSMM = readRDS("/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/subcelltype/Epithelial/monocle/monocle.rds")
deg <- read.csv("/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/subcelltype/Epithelial/monocle/monocle.DEG.xls",sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 7)
clustercolname = 'cell_type'
colorc = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
            '#e377c2')
names(colorc) = c('Epi_KRT13','Epi_KRT16','Epi_LTF','Epi_MKI67','Epi_SLIT2','Epi_TFF3','Myoepithelial')
# 重新画图
getwd()

p = plot_cell_trajectory(HSMM_myo,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)

p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)+
  scale_color_manual(values = colorc)
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)
ncols = ceiling(sqrt(length(unique(HSMM_myo@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(HSMM_myo@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)+
  facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)+
  scale_color_manual(values = colorc)
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)

p = plot_cell_trajectory(HSMM_myo,color_by=clustercolname,size=ptsize,show_backbone=TRUE)+
  scale_color_manual(values = colorc)
ggsave(paste0("monocle.trajectory.",clustercolname,".png"),p,width = 7,height = 7.5)
ggsave(paste0("monocle.trajectory.",clustercolname,".pdf"),p,width = 7,height = 7.5)

p = plot_cell_trajectory(HSMM_myo,color_by=clustercolname,size=ptsize,show_backbone=TRUE)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)+
  scale_color_manual(values = colorc)
ggsave(paste0("monocle.trajectory.",clustercolname,".faceted.png"),p,width = 4*ncols,height = 4.2*nrows)
ggsave(paste0("monocle.trajectory.",clustercolname,".faceted.pdf"),p,width = 4*ncols,height = 4.2*nrows)


pseudotime_df = HSMM_myo@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"3.monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

saveRDS(file = '3.obj_reorder.rds',HSMM_myo,compress = F)

# 逆时序山脊图
library(ggplot2)
library(ggridges)
df = HSMM_myo@phenoData@data %>% dplyr::select(cell_type,Pseudotime,group) %>% dplyr::filter(group!='NL')
p = ggplot(df, aes(x = Pseudotime, y = cell_type,fill = cell_type)) +
  geom_density_ridges(rel_min_height = 0.005,alpha = 0.6)+
  scale_fill_manual(values = colorc)+
  theme_ridges()+ theme(plot.background = element_rect('white'))
ggsave(filename = "3.monocle.pseudotime_celltype_group.png",p,h=7,w=10,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_celltype_group.pdf",p,h=7,w=10,dpi = 150)
p = ggplot(df, aes(x = Pseudotime, y = cell_type,fill = group)) +
  geom_density_ridges(rel_min_height = 0.005,alpha = 0.6)+
  # scale_y_discrete(expand = c(0.01, 0)) +
  # scale_x_continuous(expand = c(0.01, 0))+
  theme_ridges()+ theme(plot.background = element_rect('white'))
ggsave(filename = "3.monocle.pseudotime_celltype_group2.png",p,h=7,w=10,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_celltype_group2.pdf",p,h=7,w=10,dpi = 150)

# p = plot_cell_trajectory(HSMM_myo,color_by=clustercolname,size=1,show_backbone=TRUE)+
#   facet_grid(cell_type~group)

# 拟时相关的差异gene,ordergene做回归
diff_test_res <- differentialGeneTest(HSMM_myo[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >50)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = '3.obj_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 7,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "monocle.top50gene_cortime.png",p,h=8,w=6,dpi = 150)
ggsave(filename = "monocle.top50gene_cortime.pdf",p,h=8,w=6,dpi = 150)

setwd('../')
