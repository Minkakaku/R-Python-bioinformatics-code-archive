setwd('/mnt/NexData_file/home/zwl/Project/SC2510105423-01/require3')
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

# sceasy::convertFormat('../1177/celltype/anno.h5ad', from="anndata", to="seurat",outFile='1177_anno.rds')
# sceasy::convertFormat(obj, from="seurat", to="anndata",outFile='./RDS/T.h5ad')
obj <- readRDS(file = './anno_raw.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type','leiden'))

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


HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2510105423-01/require3/tdtomato_EC/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2510105423-01/require3/tdtomato_EC/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)


# 重新画图
setwd('tdtomato_EC/')
getwd()
clustercolname = 'cell_type'

HSMM_myo@phenoData@data[["tdTomato_exp"]] = factor(HSMM_myo@phenoData@data[["tdTomato_exp"]],levels = c("tdTomato+","other"))
p = plot_cell_trajectory(HSMM_myo,color_by="tdTomato_exp",size=1)
ggsave("monocle_tdtomato_pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle_tdtomato_pseudotime.png",p,width = 7,height = 7.5)

p = plot_genes_in_pseudotime(HSMM_myo['tdTomato'], color_by = "cell_type")
ggsave("0.monocle_tdtomato_pseudotime.pdf",p,width = 7,height = 4)
ggsave("0.monocle_tdtomato_pseudotime.png",p,width = 7,height = 4)

p = plot_cell_trajectory(HSMM_myo,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)

p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)
ncols = ceiling(sqrt(length(unique(HSMM_myo@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(HSMM_myo@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)

pseudotime_df = HSMM_myo@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"3.monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

saveRDS(file = '3.obj_reorder.rds',HSMM_myo,compress = F)

# 逆时序山脊图
library(ggplot2)
library(ggridges)
df = HSMM_myo@phenoData@data %>% dplyr::select(cell_type,Pseudotime,group) 
p = ggplot(df, aes(x = Pseudotime, y = cell_type,fill =  group)) +
  geom_density_ridges(rel_min_height = 0.005,alpha = 0.6)+
  theme_ridges()+ theme(plot.background = element_rect('white'))
ggsave(filename = "3.monocle.pseudotime_celltype_group.png",p,h=7,w=10,dpi = 150)
ggsave(filename = "3.monocle.pseudotime_celltype_group.pdf",p,h=7,w=10,dpi = 150)


# p = plot_cell_trajectory(HSMM_myo,color_by=clustercolname,size=1,show_backbone=TRUE)+
#   facet_grid(cell_type~group)

# 拟时相关的差异gene,ordergene做回归
diff_test_res <- differentialGeneTest(HSMM_myo[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = '3.obj_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "monocle.top50gene_cortime.png",p,h=8,w=6,dpi = 150)
ggsave(filename = "monocle.top50gene_cortime.pdf",p,h=8,w=6,dpi = 150)

setwd('../')

####################################Veen########################################
gene1177 = read.csv('./1177/3.obj_Pseudotime_ordergene.csv',header = T) %>% dplyr::filter(qval <0.05)
gene895 = read.csv('./895/3.obj_Pseudotime_ordergene.csv')%>% dplyr::filter(qval <0.05)
gene943 = read.csv('./943/3.obj_Pseudotime_ordergene.csv')%>% dplyr::filter(qval <0.05)

commen = intersect(intersect(gene1177$X, gene895$X), gene943$X)
write.table(file = './monocle_sample3_commen_gene.txt',commen,row.names = F)
library(ggvenn)
A = list(E1177 = gene1177$X, E895 = gene895$X, E943 = gene943$X)
p = ggvenn(A,fill_color= c("#E41A1C", "#1E90FF", "#FF8C00"),fill_alpha = .7,#stroke_linetype = "longdash",
           set_name_size = 8,text_size=5)+theme(plot.background = element_rect('white')) 
ggsave("./monocle_sample3_commen_gene_Veen.png", width = 6, height = 6, dpi = 300,p)
ggsave("./monocle_sample3_commen_gene_Veen.pdf", width = 6, height = 6, dpi = 300,p)

## commen gene exp in data #####
sig_gene_names = read.csv('./monocle_sample3_commen_gene.txt')
m1177 =  readRDS('./1177/3.obj_reorder.rds')
m895 = readRDS('./895/3.obj_reorder.rds')
m943 = readRDS('./943/3.obj_reorder.rds')
p = plot_pseudotime_heatmap(m1177[sig_gene_names$x,],
                            num_clusters = 2,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
p1 = plot_pseudotime_heatmap(m895[sig_gene_names$x,],
                             num_clusters = 2,
                             cores = 1,
                             show_rownames = T,return_heatmap = TRUE)
p2 = plot_pseudotime_heatmap(m943[sig_gene_names$x,],
                             num_clusters = 2,
                             cores = 1,
                             show_rownames = T,return_heatmap = TRUE)
library(cowplot)

# 组合三个热图
combined_plot <- plot_grid(
  p$gtable,  # 提取 gtable 对象
  p1$gtable,
  p2$gtable,
  labels = c("m1177", "m895", "m943"),  # 添加标签
  label_size = 12,
  nrow = 1,  # 并排显示
  align = "h",  # 水平对齐
  axis = "tb"   # 共享坐标轴
)+theme(plot.background = element_rect('white'))
ggsave(filename = './1.pseudotime_commen_gene_exp.png',combined_plot,w=10,h=8)
ggsave(filename = './1.pseudotime_commen_gene_exp.pdf',combined_plot,w=10,h=8)
#####  get commen gene ###########
g1177 <- cutree(p[["tree_row"]], k=2)
g895 <- cutree(p1[["tree_row"]], k=2)
g943 <- cutree(p2[["tree_row"]], k=2)
# 获取基因名
genes <- names(g1177)

# 拆分基因并求交集
common1 <- data.frame(gene = intersect(intersect(genes[g1177==1], genes[g895==1]), genes[g943==1])) #上调
common2 <- data.frame(gene = intersect(intersect(genes[g1177==2], genes[g895==2]), genes[g943==2])) #下调
write.csv(file = './1.pseudotime_upgene_commen.csv',common1,row.names = F)
write.csv(file = './1.pseudotime_downgene_commen.csv',common2,row.names = F)

######################################## CNV in tra ############################
cnv = read.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2507233889-01/895/celltype/CNV/metadata.csv',row.names = 1)
HSMM_myo = readRDS('./895/3.obj_reorder.rds')
HSMM_myo@phenoData@data$cnv = cnv[rownames(HSMM_myo@phenoData@data),'cnv_score_cluster']#
p = plot_cell_trajectory(HSMM_myo,x=1,y=2,color_by = 'cnv',cell_size = .5)+
  facet_wrap(~cell_type)+scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),
                                               values = c(0,0.5,0.8,1))
ggsave(filename = "0.895_monocle_pseudotime_celltype_infercnv.png",p,h=7,w=10,dpi = 150)

################################################################################
# 寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图。
Time_diff <- differentialGeneTest(HSMM_myo[ordergene,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <-Time_diff[, c(5, 3, 4)] # 把gene放前面，也可以不改
Time_genes <- Time_diff %>% filter(pval < 0.05) %>% arrange(qval) %>% distinct(gene_short_name, .keep_all = TRUE)
write.csv(Time_genes, paste0(outdir,"/Pseudotime_diff_all.csv"), row.names = F)
Time_genes <- Time_genes %>% head(10) %>% pull(gene_short_name) %>% as.character() # 提取Top10gene
# 直接查看一些基因随细胞状态等的表达变化,选择前10个top基因并将其对象取出
seuCDS_subset <- HSMM_myo[Time_genes,]
# 可视化：以state/cluster/pseudotime进行
gg=plot_genes_in_pseudotime(seuCDS_subset, color_by = "group", ncol = 5)
ggsave(paste0(outdir,"/pseudotime_genes_group.png"),limitsize = FALSE,width= 20,height= 10)
ggsave(paste0(outdir,"/pseudotime_genes_group.pdf"),limitsize = FALSE,width= 20,height= 10)
