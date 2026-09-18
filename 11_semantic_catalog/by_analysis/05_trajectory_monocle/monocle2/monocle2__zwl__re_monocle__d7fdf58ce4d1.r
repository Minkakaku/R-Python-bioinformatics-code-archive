setwd('/mnt/NexData_file/home/zwl/Project/SC2508194437/require2/monocle2/')
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
obj <- readRDS(file = './raw_anno.rds')
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


HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2508194437/require2/monocle2/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2508194437/require2/monocle2/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)


# 重新画图
setwd(sample[i])
getwd()
clustercolname = 'cell_type'

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
p = ggplot(df, aes(x = Pseudotime, y = cell_type,fill =  cell_type)) +
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

#################################clontype12####################################
BCR = read.csv(file = '../BCR_result/Result.xls',sep = '\t',row.names = 1)
BCR <- BCR %>%
  filter(airr.cc_aa_alignment == '12') %>%  # 筛选12的行
  mutate(
    barcode = str_split_fixed(rownames(.), '_', 2)[,1]  # 提取barcode（行名切分）
  )

# 2. 提取phenoData并合并
df1 <- HSMM_myo@phenoData@data %>%
  mutate(barcode = rownames(.))  # 行名转barcode列

# 3. 合并+填充NA+恢复行名（一步到位）
df2 <- merge(df1, BCR[, c('barcode','airr.cc_aa_alignment','VDJ_1_c_call')], 
             by = 'barcode', all.x = TRUE) %>%
  mutate(top_clone = ifelse(is.na(airr.cc_aa_alignment), 'other', as.character(airr.cc_aa_alignment))) %>%  # NA→other，非NA转字符
  `rownames<-`(.$barcode) 
HSMM_myo@phenoData@data = df2
p = plot_cell_trajectory(HSMM_myo,color_by = 'top_clone',cell_size = .5)#+
  # facet_wrap(~cell_type)+
  #scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'),values = c(0,0.5,0.8,1))
ggsave(filename = "0.clonetype12_monocle_pseudotime_celltype.pdf",p,h=5,w=7,dpi =200)


#############饼图
BCR = read.csv(file = '../BCR_result/Result.xls',sep = '\t',row.names = 1)

VDJC = unique(BCR$VDJ_1_c_call)
VDJC = VDJC[VDJC != ""]
colorA = brewer.pal(length(VDJC), 'Paired')
names(colorA)=VDJC

BCR <- BCR %>%
  filter(airr.cc_aa_alignment == '12')
## 计算每个类别的计数
counts <- table(BCR$VDJ_1_c_call)
total <- nrow(BCR)
percentage <- counts / total * 100
data <- as.data.frame(percentage)
colnames(data) <- c("category", "percentage")
data = data %>%filter(category != "")
print(data)
p <- ggplot(data,aes(x = "",y = percentage,fill = category)) + 
  geom_bar(stat="identity", width=1, color="white") + #绘制柱状图
  coord_polar("y", start=0) + #设置环形
  theme_void()+
  geom_text( aes (label = paste0(round(percentage, 1), " % ")), position = position_stack(vjust= 0.5 )) +
  scale_fill_manual(values = colorA)
ggsave(paste0("./2.Isotype12_pie_ALL.pdf"), p, width = 5, height = 3,dpi=200)
ggsave(paste0("./2.Isotype12_pie_ALL.png"), p, width = 5, height = 3,dpi=200)

for (si in c('OMV','PCM')) {
  BCR = read.csv(file = '../BCR_result/Result.xls',sep = '\t',row.names = 1)
  BCR <- BCR %>%
    filter(airr.cc_aa_alignment == '12')%>%
    filter(gex.group == si)
  ## 计算每个类别的计数
  counts <- table(BCR$VDJ_1_c_call)
  total <- nrow(BCR)
  percentage <- counts / total * 100
  data <- as.data.frame(percentage)
  colnames(data) <- c("category", "percentage")
  data = data %>%filter(category != "")
  print(data)
  p <- ggplot(data,aes(x = "",y = percentage,fill = category)) + 
    geom_bar(stat="identity", width=1, color="white") + #绘制柱状图
    coord_polar("y", start=0) + #设置环形
    theme_void()+
    geom_text( aes (label = paste0(round(percentage, 1), " % ")), position = position_stack(vjust= 0.5 )) +
    scale_fill_manual(values = colorA)
  ggsave(paste0("./2.Isotype12_pie_",si,".pdf"), p, width = 5, height = 3,dpi=200)
  ggsave(paste0("./2.Isotype12_pie_",si,".png"), p, width = 5, height = 3,dpi=200)
}
