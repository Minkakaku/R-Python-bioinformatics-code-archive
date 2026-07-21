rm(list = ls())
gc()
#.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# 绘图用包
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(getopt)
options(stringsAsFactors = FALSE,warn = -1)
# 可视化用到的文件：############################################################
#           seurat的umap坐标和metadata
#           auc.loom文件
# grn和ctx文件用不到
command = matrix(c(
  # 参数名称   缩写  需要的值个数 数据类型   # 参数注释
  'outpath','o',1,"character",    #输出地址
  'objumap','u',1,"character",    #obj的umap位置
  'objmeta','m',1,"character",    #obj的metadata
  'celltype','c',1,"character",   #注释后的细胞类型
  'loominput','l',2,"character"   #loom结果输入
),byrow=TRUE, ncol=4)
args=getopt(command)

#tmpdir = "/mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC"
# args$outpath = tmpdir
# args$loominput = paste0(tmpdir,"/output/aue_result.loom")
# args$objmeta = paste0(tmpdir,"/SCENIC/meta.xls")
# args$celltype = 'cell_type'
# args$objumap = paste0(tmpdir,"SCENIC/umap.xls")

#创建文件夹，读取得数据###########################################################
setwd(args$outpath)
if (!dir.exists('AUCplot')) {
  dir.create('AUCplot')}
if (!dir.exists('Scatterplot')) {
  dir.create('Scatterplot')}
if (!dir.exists('Heatmap')) {
  dir.create('Heatmap')}
packageVersion("SCENIC")

# rm(list = ls())
umap = read.csv(args$objumap,row.names = 1)
meta = read.csv(args$objmeta,row.names = 1)
meta[, args$celltype ] = gsub('/','_',meta [, args$celltype ])

loom = open_loom(args$loominput)
  # Read information from loom file:
  exprMat = get_dgem(loom)
  exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
  regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
  regulons = regulonsToGeneLists(regulons_incidMat)                     #调控子和对应的基因向量
  regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
  regulonAucThresholds = get_regulon_thresholds(loom)                   #用于二值化 regulonAUC，判断细胞是否激活某调控子
  embeddings = get_embeddings(loom)
close_loom(loom)


#每个cell type AUC热图 ###########################################################
# 划分cell_type
regulonAUC
cellsPerCluster1  =  split ( rownames ( meta ),  meta [, args$celltype ])  
cellsPerCluster <- cellsPerCluster1[sapply(cellsPerCluster1, length) >= 2]  #每个细胞类型至少有两个细胞

regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),] 
# 计算平均表达式: 
regulonActivity_byCellType  =  sapply ( cellsPerCluster , function ( cells )  rowMeans (getAUC ( regulonAUC )[, cells ])) 
write.csv(file = './Heatmap/Regulon Activity by celltype.csv',regulonActivity_byCellType)
# 尺度表达式：
regulonActivity_byCellType_Scaled  =  t ( scale ( t ( regulonActivity_byCellType ),  center  =  T ,  scale = T ))
regulonActivity_byCellType_Scaled[is.na(regulonActivity_byCellType_Scaled)] = 0
hm = Heatmap(regulonActivity_byCellType_Scaled ,  name = "Regulon Activity", cluster_rows = TRUE, show_row_names = TRUE,
             row_names_gp = grid :: gpar ( fontsize = 6 ))  # 行字体大小
pdf("Heatmap/Regulon Activity.pdf",width = 20,height = 30)
draw(hm)
dev.off()
png("Heatmap/Regulon Activity.png",width = 2000,height = 3000)
draw(hm)
dev.off()
# 获取行聚类后的基因顺序
# row_ordered_genes <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)]# 获取行顺序（索引）
# row_dend <- row_dend(hm)  # 获取行聚类树
# clusters <- cutree(as.hclust(row_dend), k = 2)  # k=2 表示分为2类
# write.csv(file = './Heatmap/Regulon Activity pdf genenames.csv',regulonActivity_byCellType)


#每个cell type AUC热图 ###########################################################
# 为了识别簇特异性调节子，调节子特异性评分（RSS）特别有用（由 Suo 等人在 2018 年为小鼠细胞图谱提出） ）。
library(plotly)
rss  <-  calcRSS ( AUC = getAUC ( regulonAUC ),
                   cellAnnotation = meta [ colnames ( regulonAUC ),  args$celltype ])
rss[is.na(rss)] = 0
rssPlot <- plotRSS(rss,cluster_columns = T)

gg = rssPlot$plot
gg = gg + 
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        axis.text.x=element_text(angle=45,size=8,hjust = 1))
ggsave("./Scatterplot/RSS.bubble.pdf",gg,height = 30)
ggsave("./Scatterplot/RSS.bubble.png",gg,height = 30)

# 散点图绘制强的转录因子
for (i in unique(colnames(rss))) {
  gg = plotRSS_oneSet(rss, setName = i)  
  ggsave(paste0("./Scatterplot/Cluster_",i,"_Scatterplot.pdf"),gg)
  ggsave(paste0("./Scatterplot/Cluster_",i,"_Scatterplot.png"),gg)
}

# 绘制转录因子的umap图
setwd('./AUCplot')
AUCell_plotTSNE(as.matrix(umap),cellsAUC = regulonAUC,plots="AUC",asPNG = T)
setwd(args$outpath)
# 活性阈值判断是否激活###########################################################
binarizeAUC <- function(auc, thresholds){
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}
revnum <- function(numthe){
  val = as.numeric(names(numthe))
  name = as.character(numthe)
  numthe = val
  names(numthe) = name
  return(numthe)
}
binaryRegulonActivity = binarizeAUC(regulonAUC,revnum(regulonAucThresholds))


# 绘制热图
setwd('../Heatmap')
gg = Heatmap(binaryRegulonActivity, name="Binarized activity", 
             col = c("white", "black"),
             cluster_rows = TRUE, cluster_columns=TRUE,
             show_column_names = F,
             row_names_gp=grid::gpar(fontsize=3),use_raster = TRUE) # row font size
pdf("Binarizedheatmap.pdf")
gg
dev.off()
png("Binarizedheatmap.png")
gg
dev.off()
