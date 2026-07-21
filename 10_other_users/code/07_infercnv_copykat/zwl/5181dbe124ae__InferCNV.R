# inferCNV用于探索肿瘤单细胞RNA-Seq数据，以确定体细胞是否存在大规模染色体拷贝数变异（Copy number variation，CNV），
# 例如整个染色体或大段染色体的增加或者减少（gains or deletions）。
# InferCNV实现的方法是，通过用一组“正常”细胞做参照，探索肿瘤基因组不同位置基因的表达强度来完成的。
# InferCNV生成的热图显示每条染色体上基因的相对表达强度，与正常细胞相比，能显而易见的看到肿瘤基因组的哪些区域过多或过少表达。
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# sudo apt install jags
# BiocManager::install("infercnv")
rm(list = ls())
gc()
library(Seurat)
library(infercnv)

args <- commandArgs(T);
options(stringsAsFactors = FALSE)
inFile = args[1];                                                #传入文件接收原始counts以及10x的matrix
CluterFile = args[2];                                            #细胞类型定义第一列CB第二列定义
GeneLocationFile = args[3];                                      #gene地址文件
# A1BG	19	58345178	58353492
# A1BG-AS1	19	58347718	58355455
# A1CF	10	50799409	50885675
# A2M	12	9067664	9116229
refGroup = args[4];                                              #ref对应的组
cutoff = as.numeric(args[5]);                                    #过滤值，df0.1

outpath =  args[6];
mincell = 3
genelist = args[8];                                             #

# 手动传参
outpath = "/lvdata/wzb/scRNA/FW2023-453/2023.10.07"
seuset = readRDS("h5_trans.rds")

allcounts = GetAssay(seuset,assay = "RNA")@counts
mincell = 10
# write.table(seuset@meta.data,"meta.txt",sep = "\t",quote = F)
# Idents(seuset) = seuset$cell_type
# data = data.frame(cellname = names(Idents(seuset)),celltype = Idents(seuset))
# unique(data[,2])
refGroup = c("B cells","DC","Endothelial cells","Macrophages","Plasma cells","T cells","Fibroblast")
GeneLocationFile = "GRch38_2020A.txt"
cutoff = 0.1
CluterFile = "Clsuter.txt"
AddTag =TRUE


print("Start")
setwd(outpath)
gc()
if(dir.exists(inFile)){
  allcounts = Read10X(data.dir=inFile)
}else{
  allcounts = read.table(inFile,sep="\t",header=T,row.names=1)
  title = read.table(inFile, he=F, nrow = 1,colClasses='character',sep="\t", row.names=1)
  colnames(allcounts) = title
}
# 细胞类型定义
data = read.table(CluterFile,sep="\t",header = T,stringsAsFactors = FALSE,check.names = F)
rownames(data) = data[,1]

# 取交集
cell = intersect(rownames(data),colnames(allcounts))
allcounts = allcounts[,cell]
data = data[cell,]

CluterFile = "Cluster.txt"
if(AddTag){
  data[,2] = paste0("Cluster",data[,2])
  refgroups = unlist(strsplit(refGroup, "[,]"))
  refgroups = paste0("Cluster",refgroups)
}else{
  refgroups = unlist(strsplit(refGroup, "[,]"))
}
cellNum = table(data[,2])
print(cellNum)

filterCluster = names(cellNum[cellNum<mincell])
finalCluster = names(cellNum[cellNum>=mincell])
refgroups = refgroups[refgroups %in% finalCluster]

print(paste0("ref group:",paste(refgroups,collapse = ",")))
if(length(filterCluster)>0){
  print(paste0("remove cluster:",paste(filterCluster,collapse = ",")))
  filter = data[,2] %in% finalCluster
  data = data[filter,]
}

write.table(data,file="Cluster.txt",sep="\t",quote=F,row.names = F,col.names = F)

if(length(args)>=8){
  gene = unique(read.table(genelist,sep="\t",header=T,stringsAsFactors = FALSE)[,1])
  gene = intersect(gene,rownames(allcounts))
  allcounts = allcounts[gene,]
}
GeneLocationFile = read.csv(GeneLocationFile,header = F,row.names = 1,sep = "\t")
CluterFile =read.csv(CluterFile,header = F,row.names = 1,sep = "\t")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=allcounts,
                                    annotations_file=CluterFile,
                                    gene_order_file=GeneLocationFile,
                                    ref_group_names=refgroups)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir=outpath, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=FALSE,
                             num_threads=64,
                             write_expr_matrix = T,
                             analysis_mode="samples"
                             )
# 绘制图片
plot_cnv(infercnv_obj, k_obs_groups = 1, cluster_by_groups = T,
         out_dir = outpath, color_safe_pal = FALSE, x.center = 1,
         title = "inferCNV",output_format = "pdf",
         obs_title = "Observations (Cells)", ref_title = "References (Cells)",
         output_filename = "infercnv")


##### use CNV data to predict cancer cells and classify clusters ----->>>>>
cnv_result= read.csv( file = "infercnv.observations.txt", stringsAsFactors=FALSE ,sep=" ")
# this "cnv_result_2" is the CNV data including the cells chosen as reference (non-malignant cells)
cnv_result_2=read.csv(file= "infercnv.references.txt",stringsAsFactors=FALSE ,sep=" ")

patients<-gsub("_.*$", "", colnames(cnv_result))
# 这里的1是均值
cnv_result = cnv_result-1
cnv_result_2 = cnv_result_2-1
cnv_result_combine = cbind(cnv_result,cnv_result_2)                      #合并恶性跟非恶性结果
patients_combine = gsub("_.*$", "", colnames(cnv_result_combine))
patients_combine = rep(1,length(patients))
patients = rep(1,length(patients))
cnv_meanSquare=apply(cnv_result_combine^2, 2, mean)
cnv_corr=c()
for(i in c(sort(unique(patients_combine))))
{
  temp_matrix_cancer=cnv_result[,patients==i]
  temp_matrix_combine=cnv_result_combine[,patients_combine==i]
  #wq: error, cnv_result_combine列数大于patients length, 那么combine的后半部分(reference)取出来的列是错的，循环patients vector。
  temp_average=apply(temp_matrix_cancer,1,mean)
  #wq: 癌细胞每个基因的均值
  temp_corr=apply(temp_matrix_combine,2,function(x) cor(x, temp_average, method = "pearson"))
  #wq: 该样本每个细胞和癌细胞均值的corr
  temp_corr[is.na(temp_corr)]=0
  cnv_corr=c(cnv_corr,temp_corr)
}
table(names(cnv_corr)==names(cnv_meanSquare))
# the name is not in the same order
# change the order of "cnv_corr" as that of names of "cnv_meanSquare"
cnv_corr = cnv_corr[order(match(names(cnv_corr),names(cnv_meanSquare)))]

mergeOne=cbind(cnv_corr,cnv_meanSquare)
write.table(mergeOne,"CNV_cells.txt",quote=F,sep="\t")
pdf("correlation_against_meanSquare.pdf",width=10,height=10,paper='special')
smoothScatter(cnv_corr, cnv_meanSquare, xlab="correlation", ylab="meanSquare", nrpoints=1000, cex=0.5)
abline(h=0.05,  lty=2,col="dodgerblue",bty="n")
abline(v=0.5,  lty=2, col="darkgreen",bty="n")
legend("left", legend=c("SNV signal 0.05","correlation with tumor cells 0.5"), col=c( "dodgerblue", "darkgreen"), lty=2, cex=1.5, bty="n")
dev.off()



group = read.delim("Cluster.txt", stringsAsFactors=F,header=F)

celllist=gsub("-", ".", group[,1])
cnv_meanSquare = cnv_meanSquare[order(match(names(cnv_meanSquare),celllist))]
cnv_result_combine = cnv_result_combine[order(match(colnames(cnv_result_combine),celllist))]

meandata = aggregate(cnv_meanSquare, list(group[,2]), FUN=mean)
cluster_mean = aggregate(t(as.matrix(cnv_result_combine)), list(group[,2]), FUN=mean)
cancer_average=apply(cnv_result,1,mean)
cluster_corr = apply(cluster_mean[,-1],1,function(x) cor(x, cancer_average, method = "pearson"))
cnv_cluster_result=cbind(meandata,cluster_corr)
# Cluster之间的相关系数
write.table(cnv_cluster_result,"CNV_clusters.txt",quote=F,sep="\t")
library(ggplot2)
gg=ggplot(cnv_cluster_result,aes(x=cnv_cluster_result[,2],y=cnv_cluster_result[,3]))+
  geom_point(shape=21,size=3)+
  geom_text(aes(y = cnv_cluster_result[,3] + .015, label = cnv_cluster_result[,1]),size = 2)+
  geom_hline(yintercept = 0.5)+geom_vline(xintercept = 0.05)

ggsave("CNVPlot.png",gg,width= 10,height= 10)


library(ggplot2)
library(dplyr)
library(stringr)
library(getopt)


corr = 0.5
meanSquare = 0.05

# 读取细胞信息，因为inferCNV的结果处理时细胞名有经过转换，这边的细胞也同样转换以匹配
# CellUmiInfo表格读取
# 检查CellUmiinfo表格是否有CellName和orig.ident两列，没有就报错
cellumiinfo = seuset@meta.data
cellumiinfo$CellName = row.names(cellumiinfo)
if (any(!(c("CellName", "sample") %in% colnames(cellumiinfo)))) {
  rm(cellumiinfo)
  stop("CellUmiInfo表格中缺少CellName或orig.ident列！")
}else {
  cellumiinfo <- cellumiinfo %>% 
    mutate(Cell_new = make.names(CellName)) %>% 
    select(Cell_new, sample)
}


# SummaryCell表格读取
summarycell <- seuset@meta.data
summarycell$Cell = row.names(summarycell)
if (any(!(c("Cell", "cell_type") %in% colnames(summarycell)))) {
  rm(summarycell)
  stop("Summarycell表格中缺少Cell或celltype列！")
} else {
  summarycell <- summarycell %>% 
    mutate(Cell_new = make.names(Cell)) %>% 
    select(Cell_new, cell_type)
}
cnv_result = read.csv("CNV_cells.txt",header = T,sep = "\t",row.names = 1)
cnv_result <- data.frame(Cell = rownames(cnv_result), cnv_result)
colnames(cnv_result) <- c("Cell", "cnv_corr", "cnv_meanSquare")
dat <- inner_join(cnv_result, cellumiinfo, by = c("Cell" = "Cell_new"))%>% 
  inner_join(summarycell, by = c("Cell" = "Cell_new"))%>% 
  mutate(TumorType = ifelse(cnv_corr > corr & cnv_meanSquare > meanSquare, "malignant", "non-malignant"))
write.table(dat, file = "AllSample_CNV_result.txt", sep = "\t", row.names = F,
            quote = F)

sample_num <- length(unique(dat$sample))
col_num <- round(sqrt(sample_num))
# violin（mean_square做）
my_plot1 <- ggplot(dat, aes(x = factor(cell_type), y = cnv_meanSquare, fill = factor(cell_type))) + 
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  geom_hline(yintercept = meanSquare, linetype = "dotted") + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "", y = "CNV Mean Square\n", fill = "cell_type")


my_plot2 <- my_plot1 + facet_wrap(vars(sample), ncol = col_num)


ggsave(my_plot1, filename = "ViolinPlot_All.png", width = 7, height = 7)
ggsave(my_plot1, filename = "ViolinPlot_All.pdf", width = 7, height = 7)
ggsave(my_plot2, filename = "ViolinPlot_EachSample.pdf", width = 7, height = 7)
ggsave(my_plot2, filename = "ViolinPlot_EachSample.png", width = 7, height = 7)

# --------------
# CNV mean square和cnv corr组成的散点图
ymax <- max(dat$cnv_corr)
xmax <- max(dat$cnv_meanSquare)

my_plot3 <- ggplot(dat, aes(x = cnv_meanSquare, y = cnv_corr, color = factor(cell_type))) + 
  geom_rect(aes(xmin = meanSquare, xmax = xmax, ymin = corr, ymax = ymax), fill = "#e6e5e6", color = NA) + 
  geom_point(alpha = .7, size = 0.6)  + 
  geom_hline(yintercept = corr, linetype = "dotted") + 
  geom_vline(xintercept = meanSquare, linetype = "dotted") +
  theme_classic() + 
  labs(x = "CNV Mean Square\n", y = "CNV Correlation", color = "cell_type")

my_plot4 <- my_plot3  + facet_wrap(vars(sample), ncol = col_num)

ggsave(my_plot3, filename = "ScatterPlot_All.png", width = 7, height = 6)
ggsave(my_plot3, filename = "ScatterPlot_All.pdf", width = 7, height = 6)
ggsave(my_plot4, filename = "ScatterPlot_EachSample.png", width = 8, height = 7)
ggsave(my_plot4, filename = "ScatterPlot_EachSample.pdf", width = 8, height = 7)

# ----------------
# UMAP/TSNE图
dat <- arrange(dat, desc(TumorType))


coord = seuset@reductions$umap@cell.embeddings
if(!is.null(coord)){
  # coord表格读取
  coord <- cbind(Cell = row.names(coord),coord)
    
  dat = merge(dat, coord)
  
  
  
  my_plot5 <- ggplot(dat, aes(x = dat[, 7], y = dat[, 8], color = TumorType)) + 
    geom_point(alpha = .7, size = 0.6) + 
    scale_color_manual(values = c("red", "grey"), limits = c("malignant", "non-malignant")) + 
    theme_classic() + labs(x = colnames(coord)[2], y = colnames(coord)[3])
  
  my_plot6 <- my_plot5 + facet_wrap(vars(sample), ncol = col_num)
  
  ggsave(my_plot5, filename = "coord_All.png", width = 7, height = 6)
  ggsave(my_plot5, filename = "coord_All.pdf", width = 7, height = 6)
  ggsave(my_plot6, filename = "coord_EachSample.png", width = 7, height = 6)
  ggsave(my_plot6, filename = "coord_EachSample.pdf", width = 7, height = 6)
  
}

