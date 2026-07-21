library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(corrplot)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/all_data_sc_doublet.rds")
# adata$group<-factor(adata$group)
# adata$sample<-factor(adata$sample)

# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = ",",head=TRUE)
geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = ",",head=FALSE)
#转换成list
gene<-as.list(geneset)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[21]<-'select_geneset_Score'




genedata<-FetchData(adata,var=c('select_geneset_Score',gene$V1))
data <- cor(genedata, method = "pearson")
write.table(data,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/cor/pearson.xls",sep="\t",quote=F,row.names=T,col.names=T)

testcor <- cor.mtest(genedata, method = "pearson",conf.level=0.95)
# testcor <- testcor$p['select_geneset_Score',]
write.table(testcor,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/cor/cortest.xls",sep="\t",quote=F,row.names=T,col.names=T)
genedata<-FetchData(adata,var=c('select_geneset_Score',gene$V1,'group'))
x_col <- 'select_geneset_Score'
y_cols <- gene$V1
# y_cols <- setdiff(gene$V1, c(x_col,'CXCR7',' ENSG00000306549','ENSG00000297914'))
# y_cols[8:10]<-c("HLA_DRA","HLA_DQA1","HLA_DPA1")
# colnames(genedata)[9:11]<-c("HLA_DRA","HLA_DQA1","HLA_DPA1")
for (y_col in y_cols) {
  # 绘制散点图
  # p <- ggplot(genedata, aes_string(x = x_col, y = y_col)) +
  #   geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
  #   geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
  #   labs(title = paste(x_col, "vs", y_col, "Cor"),
  #        x = x_col,
  #        y = y_col) +
  #   theme_bw() +  # 使用简洁主题
  #   theme(plot.title = element_text(hjust = 0.5),  # 居中标题
  #         panel.grid.major = element_blank(),  # 添加网格线
  #         panel.grid.minor = element_blank())  # 去掉次要网格线
  
  # # 打印图形
  # ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/","scatter_", y_col, ".png", sep = ""), plot = p, width = 6, height = 4)
  # ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/","scatter_", y_col, ".pdf", sep = ""), plot = p, width = 6, height = 4)
subdata<-genedata[,c("select_geneset_Score",y_col,"group")]
colnames(subdata)[2]<-'y'
cor <- subdata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(select_geneset_Score, y, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(subdata[which(subdata$group=="LAHS"),],aes(x=select_geneset_Score,y=y,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "LAHS",
         x = "select_geneset_Score",
         y = y_col) +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(subdata[which(subdata$group=="NLAHS"),],aes(x=select_geneset_Score,y=y,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "NLAHS",
         x = "select_geneset_Score",
         y = y_col) +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(subdata[which(subdata$group=="NC"),],aes(x=select_geneset_Score,y=y,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "NC",
         x = "select_geneset_Score",
         y = y_col) +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/cor/","scatter_", y_col, ".png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/cor/","scatter_", y_col, ".pdf", sep = ""), plot = a, width = 18, height = 6)
}
genedata<-FetchData(adata,var=c('group','cell_type',gene$V1))
write.table(genedata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/cor/genedata.xls",sep="\t",quote=F,row.names=T,col.names=T)


Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix="Main_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.png",width = 15,height = 10)
# data<-FetchData(adata,vars =c("annotation","select_geneset_Score"))
# p<-ggplot(data,aes(cell_type,select_geneset_Score))
# p+geom_boxplot()+theme_bw()+RotatedAxis()

# mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score"))
# a<-ggplot(mydata,aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score","group"))
a1<-ggplot(mydata[which(mydata$group=="LAHS"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="LAHS")
a2<-ggplot(mydata[which(mydata$group=="NLAHS"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NLAHS")
a3<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NC")
a<-a1+a2+a3
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main.UMAP.pdf",width = 18,height = 6)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/AddModuleScore/Main.UMAP.png",width = 18,height = 6)