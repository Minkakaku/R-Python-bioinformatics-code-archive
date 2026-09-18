# conda activate scMetabolism
library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/Main_srt_after.rds")
levels(adata)
DefaultAssay(adata)<-"RNA"
genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/script/pathway.csv",sep = "\t",head=FALSE)
genesets$V2 <- strsplit(genesets$V2, ",")
gene <- setNames(genesets$V2, genesets$V1)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[18]<-'HALLMARK_WNT_BETA_CATENIN_SIGNALING'
colnames(adata@meta.data)[19]<-'HALLMARK_TGF_BETA_SIGNALING'
colnames(adata@meta.data)[20]<-'HALLMARK_PI3K_AKT_MTOR_SIGNALING'
colnames(adata@meta.data)[21]<-'Hgf_Met'

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/"

mydata<-FetchData(adata,vars=c('Hgf','HALLMARK_WNT_BETA_CATENIN_SIGNALING','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_PI3K_AKT_MTOR_SIGNALING',"cell_type","group"))
xlim <- range(mydata$Hgf)
ylim <- range(mydata$HALLMARK_WNT_BETA_CATENIN_SIGNALING)
# xlim <- c(0,4.15)
# ylim <- c(-0.125,0.175)
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(Hgf, HALLMARK_WNT_BETA_CATENIN_SIGNALING, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="Control"),],aes(x=Hgf,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Control",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="ELA"),],aes(x=Hgf,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "ELA",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="HGF"),],aes(x=Hgf,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "HGF",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_group_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_group_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

cor <- mydata %>%
  group_by(cell_type) %>%  # 按照'group'列分组
  summarise(correlation = cor(Hgf, HALLMARK_WNT_BETA_CATENIN_SIGNALING, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = Hgf, y = HALLMARK_WNT_BETA_CATENIN_SIGNALING, colour = cell_type)) +
        geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        stat_cor(aes(label = sprintf("r == %.2f", ..r..),group = cell_type),        # 关键：按 cell_type 分组
           formula = y~x,
           output.type = "expression",
           label.x = 3.5,
					 parse = TRUE,
           size = 4) +
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "cell_type",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_cell_type_cor.png", sep = ""), plot = a1, width = 10, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_cell_type_cor.pdf", sep = ""), plot = a1, width = 10, height = 7)

cor3=round(cor$correlation[3],2)
cor4=round(cor$correlation[4],2)
a3<-ggplot(mydata[which(mydata$cell_type=="Fibroblast"),],aes(x=Hgf,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor4,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Fibroblast",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a4<-ggplot(mydata[which(mydata$cell_type=="SMCs"),],aes(x=Hgf,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor4,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "SMCs",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a3+a4
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_cell_type-fact_cor.png", sep = ""), plot = a, width = 16, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_cell_type-fact_cor.pdf", sep = ""), plot = a, width = 16, height = 7)

mydata_sub<-mydata[mydata$cell_type=="SMCs",]
a1 <- ggplot(mydata_sub,
             aes(x = Hgf, y = HALLMARK_WNT_BETA_CATENIN_SIGNALING, colour = group)) +
        geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        stat_cor(aes(label = sprintf("r == %.2f", ..r..),group = group),        # 关键：按 cell_type 分组
           formula = y~x,
           output.type = "expression",
           label.x = 3.5,
					 parse = TRUE,
           size = 4) +
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "SMCs",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_SMCs_group_cor.png", sep = ""), plot = a1, width = 8, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_SMCs_group_cor.pdf", sep = ""), plot = a1, width = 8, height = 7)

mydata_sub<-mydata[mydata$cell_type=="Fibroblast",]
a1 <- ggplot(mydata_sub,
             aes(x = Hgf, y = HALLMARK_WNT_BETA_CATENIN_SIGNALING, colour = group)) +
        geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        stat_cor(aes(label = sprintf("r == %.2f", ..r..),group = group),        # 关键：按 cell_type 分组
           formula = y~x,
           output.type = "expression",
           label.x = 3.5,
					 parse = TRUE,
           size = 4) +
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Fibroblast",
         x = "Hgf",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_Fibroblast_group_cor.png", sep = ""), plot = a1, width = 8, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_Fibroblast_group_cor.pdf", sep = ""), plot = a1, width = 8, height = 7)

cor_table1<-cor(mydata[mydata$cell_type=="Epithelials",1:4])
cor_table2<-cor(mydata[mydata$cell_type=="SMCs",1:4])
cor_table3<-cor(mydata[mydata$cell_type=="Endothelial",1:4])
cor_table4<-cor(mydata[mydata$cell_type=="Fibroblast",1:4])
data_table <- data.frame(cor_table1[-1,1],cor_table2[-1,1],cor_table3[-1,1],cor_table4[-1,1])
colnames(data_table)<-c("Epithelials","SMCs","Endothelial","Fibroblast")
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/Main-all.xls",sep="\t",row.names=TRUE)
library(pheatmap)
anno_row <- data.frame(colnames(data_table))
rownames(anno_row)<-colnames(data_table)
colnames(anno_row)<-'correlation'
ann_colors=list(group=c(Epithelials='#1f77b4',SMCs='#ff7f0e',Endothelial='#279e68',Fibroblast='#FFCC33'))
p2 <- pheatmap(t(data_table),  
              #  color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
              #  border_color = "black", 
							#  border=F, 
              #  scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = FALSE, 
               legend = TRUE, 
               legend_labels = c("low","","heigh"), 
               show_rownames = FALSE, 
               show_colnames = TRUE, 
               fontsize = 8, 
               annotation_row = anno_row,   #添加列注释信息
               annotation_legend = TRUE,    #是否显示注释信息图例
               annotation_names_row = FALSE,  #是否显示列注释的名称
               annotation_colors=ann_colors
							 )
pdf("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/Main-group-Heatmap.pdf",width = 5,height = 5)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC20251028/07.Mission/cor_AddModuleScore/Main-group-Heatmap.png",width = 5*60,height = 5*60)
p2
dev.off()