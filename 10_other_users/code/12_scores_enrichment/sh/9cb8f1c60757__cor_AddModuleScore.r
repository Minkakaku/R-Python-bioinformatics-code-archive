# conda activate scMetabolism
library(Seurat)
library(corrplot)
library(ggplot2)
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mapping/Main_mapping.rds")
levels(adata)

gene<-list(
  dsRNA  = c("IFIH1","RIGI","EIF2AK2","TLR3"),
  dsDNA = c("CGAS","STING1"),
  Signal_Transduction_Transcription = c("MAVS","TBK1","IKBKE","IRF3","IRF7"),
  Product_IFN = "IFNB1",
  ISG = c("ISG15","MX1","MX2","IFIT1","IFIT3","OAS1","RSAD2","USP18"),
  Interferon_receptor = c("IFNAR1","IFNAR2")
)

adata<-subset(adata,cell_type %in% c("F13A1 Macrophages","HLA-A Macrophages","cDC","pDC","Myofibroblasts"))
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16:21]<-c('dsRNA','dsDNA','Signal_Transduction_Transcription','Product_IFN','ISG','Interferon_receptor')

outdir <- "/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/cor_AddModuleScore/"

mydata<-FetchData(adata,vars=c('LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR','dsRNA','dsDNA','Signal_Transduction_Transcription','Product_IFN','ISG','Interferon_receptor',"cell_type","group"))

combos <- unique(mydata[, c('cell_type','group')])
combos<-combos[-4,]
setwd(outdir)
for (i in seq_len(nrow(combos))) {
  ct  <- combos$cell_type[i]
  grp <- combos$group[i]
  # 取出子集
  sub <- mydata %>% filter(cell_type == ct & group == grp)
  # 计算 Spearman 相关
  cor.mat <- cor(select(sub, -cell_type, -group),
                 method = 'spearman', use = 'pairwise.complete.obs')
  # 3. 画图 ----------------------------------------------------------------
  png(paste0(ct, '_', grp, '_corr.png'), width = 800, height = 700, res = 120)
  corrplot(cor.mat, method = 'color', type = 'upper',col = colorRampPalette(c('red','white','blue'))(100)[100:1],
           tl.cex = 0.8, tl.col = 'black', number.cex = 0.7,
           addCoef.col = 'grey50', diag = FALSE)
  dev.off()
	pdf(paste0(ct, '_', grp, '_corr.pdf'), width = 8, height = 7) 
	corrplot(cor.mat, method = 'color', type = 'upper',col = colorRampPalette(c('red','white','blue'))(100)[100:1],
           tl.cex = 0.8, tl.col = 'black', number.cex = 0.7,
           addCoef.col = 'grey50', diag = FALSE)
  dev.off()
}



#### 后面没改
xlim <- range(mydata$THBS1)
ylim <- range(mydata$HALLMARK_GLYCOLYSIS)
xlim <- c(0,4.15)
ylim <- c(-0.125,0.175)
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, HALLMARK_GLYCOLYSIS, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="LAHS"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "LAHS",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="NLAHS"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NLAHS",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NC",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/mono_group_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/mono_group_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

cor <- mydata %>%
  group_by(cell_type) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, HALLMARK_GLYCOLYSIS, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = THBS1, y = HALLMARK_GLYCOLYSIS, colour = cell_type)) +
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
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/mono_cell_type_cor.png", sep = ""), plot = a1, width = 10, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/mono_cell_type_cor.pdf", sep = ""), plot = a1, width = 10, height = 7)

cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
cor4=round(cor$correlation[4],2)
cor5=round(cor$correlation[5],2)
cor6=round(cor$correlation[6],2)
a1<-ggplot(mydata[which(mydata$cell_type=="Non-classical monocyte"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Non-classical monocyte",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$cell_type=="PTTG1 classical monocyte"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "PTTG1 classical monocyte",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$cell_type=="ZFP36L1 classical monocyte"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "ZFP36L1 classical monocyte",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a4<-ggplot(mydata[which(mydata$cell_type=="S100A8 classical monocyte"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor4,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "S100A8 classical monocyte",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a5<-ggplot(mydata[which(mydata$cell_type=="Classical monocyte"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor5,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Classical monocyte",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a6<-ggplot(mydata[which(mydata$cell_type=="Dendritic cell"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=cell_type)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor6,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Dendritic cell",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a <- (a1 | a2 | a3) /
     (a4 | a5 | a6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/mono_cell_type-fact_cor.png", sep = ""), plot = a, width = 18, height = 14)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/mono_cell_type-fact_cor.pdf", sep = ""), plot = a, width = 18, height = 14)

mydata_sub<-mydata[mydata$cell_type=="PTTG1 classical monocyte",]
cor <- mydata_sub %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, HALLMARK_GLYCOLYSIS, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata_sub[which(mydata_sub$group=="LAHS"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "LAHS",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata_sub[which(mydata_sub$group=="NLAHS"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NLAHS",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata_sub[which(mydata_sub$group=="NC"),],aes(x=THBS1,y=HALLMARK_GLYCOLYSIS,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NC",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/PTTG1_group-fact_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/PTTG1_group-fact_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

a1 <- ggplot(mydata_sub,
             aes(x = THBS1, y = HALLMARK_GLYCOLYSIS, colour = group)) +
        geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        stat_cor(aes(label = sprintf("r == %.2f", ..r..),group = group),        # 关键：按 cell_type 分组
           formula = y~x,
           output.type = "expression",
           label.x = 3.5,
					 parse = TRUE,
           size = 4) +
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "PTTG1 classical monocyte",
         x = "THBS1",
         y = "HALLMARK_GLYCOLYSIS") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/PTTG1_group_cor.png", sep = ""), plot = a1, width = 10, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_AddModuleScore/PTTG1_group_cor.pdf", sep = ""), plot = a1, width = 10, height = 7)



marker_file <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = "\t",head=FALSE)
marker <- marker_file$V1
#转换成list
gene<-list(marker)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'select_geneset'

outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/cor_select_geneset/"

mydata<-FetchData(adata,vars=c('THBS1','select_geneset',"cell_type","group"))
xlim <- range(mydata$THBS1)
ylim <- range(mydata$select_geneset)
xlim <- c(0,4.1)
ylim <- c(-0.25,1)
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, select_geneset, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="LAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "LAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="NLAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NLAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NC",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/mono_group_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/mono_group_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

cor <- mydata %>%
  group_by(cell_type) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, select_geneset, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = THBS1, y = select_geneset, colour = cell_type)) +
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
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/mono_cell_type_cor.png", sep = ""), plot = a1, width = 10, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/mono_cell_type_cor.pdf", sep = ""), plot = a1, width = 10, height = 7)

mydata_sub<-mydata[mydata$cell_type=="PTTG1 classical monocyte",]
cor <- mydata_sub %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, select_geneset, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata_sub[which(mydata_sub$group=="LAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "LAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata_sub[which(mydata_sub$group=="NLAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NLAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata_sub[which(mydata_sub$group=="NC"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NC",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/PTTG1_group-fact_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/PTTG1_group-fact_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

mydata<-FetchData(adata,vars=c('THBS1',marker,"cell_type","group"))
cor_table<-cor(mydata[,1:80])
data_table <- cor_table[,1]
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/Monocyte-all.xls",sep="\t",row.names=TRUE)

cor_table1<-cor(mydata[mydata$group=="LAHS",1:80])
cor_table2<-cor(mydata[mydata$group=="NLAHS",1:80])
cor_table3<-cor(mydata[mydata$group=="NC",1:80])
data_table <- data.frame(cor_table1[-1,1],cor_table2[-1,1],cor_table3[-1,1])
colnames(data_table)<-c("LAHS","NLAHS","NC")
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/Monocyte-all.xls",sep="\t",row.names=TRUE)
library(pheatmap)
anno_row <- data.frame(colnames(data_table))
rownames(anno_row)<-colnames(data_table)
colnames(anno_row)<-'correlation'
ann_colors=list(group=c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
p2 <- pheatmap(t(data_table),  
              #  color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
              #  border_color = "black", 
							#  border=F, 
              #  scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = TRUE, 
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
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/Monocyte-group-Heatmap.pdf",width = 20,height = 4)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/Monocyte-group-Heatmap.png",width = 20*60,height = 4*60)
p2
dev.off()

mydata_sub<-mydata[mydata$cell_type=="PTTG1 classical monocyte",]
cor_table1<-cor(mydata_sub[mydata_sub$group=="LAHS",1:80])
cor_table2<-cor(mydata_sub[mydata_sub$group=="NLAHS",1:80])
data_table <- data.frame(cor_table1[-1,1],cor_table2[-1,1])
colnames(data_table)<-c("LAHS","NLAHS")
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/PTTG1-all.xls",sep="\t",row.names=TRUE)

p2 <- pheatmap(t(data_table),  
              #  color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
              #  border_color = "black", 
							#  border=F, 
              #  scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = TRUE, 
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
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/PTTG1-group-Heatmap.pdf",width = 20,height = 3)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset/PTTG1-group-Heatmap.png",width = 20*60,height = 3*60)
p2
dev.off()


#### 分为>0 和 =0
marker_file <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = "\t",head=FALSE)
marker <- marker_file$V1
#转换成list
gene<-list(marker)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'select_geneset'

outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/cor_select_geneset_new/"

mydata<-FetchData(adata,vars=c('THBS1','select_geneset',"cell_type","group"))
mydata<-mydata[mydata$THBS1!=0,]
xlim <- range(mydata$THBS1)
ylim <- range(mydata$select_geneset)
xlim <- c(0,4.1)
ylim <- c(-0.1,0.2)
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, select_geneset, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="LAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "LAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="NLAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NLAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NC",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/mono_group_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/mono_group_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

cor <- mydata %>%
  group_by(cell_type) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, select_geneset, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = THBS1, y = select_geneset, colour = cell_type)) +
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
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/mono_cell_type_cor.png", sep = ""), plot = a1, width = 10, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/mono_cell_type_cor.pdf", sep = ""), plot = a1, width = 10, height = 7)

mydata_sub<-mydata[mydata$cell_type=="PTTG1 classical monocyte",]
cor <- mydata_sub %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(THBS1, select_geneset, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata_sub[which(mydata_sub$group=="LAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "LAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata_sub[which(mydata_sub$group=="NLAHS"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NLAHS",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata_sub[which(mydata_sub$group=="NC"),],aes(x=THBS1,y=select_geneset,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "NC",
         x = "THBS1",
         y = "select_geneset") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/PTTG1_group-fact_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/PTTG1_group-fact_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

mydata<-FetchData(adata,vars=c('THBS1',marker,"cell_type","group"))
cor_table<-cor(mydata[,1:80])
data_table <- cor_table[,1]
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/Monocyte-all.xls",sep="\t",row.names=TRUE)

cor_table1<-cor(mydata[mydata$group=="LAHS",1:80])
cor_table2<-cor(mydata[mydata$group=="NLAHS",1:80])
cor_table3<-cor(mydata[mydata$group=="NC",1:80])
data_table <- data.frame(cor_table1[-1,1],cor_table2[-1,1],cor_table3[-1,1])
colnames(data_table)<-c("LAHS","NLAHS","NC")
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/Monocyte-all.xls",sep="\t",row.names=TRUE)
library(pheatmap)
anno_row <- data.frame(colnames(data_table))
rownames(anno_row)<-colnames(data_table)
colnames(anno_row)<-'correlation'
ann_colors=list(group=c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
data_table<-na.omit(data_table)
p2 <- pheatmap(t(data_table),  
              #  color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
              #  border_color = "black", 
							#  border=F, 
              #  scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = TRUE, 
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
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/Monocyte-group-Heatmap.pdf",width = 20,height = 4)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/Monocyte-group-Heatmap.png",width = 20*60,height = 4*60)
p2
dev.off()

mydata_sub<-mydata[mydata$cell_type=="PTTG1 classical monocyte",]
cor_table1<-cor(mydata_sub[mydata_sub$group=="LAHS",1:80])
cor_table2<-cor(mydata_sub[mydata_sub$group=="NLAHS",1:80])
data_table <- data.frame(cor_table1[-1,1],cor_table2[-1,1])
colnames(data_table)<-c("LAHS","NLAHS")
write.table(data_table,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/PTTG1-all.xls",sep="\t",row.names=TRUE)
data_table<-na.omit(data_table)
p2 <- pheatmap(t(data_table),  
              #  color = colorRampPalette(c('#330066','#336699','#66CC66','#FFCC33'))(100), 
              #  border_color = "black", 
							#  border=F, 
              #  scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = TRUE, 
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
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/PTTG1-group-Heatmap.pdf",width = 20,height = 3)
p2
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cor_select_geneset_new/PTTG1-group-Heatmap.png",width = 20*60,height = 3*60)
p2
dev.off()



#### 分为>0 和 =0
marker_file1 <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = "\t",head=FALSE)
marker_file2 <- read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
marker1 <- marker_file1$V1
marker2 <- marker_file2$HALLMARK_GLYCOLYSIS
#转换成list
gene<-list(marker1,marker2)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'select_geneset'
colnames(adata@meta.data)[15]<-'HALLMARK_GLYCOLYSIS'

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/select_geneset_barplot/"

# mydata<-FetchData(adata,vars=c('THBS1','select_geneset','HALLMARK_GLYCOLYSIS',"cell_type","group"))
# mydata<-mydata[mydata$THBS1!=0,]
# xlim <- range(mydata$THBS1)
p = VlnPlot(object = adata, features = 'THBS1',pt.size = 0)+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
ggsave(paste0(outdir,'mono_cell_type-THBS1_VlnPlot.pdf'),width = 7,height = 7)
ggsave(paste0(outdir,'mono_cell_type-THBS1_VlnPlot.png'),width = 7,height = 7)
p = VlnPlot(object = adata, features = 'THBS1',pt.size = 0,group.by='thbs1_group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
ggsave(paste0(outdir,'mono-THBS1_VlnPlot.pdf'),width = 7,height = 7)
ggsave(paste0(outdir,'mono-THBS1_VlnPlot.png'),width = 7,height = 7)

thbs1.expr <- FetchData(adata, vars = "THBS1")[[1]]
adata$thbs1_group <- ifelse(thbs1.expr > 0,"THBS1_high", "THBS1_low") %>%factor(levels = c("THBS1_low","THBS1_high"))
# adata$thbs1_group <- cut(thbs1.expr, breaks = c(-Inf, 0, 1, Inf), labels = c("THBS1_zero", "THBS1_low", "THBS1_high"), right = TRUE) %>%  factor() 
marker<-marker1
outdir_sub<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/select_geneset_barplot/geneplot/"
outdir_sub<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/select_geneset_barplot/分2组/geneplot/"
for (i in marker){
		p = VlnPlot(object = adata, features = i,pt.size = 0, group.by = 'thbs1_group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
		# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
		ggsave(paste0(outdir_sub,'mono-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
		ggsave(paste0(outdir_sub,'mono-group_', i, '_VlnPlot.png'),width = 7,height = 7)

		p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)
		# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
		ggsave(paste0(outdir_sub,'mono-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
		ggsave(paste0(outdir_sub,'mono-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)
}

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/select_geneset_barplot/分2组/"
i='select_geneset'
p = VlnPlot(object = adata, features = i,pt.size = 0, group.by = 'thbs1_group')+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
# p$data %>%group_by(ident) %>%summarise(mean = mean(HIF1A, na.rm = TRUE))
ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot.pdf'),width = 7,height = 7)
ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot.png'),width = 7,height = 7)
p <- get_sig_ana(p,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",map_sig_level=T,only_sig = T,reset_y=T)#+ylim(-0.1,0.2)
# get_sig_ana(p$data,as.formula(paste0(i,'~ident')),group.by='ident',prefix="mono_group",outdir=outdir,map_sig_level=T,only_sig = T)
ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot_sig_ana.pdf'),width = 7,height = 10)
ggsave(paste0(outdir,'mono-group_', i, '_VlnPlot_sig_ana.png'),width = 7,height = 10)

