library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/Plasma cells/Plasma.rds")
levels(adata)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
#转换成list
gene<-as.list(unique(geneset))
gene$V1<-c(gene$V1,'IGHA1','IGHG1','IGHM','HLA-DRA','HLA-DQA1','HLA-DPA1','TLR7','TLR8','TLR9','TNF','IL6','IFNG','BTG3')
# genedata<-FetchData(adata,var=gene$V1)
genedata<-FetchData(adata,var=c(gene$V1,"group"))

data <- cor(genedata, method = "pearson")
# data <- data['MIF',]
write.table(data,"/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/pearson.xls",sep="\t",quote=F,row.names=T,col.names=T)

testcor <- cor.mtest(genedata, method = "pearson",conf.level=0.95)
# testcor <- testcor$p['MIF',]
write.table(testcor,"/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/cortest.xls",sep="\t",quote=F,row.names=T,col.names=T)

x_col <- 'MIF'
y_cols <- setdiff(gene$V1, c(x_col,'CXCR7',' ENSG00000306549','ENSG00000297913'))
y_cols[8:10]<-c("HLA_DRA","HLA_DQA1","HLA_DPA1")
colnames(genedata)[9:11]<-c("HLA_DRA","HLA_DQA1","HLA_DPA1")
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
subdata<-genedata[,c("MIF",y_col,"group")]
colnames(subdata)[2]<-'y'
cor <- subdata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(MIF, y, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(subdata[which(subdata$group=="Negative"),],aes(x=MIF,y=y,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Negative",
         x = "MIF",
         y = y_col) +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(subdata[which(subdata$group=="Positive"),],aes(x=MIF,y=y,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Positive",
         x = "MIF",
         y = y_col) +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(subdata[which(subdata$group=="Rheumatic Diseases"),],aes(x=MIF,y=y,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Rheumatic Diseases",
         x = "MIF",
         y = y_col) +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/","scatter_", y_col, ".png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/","scatter_", y_col, ".pdf", sep = ""), plot = a, width = 18, height = 6)
}
genedata<-FetchData(adata,var=c('group','cell_type',gene$V1))
write.table(genedata,"/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/cor/genedata.xls",sep="\t",quote=F,row.names=T,col.names=T)


MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=FALSE)
#转换成list
gene<-c(as.list(MIF),as.list(IL6))
names(gene)<-c('MIF','IL6')
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'MIF'
colnames(adata@meta.data)[17]<-'IL6'
# outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/AddModuleScore/"
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/AddModuleScore_newIL6/"
levels(adata$cell_type)
levels(adata$cell_type)[2]<-"IgG Plasma"
levels(adata$cell_type)[3]<-"IgG Plasma"
levels(adata$cell_type)[3]<-"IgM Plasma"
levels(adata$cell_type)[4]<-"IgM Plasma"
Idents(adata)<-adata$cell_type
levels(adata$group)<-c('N','P','AID')

for (pathway in c('MIF','IL6')){
	print(pathway)
# p <- VlnPlot(adata,features=pathway,pt.size = 0)+scale_fill_manual(values = c('B cells'='#1f77b4','Erythroid cells'='#ff7f0e','HSPCs'='#279e68','Megakaryocytes'='#d62728','Myeloid cells'='#aa40fc','NK cells'='#8c564b','Neutrophils'='#e377c2','Plasma cells'='#b5bd61','T cells'='#17becf'))
# ggsave(paste0(outdir,"Plasma-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
# ggsave(paste0(outdir,"Plasma-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('N'='#1f77b4','P'='#ff7f0e','AID'='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasma_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c('N'='#1f77b4','P'='#ff7f0e','AID'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasma_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

# mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2",pathway))
# a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=UMAP_1,y=UMAP_2,colour=IL6)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Negative")
# a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=UMAP_1,y=UMAP_2,colour=IL6)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Positive")
# a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=UMAP_1,y=UMAP_2,colour=IL6)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Rheumatic Diseases")
# a<-a1+a2+a3
# ggsave(paste0(outdir,"Plasma-",pathway,"-UMAP.pdf"),width = 18,height = 6)
# ggsave(paste0(outdir,"Plasma-",pathway,"-UMAP.png"),width = 18,height = 6)

}

mydata<-FetchData(adata,vars=c('MIF','IL6',"group"))
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(MIF, IL6, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Negative",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Positive",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Rheumatic Diseases",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/AddModuleScore/Plasma_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/AddModuleScore/Plasma_cor.pdf", sep = ""), plot = a, width = 18, height = 6)



adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/mapping_raw/Main_mapping.rds")
levels(adata$sub_cell_type)
levels(adata$sub_cell_type)[3]<-"IgG Plasma"
levels(adata$sub_cell_type)[4]<-"IgG Plasma"
levels(adata$sub_cell_type)[4]<-"IgM Plasma"
levels(adata$sub_cell_type)[5]<-"IgM Plasma"
Idents(adata)<-adata$sub_cell_type
adata$cell_type<-Idents(adata)
saveRDS(adata,"/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/Main_mapping.rds")

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'MIF'
colnames(adata@meta.data)[17]<-'IL6'

for (pathway in c('MIF','IL6')){
	print(pathway)
# p <- VlnPlot(adata,features=pathway,pt.size = 0)+scale_fill_manual(values = c('B cells'='#1f77b4','Erythroid cells'='#ff7f0e','HSPCs'='#279e68','Megakaryocytes'='#d62728','Myeloid cells'='#aa40fc','NK cells'='#8c564b','Neutrophils'='#e377c2','Plasma cells'='#b5bd61','T cells'='#17becf'))
# ggsave(paste0(outdir,"Plasma-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
# ggsave(paste0(outdir,"Plasma-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('Negative'='#1f77b4','Positive'='#ff7f0e','Rheumatic Diseases'='#279e68'))
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Negative='#1f77b4',Positive='#ff7f0e','Rheumatic Diseases'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2",pathway,"group"))
a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Negative")
a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Positive")
a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Rheumatic Diseases")
a<-a1+a2+a3
ggsave(paste0(outdir,"Main-",pathway,"-UMAP.pdf"),width = 18,height = 6)
ggsave(paste0(outdir,"Main-",pathway,"-UMAP.png"),width = 18,height = 6)
}

mydata<-FetchData(adata,vars=c('MIF','IL6',"group"))
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(MIF, IL6, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Negative",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Positive",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
    labs(title = "Rheumatic Diseases",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/AddModuleScore/Main_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/AddModuleScore/Main_cor.pdf", sep = ""), plot = a, width = 18, height = 6)
