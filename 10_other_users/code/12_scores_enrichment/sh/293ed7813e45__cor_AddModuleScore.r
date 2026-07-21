library(Seurat)
library(corrplot)
library(ggplot2)
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/Plasma cells/Plasma.rds")
levels(adata)
MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
#转换成list
gene<-c(as.list(MIF))
names(gene)<-c('MIF')
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'MIF'

levels(adata$cell_type)
levels(adata$cell_type)[2]<-"IgG Plasma"
levels(adata$cell_type)[3]<-"IgG Plasma"
levels(adata$cell_type)[3]<-"IgM Plasma"
levels(adata$cell_type)[4]<-"IgM Plasma"
Idents(adata)<-adata$cell_type
levels(adata$group)<-c('N','P','AID')

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
#转换成list
gene<-as.list(unique(geneset))
gene$V1<-c(gene$V1,'IGHA1','IGHG1','IGHM','HLA-DRA','HLA-DQA1','HLA-DPA1','TLR7','TLR8','TLR9','TNF','IL6','IFNG','BTG3')
# genedata<-FetchData(adata,var=gene$V1)
genedata<-FetchData(adata,var=c(gene$V1,"group"))

x_col <- 'MIF'
y_cols <- setdiff(gene$V1, c(x_col,'CXCR7',' ENSG00000306549','ENSG00000297913'))
y_cols[8:10]<-c("HLA_DRA","HLA_DQA1","HLA_DPA1")
colnames(genedata)[9:11]<-c("HLA_DRA","HLA_DQA1","HLA_DPA1")
for (y_col in y_cols) {
subdata<-genedata[,c("MIF",y_col,"group")]
colnames(subdata)[2]<-'y'
# cor <- subdata %>%
#   group_by(group) %>%  # 按照'group'列分组
#   summarise(correlation = cor(MIF, y, use = "complete.obs"))
# cor1=round(cor$correlation[1],2)
# cor2=round(cor$correlation[2],2)
# cor3=round(cor$correlation[3],2)
a1<-ggplot(subdata[which(subdata$group=="N"),],aes(x=MIF,y=y,colour=group)) + 
    # geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "solid", size = 1) +  # 添加线性拟合线
		stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE,
           size = 5) +
    scale_colour_manual(values = c("#1f77b4")) +
    scale_fill_manual(values = c("#1f77b4")) +
    labs(title = "N",
         x = "MIF",
         y = y_col) +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
		legend.position = "none"
  )
a2<-ggplot(subdata[which(subdata$group=="P"),],aes(x=MIF,y=y,colour=group)) + 
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "solid", size = 1) +  # 添加线性拟合线
		stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE,
           size = 5) +
    scale_colour_manual(values = c( "#ff7f0e")) +
    scale_fill_manual(values = c("#ff7f0e")) +
    labs(title = "P",
         x = "MIF",
         y = y_col) +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
		legend.position = "none"
  )
a3<-ggplot(subdata[which(subdata$group=="AID"),],aes(x=MIF,y=y,colour=group)) + 
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "solid", size = 1) +  # 添加线性拟合线
		stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE,
           size = 5) +
    scale_colour_manual(values = c("#279e68")) +
    scale_fill_manual(values = c("#279e68")) +
    labs(title = "AID",
         x = "MIF",
         y = y_col) +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
		legend.position = "none"
  )
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result3/","scatter_", y_col, ".png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result3/","scatter_", y_col, ".pdf", sep = ""), plot = a, width = 18, height = 6)
}
genedata<-FetchData(adata,var=c('group','cell_type',gene$V1))
write.table(genedata,"/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result3/genedata.xls",sep="\t",quote=F,row.names=T,col.names=T)



adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/Main_mapping.rds")
MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=FALSE)
#转换成list
gene<-c(as.list(MIF),as.list(IL6))
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'MIF'
colnames(adata@meta.data)[17]<-'IL6'

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
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result5/Main_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result5/Main_cor.pdf", sep = ""), plot = a, width = 18, height = 6)


a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=MIF,y=IL6,colour=group)) + 
    # geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = TRUE, alpha = 0.5, level = 0.99, linetype = "solid", size = 1) +  # 添加线性拟合线
		stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE,
           size = 5) +
    scale_colour_manual(values = c("#1f77b4")) +
    scale_fill_manual(values = c("#1f77b4")) +
    labs(title = "N",
         x = "MIF",
         y = 'IL6') +
		ylim(0, 0.3) +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
		legend.position = "none"
  )
a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_smooth(method = "lm", se = TRUE, alpha = 0.5, linetype = "solid", size = 1) +  # 添加线性拟合线
		stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE,
           size = 5) +
    scale_colour_manual(values = c( "#ff7f0e")) +
    scale_fill_manual(values = c("#ff7f0e")) +
    labs(title = "P",
         x = "MIF",
         y = 'IL6') +
		ylim(0, 0.3) +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
		legend.position = "none"
  )
a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=MIF,y=IL6,colour=group)) + 
    geom_smooth(method = "lm", se = TRUE, alpha = 0.5, linetype = "solid", size = 1) +  # 添加线性拟合线
		stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE,
           size = 5) +
    scale_colour_manual(values = c("#279e68")) +
    scale_fill_manual(values = c("#279e68")) +
    labs(title = "AID",
         x = "MIF",
         y = 'IL6') +
		ylim(0, 0.3) +
    theme_bw(base_size = 18) +  # 使用简洁主题
    theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
		legend.position = "none"
  )
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result5/Main_cor_new.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result5/Main_cor_new.pdf", sep = ""), plot = a, width = 18, height = 6)