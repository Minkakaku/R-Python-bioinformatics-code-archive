library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/Plasma cells/Plasma.rds")
levels(adata)
levels(adata$group)<-c('N','P','AID')

MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=TRUE)
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
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/cor/"
levels(adata$cell_type)

mydata<-FetchData(adata,vars=c('MIF','IL6',"group"))
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(MIF, IL6, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = MIF, y = IL6, colour = group)) +
        geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        stat_cor(aes(label = sprintf("r == %.2f", ..r..),group = group),        # 关键：按 group 分组
           formula = y~x,
           output.type = "expression",
          #  label.x = 3.5,
					 parse = TRUE,
           size = 4) +
		# coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "group",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/cor/Plasma_group_cor.png", plot = a1, width = 10, height = 7)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/cor/Plasma_group_cor.pdf", plot = a1, width = 10, height = 7)


mydata<-FetchData(adata,vars=c('MIF','IL6',"group",'sample'))
mean_by_sample <- mydata %>%
  group_by(sample, group) %>%               # 按样本+分组两列分组
  summarise(across(c(MIF, IL6), mean, na.rm = TRUE), .groups = "drop")
cor <- mean_by_sample %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(MIF, IL6, use = "complete.obs"))
a1 <- ggplot(mean_by_sample,
             aes(x = MIF, y = IL6, colour = group)) +
        geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        stat_cor(aes(label = sprintf("r == %.2f", ..r..),group = group),        # 关键：按 group 分组
           formula = y~x,
           output.type = "expression",
          #  label.x = 3.5,
					 parse = TRUE,
           size = 4) +
		# coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "group",
         x = "MIF",
         y = "IL6") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/cor/Plasma_sample_cor.png", plot = a1, width = 10, height = 7)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/cor/Plasma_sample_cor.pdf", plot = a1, width = 10, height = 7)

Idents(adata)<-'group'
for (pathway in c('MIF','IL6')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)+scale_fill_manual(values = c('N'='#1f77b4','P'='#ff7f0e','AID'='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.pdf"),width = 8,height = 10)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.png"),width = 8,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~ident')),group.by='ident',prefix=paste0("Plasma_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c('N'='#1f77b4','P'='#ff7f0e','AID'='#279e68'))
# p<-p+scale_y_continuous(limits=c(-0.2,0.8))
p<-p+scale_y_continuous(limits=c(-0.03,0.05))
# get_sig_ana(p$data,as.formula(paste0(pathway,'~ident')),group.by='ident',prefix=paste0("Plasma_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 8,height = 10)
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.png"),width = 8,height = 10)
}