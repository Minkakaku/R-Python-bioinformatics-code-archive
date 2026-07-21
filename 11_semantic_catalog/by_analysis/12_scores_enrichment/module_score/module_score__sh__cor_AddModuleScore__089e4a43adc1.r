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
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/cor/"
levels(adata$cell_type)

mydata<-FetchData(adata,vars=c('MIF','IL6',"group"))
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(MIF, IL6, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = MIF, y = IL6, colour = group, fill = group)) +
        # geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "solid", size = 1) +
        stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE, label.y = c(0.045, 0.04, 0.035),   # 三组自上而下
           size = 5) +
        scale_colour_manual(values = c("#1f77b4", "#ff7f0e", "#279e68")) +
        scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#279e68")) +
    # coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "group",
         x = "MIF",
         y = "IL6") +
  theme_bw(base_size = 18) +               # 全局字体 18 pt
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/cor/Plasma_group_cor.png", plot = a1, width = 10, height = 7)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/cor/Plasma_group_cor.pdf", plot = a1, width = 10, height = 7)

mydata1<-mydata[mydata$group=="AID",]
a1 <- ggplot(mydata1,
             aes(x = MIF, y = IL6, colour = group, fill = group)) +
        # geom_point(size = 3, alpha = 0.7, shape = 1) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "solid", size = 1) +
        stat_cor(aes(label = sprintf("r == %.2f~', '~p == %.2g", ..r.., ..p..)),
           formula = y ~ x, output.type = "expression", parse = TRUE, label.y = c(0.045),   # 三组自上而下
           size = 5) +
        scale_colour_manual(values = c("N"="#1f77b4", "P"="#ff7f0e", "AID"="#279e68")) +
        scale_fill_manual(values = c("N"="#1f77b4", "P"="#ff7f0e", "AID"="#279e68")) +
    # coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "AID",
         x = "MIF",
         y = "IL6") +
  theme_bw(base_size = 18) +               # 全局字体 18 pt
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text  = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/cor/Plasma_AID_cor.png", plot = a1, width = 10, height = 7)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/14.Mission/cor/Plasma_AID_cor.pdf", plot = a1, width = 10, height = 7)