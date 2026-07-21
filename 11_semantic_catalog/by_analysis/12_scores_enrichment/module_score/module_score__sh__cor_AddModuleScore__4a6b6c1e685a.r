# conda activate scMetabolism
library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.rds")
outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_2/"

mydata<-FetchData(adata,vars=c('LGALS9','HAVCR2','TIGIT',"cell_type","group"))
mydata<-mydata[mydata$cell_type %in% c("CD8 Cytotoxic Tcells","CD8 naive Tcells","CD8 Tcm","CD4 Treg"),]
mydata <- mydata[!(mydata$LGALS9 == 0 | mydata$HAVCR2 == 0 | mydata$TIGIT == 0), ]
xlim <- range(mydata$LGALS9)
ylim <- range(mydata$HAVCR2)
xlim <- c(0,4)
ylim <- c(0,4)
# cor <- mydata %>%
#   group_by(group) %>%  # 按照'group'列分组
#   summarise(correlation = cor(LGALS9, HAVCR2, use = "complete.obs"))
# cor1=round(cor$correlation[1],2)
# cor2=round(cor$correlation[2],2)

cor_havcr2 <- cor.test(mydata$LGALS9, mydata$HAVCR2, method = "spearman")
cor_tigit <- cor.test(mydata$LGALS9, mydata$TIGIT, method = "spearman")

a<-ggplot(mydata,aes(x=LGALS9,y=HAVCR2,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
    annotate("text", x = Inf, y = Inf, 
           label = paste0("R = ", round(cor_havcr2$estimate, 3), 
                          ", p = ", format(cor_havcr2$p.value, scientific = TRUE, digits = 3)), 
           hjust = 1.1, vjust = 1.1) +
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "control",
         x = "LGALS9",
         y = "HAVCR2") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_2/LGALS9-HAVCR2_group_cor.png", sep = ""), plot = a, width = 7, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_2/LGALS9-HAVCR2_group_cor.pdf", sep = ""), plot = a, width = 7, height = 6)

a<-ggplot(mydata,aes(x=LGALS9,y=TIGIT,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
    annotate("text", x = Inf, y = Inf, 
           label = paste0("R = ", round(cor_tigit$estimate, 3), 
                          ", p = ", format(cor_tigit$p.value, scientific = TRUE, digits = 3)), 
           hjust = 1.1, vjust = 1.1) +
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "control",
         x = "LGALS9",
         y = "TIGIT") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_2/LGALS9-TIGIT_group_cor.png", sep = ""), plot = a, width = 7, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20250829/10.Mission/result_2/LGALS9-TIGIT_group_cor.pdf", sep = ""), plot = a, width = 7, height = 6)