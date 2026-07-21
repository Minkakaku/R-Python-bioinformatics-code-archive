# conda activate scMetabolism
library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)
library(ggpubr)
library(patchwork)

load("/mnt/NexData_file/home/sh/Projects/FWSC20251028/06.Mission/1_EpiProcess_AT2s/4_srt_after_SCT.Rda")
DefaultAssay(srt_after)<-"RNA"
srt_small <- srt_after
DefaultAssay(srt_small) <- "RNA"
srt_small[['integrated']] <- NULL
srt_small@reductions <- srt_small@reductions[c("pca", "umap")]

srt_small@meta.data <- srt_small@meta.data[, !names(srt_small@meta.data) %in% c("integrated_snn_res.0.1","integrated_snn_res.0.2","integrated_snn_res.0.3","integrated_snn_res.0.4","integrated_snn_res.0.5","integrated_snn_res.0.6","integrated_snn_res.0.7","integrated_snn_res.0.8","integrated_snn_res.0.9","integrated_snn_res.1","integrated_snn_res.1.1","integrated_snn_res.1.2","integrated_snn_res.1.3","integrated_snn_res.1.4","integrated_snn_res.1.5")]

srt_small$cell_type<-srt_small$CellType
srt_small$group<-srt_small$Group
saveRDS(srt_small,"/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/Main_srt_after.rds")

df_sub = read.table("/mnt/NexData_file/home/sh/Projects/FWSC20251028/08.Mission/Markerplot_mouse/Epi_cells_features.txt", header=TRUE)
adata <- df_sub
bdata <- srt_small
# 提取两个RDS文件的meta.data行名和指定列（假设列名为'cell_type'）
adata_cell_names <- rownames(adata)
bdata_cell_names <- rownames(bdata@meta.data)
adata_cell_types <- as.character(adata$cell_type)
bdata_cell_types <- as.character(bdata@meta.data$cell_type)

# 匹配行名并替换对应列的值，保留未匹配行的原始值
mapped_indices <- match(bdata_cell_names, adata_cell_names)
matched_cell_types <- adata_cell_types[mapped_indices]
bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]

# 将新的细胞类型注释赋值给bdata，并转换回factor类型
bdata@meta.data$cell_type <- as.factor(bdata_cell_types)
Idents(bdata) <- bdata@meta.data$cell_type

adata<-subset(bdata,cell_type %in% c("Epithelials"),invert=TRUE)
levels(adata)
DefaultAssay(adata)<-"RNA"
genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/script/pathway.csv",sep = "\t",head=FALSE)
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
colnames(adata@meta.data)[21]<-'HGF_cMet_signaling'

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/cor_AddModuleScore/"

mydata<-FetchData(adata,vars=c('HGF_cMet_signaling','HALLMARK_WNT_BETA_CATENIN_SIGNALING','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_PI3K_AKT_MTOR_SIGNALING',"cell_type","group"))
xlim <- range(mydata$HGF_cMet_signaling)
ylim <- range(mydata$HALLMARK_WNT_BETA_CATENIN_SIGNALING)
# xlim <- c(0,4.15)
# ylim <- c(-0.125,0.175)
cor <- mydata %>%
  group_by(group) %>%  # 按照'group'列分组
  summarise(correlation = cor(HGF_cMet_signaling, HALLMARK_WNT_BETA_CATENIN_SIGNALING, use = "complete.obs"))
cor1=round(cor$correlation[1],2)
cor2=round(cor$correlation[2],2)
cor3=round(cor$correlation[3],2)
a1<-ggplot(mydata[which(mydata$group=="Control"),],aes(x=HGF_cMet_signaling,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor1,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "Control",
         x = "HGF_cMet_signaling",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a2<-ggplot(mydata[which(mydata$group=="ELA"),],aes(x=HGF_cMet_signaling,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor2,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "ELA",
         x = "HGF_cMet_signaling",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a3<-ggplot(mydata[which(mydata$group=="HGF"),],aes(x=HGF_cMet_signaling,y=HALLMARK_WNT_BETA_CATENIN_SIGNALING,colour=group)) + 
    geom_point(size = 3, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
		stat_poly_eq(aes(label=paste("italic(r)~`=`~",cor3,sep = "")),formula = y~x,parse=T,size=10) + 
		coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = "HGF",
         x = "HGF_cMet_signaling",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
a<-a1+a2+a3
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_group_cor.png", sep = ""), plot = a, width = 18, height = 6)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_group_cor.pdf", sep = ""), plot = a, width = 18, height = 6)

cor <- mydata %>%
  group_by(cell_type) %>%  # 按照'group'列分组
  summarise(correlation = cor(HGF_cMet_signaling, HALLMARK_WNT_BETA_CATENIN_SIGNALING, use = "complete.obs"))
a1 <- ggplot(mydata,
             aes(x = HGF_cMet_signaling, y = HALLMARK_WNT_BETA_CATENIN_SIGNALING, colour = cell_type)) +
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
         x = "HGF_cMet_signaling",
         y = "HALLMARK_WNT_BETA_CATENIN_SIGNALING") +
    theme_bw() +  # 使用简洁主题
    theme(plot.title = element_text(hjust = 0.5),  # 居中标题
          panel.grid.major = element_blank(),  # 添加网格线
          panel.grid.minor = element_blank())  # 去掉次要网格线
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_cell_type_cor.png", sep = ""), plot = a1, width = 10, height = 7)
ggsave(paste("/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/cor_AddModuleScore/HALLMARK_WNT_BETA_CATENIN_SIGNALING_cell_type_cor.pdf", sep = ""), plot = a1, width = 10, height = 7)

library(purrr)
library(Hmisc)
# 1. 读入数据 -------------------------------------------------------------
# mydata 已经是 data.frame，列名末两列为 cell_type、group
path_cols <- setdiff(names(mydata), c("cell_type", "group"))
setwd("/mnt/NexData_file/home/sh/Projects/FWSC20251028/09.Mission/cor_AddModuleScore/cor/")
levels(mydata$cell_type)[13]<-"NKs_TCs"
# 2. 按组合循环 -----------------------------------------------------------
combos <- unique(mydata[c("cell_type", "group")])
all_res <- map_dfr(1:nrow(combos), \(i) {
  ct  <- combos$cell_type[i]
  grp <- combos$group[i]
  
  sub <- mydata %>% filter(cell_type == ct & group == grp)
  if (nrow(sub) < 3) return(NULL)
  
  # 计算 spearman 相关 + P
  mat <- as.matrix(sub[path_cols])
  cc  <- rcorr(mat, type = "spearman")
  
  # 转成长表
  tmp <- data.frame(
    cell_type = ct,
    group     = grp,
    var1      = rep(path_cols, each = length(path_cols)),
    var2      = rep(path_cols, length(path_cols)),
    cor       = as.vector(cc$r),
    p.value   = as.vector(cc$P)
  ) %>% filter(var1 != var2)                # 去掉自相关
  
  # 3. 可选：画热图 ---------------------------------------------------------
  pdf(paste0(ct, "_", grp, "_corr.pdf"), width = 8, height = 7)
  corrplot(cc$r, method = "color", type = "upper",
           tl.cex = 0.8, tl.col = "black", number.cex = 0.7,
           col = colorRampPalette(c("red", "white", "blue"))(100)[100:1],
           addCoef.col = 'grey50', diag = FALSE)
  dev.off()
  png(paste0(ct, "_", grp, "_corr.png"), width = 800, height = 700)
  corrplot(cc$r, method = "color", type = "upper",
           tl.cex = 0.8, tl.col = "black", number.cex = 0.7,
           col = colorRampPalette(c("red", "white", "blue"))(100)[100:1],
           addCoef.col = 'grey50', diag = FALSE)
  dev.off()
  tmp
})

# 4. 写出总表 -------------------------------------------------------------
write.table(all_res, "all_cor_celltype_group.csv")