#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: DEG_group.R
# 功能描述: 读取一个差异基因XLS文件，提取细胞类型名称，为后续GSEA分析做准备
# 作者: zwl
# 创建时间: 2026-03-05
# 版本: v1.0
# 依赖包: 
# 使用方法: 
#   1. 直接运行: conda activate "/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/"
#                Rscript process_diff_genes.R
#   2. 需确保所有XLS文件路径正确，且文件格式为标准XLS
# 注意事项: 
#   - 路径中使用正斜杠"/"，避免反斜杠转义问题
#   - 依赖readxl包，运行前需安装：install.packages("readxl")
# ==============================================================================

# -------------------------- 加载依赖包 --------------------------
rm(list=ls())
gc()
library(dplyr)
library(Seurat)
library(tidyverse)
library(dplyr)# 用于数据处理
library(reshape2)
####vcanol####
library(ggplot2)
library(ggrepel)
library(gt) # 制作表格
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="input anno raw rds")
# argv <- add_argument(argv,"--prefix", help="group name prefix",default='F')
argv <- add_argument(argv,"--outpath", help="the output dir")
argv <- add_argument(argv,"--subcelltype", help="the cell_type used",default='all')
argv <- add_argument(argv,"--comparison", help="Group or sample for comparison eg Test:Con",default='F')
argv <- parse_args(argv)

rds <- argv$rds
outpath <- argv$outpath
if(!dir.exists(outpath)) {dir.create(outpath, recursive = T)}
subcelltype <- unlist(strsplit(argv$subcelltype, ","))
diffgroup <- argv$comparison
# diffgroup <- c("Dgal:control","EGT:Dgal","EGT:control")

setwd(outpath)
dir.create('./DiffGene')
adata<- readRDS(rds)
# Idents(adata) <- "cell_type"
levels(adata)
table(adata$cell_type)

difftype = c(levels(adata$cell_type))
diffgroup <- strsplit(diffgroup, ":")


for(n in 1:length(difftype)){
  cluster_diff_sub<-difftype[n]
  if(cluster_diff_sub=="All"){
    x = adata
  }else{
    barcodes<-rownames(adata@meta.data[which(adata$cell_type %in% cluster_diff_sub),])
    SeuratObject_splitlist<-subset(adata,cells = barcodes)
    #SeuratObject_splitlist$cell_type <- droplevels(SeuratObject_splitlist$cell_type)
    x <- SeuratObject_splitlist
  }
  
  tryCatch(
    {
      for(m in 1:length(diffgroup)){
        # 计算差异表达基因
        diff1<-diffgroup[[m]][1]
        diff2<-diffgroup[[m]][2]
        markers <- FindMarkers(x, ident.1 =  diff1, ident.2 = diff2,group.by = "group")
        markers$gene <- rownames(markers)
        # positive_fc <- markers$avg_log2FC > 0 & markers$pct.1 > 0.1
        # negative_fc <- markers$avg_log2FC < 0 & markers$pct.2 > 0.1
        # filtered_df <- markers[positive_fc | negative_fc,]
        #markers <- select(markers,6,everything())
        # 写入文件，确保文件名唯一
        cluster_diff_sub = gsub("/","_",cluster_diff_sub)
        file_path <- paste0("./DiffGene/",cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes.csv")
        write.table(data.frame(Genes= markers$gene,markers[,1:5]), file = file_path,sep = ",", row.names = FALSE)
        markers <- markers %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
        file_path <- paste0("./DiffGene/",cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes_significant.csv")
        write.table(data.frame(Genes= markers$gene,markers[,1:5]), file = file_path,sep = ",", row.names = FALSE)
        message(paste("成功处理并保存了对象", cluster_diff_sub))
      }
      return(x)
    },
    error = function(e) {
      message(paste("处理对象", cluster_diff_sub, "时出错:", e$message))
      return(NULL)
    }
  )
}

for(n in 1:length(difftype)){
  cluster_diff_sub<-difftype[n]
  cluster_diff_sub = gsub("/","_",cluster_diff_sub)
  print(cluster_diff_sub)
  for(m in 1:length(diffgroup)){
    diff1<-diffgroup[[m]][1]
    diff2<-diffgroup[[m]][2]
    print(paste0(diff1,"vs",diff2))
    if(file.exists(paste0("./DiffGene/",cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes.csv"))){
      data = read.table(paste0("./DiffGene/",cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes.csv"),sep = ",",head=TRUE)
      if(nrow(data)==0){
        next
      }
    }else{
      print(paste0(cluster_diff_sub,":",diff1,"vs",diff2," not exists"))
      next
    }
    # markers <- data %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
    # file_path <- paste0(cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes_significant.csv")
    # write.table(markers[,1:6], file = file_path, row.names = FALSE)
    # 
    
    # colnames(data) <- c("gene_name", "PValue", "log2FC", "pct_nz_group" , "pct_nz_reference", "FDR")
    data$log10p_val_adj = -log10(data$p_val_adj)
    
    #data$log10p_val_adj[is.infinite(data$log10p_val_adj)] = 300
    cut_off_PValue =0.05 #设置FDR的阈值
    cut_off_log2FC =0.25 #设置log2FC的阈值
    data$Sig = ifelse(data$p_val_adj < cut_off_PValue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                        abs(data$avg_log2FC) >= cut_off_log2FC,  #abs绝对值
                      ifelse(data$avg_log2FC > cut_off_log2FC ,'Up','Down'),'no')
    data = data.frame(data)
    table(data$Sig)
    # data = data %>%
    #   filter(pct_nz_group >= 0.2) %>%
    #   filter(pct_nz_reference <= 0.6) %>%
    #   arrange(desc(abs(log2FC)),FDR)
    xlim_min<-min(data$avg_log2FC)-0.3
    xlim_max<-max(data$avg_log2FC)+0.3
    ylim_max <- min(ceiling(max(data$log10p_val_adj) / 50) * 50,300)
    ###绘图——基础火山图###
    colmap = data.frame(type = c("Down","no","Up"),color=c("#546de5","#d2dae2","#ff4757"))
    coloruse = colmap[colmap$type %in% data$Sig,]
    
    p1 <- ggplot(data, aes(x =avg_log2FC, y=log10p_val_adj, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
      geom_point(alpha=0.65, size=3) +  #点的透明度、大小
      scale_color_manual(values=coloruse$color) + 
      # xlim(-7,7) +  #调整点的颜色和x轴的取值范围
      geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
      geom_hline(yintercept = -log10(cut_off_PValue), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
      labs(x="log2FC", y="-log10p_val_adj") +  #x、y轴标签
      ggtitle("Volcano Plot") + #标题
      scale_x_continuous(limits = c(xlim_min, xlim_max)) +
      scale_y_continuous(limits = c(0, ylim_max), breaks = seq(0, ylim_max, 50)) +
      theme_bw() + # 主题，help(theme)查找其他个性化设置
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 30),
            legend.text = element_text(size = 15),
            legend.position="right", 
            legend.title = element_blank()
      ) 
    # p1
    
    Up_top_10 =(     #筛选差异显著上调的前10个Gene
      data %>%
        filter(Sig == 'Up') %>%
        arrange(desc(abs(avg_log2FC)),p_val_adj) %>%
        head(5)
    )
    Up_top_10 %>% gt() #数据制成表
    Down_top_10 = (       #筛选差异显著下调的前10个Gene
      data %>%
        filter(Sig == 'Down') %>%
        arrange( desc(abs(avg_log2FC)),p_val_adj) %>%
        head(5)
    )
    Down_top_10 %>% gt() #数据制成表
    
    p2 <- p1 +
      geom_label_repel(data = Up_top_10,
                       aes(avg_log2FC, log10p_val_adj, label = Genes),
                       size = 3, fill="#CCFFFF",
                       alpha = 0.65, color = "black")+
      geom_label_repel(data = Down_top_10,
                       aes(avg_log2FC, log10p_val_adj, label = Genes),
                       size = 3, fill="#FFCCCC",
                       alpha = 0.65, color = "black",max.overlaps=10)
    
    p2 #出图
    ggsave(paste0("./DiffGene/",cluster_diff_sub,"_",diff1,"-",diff2, ".Volcanoplot.pdf"),width = 10,height = 10)
    ggsave(paste0("./DiffGene/",cluster_diff_sub,"_",diff1,"-",diff2, ".Volcanoplot.png"),width = 10,height = 10)
  }
}




