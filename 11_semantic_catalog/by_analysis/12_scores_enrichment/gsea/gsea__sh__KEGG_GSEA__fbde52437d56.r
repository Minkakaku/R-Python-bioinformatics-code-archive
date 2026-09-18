setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/Figure1E/")
library(clusterProfiler)

# 1. 读取 GMT 文件（clusterProfiler 标准函数）
# gmt_df <- read.gmt("/mnt/NexData_file/home/zwl/Project/FWAE20242338/mmu.pathway.gmt")

# # 2. 你要保留的通路
# selected_pathways <- c(
#   "IL-17 signaling pathway",
#   "AGE-RAGE signaling pathway in diabetic complications",
#   "Inflammatory bowel disease",
#   "Thyroid cancer",
#   "HIF-1 signaling pathway",
#   "Fluid shear stress and atherosclerosis",
#   "Platinum drug resistance",
#   "Allograft rejection",
#   "Prolactin signaling pathway",
#   "Human cytomegalovirus infection",
#   "ErbB signaling pathway",
#   "Inflammatory mediator regulation of TRP channels",
#   "MAPK signaling pathway",
#   "Central carbon metabolism in cancer",
#   "Mitophagy - animal",
#   "Pathways in cancer",
#   "PI3K-Akt signaling pathway",
#   "Transcriptional misregulation in cancer",
#   "TGF-beta signaling pathway",
#   "Oocyte meiosis"
# )

# # 3. 筛选通路
# filtered_df <- subset(gmt_df, term %in% selected_pathways)

# DEG = read.csv(file = "./clusterH0/GO_Biological_Process_2023.Mouse.enrichr.reports.txt",header = T,sep='\t')
# DEG<-DEG[1:5,]
# DEG$group<-'H0'

groups <- c("H0", "H2", "H8", "H24", "H72")
# 初始化空列表存储结果
all_top5 <- list()
# 循环读取每个分组的数据
for(group in groups) {
  # 构建文件路径
  file_path <- paste0("./cluster", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  # 读取数据
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  # 提取前5行
  deg_top5 <- deg[1:5, ]
  # 添加分组列
  deg_top5$group <- group
  # 存储到列表
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
# KEGG_res_df <- KEGG_res_df %>% arrange(GeneRatio_num)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
library(RColorBrewer)
library(tidyverse)
# 绘图

p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = rev(brewer.pal(7, "RdYlBu"))) +
  scale_size_continuous(
    range = c(0, 15),
    breaks = c(5, 10, 15, 20, 25),
    labels = c("5", "10", "15", "20", "25"),
    name = "Count"
  ) +
  theme_bw() +
  labs(x = "group", fill = "Adjusted.P.value")
  
  # ✅ X 轴只显示 4 个标签，刻度保留
  # scale_x_continuous(
  #   breaks = scales::breaks_extended(n = 4),  # 自动生成 4 个刻度
  #   labels = function(x) round(x, 3)          # 保留 3 位小数
  # )

ggsave("GOBP_dotplot.pdf",width = 10,height = 10,dpi=300,p)
ggsave("GOBP_dotplot.png",width = 10,height = 10,dpi=300,p)


setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/Proximal tubule cells/")
groups <- c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4","cluster5", "cluster6", "cluster7", "cluster8", "cluster9")
all_top5 <- list()
for(group in groups) {
  file_path <- paste0("./", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  deg_top5 <- deg[1:5, ]
  deg_top5$group <- group
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
# 绘图
p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = rev(brewer.pal(7, "RdYlBu"))) +
  scale_size_continuous(
    range = c(0, 15),
    breaks = c(50, 100, 150, 200, 250, 300),
    labels = c("50", "100", "150", "200", "250", "300"),
    name = "Count"
  ) +
  theme_bw() +
  labs(x = "group", fill = "Adjusted.P.value")

ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/GOBP_dotplot.pdf",width = 15,height = 20,dpi=300,p)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/GOBP_dotplot.png",width = 15,height = 20,dpi=300,p)


################################################################################
#######GSEA 后面不用做
################################################################################
data_all_sort <- DEG %>% arrange(desc(avg_log2FC))
geneList <- data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$gene
# names(geneList) %in% gmt$gene

res <- GSEA(geneList,TERM2GENE = filtered_df,  minGSSize = 1,
            pvalueCutoff = 0.9,pAdjustMethod = 'BH',
            maxGSSize = 500)
write.csv(file ="/mnt/NexData_file/home/zwl/Project/FWAE20242338/5.KEGG_GSEA.csv",res@result)
num = length(res@result$ID)
print(num)
for (i in 1:num) {
  p <- enrichplot::gseaplot2(res, geneSetID = i,
                             color = 'green2',
                             title = res$Description[i])
  
  # 保存 PDF
  cairo_pdf(filename = paste0('/mnt/NexData_file/home/zwl/Project/FWAE20242338/5.',res@result[["ID"]][i],'_GSEA_plot.pdf'),
            width=8,height=6)
  print(p)
  grid.text(label = paste0("NES = ", round(res$NES[i],2), "\n", "p value = ", signif(res$pvalue[i], 2)),
            x = unit(0.95,"npc"), y = unit(0.94,"npc"),
            just = c("right","top"), gp = gpar(fontsize=13,fontface="bold"))
  dev.off()
    
}

