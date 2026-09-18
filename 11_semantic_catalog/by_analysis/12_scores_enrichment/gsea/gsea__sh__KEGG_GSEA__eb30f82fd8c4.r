setwd('/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/')
library(clusterProfiler)
library(grid)
# 1. 读取 GMT 文件（clusterProfiler 标准函数）
gmt_df <- read.gmt("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/mmu.pathway.gmt")

# 2. 你要保留的通路
selected_pathways <- c(
  "IL-17 signaling pathway",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Inflammatory bowel disease",
  "Thyroid cancer",
  "HIF-1 signaling pathway",
  "Fluid shear stress and atherosclerosis",
  "Platinum drug resistance",
  "Allograft rejection",
  "Prolactin signaling pathway",
  "Human cytomegalovirus infection",
  "ErbB signaling pathway",
  "Inflammatory mediator regulation of TRP channels",
  "MAPK signaling pathway",
  "Central carbon metabolism in cancer",
  "Mitophagy - animal",
  "Pathways in cancer",
  "PI3K-Akt signaling pathway",
  "Transcriptional misregulation in cancer",
  "TGF-beta signaling pathway",
  "Oocyte meiosis"
)

# 3. 筛选通路
filtered_df <- subset(gmt_df, term %in% selected_pathways)

DEG = read.csv(file = "/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/Comp Fibroblasts_MI_DZ-MI_diffgenes.csv",header = T)
head(DEG$gene)
head(filtered_df)

KEGG_res <- enricher(DEG$gene, TERM2GENE=filtered_df)
write.csv(file ="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/4.KEGG_enrichment.csv",KEGG_res@result)

KEGG_res_df = KEGG_res@result
KEGG_res_df$LogP = -log2(KEGG_res_df$pvalue)

KEGG_res_df <- KEGG_res_df %>%
  mutate(
    Num = as.numeric(sapply(strsplit(GeneRatio, "/"), "[[", 1)),
    Den = as.numeric(sapply(strsplit(GeneRatio, "/"), "[[", 2)),
    GeneRatio_num = Num / Den  # 最终数值
  )

KEGG_res_df <- KEGG_res_df %>% arrange(GeneRatio_num)
KEGG_res_df$Description <- factor(KEGG_res_df$Description, levels = unique(KEGG_res_df$Description))
library(RColorBrewer)
library(tidyverse)
# 绘图

p = ggplot(KEGG_res_df, aes(x = Count, y = Description, fill = pvalue, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(
    colors = brewer.pal(7, "RdYlBu"),
		limits=c(0,0.5),
    label = c( "0.01", "0.05", "0.1"," 0.25"),                          # 0-1 范围
    breaks = c( 0.01, 0.05, 0.1, 0.25),
    name = "p.value"
  ) +
	ggtitle('Comp Fibroblasts MI_DZ VS MI DOWN KEGG enrichment plot') +
  
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),  # 居中加粗
    plot.margin = margin(5, 5, 5, 5),                                   # 调整边距
    # 调整图例样式
    legend.key.height = unit(1, "cm"),      # 图例色块高度
    legend.key.width = unit(0.5, "cm"),       # 图例色块宽度
    legend.text = element_text(size = 8),     # 图例文字大小
    legend.title = element_text(size = 9),    # 图例标题大小
    legend.spacing.y = unit(0.1, "cm")        # 图例项间距
  ) +
  labs(x = "Count", fill = "p.value") +
  
  # ✅ X 轴只显示 4 个标签，刻度保留
  scale_x_continuous(
    breaks = scales::breaks_extended(n = 4),  # 自动生成 4 个刻度
    labels = function(x) round(x, 3)          # 保留 3 位小数
  )

ggsave("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/4.KEGG_dotplot.pdf",width = 8,height = 6,dpi=300,p)
ggsave("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/4.KEGG_dotplot.png",width = 8,height = 6,dpi=300,p)

################################################################################
#######GSEA
################################################################################
data_all_sort <- DEG %>% arrange(desc(avg_log2FC))
geneList <- data_all_sort$avg_log2FC
names(geneList) <- data_all_sort$gene
# names(geneList) %in% gmt$gene

res <- GSEA(geneList,TERM2GENE = filtered_df,  minGSSize = 1,
            pvalueCutoff = 0.9,pAdjustMethod = 'BH',
            maxGSSize = 500)
write.csv(file ="/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/5.KEGG_GSEA.csv",res@result)
num = length(res@result$ID)
print(num)
for (i in 1:num) {
  p <- enrichplot::gseaplot2(res, geneSetID = i,
                             color = 'green2',
                             title = res$Description[i])
  
  # 保存 PDF
  cairo_pdf(filename = paste0('/mnt/NexData_file/home/sh/Projects/AE2511106372-01/10.Mission/5.',res@result[["ID"]][i],'_GSEA_plot.pdf'),
            width=8,height=6)
  print(p)
  grid.text(label = paste0("NES = ", round(res$NES[i],2), "\n", "p value = ", signif(res$pvalue[i], 2)),
            x = unit(0.95,"npc"), y = unit(0.94,"npc"),
            just = c("right","top"), gp = gpar(fontsize=13,fontface="bold"))
  dev.off()
    
}

