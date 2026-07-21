# ==========================================
# 1. 环境准备与加载包
# ==========================================
# 如果未安装，请取消下一行的注释运行一次
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Rn.eg.db")

library(clusterProfiler)
library(org.Rn.eg.db)  # 大鼠官方基因注释包
library(data.table)    # 极其高效的 fread 读取包
library(dplyr)
library(stringr)
library(ggplot2)

# ==========================================
# 2. 数据读取与交集基因提取
# ==========================================
data_BA <- fread(file.choose()) %>% 
  filter(regulate == "Up")

data_DB <- fread(file.choose()) %>% 
  filter(regulate == "Down")

target_gene <- intersect(data_BA$gene, data_DB$gene)

if(length(target_gene) == 0) {
  stop("🚨【错误】交集基因数量为 0！请检查两张表中的 gene 列名字是否匹配。")
}
message(">>> 成功提取交集核心基因共: ", length(target_gene), " 个")

gene_mapped <- bitr(target_gene, 
                    fromType = "ENSEMBL", 
                    toType   = "ENTREZID", 
                    OrgDb    = org.Rn.eg.db)

# ==========================================
# 3. 在线KEGG 富集分析
# ==========================================

kegg_res <- enrichKEGG(gene         = gene_mapped$ENTREZID,
                       organism     = "rno",   # 大鼠在KEGG官方的物种代号
                       keyType      = "kegg",  # 默认的标准 Entrez ID 模式
                       pvalueCutoff = 1,       # 设为1全放开，交给你在画图时筛选
                       qvalueCutoff = 1)

# 转换为标准数据框
res_df <- as.data.frame(kegg_res)

# 检查是否成功捞到通路
if(nrow(res_df) == 0) {
  stop("🚨 很遗憾，即使调用官方总表也没有富集到任何 KEGG 通路，说明这些基因在 KEGG 里全挂空挡了。")
}

# ==========================================
# ==========================================
# 4. negLogP排序
# ==========================================
plot_df <- res_df %>%
  arrange(pvalue) %>%         
  head(15) %>%                 
  mutate(
    Gene_hit = as.numeric(sub("/.*", "", GeneRatio)),
    Gene_total = as.numeric(sub(".*/", "", GeneRatio)),
    GeneRatio_num = Gene_hit / Gene_total,
    
    negLogP = -log10(pvalue),
    Terms = reorder(Description, negLogP) 
  )

# ==========================================
# 5. 绘制
# ==========================================
p_bar <- ggplot(plot_df, aes(x = GeneRatio_num, y = Terms, fill = pvalue)) +  geom_col(color = "black", width = 0.8, linewidth = 0.3) +
  
  scale_fill_gradientn(colors = c("#E64B35FF", "#F39B7FFF", "#4DBBD5FF")) +
  
  theme_bw() +
  labs(
    title = "",
    x = "GeneRatio", 
    y = "",
    fill = "p-value"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

# 显示图表
print(p_bar)
ggsave("KEGG_BA(Up)_DB(Down).pdf", p_bar, width = 8, height = 6, device = "pdf")
ggsave("KEGG_BA(Up)_DB(Down).png", p_bar, width = 8, height = 6, dpi = 300)

