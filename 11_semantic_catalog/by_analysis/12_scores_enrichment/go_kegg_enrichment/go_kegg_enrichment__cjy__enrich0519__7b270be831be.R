library(clusterProfiler)
library(org.Rn.eg.db)  
library(data.table)    
library(dplyr)
library(stringr)
library(ggplot2)

p_cutoff      <- 0.05
log2FC_cutoff <- log2(1.5)  

#-----------------Up-------
dataraw_BA <- fread(file.choose()) %>% 
  filter(pvalue < p_cutoff & log2FoldChange > log2FC_cutoff) %>% 
  mutate(regulate = "Up")
nrow(dataraw_BA)

dataraw_DB <- fread(file.choose()) %>% 
  filter(pvalue < p_cutoff & log2FoldChange < -log2FC_cutoff) %>% 
  mutate(regulate = "Down")
nrow(dataraw_DB)

dataraw_DC <- fread(file.choose()) %>% 
  filter(pvalue < p_cutoff & log2FoldChange < -log2FC_cutoff) %>% 
  mutate(regulate = "Down")
nrow(dataraw_DC)

target_gene <- intersect(dataraw_BA$gene, dataraw_DC$gene)

length(target_gene)


gene_mapped <- bitr(target_gene, 
                    fromType = "ENSEMBL", 
                    toType   = "ENTREZID", 
                    OrgDb    = org.Rn.eg.db)

kegg_res <- enrichKEGG(gene         = gene_mapped$ENTREZID,
                       organism     = "rno",   
                       keyType      = "kegg",  
                       pvalueCutoff = 1,       
                       qvalueCutoff = 1)
res_df <- as.data.frame(kegg_res)

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
ggsave("KEGG_BA(Up)_DC(Down).pdf", p_bar, width = 8, height = 6, device = "pdf")
ggsave("KEGG_BA(Up)_DC(Down).png", p_bar, width = 8, height = 6, dpi = 300)
