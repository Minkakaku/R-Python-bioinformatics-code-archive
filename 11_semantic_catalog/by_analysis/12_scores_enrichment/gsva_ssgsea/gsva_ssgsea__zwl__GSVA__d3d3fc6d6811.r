library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)
# library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyr)
library(stringr)
library(GSVA)
library(patchwork)
library(limma)

adata<-readRDS("/mnt/NexData_file/home/zwl/Project/AE2507033447-08/require3/mapping_raw/Main_mapping.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/Enrich/genelist_enrich.go.csv",sep = ",",head=TRUE)
# # geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = "\t",head=TRUE)
# geneset<-geneset[,'Description']
# # geneset<-toupper(paste0('GOBP_',gsub(" ","_",geneset)))
# geneset<-gsub(" ","_",geneset)

# GO_geneset <- clusterProfiler::read.gmt("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/c5.go.v2025.1.Hs.symbols.gmt")
# result <- GO_geneset %>%
#   group_by(term) %>%
#   mutate(row = row_number()) %>%  # 修正 row_number 的用法
#   pivot_wider(names_from = row, values_from = gene)

# combined_pattern <- paste(geneset, collapse = "|")
# matched_rows <- grepl(combined_pattern, result$term, ignore.case = TRUE)
# result_geneset <- result[matched_rows,]

# matched_rows <- which(result$term %in% geneset) 
# result_geneset <- result[matched_rows,]
# # gene_ids <- select(org.Hs.eg.db, keys = 'GO:0001819', columns = c("ACCNUM","ALIAS","ENSEMBL","ENSEMBLPROT","ENSEMBLTRANS","ENTREZID","ENZYME","EVIDENCE","EVIDENCEALL","GENENAME","GENETYPE","GO","GOALL","IPI","MAP","OMIM","ONTOLOGY","ONTOLOGYALL","PATH","PFAM","PMID","PROSITE","REFSEQ","SYMBOL","UCSCKG","UNIPROT"), keytype = "GO")
# write.table(result_geneset,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/pathway.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
# write.table(geneset,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/geneset_Description.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)

#### 做的分组的GSVA
# genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/pathway.csv",sep = "\t",head=FALSE)
# #转换成list
# genesets$V2 <- strsplit(genesets$V2, ",")
# # 创建一个新的列表，以 V1 为键，V2 为值
# gene <- setNames(genesets$V2, genesets$V1)

# # Idents(adata) <- "annotation"
# Idents(adata) <- "group"
# expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
# expr <- expr[rowSums(expr)>0,] #选取非零基因
# expr <- as.matrix(expr)
# head(expr)

# outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_group/"
# dir.create(outdir)
# # gsva默认开启全部线程计算
# gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
# gsva_es <- gsva(gsvapar)

# head(gsva_es)
# gsva.res <- gsva_es
# # gsva.res <- gsva(expr, genesets, method="gsva")
# gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
# write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_group/gsva.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

# gsva.df$max <- apply(gsva.df[, -1], 1, max)
# gsva.df <- gsva.df[gsva.df$LAHS == gsva.df$max, ]
# gsva.df <- head(gsva.df[order(-gsva.df$LAHS), ], 20)
# gsva.res <- gsva.res[rownames(gsva.df),]

# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 10,height = 6)
# pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
#                    scale = "none",angle_col = "45",
#                    color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# dev.off()
# png(paste0(outdir,"pathway-heatmap.png"),width = 700,height = 600)
# pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
#                    scale = "none",angle_col = "45",
#                    color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# dev.off()

# write.table(gsva.res,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_group/gsva_heatmap.csv", sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)



genesets<-read.table("/mnt/NexData_file/home/zwl/Project/AE2507033447-08/require3/k4_top_genes_by_score.csv",sep = ",",head=T)
#转换成list
gene_list <- split(genesets$gene, genesets$module)  
# genesets$gene <- strsplit(genesets$gene, ",")
# 创建一个新的列表，以 V1 为键，V2 为值
# gene <- setNames(genesets$gene, genesets$V1)
# adata <- subset(adata,subset=(annotation %in% c("Monocyte")))
# Idents(adata) <- "annotation"
Idents(adata) <- "group"
expr <- AverageExpression(adata,group.by = 'sample', assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/zwl/Project/AE2507033447-08/require3/GSVA/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene_list, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,"/mnt/NexData_file/home/zwl/Project/AE2507033447-08/require3/GSVA/gsva.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

attr(gsva.res, "geneSets") <- NULL
colnames(gsva.res) <- colnames(gsva.res) %>%
  gsub("^g", "", .) %>%      # 去掉开头的 g
  gsub("-", "_", .)   
 
annot_df <- unique(score_res[,c('sample','group')]) %>%
  mutate(
    group = factor(group, levels = c("non_pcr_A","non_pcr_B","pcr_A","pcr_B"))
  ) %>%
  arrange(group) %>%  # 按分组排序
  as.data.frame()

rownames(annot_df) <- annot_df$sample
annot_df <- annot_df["group"] # 只留分组

# 2. 矩阵按排序后的样本顺序对齐
gsva.res <- gsva.res[, rownames(annot_df)]

pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 15,height = 6)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=T,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50) )
dev.off()
png(paste0(outdir,"pathway-heatmap.png"),width = 1500,height = 500)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=T,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   margins = c(5,12) )
dev.off()

gsva_df <- as.data.frame(t(gsva.res))  # 转置：行为样本，列为M1-M4

# 合并分组信息
final_df <- cbind(
  sample = rownames(gsva_df),  # 样本名
  annot_df,                    # 分组 group
  gsva_df                      # M1 M2 M3 M4
)

# 重设行名
rownames(final_df) <- NULL

for (i in colnames(score_res)[3:6]) {
  
  # 动态公式：用 !!sym(i) 把字符串变成列变量
  p_value <- df_plot %>% 
    rstatix::wilcox_test(as.formula(paste(i, "~ group"))) %>%
    add_significance() %>% 
    add_xy_position()
  
  # 绘图（y轴、标题全部自动变）
  p <- ggplot(df_plot, aes(x = group, y = !!sym(i), color = group)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2) +
    stat_pvalue_manual(p_value, label = "p.adj.signif", hide.ns = TRUE) +
    labs(
      x = "Group",
      y = paste0(i),  # 自动变 M1/M2/M3/M4
      color = "Group"
    ) +
    theme_classic()
  
  # 自动保存不同文件名（关键！每张图不覆盖）
  ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/AE2507033447-08/require3/GSVA/1.',i,'_group.png'),
         p, width=6, height=4, dpi=150)
  ggsave(filename = paste0('/mnt/NexData_file/home/zwl/Project/AE2507033447-08/require3/GSVA/1.',i,'_group.pdf'),
         p, width=6, height=4, dpi=150)
}






# write.table(gsva.res,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/gsva_heatmap.csv", sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)
# 
# 
# group="NC"
# adata <- subset(adata_raw,subset=(group %in% c("NC")))
# Idents(adata) <- "annotation"
# expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
# expr <- expr[rowSums(expr)>0,] #选取非零基因
# expr <- as.matrix(expr)
# head(expr)
# 
# outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/NC/"
# dir.create(outdir)
# # gsva默认开启全部线程计算
# gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
# gsva_es <- gsva(gsvapar)
# 
# head(gsva_es)
# gsva.res <- gsva_es
# # gsva.res <- gsva(expr, genesets, method="gsva")
# gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
# write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/NC/gsva.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
# 
# gsva.df$max <- apply(gsva.df[, -1], 1, max)
# gsva.df <- gsva.df[gsva.df$Monocyte == gsva.df$max, ]
# gsva.df <- head(gsva.df[order(-gsva.df$Monocyte), ], 20)
# gsva.res <- gsva.res[rownames(gsva.df),]
# 
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 13,height = 7)
# pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
#                    scale = "none",angle_col = "45",
#                    color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# dev.off()
# png(paste0(outdir,"pathway-heatmap.png"),width = 1000,height = 700)
# pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
#                    scale = "none",angle_col = "45",
#                    color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
# dev.off()
# 
# write.table(gsva.res,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/NC/gsva_heatmap.csv", sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)