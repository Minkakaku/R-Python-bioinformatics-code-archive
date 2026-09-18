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
library(RColorBrewer)

# adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/Plasma cells/Plasma.rds")
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/mapping_raw/Main_mapping.rds")
# adata <- JoinLayers(adata)
adata <- subset(adata,subset=(sub_cell_type %in% c("CCL2 Classical monocyte","DUSP1 Classical monocyte","Non-classical monocyte","VCAN Classical monocyte","cDCs","pDCs")))
# adata <- subset(adata,subset=(sub_cell_type %in% c("Memory B C1","Memory B C2","Naive B")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/pathway.csv",sep = "\t",head=FALSE)
genesets$V1 <- substring(genesets$V1, 6) |>       # 去掉前 5 个字符
              str_to_title()
#转换成list
genesets$V2 <- strsplit(genesets$V2, ",")
# 创建一个新的列表，以 V1 为键，V2 为值
gene <- setNames(genesets$V2, genesets$V1)
# Idents(adata) <- "annotation"
# Idents(adata) <- "group"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig1_D_GSVA/Myeloid/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,paste0(outdir,"gsva.csv"), sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

# gsva.df$max <- apply(gsva.df[, -1], 1, max)
# gsva.df <- gsva.df[gsva.df$LAHS == gsva.df$max, ]
# gsva.df <- head(gsva.df[order(-gsva.df$LAHS), ], 20)
# gsva.res <- gsva.res[rownames(gsva.df),]

mat   <- as.matrix(gsva.df[,-1])      # df 是你的原始数据框
picked <- c()
for (j in seq_len(ncol(mat))) {
  remain <- setdiff(rownames(mat), picked)
  colVec <- mat[remain, j]                     # 关键：别漏 drop
  picked <- c(picked, remain[order(-colVec)][1:3])
}
gsva.res <- mat[picked, , drop = FALSE]
gsva.res <- gsva.res[c(3, 1:2, 4:nrow(gsva.res)), , drop = FALSE]

pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 8,height = 7)
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 19,height = 10)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"pathway-heatmap.png"),width = 800,height = 700)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()

write.table(gsva.res,paste0(outdir,"gsva_heatmap.csv"), sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)

# sed 's/,/\t/' pathway.csv > out.csv
genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/script/out.csv",sep = "\t",head=FALSE)
genesets$V2 <- strsplit(genesets$V2, ",")
gene <- setNames(genesets$V2, genesets$V1)

expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

# outdir<-"/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig1_D_GSVA/"
# dir.create(outdir)
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,paste0(outdir,"gsva_KEGG.csv"), sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

mat   <- as.matrix(gsva.df[,-1])      # df 是你的原始数据框
picked <- c()
for (j in seq_len(ncol(mat))) {
  remain <- setdiff(rownames(mat), picked)
  colVec <- mat[remain, j]                     # 关键：别漏 drop
  picked <- c(picked, remain[order(-colVec)][1:3])
}
gsva.res <- mat[picked, , drop = FALSE]

# pdf(paste0(outdir,"pathway-heatmap_KEGG.pdf"),width = 18,height = 10)
pdf(paste0(outdir,"pathway-heatmap_KEGG.pdf"),width = 9,height = 8)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"pathway-heatmap_KEGG.png"),width = 800,height = 700)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()

write.table(gsva.res,paste0(outdir,"gsva_heatmap_KEGG.csv"), sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)