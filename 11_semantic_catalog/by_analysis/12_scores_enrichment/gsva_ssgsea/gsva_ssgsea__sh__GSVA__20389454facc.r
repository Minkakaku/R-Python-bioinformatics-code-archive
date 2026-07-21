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

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2512227671-05/01.Mission/PTC/PTC.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
#### 做的分组的GSVA
lipid_geneset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/script/lipid_pathway.csv",sep = ",",head=FALSE)
#转换成list
lipid_geneset$V2 <- strsplit(lipid_geneset$V2, ",")
# 创建一个新的列表，以 V1 为键，V2 为值
gene <- setNames(lipid_geneset$V2, lipid_geneset$V1)

# Idents(adata) <- "annotation"
Idents(adata) <- "group"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/GSVA/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/GSVA/gsva_group.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

gsva.res <- data.frame(gsva.res, check.names = F)
pdf(paste0(outdir,"pathway-heatmap_group.pdf"),width = 10,height = 6)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,"pathway-heatmap_group.png"),width = 700,height = 600)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()


# Idents(adata) <- "annotation"
Idents(adata) <- "cell_type"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/GSVA/gsva_cell_type.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

gsva.res <- data.frame(gsva.res, check.names = F)
pdf(paste0(outdir,"pathway-heatmap_cell_type.pdf"),width = 10,height = 6)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,"pathway-heatmap_cell_type.png"),width = 700,height = 600)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()


for(n in 1:length(levels(adata))){
	cluster_sub<-levels(adata)[[n]]
	adata_sub<-subset(adata,idents=cluster_sub)
	# Idents(adata) <- "annotation"
Idents(adata_sub) <- "group"
expr <- AverageExpression(adata_sub, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,paste0("/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/GSVA/",cluster_sub,"_gsva_group.csv"), sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

gsva.res <- data.frame(gsva.res, check.names = F)
pdf(paste0(outdir,cluster_sub,"_pathway-heatmap_group.pdf"),width = 10,height = 6)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,cluster_sub,"_pathway-heatmap_group.png"),width = 700,height = 600)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
}