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

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/all_data_sc_doublet.rds")
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



genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/pathway.csv",sep = "\t",head=FALSE)
#转换成list
genesets$V2 <- strsplit(genesets$V2, ",")
# 创建一个新的列表，以 V1 为键，V2 为值
gene <- setNames(genesets$V2, genesets$V1)
adata <- subset(adata,subset=(annotation %in% c("Monocyte")))
# Idents(adata) <- "annotation"
Idents(adata) <- "group"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/gsva.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

gsva.df$max <- apply(gsva.df[, -1], 1, max)
gsva.df <- gsva.df[gsva.df$LAHS == gsva.df$max, ]
gsva.df <- head(gsva.df[order(-gsva.df$LAHS), ], 20)
gsva.res <- gsva.res[rownames(gsva.df),]

pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 10,height = 6)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,"pathway-heatmap.png"),width = 700,height = 600)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

write.table(gsva.res,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/gsva_heatmap.csv", sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)


group="NC"
adata <- subset(adata_raw,subset=(group %in% c("NC")))
Idents(adata) <- "annotation"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/NC/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/NC/gsva.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

gsva.df$max <- apply(gsva.df[, -1], 1, max)
gsva.df <- gsva.df[gsva.df$Monocyte == gsva.df$max, ]
gsva.df <- head(gsva.df[order(-gsva.df$Monocyte), ], 20)
gsva.res <- gsva.res[rownames(gsva.df),]

pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 13,height = 7)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,"pathway-heatmap.png"),width = 1000,height = 700)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

write.table(gsva.res,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA_Monocyte/NC/gsva_heatmap.csv", sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)