library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
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

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/Enrich/genelist_enrich.go.csv",sep = ",",head=TRUE)
# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = "\t",head=TRUE)
geneset<-geneset[,'Description']
# geneset<-toupper(paste0('GOBP_',gsub(" ","_",geneset)))
geneset<-gsub(" ","_",geneset)

GO_geneset <- clusterProfiler::read.gmt("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/c5.go.v2025.1.Hs.symbols.gmt")
result <- GO_geneset %>%
  group_by(term) %>%
  mutate(row = row_number()) %>%  # 修正 row_number 的用法
  pivot_wider(names_from = row, values_from = gene)

combined_pattern <- paste(geneset, collapse = "|")
matched_rows <- grepl(combined_pattern, result$term, ignore.case = TRUE)
result_geneset <- result[matched_rows,]

matched_rows <- which(result$term %in% geneset) 
result_geneset <- result[matched_rows,]
# gene_ids <- select(org.Hs.eg.db, keys = 'GO:0001819', columns = c("ACCNUM","ALIAS","ENSEMBL","ENSEMBLPROT","ENSEMBLTRANS","ENTREZID","ENZYME","EVIDENCE","EVIDENCEALL","GENENAME","GENETYPE","GO","GOALL","IPI","MAP","OMIM","ONTOLOGY","ONTOLOGYALL","PATH","PFAM","PMID","PROSITE","REFSEQ","SYMBOL","UCSCKG","UNIPROT"), keytype = "GO")
write.table(result_geneset,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/pathway.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
write.table(geneset,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/geneset_Description.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)

genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/script/pathway.csv",sep = "\t",head=FALSE)
#转换成list
genesets$V2 <- strsplit(genesets$V2, ",")
# 创建一个新的列表，以 V1 为键，V2 为值
gene <- setNames(genesets$V2, genesets$V1)

Idents(adata) <- "annotation"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/GSVA/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/GSVA/gsva.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

gsva.df$max <- apply(gsva.df[, -1], 1, max)
gsva.df <- gsva.df[gsva.df$Monocyte == gsva.df$max, ]
gsva.df <- head(gsva.df[order(-gsva.df$Monocyte), ], 20)
gsva.res <- gsva.res[rownames(gsva.df),]

pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 10,height = 6)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
png(paste0(outdir,"pathway-heatmap.png"),width = 1000,height = 600)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

write.table(gsva.res,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/15.Mission/GSVA/gsva_heatmap.csv", sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)








adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'HALLMARK_INFLAMMATORY_RESPONSE'
colnames(adata@meta.data)[15]<-'HALLMARK_GLYCOLYSIS'
colnames(adata@meta.data)[16]<-'CYTOKINE_RECEPTOR_INTERACTION'
colnames(adata@meta.data)[14]<-'CYTOKINE'
mydata<-FetchData(adata,vars=c("cell_type","group","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_GLYCOLYSIS","CYTOKINE_RECEPTOR_INTERACTION"))
mean_scores <- mydata %>%
  group_by(cell_type, group) %>%
  summarise(across(everything(), mean))

score_select<-mean_scores$CYTOKINE
score_select <- matrix(score_select, nrow = 6, ncol = 3, byrow = TRUE)
row.names(score_select)<-c("Non-classical monocyte","PTTG1 classical monocyte","ZFP36L1 classical monocyte","S100A8 classical monocyte","Classical monocyte","Dendritic cell")
colnames(score_select)<-c("LAHS","NLAHS","NC")
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/CYTOKINE.pdf",width = 7,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/CYTOKINE.png",width = 480,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

write.table(mean_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/CYTOKINE_AddModuleScore.csv", sep = ",", append = TRUE, row.names = FALSE)