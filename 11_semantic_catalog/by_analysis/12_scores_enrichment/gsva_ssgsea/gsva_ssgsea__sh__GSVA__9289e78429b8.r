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

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/Main_mapping.rds")
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# levels(adata$cell_type)<-c("Myeloid cells","Myeloid cells","Myeloid cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Megakaryocytes","B cells","B cells","B cells","NK cells","NKT cells","T cells","Myeloid cells")
# DimPlot(adata,reduction='umap',group.by='cell_type',raster=FALSE)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/fig4_F_GSVA/umap.png",width =7,height =6)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/fig4_F_GSVA/umap.pdf",width =7,height =6)

genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/script/pathway.csv",sep = "\t",head=FALSE)
genesets$V2 <- strsplit(genesets$V2, ",")
gene <- setNames(genesets$V2, genesets$V1)

adata$raw_cell_type <- factor(adata$raw_cell_type,levels=c("pDCs","Myeloid cells","B cells","Plasma cells","Megakaryocyte","NK cells","NKT cells","T cells"))
adata$cell_type_group <- factor(
  paste0(adata$raw_cell_type, "_", adata$group),
  levels = paste0(
    rep(levels(adata$raw_cell_type), each = length(levels(adata$group))),
    "_",
    rep(levels(adata$group), times = length(levels(adata$raw_cell_type)))
  )
)
Idents(adata) <- "cell_type_group"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/fig4_F_GSVA/"
dir.create(outdir)
# gsva默认开启全部线程计算
gsvapar <- gsvaParam(expr, gene, maxDiff=TRUE)
gsva_es <- gsva(gsvapar)

head(gsva_es)
gsva.res <- gsva_es
# gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.table(gsva.df,paste0(outdir,"gsva.csv"), sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)

mat   <- as.matrix(gsva.df[1,-1])      # df 是你的原始数据框
gsva.res <- matrix(mat, ncol = 3, byrow = TRUE)
colnames(gsva.res)<-levels(adata$group)
rownames(gsva.res)<-levels(adata$raw_cell_type)
pdf(paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING-heatmap.pdf"),width = 6,height = 8)
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 19,height = 10)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING-heatmap.png"),width = 600,height = 800)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
write.table(gsva.res,paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING_heatmap.csv"), sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)

mat   <- as.matrix(gsva.df[2,-1])      # df 是你的原始数据框
gsva.res <- matrix(mat, ncol = 3, byrow = TRUE)
colnames(gsva.res)<-levels(adata$group)
rownames(gsva.res)<-levels(adata$raw_cell_type)
pdf(paste0(outdir,"mTOR_pathway-heatmap.pdf"),width = 6,height = 8)
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 19,height = 10)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"mTOR_pathway-heatmap.png"),width = 600,height = 800)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
write.table(gsva.res,paste0(outdir,"mTOR_pathway.csv"), sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)

gene<-'MTOR'
df.plot <- FetchData(adata,vars = c(gene, "raw_cell_type", "group"), assay  = "RNA") %>%
  group_by(raw_cell_type, group) %>%
  summarise(
    pct   = sum(.data[[gene]] > 0) / n() * 100,   # 表达百分比
    mean  = mean(.data[[gene]])                   # 平均表达量
  )
# 2. 气泡图
df.plot$group<-factor(df.plot$group,levels=c("P_after","P","HD"))
ggplot(df.plot, aes(x = raw_cell_type, y = group)) +
  geom_point(aes(size = pct, color = mean), alpha = 0.8) +
  scale_size_continuous(name   = "Percent Expression",
                        range  = c(1, 10)) +
  scale_color_gradientn(colours = c("navy", "white", "firebrick3"),
                        name    = "Average Expression") +
  theme_classic(base_size = 14) +
  labs(x = "Cell Type", y = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outdir,'DotPlot_MTOR.pdf'),width = 7,height = 7,bg = 'white')
ggsave(paste0(outdir,'DotPlot_MTOR.png'),width = 7,height = 7,bg = 'white')

gene<-'PI3'
df.plot <- FetchData(adata,vars = c(gene, "raw_cell_type", "group"), assay  = "RNA") %>%
  group_by(raw_cell_type, group) %>%
  summarise(
    pct   = sum(.data[[gene]] > 0) / n() * 100,   # 表达百分比
    mean  = mean(.data[[gene]])                   # 平均表达量
  )
# 2. 气泡图
df.plot$group<-factor(df.plot$group,levels=c("P_after","P","HD"))
ggplot(df.plot, aes(x = raw_cell_type, y = group)) +
  geom_point(aes(size = pct, color = mean), alpha = 0.8) +
  scale_size_continuous(name   = "Percent Expression",
                        range  = c(1, 10)) +
  scale_color_gradientn(colours = c("navy", "white", "firebrick3"),
                        name    = "Average Expression") +
  theme_classic(base_size = 14) +
  labs(x = "Cell Type", y = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(outdir,'DotPlot_PI3.pdf'),width = 7,height = 7,bg = 'white')
ggsave(paste0(outdir,'DotPlot_PI3.png'),width = 7,height = 7,bg = 'white')