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
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/script/Main_mapping_clean.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sub_cell_type %in% c("CCL2 Classical monocyte","DUSP1 Classical monocyte","Non-classical monocyte","VCAN Classical monocyte","cDCs","pDCs")))
# adata <- subset(adata,subset=(sub_cell_type %in% c("Memory B C1","Memory B C2","Naive B")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
levels(adata$cell_type)<-c("Myeloid cells","Myeloid cells","Myeloid cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Megakaryocytes","B cells","B cells","B cells","NK cells","NK cells","Myeloid cells","pDCs","T cells","Myeloid cells")
adata$cell_type <- factor(adata$cell_type,levels=c("pDCs","Myeloid cells","B cells","Plasma cells","Megakaryocytes","NK cells","T cells"))

adata<-subset(adata,group %in% c('N','P'))
adata$group<-droplevels(adata$group)

DimPlot(adata,reduction='umap',group.by='cell_type',raster=FALSE)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/22.Mission/mTOR_heatmap_dotplot/umap.png",width =7,height =6)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251216/22.Mission/mTOR_heatmap_dotplot/umap.pdf",width =7,height =6)

genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/script/pathway.csv",sep = "\t",head=FALSE)
genesets$V2 <- strsplit(genesets$V2, ",")
gene <- setNames(genesets$V2, genesets$V1)

adata$cell_type_group <- factor(
  paste0(adata$cell_type, "_", adata$group),
  levels = paste0(
    rep(levels(adata$cell_type), each = length(levels(adata$group))),
    "_",
    rep(levels(adata$group), times = length(levels(adata$cell_type)))
  )
)
Idents(adata) <- "cell_type_group"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWAE20251216/22.Mission/mTOR_heatmap_dotplot/"
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
gsva.res <- matrix(mat, ncol = 2, byrow = TRUE)
colnames(gsva.res)<-levels(adata$group)
rownames(gsva.res)<-levels(adata$cell_type)
pdf(paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING-heatmap.pdf"),width = 4,height = 8)
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 19,height = 10)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING-heatmap.png"),width = 400,height = 800)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
write.table(gsva.res,paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING_heatmap.csv"), sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)

gsva.res<-gsva.res[c(1,2,3),]
pdf(paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING-heatmap_select3.pdf"),width = 4,height = 4)
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 19,height = 10)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"HALLMARK_PI3K_AKT_MTOR_SIGNALING-heatmap_select3.png"),width = 400,height = 400)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()


mat   <- as.matrix(gsva.df[2,-1])      # df 是你的原始数据框
gsva.res <- matrix(mat, ncol = 2, byrow = TRUE)
colnames(gsva.res)<-levels(adata$group)
rownames(gsva.res)<-levels(adata$cell_type)
pdf(paste0(outdir,"mTOR_pathway-heatmap.pdf"),width = 4,height = 8)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"mTOR_pathway-heatmap.png"),width = 400,height = 800)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
write.table(gsva.res,paste0(outdir,"mTOR_pathway_heatmap.csv"), sep = ",", append = TRUE, row.names = TRUE, col.names = TRUE)

gsva.res<-gsva.res[c(1,2,3),]
pdf(paste0(outdir,"mTOR_pathway_heatmap-heatmap_select3.pdf"),width = 4,height = 4)
# pdf(paste0(outdir,"pathway-heatmap.pdf"),width = 19,height = 10)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()
png(paste0(outdir,"mTOR_pathway_heatmap-heatmap_select3.png"),width = 400,height = 400)
pheatmap::pheatmap(gsva.res, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",fontsize_row  = 14, fontsize_col  = 14,
                   color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100))
dev.off()


bdata<-subset(adata,cell_type %in% c("Myeloid cells","B cells","pDCs"))
bdata$cell_type_group<-factor(bdata$cell_type_group,levels=c("Myeloid cells_N", "Myeloid cells_P", "B cells_N", "B cells_P", "pDCs_N", "pDCs_P"))
p<-DotPlot(bdata, features = gene$mTOR_pathway,group.by = 'cell_type_group',dot.scale = 15,scale =TRUE) + 
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#375093','#e8edf1','#831a21')) #颜色
		ggsave(paste0(outdir,'mTOR_pathway_DotPlot.pdf'),width = 12,height = 7,bg = 'white')
		ggsave(paste0(outdir,'mTOR_pathway_DotPlot.png'),width = 12,height = 7,bg = 'white')

gene<-'MTOR'
df.plot <- FetchData(bdata,vars = c(gene, "cell_type", "group"), assay  = "RNA") %>%
  group_by(cell_type, group) %>%
  summarise(
    pct   = sum(.data[[gene]] > 0) / n() * 100,   # 表达百分比
    mean  = mean(.data[[gene]])                   # 平均表达量
  )
# 2. 气泡图
df.plot$group<-factor(df.plot$group,levels=c("N","P","AID"))
ggplot(df.plot, aes(x = cell_type, y = group)) +
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