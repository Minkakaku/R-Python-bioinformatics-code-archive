setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require3/')
rm(list = ls())
gc()
library(clusterProfiler)
library(org.Mm.eg.db)
DEG = read.csv(file = '../require2/top10clonotypes_others/all_cluster_marker.xls',sep = '\t')
DEG = read.csv(file = './top10clontype_others/EDA1/all_cluster_marker.xls',sep = '\t')
DEG <- DEG %>% mutate(logfoldchanges = ifelse(cluster == "nan", -logfoldchanges, logfoldchanges))

eg = bitr(DEG$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

DEG_ALL = merge(DEG,eg,by.x='Symbol',by.y='SYMBOL')%>%arrange(.,desc(logfoldchanges))
geneList <- DEG_ALL$logfoldchanges
names(geneList) <- DEG_ALL$ENTREZID


ego_BP <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego_BP <- setReadable(ego_BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
res = ego_BP@result
write.csv(file = paste0('./top10clontype_others/GSEA/1.GO_BP_pathway_res.csv'),res)

ego_CC <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego_CC <- setReadable(ego_CC, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
res = ego_CC@result
write.csv(file = paste0('./top10clontype_others/GSEA/1.GO_CC_pathway_res.csv'),res)

ego_MF <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "MF",
              pvalueCutoff = 0.05,
              verbose      = FALSE)
ego_MF <- setReadable(ego_MF, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
res = ego_MF@result
write.csv(file = paste0('./top10clontype_others/GSEA/1.GO_MF_pathway_res.csv'),res)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mouse',
               pvalueCutoff = 0.5,
               verbose      = FALSE)
kk2 <- setReadable(kk2, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(file = paste0('./top10clontype_others/GSEA/1.KEGG_pathway_res.csv'),kk2@result)


library(DOSE)
library(enrichplot)
library(scales)
library(RColorBrewer)
p = gseaplot2(ego_BP, geneSetID = 1:5)
ggsave(filename = './top10clontype_others/GSEA/2.GO_BPpathway_group_top10_vs_other.png',w=6,h=5,dpi = 200,p)
ggsave(filename = './top10clontype_others/GSEA/2.GO_BPpathway_group_top10_vs_other.png',w=6,h=5,dpi = 200,p)
p2 = gseaplot2(ego_CC, geneSetID = 1:2)
ggsave(filename = './top10clontype_others/GSEA/2.GO_CCpathway_group_top10_vs_other.png',w=6,h=4,dpi = 200,p2)
ggsave(filename = './top10clontype_others/GSEA/2.GO_CCpathway_group_top10_vs_other.png',w=6,h=4,dpi = 200,p2)
# p3 = gseaplot2(ego_MF, geneSetID = 1:5)
# ggsave(filename = './top10clontype_others/GSEA/2.GO_MFpathway_group_top10_vs_other.png',w=6,h=4,dpi = 200,p3)
# ggsave(filename = './top10clontype_others/GSEA/2.GO_MFpathway_group_top10_vs_other.png',w=6,h=4,dpi = 200,p3)
p4 = gseaplot2(kk2, geneSetID = 1)
ggsave(filename = './top10clontype_others/GSEA/2.KEGG_pathway_group_top10_vs_other.png',w=6,h=4,dpi = 200,p4)
ggsave(filename = './top10clontype_others/GSEA/2.KEGG_pathway_group_top10_vs_other.png',w=6,h=4,dpi = 200,p4)


# cellchat data##################################################
library(sceasy)
library(Seurat)
library(reticulate)
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require3/cellchat')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require3/cellchat/anno_T_Mye.h5ad', from="anndata", to="seurat",
                      outFile='anno_T_Mye.rds')
obj <- readRDS(file = './anno_T_Mye.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type'))
