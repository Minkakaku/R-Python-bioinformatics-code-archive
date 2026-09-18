
#差异基因的确定
library(ggvenn)

diffSig.TCGA <- read.xlsx("../15-TCGA差异分析/diffSig.xlsx", sheet = 1)
diffSig.GSE68465 <- read.xlsx("../18-GSE68465数据的差异分析/diffSig.xlsx", sheet = 1)

venn_list <- list("TCGA UP" = diffSig.TCGA$row.names[diffSig.TCGA$group %in% "UP"], 
                  "TCGA DOWN" = diffSig.TCGA$row.names[diffSig.TCGA$group %in% "DOWN"],
                  "GSE68465 UP" = diffSig.GSE68465$row.names[diffSig.GSE68465$group %in% "UP"],
                  "GSE68465 DOWN" = diffSig.GSE68465$row.names[diffSig.GSE68465$group %in% "DOWN"])
p <- ggvenn(data = venn_list,
            columns = c("TCGA UP", "GSE68465 UP", "TCGA DOWN", "GSE68465 DOWN"),
            show_percentage = T, stroke_color = "white", 
            fill_color = c("#E98DAFCC", "#E9CF91CC", "#49B1DDCC", "#95C6C6CC"),
            set_name_color = c("#E98DAF", "#E9CF91", "#49B1DD", "#95C6C6"),
            set_name_size = 6)
p
ggsave(plot = p, filename = "Venn.pdf", dpi = 300, width = 6, height = 6)
#----------------------------------------------------------------
library(openxlsx)
library(clusterProfiler)
library(R.utils)
setOption("clusterProfiler.download.method", "auto")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(VennDiagram)

inter <- get.venn.partitions(venn_list)

diffSig <- as.vector(inter[inter$..count.. == 59,]$..values..[[1]])

gene.sel <- "DOWN"
write.xlsx(as.data.frame(diffSig), file = paste0(gene.sel, "_DEG.xlsx"), asTable = T, overwrite = T, rowNames = T)

#symbol ID转换为ENTREZ ID
DEG.gene_symbol <- as.character(diffSig)
DEG.entrez_id <- mapIds(x = org.Hs.eg.db,
                        keys = DEG.gene_symbol,
                        keytype = "SYMBOL",
                        column = "ENTREZID")
DEG.entrez_id <- na.omit(DEG.entrez_id)

#GO分析
erich.go.BP <- enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable= TRUE)
erich.go.BP <- as.data.frame(erich.go.BP)
erich.go.BP <- erich.go.BP[erich.go.BP$pvalue < 0.05,]
write.xlsx(erich.go.BP, paste0("GO.bp.Sig.xlsx"), overwrite = T, asTable = T)

erich.go.CC <- enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "CC",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable= TRUE)
erich.go.CC <- as.data.frame(erich.go.CC)
erich.go.CC <- erich.go.CC[erich.go.CC$pvalue < 0.05,]
write.xlsx(erich.go.CC, paste0("GO.cc.Sig.xlsx"), overwrite = T, asTable = T)

erich.go.MF <- enrichGO(gene = DEG.entrez_id,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable= TRUE)
erich.go.MF <- as.data.frame(erich.go.MF)
erich.go.MF <- erich.go.MF[erich.go.MF$pvalue < 0.05,]
write.xlsx(erich.go.MF, paste0("GO.mf.Sig.xlsx"), overwrite = T, asTable = T)

#KEGG分析
egoPathway <- enrichKEGG(gene = DEG.entrez_id,
                         organism = "hsa",
                         keyType = "kegg",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)
#将entrez重新转换为symbol：
egoPathway2 <- setReadable(egoPathway, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#转换前：
egoPathway@result$geneID[1]
#转换后：
egoPathway2@result$geneID[1]
egoPathway2 <- as.data.frame(egoPathway2)
egoPathway2 <- egoPathway2[egoPathway2$pvalue < 0.05,]
egoPathway2$Description <- gsub(pattern = " \\- Mus musculus \\(house mouse\\)", replacement = "", x = egoPathway2$Description)
write.xlsx(egoPathway2, paste0("KEGG.Sig.xlsx"), overwrite = T, asTable = T)
#----------------------------------------------------------------
#将富集结果整合做长条形图
library(ggplot2)
library(stringr)
library(openxlsx)
library(dplyr)

GO.bp.sel <- read.xlsx("GO.bp.Sig.xlsx", sheet = 2)
GO.bp.sel$Category <- "GO BP"
GO.cc.sel <- read.xlsx("GO.cc.Sig.xlsx", sheet = 2)
GO.cc.sel$Category <- "GO CC"
GO.mf.sel <- read.xlsx("GO.mf.Sig.xlsx", sheet = 2)
GO.mf.sel$Category <- "GO MF"
kegg.sel <- read.xlsx("KEGG.Sig.xlsx", sheet = 2)
kegg.sel$Category <- "KEGG"

enrich <- rbind(GO.bp.sel, GO.cc.sel)
enrich <- rbind(enrich, GO.mf.sel)
enrich <- rbind(enrich, kegg.sel)
##计算Rich Factor（富集因子）：
enrich <- mutate(enrich, RichFactor = Count/as.numeric(sub("/\\d+", "", BgRatio)))
##计算Fold Enrichment（富集倍数）：
enrich$FoldEnrichment <- apply(enrich,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  foldEnrichment <- round(GeneRatio/BgRatio, 2)
  foldEnrichment})

#指定绘图顺序（转换为因子）：
enrich$Description <- factor(paste(enrich$Description), levels = rev(paste(enrich$Description)))
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

#先在自定义主题中隐去y轴文本标签显示
mytheme2 <- mytheme + theme(axis.text.y = element_blank())
enrich$text_x <- rep(0.005, nrow(enrich)) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
p2 <- ggplot(data = enrich, aes(x = RichFactor, y = Description))
p2 <- p2 + geom_bar(aes(fill = Category), stat = "identity", width = 0.8, alpha = 0.7)
p2 <- p2 + scale_fill_manual(values = c("GO BP" = "#067467", "GO CC" = "#F5C35E",
                                        "GO MF" = "#E98DAF", "KEGG" = "#49B1DD"))
p2 <- p2 + labs(x = "RichFactor", y = "pathway", title = "Top GO & KEGG enrichment barplot")
p2 <- p2 + geom_text(aes(x = text_x, label = Description), hjust = 0) #hjust=0，左对齐
p2 <- p2 + theme_classic() + mytheme2
p2
ggsave(p2, filename = "enrich_barplot.pdf", width = 5, height = 14)
ggsave(p2, filename = "enrich_barplot.png", width = 5, height = 14)

