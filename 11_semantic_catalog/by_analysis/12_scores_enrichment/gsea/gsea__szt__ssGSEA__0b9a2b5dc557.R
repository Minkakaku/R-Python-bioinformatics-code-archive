
#ssGSEA分析
options(stringsAsFactors = F)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(openxlsx)
library(ComplexHeatmap)

#读取signatures
genesets <- getGmt("h.all.v2023.1.Hs.symbols.gmt")

eSet_mRNA <- readRDS("../../17-验证集数据的数据处理/GSE68465/eSet_use.rds")
pData <- read.xlsx("../../16-验证集数据的原始数据下载/GSE68465/Samples.xlsx", sheet = 1)

eSet_tumor <- eSet_mRNA[,pData$Sample_geo_accession[pData$Group %in% "Tumor"]]

#如果是log2(TPM+1)或者是DESeq2 normalized counts+vst()，设置为"Gaussian"
#如果是integer count，设置为"Poisson"
resgsva <- gsva(as.matrix(eSet_tumor), genesets, kcdf="Gaussian", method = "ssgsea")
resgsva <- as.data.frame(resgsva)
write.xlsx(resgsva, file = "ssGSEA_res.xlsx", rowNames = T, overwrite = T)
#-------------------------------------------------------------------------
nmf.group <- read.xlsx("../../25-三级淋巴结构下的生存分析/GSE68465/sur.cat.xlsx", sheet = 1)

resgsva <- resgsva[,nmf.group$Samples]

#对功能进行差异分析
library(limma)

gsva_design <- model.matrix(~0+nmf.group$Cell)
colnames(gsva_design) <- levels(factor(nmf.group$Cell))

gsva_h_fit <- lmFit(resgsva, gsva_design) #exp matrix with samples in columns
gsva_contrast <- makeContrasts(coef1 = high-low, levels = gsva_design)
gsva_h_contrastfit <- contrasts.fit(gsva_h_fit, gsva_contrast)
gsva_h_contrastfit <- eBayes(gsva_h_contrastfit)
gsva_h_contrastfitTable <- topTable(gsva_h_contrastfit, n= 'Inf', coef = 1)

gsva_h_contrastfitTable_sig <- gsva_h_contrastfitTable[which(gsva_h_contrastfitTable$adj.P.Val<0.01),]
write.xlsx(gsva_h_contrastfitTable_sig, file = "gsva_h_contrastfitTable_sig.xlsx", asTable = T, overwrite = T, rowNames = T)
#-------------------------------------------------------------------------
library(ggplot2)
#order the significant hallmarks by logFC high to low
gsva_h_contrastfitTable_sigord <- gsva_h_contrastfitTable_sig[order(gsva_h_contrastfitTable_sig$logFC, decreasing = F),]

p <- ggplot(gsva_h_contrastfitTable_sigord, aes(x = 1, y = rownames(gsva_h_contrastfitTable_sigord))) 
p <- p + geom_point(aes(color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) 
p <- p + scale_y_discrete(limits = rownames(gsva_h_contrastfitTable_sigord), name = "Hallmarks", expand=c(0,2))
p <- p + scale_color_gradient(low="#219EBC", high="#E73F74")
p <- p + scale_size(range = c(2, 10)) 
p <- p + theme(axis.text=element_text(size = 8, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1, color="grey" ))
p
ggsave(plot = p, filename = "HALLMARK_ssGSEA.pdf", dpi = 300, width = 5.5, height = 10)
#-------------------------------------------------------------------------
col_fun <- colorRampPalette(c("#4DBBD5FF", "white", "#E64B35FF"))(50)

all(colnames(resgsva)==nmf.group$Samples)
resgsva_plot <- resgsva[rownames(gsva_h_contrastfitTable_sigord),]
gsva_h_contrastfitTable_sigord$Group <- ifelse(gsva_h_contrastfitTable_sigord$logFC > 0, yes = "UP", no = "DOWN")

ht <- Heatmap(t(scale(t(resgsva_plot))), col = col_fun,
              row_labels = gsub(x = rownames(resgsva_plot), pattern = "HALLMARK_", replacement = ""),
              border = T, border_gp = gpar(col = c("black"), lty = 2),
              column_split = nmf.group$Cell, row_split = gsva_h_contrastfitTable_sigord$Group, 
              cluster_rows = F,
              show_column_names = F, column_names_gp = gpar(fontsize = 4),
              show_row_names = T, row_names_gp = gpar(fontsize = 8),
              top_annotation = HeatmapAnnotation(border = T, 
                                                 subGroup = as.factor(nmf.group$Cell), 
                                                 col = list(subGroup = c("high"="#E73F74", "low"="#219EBC"))))

pdf("heatmap.pdf", onefile = F, width = 8, height = 12)
print(ht)
dev.off()
########################################################################

#C2

########################################################################
library(msigdbr)

hs_msigdbgenesets_c2 <- msigdbr(species = "Homo sapiens", category = "C2")

hs_msigdbgenesets_c2 <- split(x = hs_msigdbgenesets_c2$gene_symbol, f = hs_msigdbgenesets_c2$gs_name)

#Msigdb C2
gsva_msigdb_c2 <- gsva(as.matrix(eSet_tumor), hs_msigdbgenesets_c2, kcdf="Gaussian", method = "ssgsea")

gsva_msigdb_c2 <- gsva_msigdb_c2[,nmf.group$Samples]
#-------------------------------------------------------------------------
#对功能进行差异分析
library(limma)

gsva_design_c2 <- model.matrix(~0+nmf.group$Cell)
colnames(gsva_design_c2) <- levels(factor(nmf.group$Cell))

gsva_c2_fit <- lmFit(gsva_msigdb_c2, gsva_design) #exp matrix with samples in columns
gsva_contrast <- makeContrasts(coef1 = high-low, levels = gsva_design_c2)
gsva_c2_contrastfit <- contrasts.fit(gsva_c2_fit, gsva_contrast)
gsva_c2_contrastfit <- eBayes(gsva_c2_contrastfit)
gsva_c2_contrastfitTable <- topTable(gsva_c2_contrastfit, n= 'Inf', coef = 1)

gsva_c2_contrastfitTable_sig <- gsva_c2_contrastfitTable[which(gsva_c2_contrastfitTable$adj.P.Val < 0.01),]
write.xlsx(gsva_c2_contrastfitTable_sig, file = "gsva_c2_contrastfitTable_sig.xlsx", asTable = T, overwrite = T, rowNames = T)
#-------------------------------------------------------------------------
#there were a lot of significant hits in the C2 comparison - need to manually choose which to include
TCRrows <- grep("TCR", rownames(gsva_c2_contrastfitTable_sig))
Tcellrows <- grep("TCELL*", rownames(gsva_c2_contrastfitTable_sig))
CD8rows <- grep("CD8", rownames(gsva_c2_contrastfitTable_sig))
CYTOKINErows <- grep("CYTOKINE*", rownames(gsva_c2_contrastfitTable_sig))
chemorows <- grep("CHEMOKINE*", rownames(gsva_c2_contrastfitTable_sig))
#slrprows <- grep("SLRP*", rownames(gsva_c2_contrastfitTable_sig))
#proteoglycanrows <- grep("PROTEOGLYCAN*", rownames(gsva_c2_contrastfitTable_sig))

gsva_c2_contrastfitTable_CD8 <- gsva_c2_contrastfitTable_sig[unique(c(CD8rows, TCRrows, Tcellrows, CYTOKINErows, chemorows)),]

#these all look significant so don't need to filter further but do need to order them
gsva_c2_contrastfitTable_CD8 <- gsva_c2_contrastfitTable_CD8[order(gsva_c2_contrastfitTable_CD8$logFC, decreasing = F),]
gsva_c2_contrastfitTable_CD8

p <- ggplot(gsva_c2_contrastfitTable_CD8, aes(x = 1, y = rownames(gsva_c2_contrastfitTable_CD8))) 
p <- p + geom_point(aes(color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) 
p <- p + scale_y_discrete(limits = str_wrap(rownames(gsva_c2_contrastfitTable_CD8), width = 40), name = "ECM related signatures", expand=c(0,2))
p <- p + scale_color_gradient(low="#E98DAF", high="#E73F74")
p <- p + scale_size(range = c(4, 10)) 
p <- p + theme(axis.text=element_text(size = 8, face="bold"),
               axis.title=element_text(size=10, face="bold"),
               plot.title=element_text(size = 12, face="bold", hjust= 0.5),
               panel.grid.major.x = element_blank(),
               panel.grid.minor= element_blank(),
               legend.position = "right",
               axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.ticks.y=element_blank(),
               panel.background = element_blank(),
               panel.grid.major.y = element_line(linewidth = 0.1, color="grey" ))
p
ggsave(plot = p, filename = "Tcell_ssGSEA.pdf", dpi = 300, width = 8, height = 9)
#-------------------------------------------------------------------------
#cisplati
gsva_c2_contrastfitTable <- gsva_c2_contrastfitTable_sig[grep("*CISPLATIN*", rownames(gsva_c2_contrastfitTable_sig)),]
gsva_c2_contrastfitTable_up <- gsva_c2_contrastfitTable[grep("*_UP",rownames(gsva_c2_contrastfitTable)),]

p <- ggplot(gsva_c2_contrastfitTable_up, aes(x = 1, y = rownames(gsva_c2_contrastfitTable_up))) 
p <- p + geom_point(aes(color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) 
p <- p + scale_y_discrete(limits = rownames(gsva_c2_contrastfitTable_up), name = "Cisplatin Signatures", expand=c(0,2))
p <- p + scale_color_gradient(low="#219EBC", high="#E73F74")
p <- p + scale_size(range = c(4, 10)) 
p <- p + theme(axis.text=element_text(size = 8, face="bold"),
               axis.title=element_text(size=10, face="bold"),
               plot.title=element_text(size = 12, face="bold", hjust= 0.5),
               panel.grid.major.x = element_blank(),
               panel.grid.minor= element_blank(),
               legend.position = "right",
               axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.ticks.y=element_blank(),
               panel.background = element_blank(),
               panel.grid.major.y = element_line(linewidth = 0.1, color="grey" ))
p
ggsave(plot = p, filename = "Cisplatin_ssGSEA.pdf", dpi = 300, width = 5.5, height = 4)
