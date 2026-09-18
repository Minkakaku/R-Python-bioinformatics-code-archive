
#ssGSEA分析
options(stringsAsFactors = F)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(openxlsx)
library(ComplexHeatmap)
#--------------------------------------------------------------------------
#做gmt文件
cellMarker <- read.xlsx("../../22-AUCell分析/三级淋巴结构的AUCell分析/signature.xlsx", sheet = 2)
cellMarker <- data.frame("TSL" = "TSL", "Marker" = unique(cellMarker$Marker))

GeneID2Cell_list <- tapply(cellMarker[,1],as.factor(cellMarker[,2]),function(x) x)
Cell2GeneID_list <- tapply(cellMarker[,2],as.factor(cellMarker[,1]),function(x) x)

geneSet <- Cell2GeneID_list
gmt_file <- 'Cell2GeneID.gmt'
sink(gmt_file)
for (i in 1:length(geneSet)){
  cat(names(geneSet)[i])
  cat('\tNA\t')
  cat(paste(geneSet[[i]],collapse = '\t'))
  cat('\n')
}
sink()
#--------------------------------------------------------------------------
#读取signatures进行ssGSEA分析
genesets <- getGmt("Cell2GeneID.gmt")

eSet_mRNA <- readRDS("../../17-验证集数据的数据处理/GSE50081/eSet_use.rds")
pData <- read.xlsx("../../16-验证集数据的原始数据下载/GSE50081/Samples.xlsx", sheet = 1)

#如果是log2(TPM+1)或者是DESeq2 normalized counts+vst()，设置为"Gaussian"
#如果是integer count，设置为"Poisson"
resgsva <- gsva(as.matrix(eSet_mRNA), genesets, kcdf="Gaussian", method = "ssgsea")
resgsva <- as.data.frame(resgsva)
write.xlsx(resgsva, file = "ssGSEA_res.xlsx", rowNames = T, overwrite = T)
#-------------------------------------------------------------------------
nmf.group <- read.xlsx("../15-NMF计算Matreotype分型/group.xlsx", sheet = 1)
nmf.group <- nmf.group[!nmf.group$Group %in% "Normal",]

#对功能进行差异分析
library(limma)

gsva_design <- model.matrix(~0+nmf.group$Group)
colnames(gsva_design) <- c("K1", "K2")

gsva_h_fit <- lmFit(resgsva, gsva_design) #exp matrix with samples in columns
gsva_contrast <- makeContrasts(coef1 = K1-K2, levels = gsva_design)
gsva_h_contrastfit <- contrasts.fit(gsva_h_fit, gsva_contrast)
gsva_h_contrastfit <- eBayes(gsva_h_contrastfit)
gsva_h_contrastfitTable <- topTable(gsva_h_contrastfit, n= 'Inf', coef = 1)

gsva_h_contrastfitTable_sig <- gsva_h_contrastfitTable[which(gsva_h_contrastfitTable$adj.P.Val<0.05),]
write.xlsx(gsva_h_contrastfitTable_sig, file = "gsva_h_contrastfitTable_sig.xlsx", asTable = T, overwrite = T)
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
              column_split = nmf.group$Group, row_split = gsva_h_contrastfitTable_sigord$Group, 
              cluster_rows = F,
              show_column_names = F, column_names_gp = gpar(fontsize = 4),
              show_row_names = T, row_names_gp = gpar(fontsize = 8),
              top_annotation = HeatmapAnnotation(border = T, 
                                                 subGroup = as.factor(nmf.group$Group), 
                                                 col = list(subGroup = c("1"="#FB7F0D", "2"="#1A891E"))))

pdf("heatmap.pdf", onefile = F, width = 8, height = 12)
print(ht)
dev.off()
########################################################################

#C2

########################################################################
library(msigdbr)

hs_msigdbgenesets_c2 <- msigdbr(species = "Homo sapiens", category = "C2")

hs_msigdbgenesets_c2 <- split(x = hs_msigdbgenesets_c2$gene_symbol, f = hs_msigdbgenesets_c2$gs_name)
msigdb_c2 <- hs_msigdbgenesets_c2

#Msigdb C2
gsva_msigdb_c2 <- gsva(as.matrix(eSet_tumor), msigdb_c2, kcdf="Gaussian", method = "ssgsea")
#-------------------------------------------------------------------------
#对功能进行差异分析
library(limma)

gsva_design <- model.matrix(~0+nmf.group$Group)
colnames(gsva_design) <- c("K1", "K2")

gsva_c2_fit <- lmFit(gsva_msigdb_c2, gsva_design) #exp matrix with samples in columns
gsva_contrast <- makeContrasts(coef1 = K1-K2, levels = gsva_msigdb_c2)
gsva_c2_contrastfit <- contrasts.fit(gsva_c2_fit, gsva_contrast)
gsva_c2_contrastfit <- eBayes(gsva_c2_contrastfit)
gsva_c2_contrastfitTable <- topTable(gsva_c2_contrastfit, n= 'Inf', coef = 1)

gsva_c2_contrastfitTable_sig <- gsva_c2_contrastfitTable[which(gsva_c2_contrastfitTable$adj.P.Val < 0.05),]
write.xlsx(gsva_c2_contrastfitTable_sig, file = "gsva_c2_contrastfitTable_sig.xlsx", asTable = T, overwrite = T)
#-------------------------------------------------------------------------
#there were a lot of significant hits in the C2 comparison - need to manually choose which to include
ECMrows <- grep("ECM", rownames(gsva_c2_contrastfitTable_sig))
colrows <- grep("COLLAGEN*", rownames(gsva_c2_contrastfitTable_sig))
fibrows <- grep("FIBRONECTIN", rownames(gsva_c2_contrastfitTable_sig))
integrinrows <- grep("INTEGRIN*", rownames(gsva_c2_contrastfitTable_sig))
matrisomerows <- grep("MATRISOM*", rownames(gsva_c2_contrastfitTable_sig))
slrprows <- grep("SLRP*", rownames(gsva_c2_contrastfitTable_sig))
proteoglycanrows <- grep("PROTEOGLYCAN*", rownames(gsva_c2_contrastfitTable_sig))

gsva_c2_contrastfitTable_ecm <- gsva_c2_contrastfitTable_sig[c(ECMrows, colrows, fibrows, 
                                                               integrinrows, matrisomerows, slrprows, 
                                                               proteoglycanrows),]

#these all look significant so don't need to filter further but do need to order them
gsva_c2_contrastfitTable_ecm <- gsva_c2_contrastfitTable_ecm[order(gsva_c2_contrastfitTable_ecm$logFC, decreasing = F),]
gsva_c2_contrastfitTable_ecm

p <- ggplot(gsva_c2_contrastfitTable_ecm, aes(x = 1, y = rownames(gsva_c2_contrastfitTable_ecm))) 
p <- p + geom_point(aes(color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) 
p <- p + scale_y_discrete(limits = rownames(gsva_c2_contrastfitTable_ecm), name = "ECM related signatures", expand=c(0,2))
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
ggsave(plot = p, filename = "ECM_ssGSEA.pdf", dpi = 300, width = 5.5, height = 4.5)
#-------------------------------------------------------------------------
#cisplatin
gsva_c2_contrastfitTable_cisplatin <- gsva_c2_contrastfitTable_sig[grep("*CISPLATIN*", rownames(gsva_c2_contrastfitTable_sig)),]
gsva_c2_contrastfitTable_cisplatin_up <- gsva_c2_contrastfitTable_cisplatin[grep("*_UP",rownames(gsva_c2_contrastfitTable_cisplatin)),]

p <- ggplot(gsva_c2_contrastfitTable_cisplatin_up, aes(x = 1, y = rownames(gsva_c2_contrastfitTable_cisplatin_up))) 
p <- p + geom_point(aes(color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) 
p <- p + scale_y_discrete(limits = rownames(gsva_c2_contrastfitTable_cisplatin_up), name = "Cisplatin Signatures", expand=c(0,2))
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
ggsave(plot = p, filename = "Cisplatin_ssGSEA.pdf", dpi = 300, width = 5.5, height = 4)
########################################################################
#顺铂耐药基因相关的表达
library(ggplot2)
library(ggpubr)

my_comparisons <- list(c("1", "2"))

#CHEK2 YBX1
gene.sel <- "CHEK2"

all(colnames(eSet_tumor)==nmf.group$Samples)
geneSet.plot <- data.frame(Gene = as.numeric(eSet_tumor[gene.sel,]), Group = nmf.group$Group)
colnames(geneSet.plot)[1] <- gene.sel

q <- ggplot(data = geneSet.plot, aes(x = Group, #分组列名
                                     y = get(gene.sel), #连续变量列名
                                     fill = Group))
q <- q + scale_fill_manual(values = c("1"="#FB7F0D", "2"="#1A891E"))
q <- q + geom_violin(alpha = 0.4, position = position_dodge(width = .75), linewidth = 0.8, color="black") 
q <- q + geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7) 
q <- q + geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha = 1) 
q <- q + theme_classic() + ylab(paste0(gene.sel)) + xlab("subGroup") 
q <- q + theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
               axis.ticks = element_line(linewidth=0.2, color="black"),
               axis.ticks.length = unit(0.2, "cm"),
               legend.position = "none",
               axis.title = element_text(size = 12),
               axis.text = element_text(size = 12)) 
q <- q + stat_compare_means(comparisons = my_comparisons, method = "t.test", hide.ns = F, label = "p.signif") 
q <- q + stat_compare_means(method = "t.test", label.y = min(geneSet.plot[,gene.sel])-0.5)
q <- q + theme_bw() + theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
                            axis.text.y = element_text(hjust = 0.5, size = 10, color = "black"),
                            axis.ticks = element_line(linewidth = 0.2, color="black"),
                            axis.ticks.length = unit(0.2, "cm"),
                            legend.position = "none",
                            axis.title = element_text(size = 12, color = "black"),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q
ggsave(plot = q, filename = paste0(gene.sel, "_Group.pdf"), width = 3.5, height = 5)

