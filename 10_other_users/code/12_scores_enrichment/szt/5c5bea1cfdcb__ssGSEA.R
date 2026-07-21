
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

eSet_mRNA <- read.delim("../../14-TCGA LUAD数据处理/eSet_mRNA.txt", header = T, sep = "\t", 
                        row.names = 1, check.names = F, stringsAsFactors = F)
pData <- read.xlsx("../../14-TCGA LUAD数据处理/LUAD_Group.xlsx", sheet = 1)

#如果是log2(TPM+1)或者是DESeq2 normalized counts+vst()，设置为"Gaussian"
#如果是integer count，设置为"Poisson"
resgsva <- gsva(as.matrix(eSet_mRNA), genesets, kcdf="Gaussian", method = "ssgsea")
resgsva <- as.data.frame(resgsva)
write.xlsx(resgsva, file = "ssGSEA_res.xlsx", rowNames = T, overwrite = T)
#--------------------------------------------------------------------------
#分组间的差异分析
library(openxlsx)
library(reshape2)
library(ggpubr)

cibersort <- read.xlsx("ssGSEA_res.xlsx", sheet = 1)
plot.info <- melt(cibersort[,c("Samples", "TSL")])
colnames(plot.info) <- c("Samples", "CellType", "Composition")
plot.info$CellType <- "TLS markers in TCGA Data"

immune_group <- read.xlsx("../../14-TCGA LUAD数据处理/LUAD_Group.xlsx", sheet = 1)

plot.info <- merge(plot.info, immune_group[,c("Samples", "Group")], by.x = "Samples", by.y = "Samples")
plot.info <- plot.info[!is.na(plot.info$Group),]
plot.info$Group <- as.factor(plot.info$Group)

my_comparisons <- list(c("Tumor", "Normal"))

p <- ggplot(data = plot.info, aes(x = Group, y = Composition, fill = Group, group = Group)) 
p <- p + scale_fill_manual(values = c("Tumor"="#E73F74", "Normal" = "#219EBC")) 

p <- p + geom_violin(alpha = 0.4, position = position_dodge(width = .75), linewidth = 0.5, color="black")
p <- p + geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.5, alpha = 0.7)
p <- p + geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha = 1)

library(ggh4x)
p <- p + facet_wrap2(.~CellType, scales = "free_y", nrow = 3, strip = strip_themed(background_x = elem_list_rect(fill = "#FEB705"))) #分面，此命令可以给分面bar改变颜色

#添加统计结果
p <- p + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", vjust = 0.5, size = 5)
p <- p + stat_compare_means(method = "t.test", label.y = min(plot.info$Composition)-0.05)

p <- p + ggtitle("") + xlab("") + ylab("ssGSEA score of TLS markers") 

p <- p + theme_bw() + theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
                            axis.text.y = element_text(hjust = 0.5, size = 10, color = "black"),
                            axis.ticks = element_line(linewidth = 0.2, color="black"),
                            axis.ticks.length = unit(0.2, "cm"),
                            legend.position = "none",
                            axis.title = element_text(size = 12, color = "black"),
                            plot.title = element_text(size = 12, color = "black"),
                            strip.text = element_text(size = 12, color = "white"),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave(plot = p, "TLS of tcga.pdf", width = 3.5, height = 5, dpi = 300)
ggsave(plot = p, "TLS of tcga.tif", width = 3.5, height = 5, dpi = 300, device = "tiff")

