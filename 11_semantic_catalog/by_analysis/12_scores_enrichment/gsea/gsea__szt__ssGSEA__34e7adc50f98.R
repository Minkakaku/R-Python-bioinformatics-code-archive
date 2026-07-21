
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

eSet_mRNA <- readRDS("../../17-验证集数据的数据处理/GSE72094/eSet_use.rds")
pData <- read.xlsx("../../16-验证集数据的原始数据下载/GSE72094/Samples.xlsx", sheet = 1)

#如果是log2(TPM+1)或者是DESeq2 normalized counts+vst()，设置为"Gaussian"
#如果是integer count，设置为"Poisson"
resgsva <- gsva(as.matrix(eSet_mRNA), genesets, kcdf="Gaussian", method = "ssgsea")
resgsva <- as.data.frame(resgsva)
write.xlsx(resgsva, file = "ssGSEA_res.xlsx", rowNames = T, overwrite = T)

