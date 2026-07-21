
library(Seurat)
library(dplyr)
library(Matrix)
library(scales)
library(ggplot2)
library(future)
library(copykat)


all_data <- readRDS("/mnt/NexData_file/home/szt/20250605_FWSP20241543_brain_hd/Astrocyte_subcluster/Astrocyte.rds")
Idents(all_data) <- "annotation_sub"

# all_data <- subset(all_data, subset = RCTD_subtype2 != "Unknown")


# table(all_data@meta.data$RCTD_subtype2)
# allCells = names(Idents(all_data))
# allType = levels(Idents(all_data))

# choose_Cells = unlist(lapply(allType, function(x){
#   cgCells = allCells[Idents(all_data)== x ]
#   num = round(length(cgCells)*0.1)
#   cg=sample(cgCells,num)
#   cg
# }))

# all_data <- subset(all_data, cells = choose_Cells)
table(all_data@meta.data$annotation_sub)


counts <- GetAssayData(object = JoinLayers(all_data)[["Spatial.Polygons"]], layer = "counts")


cnv <- copykat(rawmat = counts,
               ngene.chr = 1,
               sam.name = "Luminal",
               n.cores = 32,
               LOW.DR = 0.01,
               UP.DR = 0.3,
               min.gene.per.cell = 20
               )

saveRDS(cnv, "/mnt/NexData_file/home/szt/20250605_FWSP20241543_brain_hd/Astrocyte_subcluster/copykat/cnv.rds")


mallignant <- read.delim("/mnt/NexData_file/home/szt/20250605_FWSP20241543_brain_hd/Astrocyte_subcluster/copykat/Luminal_copykat_prediction.txt")
mallignant <- mallignant[mallignant$copykat.pred != "not.defined",]
rownames(mallignant) <- mallignant[,1]


# 把细胞的良恶性信息加入metadata
all_data <- AddMetaData(all_data, metadata = mallignant)
all_data@meta.data$copykat.pred[is.na(all_data@meta.data$copykat.pred)] <- "undefined"
p1 <- DimPlot(all_data, group.by = "annotation_sub", label = T)
p2 <- DimPlot(all_data, group.by = "copykat.pred") + scale_color_manual(values = c("red", "blue", "gray"))
pc <- p1 + p2
ggsave("/mnt/NexData_file/home/szt/20250605_FWSP20241543_brain_hd/Astrocyte_subcluster/copykat/pred_mallignant.pdf", pc, width = 12, height = 5)

p3 <- DimPlot(all_data, group.by = "copykat.pred", split.by = "copykat.pred") + scale_color_manual(values = c("red", "blue", "gray"))
ggsave("/mnt/NexData_file/home/szt/20250605_FWSP20241543_brain_hd/Astrocyte_subcluster/copykat/pred_mallignant_split.pdf", p3, width = 15, height = 5)


write.table(all_data@meta.data, '/mnt/NexData_file/home/szt/20250605_FWSP20241543_brain_hd/Astrocyte_subcluster/copykat/meta.txt', sep='\t', quote=F, row.names=T)



