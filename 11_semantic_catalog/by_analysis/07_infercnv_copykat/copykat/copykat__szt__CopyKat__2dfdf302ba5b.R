script_path <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(script_path)
library(Seurat) 
library(copykat) 
library(tidyverse) 
library(tidydr)
library(openxlsx)
rm(list = ls()) 


my_theme <- theme_void()+theme(plot.margin = margin(5.5, 15, 5.5, 5.5))+
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, units = "inches"), ends = "last", type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 15, color = "black"))
#----------------------------------------------------------------------
#准备数据
cell_anno <- readRDS("../7-细胞注释/alldata_anno.rds") 
# scRNA <- UpdateSeuratObject(scRNA) 

#做CNV输入的是原始的count数据
#有文献提示，使用count数据和log后的data数据来做分析没有什么差别（copyKAT和inferCNV都是）
#但还是建议使用count数据，因为它自己会进行log标准化
#因为运行时间过长，在服务器上运行
copykat.test <- copykat(rawmat=as.matrix(cell_anno@assays$RNA@counts),
                        id.type="S", 
                        ngene.chr=5, 
                        win.size=25, 
                        KS.cut=0.1, 
                        sam.name="LUAD",
                        distance="euclidean")
pred.test <- data.frame(copykat.test$prediction)

save.image("copykat.RData")
##################################################################
load("copykat.RData")

pred.test <- read.delim("LUAD_copykat_prediction.txt", header = T, sep = "\t", check.names = F, stringsAsFactors = F)
  
dim(cell_anno)
cell_anno@meta.data$copykat <- pred.test[match(colnames(cell_anno), pred.test$cell.names), 'copykat.pred']
cell_anno$copykat <- ifelse(test = cell_anno$copykat %in% "not.defined", yes = "not.defined", 
                            no = ifelse(test = cell_anno$copykat %in% "diploid", yes = "aneuploid", no = "diploid"))
table(cell_anno@meta.data$copykat)

my_cols <- c("diploid" = "#CFCFCFCC", "aneuploid" = "#E73F74CC", "not.defined" = "#F5F5F5CC")
DimPlot(cell_anno, reduction = "tsne", #label=T, label.box = T,
        group.by = 'copykat', cols = my_cols, pt.size = 0.75, 
        shuffle = T, raster = F#, order = c("aneuploid", "diploid", "not.defined")
        
        ) + my_theme
ggsave(filename = paste0("copykat",'_tSNE.png'), width = 8.5, height = 8)
ggsave(filename = paste0("copykat",'_tSNE.pdf'), width = 8.5, height = 8)

DimPlot(cell_anno, reduction = "umap", #label=T, label.box = T,
        group.by = 'copykat', cols = my_cols, pt.size = 0.75, 
        shuffle = T, raster = F#, order = c("aneuploid", "diploid", "not.defined")
        ) + my_theme
ggsave(filename = paste0("copykat",'_umap.png'), width = 8.5, height = 8)
ggsave(filename = paste0("copykat",'_umap.pdf'), width = 8.5, height = 8)

write.xlsx(cell_anno@meta.data, file = "copykat_metadata.xlsx", asTable = T, overwrite = T, rowNames = T)
##################################################################
#绘制比例
library(GGally, quietly = TRUE)

#prop <- as.data.frame(table(cell_anno@meta.data$scissor, cell_anno@meta.data$copykat))
prop <- as.data.frame(table(cell_anno@meta.data$S.Score, cell_anno@meta.data$copykat))
prop$Var2 <- factor(prop$Var2, levels = rev(c("not.defined", "diploid", "aneuploid")))
table(cell_anno@meta.data$copykat)

p <- ggplot(data = prop)
p <- p + aes(x = Var1, weight = Freq, fill = Var2, by = Var1)
p <- p + geom_bar(position = "fill")
p <- p + geom_text(stat = "prop", position = position_fill(.5))
p <- p + scale_fill_manual(values = my_cols)
p <- p + coord_flip() + theme_light()
p
ggsave("prop.pdf", plot = p, width = 15, height = 5)
##################################################################
saveRDS(object = cell_anno, file = "cell_copyKAT.rds")



