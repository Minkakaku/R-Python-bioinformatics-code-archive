setwd('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6')
rm(list = ls())
gc()
#.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# 绘图用包
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(getopt)
options(stringsAsFactors = FALSE,warn = -1)
# 可视化用到的文件：############################################################
#           seurat的umap坐标和metadata
#           auc.loom文件
# grn和ctx文件用不到


  
# rm(list = ls())
umap = read.csv('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L2-3IT-MT/umap_coordinates.csv',row.names = 1)
meta = read.csv('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L2-3IT-MT/metadata.csv',row.names = 1)
meta[, 'cell_type' ] = gsub('/','_',meta [,'cell_type'])

loom = open_loom('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L2-3IT-MT/aue_result.loom')
# Read information from loom file:
exprMat = get_dgem(loom)
exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
regulons = regulonsToGeneLists(regulons_incidMat)                     #调控子和对应的基因向量
regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
regulonAucThresholds = get_regulon_thresholds(loom)                   #用于二值化 regulonAUC，判断细胞是否激活某调控子
embeddings = get_embeddings(loom)
close_loom(loom)


library(jsonlite)
write_json(path = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L2-3IT-MT/res/regulon_gene.json',regulons)
regulonAUC
regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),] 
AUC_res = cbind(meta[,c('sample','group','cell_type','leiden')],t(getAUC ( regulonAUC ))[rownames(meta),])
write.csv(file = './L2-3IT-MT/res/Regulon_AUC_metadata.csv',AUC_res)

metadata_long = AUC_res %>% pivot_longer( cols = -c(sample, group, cell_type, leiden),
                                           names_to = "regulon",
                                           values_to = "AUC")
p_value = metadata_long%>%
  group_by(regulon) %>% 
  rstatix::wilcox_test(formula = AUC~group)%>%
  add_significance() %>% add_xy_position(x = 'group')
write.csv(file = './L2-3IT-MT/res/Regulon_AUC_group_pvalue.csv',p_value%>%dplyr::select(-c('groups')))

# 划分cell_type
cellsPerCluster1  =  split ( rownames ( meta ),  meta [, 'sample' ])  
cellsPerCluster <- cellsPerCluster1[sapply(cellsPerCluster1, length) >= 2]  #每个细胞类型至少有两个细胞
# # 计算平均表达式: 
# regulonAuc_byCellType  =  sapply ( cellsPerCluster , 
#                                    function ( cells )  rowMeans (getAUC ( regulonAUC )[, cells ]))%>%
#   data.frame()
# regulonAuc_byCellType$regulon = rownames(regulonAuc_byCellType)
# regulonAucsample_long = regulonAuc_byCellType %>% 
#   pivot_longer( cols = CASP3_1:CT_3,
#                 names_to = "sample",
#                 values_to = "mean_AUC")
# regulonAucsample_long$group = do::Replace(regulonAucsample_long$sample,
#                                           pattern = c("CASP3_1:CASP3","CASP3_2:CASP3","CASP3_3:CASP3",
#                                                       "CT_1:CT","CT_2:CT","CT_3:CT"))
# p_value_sample = regulonAucsample_long%>%
#   group_by(regulon) %>% 
#   rstatix::wilcox_test(formula = mean_AUC~group)%>%
#   add_significance() %>% add_xy_position(x = 'group')

# rm(list = ls())
umap = read.csv('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/L5-6-NP2/umap_coordinates.csv',row.names = 1)
meta = read.csv('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/L5-6-NP2/metadata.csv',row.names = 1)
meta[, 'cell_type' ] = gsub('/','_',meta [,'cell_type'])

loom = open_loom('/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/L5-6-NP2/aue_result.loom')
# Read information from loom file:
exprMat = get_dgem(loom)
exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
regulons = regulonsToGeneLists(regulons_incidMat)                     #调控子和对应的基因向量
regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
regulonAucThresholds = get_regulon_thresholds(loom)                   #用于二值化 regulonAUC，判断细胞是否激活某调控子
embeddings = get_embeddings(loom)
close_loom(loom)

library(jsonlite)
write_json(path = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/scenic/L5-6-NP2/res/regulon_gene.json',regulons)
regulonAUC
regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),] 
AUC_res = cbind(meta[,c('sample','group','cell_type','leiden')],t(getAUC ( regulonAUC ))[rownames(meta),])
write.csv(file = './L5-6-NP2/res/Regulon_AUC_metadata.csv',AUC_res)

metadata_long = AUC_res %>% pivot_longer( cols = -c(sample, group, cell_type, leiden),
                                          names_to = "regulon",
                                          values_to = "AUC")
p_value = metadata_long%>%
  group_by(regulon) %>% 
  rstatix::wilcox_test(formula = AUC~group)%>%
  add_significance() %>% add_xy_position(x = 'group')
write.csv(file = './L5-6-NP2/res/Regulon_AUC_group_pvalue.csv',p_value%>%dplyr::select(-c('groups')))

comen = intersect(p_value%>%dplyr::filter(p<0.05)%>%pull(regulon),
                  p_value1%>%dplyr::filter(p<0.05)%>%pull(regulon))
write.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require6/L2-3IT-MT_L5-6NP2_commen_Regulon_name.csv',comen)
