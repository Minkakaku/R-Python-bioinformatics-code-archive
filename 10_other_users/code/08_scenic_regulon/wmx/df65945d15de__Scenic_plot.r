setwd('/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/')
rm(list = ls())
gc()
#.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(tidyverse)
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

# 读入坐标信息和loom文件,原始细胞数
# L2-3 IT Cmss1+     L2-3 IT-MT 
#          30358           2435 
# L2_3_IT_MT 细胞类型
umap = read.csv("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/umap_coordinates.csv",row.names = 1)
meta = read.csv("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/metadata.csv",row.names = 1)
meta[, 'cell_type' ] = gsub('/','_',meta [,'cell_type'])

loom = open_loom("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/aue_result.loom")
# Read information from loom file:
exprMat = get_dgem(loom)
exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
regulons = regulonsToGeneLists(regulons_incidMat)                     #调控子和对应的基因向量
regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
regulonAucThresholds = get_regulon_thresholds(loom)                   #用于二值化 regulonAUC，判断细胞是否激活某调控子
embeddings = get_embeddings(loom)
close_loom(loom)

# 这里先保存regulon和对应的基因列表
library(jsonlite)
write_json(path = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/regulon_gene.json',regulons) 
regulonAUC
dim(regulonAUC)
regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),]  # query the regulons that are not duplicated 
dim(regulonAUC)
AUC_res = cbind(meta[,c('sample','group','cell_type','leiden')],t(getAUC ( regulonAUC ))[rownames(meta),])
AUC_res %>% head()
write.csv(file = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/Regulon_AUC_metadata.csv',AUC_res)

# 换成长表形式并计算pvalue（wilcox_test）
metadata_long = AUC_res %>% pivot_longer( cols = -c(sample, group, cell_type, leiden),
                                           names_to = "regulon",
                                           values_to = "AUC")

library(dplyr,lib.loc =  "/home/szt/miniforge3/envs/rpy/lib/R/library")
library(rstatix,lib.loc =  "/home/szt/miniforge3/envs/rpy/lib/R/library")
library(ggpubr,lib.loc =  "/home/szt/miniforge3/envs/rpy/lib/R/library")

my_comparisons <- list(
  c("CASP3_P80", "CT_P80"),
  c("CASP3_P42", "CT_P42")
)

p_value <- metadata_long %>%
  group_by(regulon) %>%
  rstatix::wilcox_test(
    formula = AUC ~ group,
    comparisons = my_comparisons
  ) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "group")

p_value

write.csv(file = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/Regulon_AUC_group_pvalue.csv',p_value%>%dplyr::select(-c('groups')))
# 分别获取2组比较间都有差异的regulons，及同一个细胞类型，CPAS80vsCT_P80和CASP3_P42vsCT_P42都有差异的上下调regulons
# 获取的即是在早晚期都存在差异的regulons
p_value %>% filter(p.adj < 0.05) %>% group_by(regulon) %>% summarise(n = n()) %>% filter(n == 2) %>% pull(regulon) -> regulons_diff
write(regulons_diff,file = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_MT/res/regulons_diff.txt',sep = '\n')



# L2_3_IT_Cmss1 细胞类型
umap = read.csv("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/umap_coordinates.csv",row.names = 1)
meta = read.csv("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/metadata.csv",row.names = 1)
meta[, 'cell_type' ] = gsub('/','_',meta [,'cell_type'])

loom = open_loom("/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/aue_result.loom")
# Read information from loom file:
exprMat = get_dgem(loom)
exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")    
regulons = regulonsToGeneLists(regulons_incidMat)                     #调控子和对应的基因向量
regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
regulonAucThresholds = get_regulon_thresholds(loom)                   #用于二值化 regulonAUC，判断细胞是否激活某调控子
embeddings = get_embeddings(loom)
close_loom(loom)

# 这里先保存regulon和对应的基因列表
library(jsonlite)
write_json(path = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/regulon_gene.json',regulons) 
regulonAUC
dim(regulonAUC)
regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),]  # query the regulons that are not duplicated 
dim(regulonAUC)
AUC_res = cbind(meta[,c('sample','group','cell_type','leiden')],t(getAUC ( regulonAUC ))[rownames(meta),])
AUC_res %>% head()
write.csv(file = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/Regulon_AUC_metadata.csv',AUC_res)

# 换成长表形式并计算pvalue（wilcox_test）
metadata_long = AUC_res %>% pivot_longer( cols = -c(sample, group, cell_type, leiden),
                                           names_to = "regulon",
                                           values_to = "AUC")

library(dplyr,lib.loc =  "/home/szt/miniforge3/envs/rpy/lib/R/library")
library(rstatix,lib.loc =  "/home/szt/miniforge3/envs/rpy/lib/R/library")
library(ggpubr,lib.loc =  "/home/szt/miniforge3/envs/rpy/lib/R/library")

my_comparisons <- list(
  c("CASP3_P80", "CT_P80"),
  c("CASP3_P42", "CT_P42")
)

p_value <- metadata_long %>%
  group_by(regulon) %>%
  rstatix::wilcox_test(
    formula = AUC ~ group,
    comparisons = my_comparisons
  ) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "group")

p_value

write.csv(file = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/Regulon_AUC_group_pvalue.csv',p_value%>%dplyr::select(-c('groups')))

# 分别获取2组比较间都有差异的regulons，及同一个细胞类型，CPAS80vsCT_P80和CASP3_P42vsCT_P42都有差异的上下调regulons
# 获取的即是在早晚期都存在差异的regulons
p_value %>% filter(p.adj < 0.05) %>% group_by(regulon) %>% summarise(n = n()) %>% filter(n == 2) %>% pull(regulon) -> regulons_diff
write(regulons_diff,file = '/mnt/NexData_file/home/wmx/Project/SC2604092832-01/shouhou4/scenic/L2_3_IT_Cmss1/res/regulons_diff.txt',sep = '\n')
