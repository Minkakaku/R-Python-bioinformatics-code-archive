library(dplyr)
library(ggplot2)
library(Seurat)
library(monocle)
library(reshape2)
library(ggpubr)
library(ggsci)
library(cowplot)

setwd('G:\\FW2024-302\\售后3')
scdata = dior::read_h5(file='F:/FW2023-629/售后8/adata_final.h5', target.object = 'seurat')
#seuset=readRDS('F:/FW2023-629/售后2/h5_trans.rds')
seuset=scdata
Idents(seuset)='cell_type'
#seuset= subset(seuset, idents = c("Hepatic stellate cells", "Fibroblasts","Hepatocytes","Cholangiocyte"))
seuset= subset(seuset, idents = c("HSCs","LECs" ,"Fibroblasts","Hepatocytes","Cholangiocytes"))
Idents(seuset)='group'
#seuset= subset(seuset, idents = c("Hepatic stellate cells", "Fibroblasts","Hepatocytes","Cholangiocytes"))
scdata=subset(seuset, idents = c("WT-ccl4-6"))# "KO-ccl4-010" "KO-ccl4-488" "WT-ccl4-4"   "WT-ccl4-6"  
scdata=seuset
Idents(scdata)='cell_type'
scdata=FindVariableFeatures(scdata)
genes=VariableFeatures(scdata)
# Pseudotime analysis
data <- as(as.matrix(scdata@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scdata@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# 1.CellDataSet
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

# 2.Estimate size factors and dispersions 
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# 3.Filtering low-quality cells
#HSMM <- detectGenes(HSMM, min_expr = 3 )
HSMM <- detectGenes(HSMM, min_expr = 1.5 )
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
#expressed_genes <- row.names(subset(fData(HSMM)))
# Trajectory step 1: choose genes that define a cell's progress
#expressed_genes=genes
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = "~ cell_type")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

# Trajectory step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, 
                        max_components = 2,
                        norm_method = 'log',
                        num_dim = 20,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)

a1 <- plot_cell_trajectory(HSMM, color_by = "State") + scale_color_npg()
a2 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
a3 <- plot_cell_trajectory(HSMM, color_by = "cell_type") + scale_color_npg()
a4<- plot_cell_trajectory(HSMM, color_by = "group")

p <- (a1+a2)/(a3+a4)

ggsave(a3, filename = 'pseudotime_celltype.png', width = 8, height = 8)
ggsave(a3, filename = 'pseudotime_celltype.pdf', width = 8, height = 8)
p <- plot_cell_trajectory(HSMM, color_by = "cell_type") +
  facet_wrap(~cell_type, nrow = 2)
ggsave(p, filename = 'pseudotime_split_celltype.png', width = 8, height = 8)
ggsave(p, filename = 'pseudotime_split_celltype.pdf', width = 8, height = 8)


p <- plot_cell_trajectory(HSMM, color_by = "group") +
  facet_wrap(~group, nrow = 2)
ggsave(p, filename = 'pseudotime_split_group.png', width = 8, height = 8)
ggsave(p, filename = 'pseudotime_split_group.pdf', width = 8, height = 8)

p <- plot_cell_trajectory(HSMM, color_by = "sample") +
  facet_wrap(~group, nrow = 2)
ggsave(p, filename = 'pseudotime_split_sample.png', width = 8, height = 8)
ggsave(p, filename = 'pseudotime_split_sample.pdf', width = 8, height = 8)
# modify root state
test=HSMM
test=orderCells(test,root_state = 7) 
a1 <- plot_cell_trajectory(test, color_by = "State") + scale_color_npg()
a2 <- plot_cell_trajectory(test, color_by = "Pseudotime")
a3 <- plot_cell_trajectory(test, color_by = "cell_type") + scale_color_npg()
a4<- plot_cell_trajectory(test, color_by = "group")



plot_cell_trajectory(test,color_by = "Pseudotime")

# density plot
df <- pData(test)

p <- ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2() + scale_fill_npg()
ggsave(p, filename = 'pseudotime_density.png', width = 10, height = 8)
ggsave(p, filename = 'pseudotime_density.pdf', width = 10, height = 8)
# gene expression view
HSMM_subset <- test[head(expressed_genes,5),]
p=plot_genes_in_pseudotime(HSMM_subset, color_by = "cell_type") + scale_color_npg()
ggsave(p, filename = 'genes_in_pseudotime.png', width = 10, height = 8)
ggsave(p, filename = 'genes_in_pseudotimey.pdf', width = 10, height = 8)
# save(HSMM, expressed_genes, Tim_gene, file='tmp.Rdata')
save(test, file='tmp.rds')
# heatmap
Tim_diff <- differentialGeneTest(test[expressed_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
Tim_gene <- top_n(Tim_diff, n=20, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(test[Tim_gene,],
                             num_clusters = 3,
                             show_rownames = T,
                             return_heatmap=T)
ggsave(p, filename = 'pseudotime_heatmap.png', width=14, height = 10)
ggsave(p, filename = 'pseudotime_heatmap.pdf', width=14, height = 10)
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.table(clustering, file = 'pseudotime_cluster.xls', sep = '\t')


HSMM=readRDS('G:\\FW2024-302\\result\\pseudotime\\Neutrophils\\HSMM.rds')
p <- (a1+a2)/(a3+a4)

keygenes=c('FOXO3','SOD2','CYCs','GPX4','PPARGC1A')
p1=plot_genes_in_pseudotime(HSMM["FOXO3", ], color_by = "Pseudotime")
p2=plot_genes_in_pseudotime(HSMM["FOXO3", ], color_by = "State")
p3=plot_genes_in_pseudotime(HSMM["FOXO3", ], color_by = "group")
p4=plot_genes_in_pseudotime(HSMM["FOXO3", ], color_by = "celltypeNew")

P=p1|p2|p3|p4
P
ggsave(P, filename = 'FOXO3_pseudotime.png', width = 16, height = 4)
ggsave(P, filename = 'FOXO3_pseudotime.pdf', width = 16, height = 4)


library(dplyr)
library(ggplot2)
library(Seurat)
library(monocle)
library(reshape2)
library(ggpubr)
library(ggsci)
library(cowplot)
library(patchwork)  # 用于拼接多个 ggplot 图

keygenes = c('FOXO3', 'SOD2', 'CYCS', 'GPX4', 
             'PPARGC1A','PPARGC1B','TFAM','TOM20',
             'NDUFA6','AK2','FH','BCL2L11','BBC3',
             'ATG7','MAP1LC3A','PDK4','GLS','GLS2',
             'IL10','SLC7A11')

# 获取 HSMM 中所有基因名
all_genes = rownames(HSMM)

# 遍历 keygenes 生成所有图
for (gene in keygenes) {
  
  # **检查该基因是否在 HSMM 中**
  if (!(gene %in% all_genes)) {
    message(paste("Warning:", gene, "not found in HSMM! Skipping..."))
    next  # 跳过该基因，继续下一个
  }
  
  p1 = plot_genes_in_pseudotime(HSMM[gene, ], color_by = "Pseudotime")
  p2 = plot_genes_in_pseudotime(HSMM[gene, ], color_by = "State")
  p3 = plot_genes_in_pseudotime(HSMM[gene, ], color_by = "group")
  p4 = plot_genes_in_pseudotime(HSMM[gene, ], color_by = "celltypeNew")
  
  # 组合图像
  P = p1 | p2 | p3 | p4
  
  # 显示图像
  print(P)
  
  # 保存图像
  ggsave(P, filename = paste0(gene, "_pseudotime.png"), width = 16, height = 4)
  ggsave(P, filename = paste0(gene, "_pseudotime.pdf"), width = 16, height = 4)
}



"CYCS" %in% rownames(fData(HSMM))

















pdata <- Biobase::pData(HSMM)
write.csv(pdata,'pdata.csv')


# 查看 CellDataSet 中有哪些基因
gene_list <- rownames(fData(HSMM))
head(gene_list)  # 查看前几个基因名

expression_matrix <- exprs(HSMM)
dim(expression_matrix)  # 确认行（基因）和列（细胞）数

# 查看某个基因的表达量
expression_matrix["FOXO3", ]


