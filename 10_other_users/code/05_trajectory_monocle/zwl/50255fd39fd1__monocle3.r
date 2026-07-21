##############################CD8Tmonocle3###############################
#############h5ad转变成seurat对象####################
setwd('./CD8 Tex_monocle')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/FWSC20242216/require2/CD8 Tex/result.h5ad', from="anndata", to="seurat",
                      outFile='CD8Tex_result.rds')

################导入注释信息####################

obj <- readRDS(file = '/mnt/NexData_file/home/zwl/Project/FWSC20242216/require2/CD8 Tex_monocle/CD8Tex_result.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('cell_type','leiden'))

library(monocle3)
f_data = data.frame(gene_short_name = rownames(obj),
                    row.names = rownames(obj))
expression_matrix = as.matrix(LayerData(obj, assay="RNA", layer='counts'))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = obj@meta.data,
                         gene_metadata =  f_data)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, color_cells_by="cell_type")

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds,color_cells_by = "cell_type", show_trajectory_graph = FALSE) + 
  ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
wrap_plots(p1, p2)
#用注释后的cluster
library(do)
obj$cluster1 = 1
cds@clusters@listData[["UMAP"]][["partitions"]] = obj$
  plot_cells(cds,color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by newclusterID")

## 识别轨迹
cds <- learn_graph(cds,close_loop=FALSE, #分支不画圈
                   use_partition = FALSE)# 所有细胞看称整体画分支
#learn_graph_control=list(minimal_branch_len = 6))#保留一些比较短的分支
p = plot_cells(cds, color_cells_by = "leiden",label_groups_by_cluster = FALSE,
               label_leaves = FALSE,group_label_size  = 5, 
               label_branch_points = FALSE)

cds <- order_cells(cds)
p2 = plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE,
                graph_label_size=3)
ggsave(filename = 'monocle3_trajectory_pseudotime_umap.png',p+p2,w=10,h=5)
save(file = 'CD8T_order.Rdata',cds)