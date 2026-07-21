setwd('/mnt/NexData_file/home/zwl/Project/FWAE20250109/require5/')
library(monocle3)
rm(list = ls())
gc()
library(Seurat)
HSPC = readRDS('/mnt/NexData_file/home/zwl/Project/FWAE20250109/rds_h5ad/HSCs/P20121501.diff_PRO.rds')
HSPC_meta = read.csv('/mnt/NexData_file/home/zwl/Project/FWAE20250109/rds_h5ad/HSCs/HSCs_anno_metadata.csv',row.names = 1)
HSPC@meta.data = HSPC_meta
group = unique(paste0(HSPC$type,'_',HSPC$age))

Seurat::DimPlot(HSPC,reduction = 'umap',group.by = c('cell_type','type','age'))

# color cell_type
color_cell = RColorBrewer::brewer.pal(name = 'Set2',n=8)
names(color_cell) <- unique(HSPC$cell_type)
library(monocle3)
library(ggplot2)
for (ni in group) {
  setwd('/mnt/NexData_file/home/zwl/Project/FWAE20250109/require5/')
  if (!dir.exists(ni)) {
    dir.create(ni)}
  
  obj = subset(HSPC,subset = (age == unlist(strsplit(ni[1],'_'))[2])&(type == unlist(strsplit(ni[1],'_'))[1]))
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
  p1 <- plot_cells(cds,color_cells_by = "cell_type", show_trajectory_graph = FALSE,cell_size = 1,
                   group_cells_by = 'cell_type',group_label_size  = 5 )
  p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,cell_size = 1,group_label_size  = 5) + 
    ggtitle("label by partitionID")
  library(patchwork)
  ggsave(filename = paste0('./',ni,'/1.monocle3_cell_type_umap.png'),wrap_plots(p1, p2),w=10,h=5)
  ggsave(filename = paste0('./',ni,'/1.monocle3_cell_type_umap.pdf'),wrap_plots(p1, p2),w=10,h=5)
  
  
  # #用注释后的cluster
  # library(do)
  # obj$cluster1 = 1
  # monocle3::partitions(cds) = obj$cluster1 #cell_type

  
  
  ## 识别轨迹
  cds <- learn_graph(cds,close_loop=FALSE, #分支不画圈
                     use_partition = FALSE,
                     learn_graph_control=list(minimal_branch_len = 6))# 所有细胞看称整体画分支
  #learn_graph_control=list(minimal_branch_len = 6))#保留一些比较短的分支
  p = plot_cells(cds, color_cells_by = "sample",label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,group_label_size  = 5, cell_size = 1,
                 label_branch_points = FALSE)
  ggsave(filename = paste0('./',ni,'/2.monocle3_sample_umap.png'),p,w=5,h=4)
  ggsave(filename = paste0('./',ni,'/2.monocle3_sample_umap.pdf'),p,w=5,h=4)
  
  
  
  ##选择起点
  get_earliest_principal_node  <- function(cds, time_bin="HSC"){
    #function中的参数time_bin选择自己想设置为root的细胞类型，比如我这里选的是#Monocyte  
    cell_ids <- which(colData(cds)[, "cell_type"] == time_bin)    
    closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex 
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])  
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]    
    root_pr_nodes}
  cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
  
  p2 = plot_cells(cds,
                  color_cells_by = "pseudotime",cell_size = 1,
                  label_cell_groups=FALSE,label_leaves=FALSE,
                  label_branch_points=FALSE,graph_label_size=3)
  ggsave(filename = paste0('./',ni,'/3.monocle3_trajectory_pseudotime_umap.png'),p2,w=5,h=4)
  ggsave(filename = paste0('./',ni,'/3.monocle3_trajectory_pseudotime_umap.pdf'),p2,w=5,h=4)
  save(file = paste0('./',ni,'/monocle3_order.Rdata'),cds)
  
  pse = data.frame(monocle3::pseudotime(cds))
  obj$pseudotime = pse[colnames(obj),]
  
  obj$cell_type = factor(obj$cell_type,levels = c('HSC','MPP/LMPP','GMP','CLP/proB','MEP','EBMP','MEP/MKP','MEP/EP'))
  p2 = ggplot(obj@meta.data,aes(x = cell_type,y=pseudotime,fill = cell_type))+
    geom_boxplot()+theme_classic()+
    scale_fill_manual(values = color_cell)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(title = ni,x=NULL)
  ggsave(filename = paste0('./',ni,'/4.monocle3_celltype_pseudotime_umap.png'),p2,w=7,h=4)
  ggsave(filename = paste0('./',ni,'/4.monocle3_celltype_pseudotime_umap.pdf'),p2,w=7,h=4)
  
}

