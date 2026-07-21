# conda activate "/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle3/"
# install.packages("devtools")
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')
# conda install bioconductor-rhdf5filtersn
# conda install bioconductor-lme4

library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)
library(Seurat)

# library(dior)
# adata = dior::read_h5(file='G:/FW2024-629/售后19/monocle3/Fibroblast_0.2/result.h5')
# clo=colnames(adata@meta.data)
# Idents(adata) = adata$leiden
# scRNA =adata
scRNA =readRDS("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/CD8T/CD8T.rds")
# scRNA<-subset(scRNA,cell_type %in% c("Immature Mac C1","Immature Mac C2","Pro-inflammatory Mac","Pro-repair Mac","Proliferative Mac","Remodeling Mac-Svep1"))
scRNA[["RNA"]] <- as(object = scRNA[["RNA"]], Class = "Assay") #V5版本替换
DefaultAssay(scRNA) <- 'RNA'
#meta=read.csv('G:\\FW2024-348\\FW2024-348_10\\售后4\\monocle3\\all_features_BP.xls',sep = '\t',header = 1,row.names = 1)
#scRNA@meta.data=meta
# 构建cds矩阵
#scRNA$celltype1=gsub('proliferating Pre B','B cells',scRNA$cell_type)
#scRNA$celltype1=gsub('Pre_immature B','B cells',scRNA$celltype1)
setwd("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/monocle3/")
rds2cds <- function(seuset,umap2dfile){
  data = seuset@assays$RNA@counts
  gene_metadata = data.frame(gene_short_name = rownames(data),check.names = F)
  rownames(gene_metadata) = gene_metadata[,1]
  cell_metadata = as.data.frame(cbind(seuset@meta.data,Clusters = as.character(seuset@active.ident)))
  cds <- new_cell_data_set(
    data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  cds <- estimate_size_factors(cds)
  ########################## 添加UMAP降维坐标 ##########################
  #首先统一umap文件顺序
  if(is.null(umap2dfile)){
    umap1 = as.data.frame(seuset@reductions$umap@cell.embeddings)
  }else{
    umap1 = read.table(umap2dfile,sep = "\t",row.names = 1,header= T)
    umap1 = umap1[is.finite(match(rownames(umap1),colnames(cds))),]
  }
  reducedDims(cds)$UMAP <- umap1
  ################ 走个流程得到cds的降维格式
  cds <- cluster_cells(cds,k = 30,reduction_method = "UMAP",cluster_method = "louvain")
  # 把rds的clsuter替换进去 但是不替换partition(monocle3给umap区分的区域 即曲线条数)
  # 这样的好处是learn_graph步骤会在每个曲线区域根据cluster分布去生成轨迹 
  # 这里明确 直接生成的names(cds@clusters$UMAP$clusters) 和 colnames(cds) 和 rownames(cds@colData)是相同的 所以这里同步替换掉
  names(cds@colData$Clusters) = rownames(cds@colData)
  cds@colData$Clusters = as.factor(cds@colData$Clusters)
  cds@clusters$UMAP$clusters = cds@colData$Clusters
  return(cds)
}

cds = rds2cds(scRNA,NULL)


# cds@clusters@listData[["UMAP"]][["partitions"]] = scRNA$RNA3_snn_res.0.4
gg = plot_cells(cds,color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by newclusterID")
ggsave("Trajectory_partition_Monocle3_Plot.png",gg,width = 7,height = 7)
ggsave("Trajectory_partition_Monocle3_Plot.pdf",gg,width = 7,height = 7)
# 最小分枝长度
minimal_branch_len1 = 25
# 生成轨迹最近邻K
nn_k1 = 30
param = list()
param[["minimal_branch_len"]] = minimal_branch_len1
param[["nn.k"]] = nn_k1
use_partition1 = T 

#cds=readRDS('G:\\FW2024-348\\FW2024-348_10\\售后4\\monocle3\\result1\\allClusterCDS.rds')
cds <- learn_graph(cds,close_loop = F,learn_graph_control = param,verbose = T,use_partition =FALSE)
# dir.create('TrajectoryResult')
gg = plot_cells(cds,
                reduction_method = "UMAP",
                color_cells_by = "Clusters",
                show_trajectory_graph = F,
                label_cell_groups=T,#在图上表示cluster等标签
                label_leaves=F,#不显示曲线演化终点的白圈
                label_branch_points=F, #不显示曲线分支点的黑圈
                group_label_size = 5, # 在当前参数下无用 
                graph_label_size = 2.5,# 曲线上黑白圈节点字号大小 可能会导致重叠 
                cell_size = 1) #加大细胞点大小
#trajectory_graph_segment_size = 1.5)# 曲线粗细 
ggsave("./Trajectory_cluster_Monocle3_Plot.png",gg,width = 7,height = 7)
ggsave("./Trajectory_cluster_Monocle3_Plot.pdf",gg,width = 7,height = 7)

gg = plot_cells(cds,
                reduction_method = "UMAP",  
                color_cells_by = "cluster",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用 
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠 
                trajectory_graph_segment_size = 2,# 曲线粗细 
                cell_size = 1) #加大细胞点大小
ggsave("./Trajectory_branch_Monocle3_Plot.png",gg,width = 8,height = 7)
ggsave("./Trajectory_branch_Monocle3_Plot.pdf",gg,width = 8,height = 7)

saveRDS(cds,"./allClusterCDS.rds")
cds=readRDS('./allClusterCDS.rds')

get_earliest_principal_node <- function(cds, time_bin="2"){
  cell_ids <- which(colData(cds)[, "cell_type"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
get_earliest_principal_node(cds,time_bin="Top2a+")
cds <- order_cells(cds, 
                   root_pr_nodes="Y_28")

# allorder = function(cds,allnode = NULL,cell = NULL,outpath = NULL){
#   if(!is.null(allnode)){
#     cds <- order_cells(cds,reduction_method = "UMAP",root_pr_nodes = allnode)
#   }else if(!is.null(cell)){
#     cds <- order_cells(cds,reduction_method = "UMAP",root_cells = cell)
#   }
# }
#cds = allorder(cds)
gg = plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups=FALSE,
  label_leaves=FALSE,
  label_branch_points=FALSE,
  trajectory_graph_segment_size = 1.5,
  graph_label_size=3,
  cell_size = 1) #加大细胞点大小
ggsave("./Trajectory_pseudotime_Monocle3_Plot.png",gg,width = 8,height = 7)
ggsave("./Trajectory_pseudotime_Monocle3_Plot.pdf",gg,width = 8,height = 7)

gg = plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "cluster",
  label_cell_groups=TRUE,
	group_label_size = 7,
  label_leaves=FALSE,
  label_branch_points=FALSE,
  trajectory_graph_segment_size = 1.5,
  graph_label_size=3,
  cell_size = 1) #加大细胞点大小
ggsave("./Trajectory_pseudotime_Monocle3_Plot_cell_type.png",gg,width = 8,height = 7)
ggsave("./Trajectory_pseudotime_Monocle3_Plot_cell_type.pdf",gg,width = 8,height = 7)


# 计算拟时序差异基因？
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- as.matrix(scRNA@reductions[["umap"]]@cell.embeddings)
Track_genes <- graph_test(cds,neighbor_graph="principal_graph",cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value<1e-3)
write.csv(Track_genes,"Trajectory_genes.csv",row.names=F)

Track_genes_sig<-Track_genes %>% top_n(n=10,morans_I) %>% pull(gene_short_name)%>%as.character()
p <- plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="cell_type", min_expr=0.5, ncol=2)
ggsave("./Trajectory_Track_genes_Monocle3_Plot.png",p,width = 8,height = 8)
ggsave("./Trajectory_Track_genes_Monocle3_Plot.pdf",p,width = 8,height = 8)
p <- plot_cells(cds,genes=Track_genes_sig,show_trajectory_graph=FALSE,label_cell_groups=FALSE,label_leaves=FALSE)
p$facet$paramssncol<-5
ggsave("./Trajectory_Track_genes_Featureplot_Monocle3_Plot.png",p,width = 20,height = 15)
ggsave("./Trajectory_Track_genes_Featureplot_Monocle3_Plot.pdf",p,width = 20,height = 15)