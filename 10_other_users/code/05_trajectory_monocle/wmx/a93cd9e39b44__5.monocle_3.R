.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle3/lib/R/library/")
library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)
library(Seurat)
scRNA = readRDS("/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/data/Neu_fixed.rds")
scRNA@meta.data <- scRNA@meta.data[, !grepl("cluster", colnames(scRNA@meta.data), ignore.case = TRUE)]
Idents(scRNA) <- 'cell_type'
# 构建cds矩阵
od = "/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/Neu_Pseudotime/Monocle3"
dir.create(od)
setwd(od)

# 构建cds矩阵
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
# 最小分枝长度
minimal_branch_len1 = 20
# 生成轨迹最近邻K
nn_k1 = 30
param = list()
param[["minimal_branch_len"]] = minimal_branch_len1
param[["nn.k"]] = nn_k1
use_partition1 = T 
cds <- learn_graph(cds,use_partition = use_partition1,close_loop = F,learn_graph_control = param,verbose = T)

library(ggplot2)
dir.create('TrajectoryResult')
gg = plot_cells(cds,
                reduction_method = "UMAP",  
                color_cells_by = "cell_type",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用 
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠 
                trajectory_graph_segment_size = 2,# 曲线粗细 
                cell_size = 0.5) #加大细胞点大小
#trajectory_graph_segment_size = 1.5)# 曲线粗细 
ggsave("./TrajectoryResult/Trajectory_celltype_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_celltype_Monocle3_Plot.pdf",gg,width = 7,height = 5)


saveRDS(cds,"./TrajectoryResult/allClusterCDS.rds")
# cds=readRDS('./allClusterCDS.rds')
#get_earliest_principal_node <- function(cds, time_bin="Plasma"){

# 根据cytotrace的结果 
get_earliest_principal_node <- function(cds, time_bin="immature_N2-like_TANs"){
  cell_ids <- which(colData(cds)[, "cell_type"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds,
                   root_pr_nodes=get_earliest_principal_node(cds))
# allorder = function(cds,allnode = NULL,cell = NULL,outpath = NULL){
#   if(!is.null(allnode)){
#     cds <- order_cells(cds,reduction_method = "UMAP",root_pr_nodes = allnode)
#   }else if(!is.null(cell)){
#     cds <- order_cells(cds,reduction_method = "UMAP",root_cells = cell)
#   }
# }
# cds = allorder(cds)
gg = plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups=FALSE,
  label_leaves=FALSE,
  label_branch_points=FALSE,
  trajectory_graph_segment_size = 1.5,
  graph_label_size=3,
  cell_size = 0.5) #加大细胞点大小
ggsave("./TrajectoryResult/Trajectory_pseudotime_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_pseudotime_Monocle3_Plot.pdf",gg,width = 7,height = 5)



gg = plot_cells(cds,
                reduction_method = "UMAP",
                color_cells_by = "group",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠
                trajectory_graph_segment_size = 2,# 曲线粗细
                cell_size = 1) #加大细胞点大小

ggsave("./TrajectoryResult/Trajectory_group_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_group_Monocle3_Plot.pdf",gg,width = 7,height = 5)


gg = plot_cells(cds,
                reduction_method = "UMAP",
                color_cells_by = "sample",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠
                trajectory_graph_segment_size = 2,# 曲线粗细
                cell_size = 1) #加大细胞点大小

ggsave("./TrajectoryResult/Trajectory_sample_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_sample_Monocle3_Plot.pdf",gg,width = 7,height = 5)

# cytotrace 
cyto = read.delim2("/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/Neu_Pseudotime/Cytotrace/CytoTRACE_value.txt")
# > identical(cyto$Cell,rownames(cds@colData))
# [1] TRUE
cds@colData$cytotrace_value = as.numeric(cyto$CytoTRACE_value)

gg = plot_cells(cds,
                reduction_method = "UMAP",
                color_cells_by = "cytotrace_value",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠
                trajectory_graph_segment_size = 2,# 曲线粗细
                cell_size = 1) #加大细胞点大小

ggsave("./TrajectoryResult/Trajectory_cytotrace_value_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_cytotrace_value_Monocle3_Plot.pdf",gg,width = 7,height = 5)



# 增加对monocle2的结果绘制一个cytotrace的结果 -------
mono2 <- readRDS("/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/Neu_Pseudotime/Monocle2/Pseudotime/Pseudotime.monocle.rds")
# > identical(cyto$Cell,rownames(mono2@phenoData))
# [1] TRUE

.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle2/lib/R/library/")
library(monocle)
mono2$cytotrace_value = as.numeric(cyto$CytoTRACE_value)
plot_theme <- function() {
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )
}
gg <- plot_cell_trajectory(mono2, color_by = 'cytotrace_value') + 
  plot_theme()
save_plot_pdf_png <- function(gg, prefix, width = 8, height = 6) {
  ggsave(paste0(prefix, ".png"), gg, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(prefix, ".pdf"), gg, width = width, height = height, limitsize = FALSE)
}

save_plot_pdf_png(
  gg = gg,
  prefix = paste0("Pseudotime.", 'cytotrace_value'),
  width = 8,
  height = 6
)