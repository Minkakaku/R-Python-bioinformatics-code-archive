.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle3/lib/R/library/")
library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)
library(Seurat)
scRNA = readRDS('/mnt/NexData_file/home/wmx/Project/SC2601080130/shouhou6/anno/adata_anno_count.rds')
scRNA@meta.data <- scRNA@meta.data[, !grepl("cluster", colnames(scRNA@meta.data), ignore.case = TRUE)]
Idents(scRNA) <- 'cell_type'
# 构建cds矩阵
setwd('/mnt/NexData_file/home/wmx/Project/SC2601080130/shouhou6')
#meta=read.csv('G:\\FW2024-348\\FW2024-348_10\\售后4\\monocle3\\all_features_BP.xls',sep = '\t',header = 1,row.names = 1)
#scRNA@meta.data=meta
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
minimal_branch_len1 = 25
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
                color_cells_by = "Clusters",
                show_trajectory_graph = F,
                label_cell_groups=T,#在图上表示cluster等标签
                label_leaves=F,#不显示曲线演化终点的白圈
                label_branch_points=F, #不显示曲线分支点的黑圈
                group_label_size = 5, # 在当前参数下无用 
                graph_label_size = 2.5,# 曲线上黑白圈节点字号大小 可能会导致重叠 
                cell_size = 0.5) #加大细胞点大小
#trajectory_graph_segment_size = 1.5)# 曲线粗细 
ggsave("./TrajectoryResult/Trajectory_cluster_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_cluster_Monocle3_Plot.pdf",gg,width = 7,height = 5)



gg = plot_cells(cds,
                reduction_method = "UMAP",  
                color_cells_by = "leiden",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用 
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠 
                trajectory_graph_segment_size = 2,# 曲线粗细 
                cell_size = 0.5) #加大细胞点大小

ggsave("./TrajectoryResult/Trajectory_branch_Monocle3_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_branch_Monocle3_Plot.pdf",gg,width = 7,height = 5)
saveRDS(cds,"./allClusterCDS.rds")
# cds=readRDS('./allClusterCDS.rds')
#get_earliest_principal_node <- function(cds, time_bin="Plasma"){
get_earliest_principal_node <- function(cds, time_bin="Pre-GC"){
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
ggsave("./TrajectoryResult/Trajectory_pseudotime_Monocle3_GC_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_pseudotime_Monocle3_GC_Plot.pdf",gg,width = 7,height = 5)

gg = plot_cells(cds,
                reduction_method = "UMAP",
                color_cells_by = "cell_type",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠
                trajectory_graph_segment_size = 2,# 曲线粗细
                cell_size = 0.5)

ggsave("./TrajectoryResult/Trajectory_celltype_Monocle3_GC_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_celltype_Monocle3_GC_Plot.pdf",gg,width = 7,height = 5)


get_earliest_principal_node <- function(cds, time_bin="Pre-MBCs"){
  cell_ids <- which(colData(cds)[, "cell_type"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds,
                   root_pr_nodes=get_earliest_principal_node(cds))

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
ggsave("./TrajectoryResult/Trajectory_pseudotime_Monocle3_B_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_pseudotime_Monocle3_B_Plot.pdf",gg,width = 7,height = 5)
# 
gg = plot_cells(cds,
                reduction_method = "UMAP",
                color_cells_by = "cell_type",
                label_cell_groups=FALSE,#不在图上表示cluster等标签
                label_leaves=TRUE,#显示曲线演化终点的白圈
                label_branch_points=TRUE, #显示曲线分支点的黑圈
                # group_label_size = 5, # 在当前参数下无用
                graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠
                trajectory_graph_segment_size = 2,# 曲线粗细
                cell_size = 0.5)

ggsave("./TrajectoryResult/Trajectory_celltype_Monocle3_B_Plot.png",gg,width = 7,height = 5)
ggsave("./TrajectoryResult/Trajectory_celltype_Monocle3_B_Plot.pdf",gg,width = 7,height = 5)


# plot_trajectory <- function(cds, column_name) {
#   gg <- plot_cells(cds,
#                    reduction_method = "UMAP",
#                    color_cells_by = column_name,
#                    label_cell_groups = FALSE,
#                    label_leaves = TRUE,
#                    label_branch_points = TRUE,
#                    graph_label_size = 2.5,
#                    trajectory_graph_segment_size = 2,
#                    cell_size = 0.5)
#   
#   return(gg)
# }
# output_dir='G:\\FW2024-348\\FW2024-348_10\\售后4\\monocle3\\result1\\TrajectoryResult\\TotalSeqC\\'
# 循环绘制每个列名
# for (col in clo[17:27]) {
#   tryCatch({
#     # 生成图形
#     gg <- plot_trajectory(cds, col)
#     
#     # 保存图形
#     png_file <- paste0(output_dir, col ,"_Monocle3_Plot.png")
#     pdf_file <- paste0(output_dir, col, "_Monocle3_Plot.pdf")
#     
#     ggsave(png_file, gg, width = 7, height = 5)
#     ggsave(pdf_file, gg, width = 7, height = 5)
#     
#     message("成功绘制并保存: ", col)
#   }, error = function(e) {
#     message("绘制 ", col, " 时出错: ", e$message)
#   })
# }

# gg = plot_cells(cds,
#                 reduction_method = "UMAP",  
#                 color_cells_by = "CytoTRACE_value",
#                 label_cell_groups=FALSE,#不在图上表示cluster等标签
#                 label_leaves=TRUE,#显示曲线演化终点的白圈
#                 label_branch_points=TRUE, #显示曲线分支点的黑圈
#                 # group_label_size = 5, # 在当前参数下无用 
#                 graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠 
#                 trajectory_graph_segment_size = 2,# 曲线粗细 
#                 cell_size = 1) #加大细胞点大小
# 
# ggsave("./TrajectoryResult/Trajectory_CytoTRACE_value_Monocle3_Plot.png",gg,width = 7,height = 5)
# ggsave("./TrajectoryResult/Trajectory_CytoTRACE_value_Monocle3_Plot.pdf",gg,width = 7,height = 5)
# 
# 
# 
# # 创建绘图
# p <- plot_cells(cds,
#                 reduction_method = "UMAP",  
#                 color_cells_by = "cnv_score",
#                 label_cell_groups = FALSE,
#                 label_leaves = TRUE,
#                 label_branch_points = TRUE,
#                 graph_label_size = 2.5,
#                 trajectory_graph_segment_size = 2,
#                 cell_size = 1) +
#   scale_color_gradientn(
#     colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),  # 蓝-白红11色阶
#     values = scales::rescale(c(0, 0.004, 0.012)), 
#     limits = c(0, 0.012),
#     name = "CNV Score")
# 
# ggsave("./TrajectoryResult/Trajectory_cnvscore_rb_Monocle3_Plot.png",p,width = 7,height = 5)
# ggsave("./TrajectoryResult/Trajectory_cnvscore_rb_Monocle3_Plot.pdf",p,width = 7,height = 5)
# 
# # +
# #ggtitle("CNV Score along Trajectory") +
# #theme(plot.title = element_text(hjust = 0.5))
# 
# gg = plot_cells(cds,
#                 reduction_method = "UMAP",  
#                 color_cells_by = "CytoTRACE_value",
#                 label_cell_groups=FALSE,#不在图上表示cluster等标签
#                 label_leaves=TRUE,#显示曲线演化终点的白圈
#                 label_branch_points=TRUE, #显示曲线分支点的黑圈
#                 # group_label_size = 5, # 在当前参数下无用 
#                 graph_label_size = 2.5, # 曲线上黑白圈节点字号大小 可能会导致重叠 
#                 trajectory_graph_segment_size = 2,# 曲线粗细 
#                 cell_size = 1)+
#   scale_color_gradientn(colors = viridisLite::turbo(256),
#                         #limits = c(0, 1),
#                         name = "CytoTRACE_value")
# 
# ggsave("./TrajectoryResult/Trajectory_CytoTRACE_value_turbo_Monocle3_Plot.png",gg,width = 7,height = 5)
# ggsave("./TrajectoryResult/Trajectory_CytoTRACE_value_turbo_Monocle3_Plot.pdf",gg,width = 7,height = 5)






