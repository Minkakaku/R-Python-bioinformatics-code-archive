setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/10.Macro_monocle/')
library(monocle3)
rm(list = ls())
gc()

library(sceasy)
library(Seurat)
library(reticulate)
# sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/monocle3/anno.h5ad', from="anndata", to="seurat",
#                       outFile='anno.rds')


ALLobj = readRDS('/mnt/NexData_file/home/zwl/Project/SC2511076348/require7/Myeloids/annoraw.rds')
Seurat::DimPlot(ALLobj,reduction = 'umap',group.by = c('cell_type','group'))
subobj = subset(ALLobj,cell_type %in% c("Macro1","Osteoclast","Macro2"))
# ALLobj$group = paste0(ALLobj$group_plot,'_',ALLobj$group_time)
unique(ALLobj$group)
library(monocle3)
library(ggplot2)

obj = ALLobj
f_data = data.frame(gene_short_name = rownames(obj),
                    row.names = rownames(obj))
expression_matrix = as.matrix(LayerData(obj, assay="RNA", layer='counts'))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = obj@meta.data,
                         gene_metadata =  f_data)
cds <- preprocess_cds(cds, num_dim = 100)
#plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, color_cells_by="cell_type",cell_size = 1,group_label_size  = 5)
plot_cells(cds, color_cells_by="group")

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds,color_cells_by = "cell_type", 
                 show_trajectory_graph = FALSE,
                 cell_size = 1,
                 group_cells_by = 'cell_type',group_label_size  = 5 )
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,cell_size = 1,group_label_size  = 5) + 
  ggtitle("label by partitionID")
library(patchwork)
ggsave(filename = paste0('./1.monocle3_cell_type_umap.png'),wrap_plots(p1, p2),w=10,h=5)
ggsave(filename = paste0('./1.monocle3_cell_type_umap.pdf'),wrap_plots(p1, p2),w=10,h=5)


# #用注释后的cluster
# library(do)
# obj$cluster1 = 1
monocle3::partitions(cds, "UMAP") = obj$cell_type #cell_type



## 识别轨迹
cds <- learn_graph(cds,close_loop=FALSE, #分支不画圈
                   use_partition = FALSE)
# learn_graph_control=list(minimal_branch_len = 6))# 所有细胞看称整体画分支
#learn_graph_control=list(minimal_branch_len = 6))#保留一些比较短的分支
p = plot_cells(cds, color_cells_by = "group",
               label_groups_by_cluster = FALSE,
               label_leaves = FALSE,group_label_size  = 5, 
               cell_size = 1,
               label_cell_groups = F,
               label_branch_points = FALSE)
ggsave(filename = paste0('./2.monocle3_group_umap.png'),p,w=5,h=4)
ggsave(filename = paste0('./2.monocle3_group_umap.pdf'),p,w=5,h=4)



##选择起点
get_earliest_principal_node  <- function(cds, time_bin="GMP"){
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
ggsave(filename = paste0('./3.monocle3_trajectory_pseudotime_umap.png'),p2,w=5,h=4)
ggsave(filename = paste0('./3.monocle3_trajectory_pseudotime_umap.pdf'),p2,w=5,h=4)
save(file = paste0('./monocle3_order.Rdata'),cds)

Lip_genes <- c('CD36','PLIN2','LPL','PPARG','APOE','APOC1','CD9',
               'LGALS3','FABP4','FABP5','ABCA1','ABCG1')
Lip_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Lip_genes,]

plot_genes_in_pseudotime(Lip_lineage_cds,
                         color_cells_by="cell_type",
                         min_expr=0.5)

plot_genes_violin(Lip_lineage_cds, group_cells_by="cell_type", ncol=2) 
p2 = plot_cells(cds, genes=Lip_genes,
                show_trajectory_graph=T,
                label_cell_groups=FALSE,
                label_leaves=FALSE)+
  scale_color_gradientn(colors = c("#39bbec","#f38466","#b81f25"))
ggsave(filename = paste0('./4.monocle3_trajectory_gene_umap.png'),p2,w=12,h=9)
ggsave(filename = paste0('./4.monocle3_trajectory_gene_umap.pdf'),p2,w=12,h=9)
## 提取具有显著时间成分的基因
gene_fits <- fit_models(cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "pseudotime")%>% dplyr::select(-c("model","model_summary","status" ))
write.csv(file = "/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/monocle3/cor_pseudotime_gene.csv",emb_time_terms)
# pse = monocle3::pseudotime(cds)
# write.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/monocle3/pseudotime.csv',pse)

# pse = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/monocle3/pseudotime.csv',row.names = 1)
# emb_time_terms = read.csv('/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/monocle3/cor_pseudotime_gene.csv',row.names = 1)
# 
# library(dplyr)
# emb_time_terms = emb_time_terms %>% filter (q_value < 0.05,num_cells_expressed>3) 
# mat = expression_matrix[emb_time_terms$gene_short_name,rownames(pse)[order(pse$x)]]
# write.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require6/monocle3/pseudotime_count_gene_exp.csv',mat)
emb_time_terms1 = emb_time_terms %>% filter (q_value < 0.01,num_cells_expressed>50) %>% arrange(q_value) %>% head(500)
mat <- ClusterGVis::plot_pseudotime_heatmap2(cds_subset = cds[,emb_time_terms1])
# 先加载igraph，手动暴露.nei对象
library(igraph)
library(ClusterGVis)#要把所有R文件的.nei换成neighbors
# kmeans
km <- kmeans(mat, centers = 5)

df <- as.data.frame(mat)
df$cluster <- as.factor(km$cluster)

ck <- list(
  data = df,                  # ⭐这里必须带 cluster 列
  cluster = km$cluster,
  cluster.num = 5,
  sample = colnames(mat),
  gene = rownames(mat),
  type = "gene"
)

visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = FALSE)
# path <- "/mnt/NexData_file/home/zwl/data/ClusterGVis-main"


