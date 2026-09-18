setwd('/mnt/NexData_file/home/zwl/Project/other/SC2603242447')
library(monocle3)
library(tradeSeq)
library(magrittr)
library(Seurat)
load('/mnt/NexData_file/home/zwl/Project/other/SC2603242447/monocle3.Rdata')
# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

# Get the root vertices
# It is the same node as above
root <- cds@principal_graph_aux$UMAP$root_pr_nodes

# Get the other endpoints
mst <- principal_graph(cds)$UMAP
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

# For each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  # We find the cells that map along that path
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()
rownames(cellWeights) <- colnames(cds)

pt <- pseudotime(cds)
pseudotime <- matrix(
  pt,
  nrow = length(pt),
  ncol = ncol(cellWeights)
)
rownames(pseudotime) <- colnames(cds)
colnames(pseudotime) <- colnames(cellWeights)

## 检查
identical(rownames(pseudotime), rownames(cellWeights))
identical(colnames(pseudotime), colnames(cellWeights))

counts <- SingleCellExperiment::counts(cds)

sc <- FindVariableFeatures(
  sc,
  selection.method = "vst",
  nfeatures = 2000
)
hvg <- VariableFeatures(sc)

counts_hvg <- counts[hvg, ]



#### 抽样
set.seed(123)

cluster <- sc$cell_type
cell_ids <- rownames(colData(cds))

target_max <- 2000   # 超过这个就抽样
target_min_keep <- 1 # 小群全部保留

cells.keep <- character()

for (cl in unique(cluster)) {
  
  cells <- cell_ids[cluster == cl]
  
  n <- length(cells)
  
  # 如果有 pseudotime，建议过滤 NA
  pt.sub <- pseudotime(cds)[cells]
  keep_valid <- !is.na(pt.sub)
  cells <- cells[keep_valid]
  pt.sub <- pt.sub[keep_valid]
  
  n <- length(cells)
  
  if (n == 0) next
  
  # 排序（保证轨迹连续性）
  cells <- cells[order(pt.sub)]
  
  # ===== 核心规则 =====
  if (n > target_max) {
    
    target <- target_max
    
    idx <- unique(round(seq(
      1, n,
      length.out = target
    )))
    
    cells.sel <- cells[idx]
    
  } else {
    
    # 小 cluster 全保留
    cells.sel <- cells
  }
  
  cells.keep <- c(cells.keep, cells.sel)
}

length(cells.keep)


counts_hvg_sub <- counts_hvg[, cells.keep]
pseudotime_sub <- pseudotime[cells.keep, , drop = FALSE]
cellWeights_sub <- cellWeights[cells.keep, , drop = FALSE]

pseudotime_sub <- pseudotime_sub[,c('Y_49','Y_87','Y_79','Y_112')]
cellWeights_sub <- cellWeights_sub[,c('Y_49','Y_87','Y_79','Y_112')]

keep <- rowSums(cellWeights_sub) > 0
counts_hvg_sub <- counts_hvg_sub[, keep]
pseudotime_sub <- pseudotime_sub[keep, , drop = FALSE]
cellWeights_sub <- cellWeights_sub[keep, , drop = FALSE]

save(file = '/mnt/NexData_file/home/zwl/Project/other/SC2603242447/res.RData',counts_hvg_sub,pseudotime_sub,cellWeights_sub)
