detach("package:monocle")
# detach("package:monocle3")


library(monocle3)
library(Seurat)
library(tidyverse)

# outdir_monocle3 <- "monocle3 - annov6"
# dir.create(outdir_monocle3)


all_data_full <- readRDS('all_data_full_v6.rds')

sc_epi <- all_data_full[, all_data_full$sub_celltype2 %in% c(
  'chief_cells',
  # 'endocrine_cells',
  # 'enteroendocrine_cells',
  'gmc',
  'parietal_cells',
  'pmc',
  'pmc_infla',
  'enterocytes',
  'goblet_cells',
  'tumor_cells'
)]
DefaultAssay(sc_epi) <- "Spatial.008um" 

# sc_epi$site %>% unique
# [1] distal_stomach Gastric_antrum body

# QC
sc_epi2 <- sc_epi[, 
                  sc_epi$nCount_Spatial.008um > 100 & 
                    sc_epi$site %in% c('body') & 
                    (! sc_epi$group %in% c("TLS", "other"))
]

# sc_epi2$annotation_full <- sc_epi2$sub_celltype2 %>% fct_drop()

# subset
set.seed(1234)
# meta_sc_epi2_subset <-
#   sc_epi2[[]] %>%
#   rownames_to_column("cell_id") %>%
#   # group_by( sub_celltype , project ) %>%
#   # slice_sample(n = 200)
#   slice_sample(n = 5e4)

meta_sc_epi2_subset <- sc_epi2[['sub_celltype']] %>% rownames_to_column("cell_id")
sc_epi2_subset <- sc_epi2[, meta_sc_epi2_subset$cell_id]

# factor 
sc_epi2_subset$sub_celltype2  <- sc_epi2_subset$sub_celltype2 %>% fct_drop()

# set default ------- 
# outdir_monocle3 <- "monocle3 - v2 - gastric body" 
outdir_monocle3 <- "monocle3 - v3 - gastric body - all cells" # 20260115 using all cells
dir.create(outdir_monocle3) 


sc_epi2_subset <- FindVariableFeatures(sc_epi2_subset, nfeatures = 3000)
sc_epi2_subset <- ScaleData(
  sc_epi2_subset,
  vars.to.regress = c(
    'nCount_Spatial.008um'
  )
)

sc_epi2_subset <- RunPCA(sc_epi2_subset, features = VariableFeatures(sc_epi2_subset))

# pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(sc_epi2_subset, ndims = 30)
# dev.off()

dims <- 10
sc_epi2_subset <- RunUMAP(sc_epi2_subset, dims = 1:dims, return.model = TRUE) # , assay  = 'sketch')


# pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_epi2_subset, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi2_subset, group.by = "sub_celltype2", label = TRUE, reduction = "umap") + NoLegend()
# dev.off()


sc_epi2_subset <- sc_epi2_subset %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

sc_epi2_subset <- RunUMAP(
  sc_epi2_subset,
  dims = 1:dims,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)
DimPlot(sc_epi2_subset, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi2_subset, group.by = "sub_celltype", label = TRUE, reduction = "harmonyumap") + NoLegend()


saveRDS(sc_epi2_subset, file = file.path(outdir_monocle3, "sc_epi2_subset.rds"))

# cds --------

# function `add_minicoord` is in dat_analysis_sc_v4.R

# new monocle3 object ------------------
expr_dat_monocle3 <- as(LayerData(JoinLayers(sc_epi2_subset), "counts"), "sparseMatrix")

p_data_monocle3 <- sc_epi2_subset[[]] 
f_data_monocle3 <- data.frame(
  gene_short_name = rownames(sc_epi2_subset),
  row.names = rownames(sc_epi2_subset)
)

cds3 <- new_cell_data_set(
  expr_dat_monocle3,
  cell_metadata = p_data_monocle3,
  gene_metadata = f_data_monocle3
)

cds3 <- estimate_size_factors(cds3)

# add umap
reducedDims(cds3)$UMAP <- Embeddings(sc_epi2_subset[['harmonyumap']])

# cluster_cells using louvain/leiden to compute cluster
# also, using kNN method to compute partitions
# partitions like supercluster
cds3 <- cluster_cells(
  cds3,
  k = 100,
  partition_qval = 0.9,
  reduction_method = "UMAP",
  cluster_method = "louvain"
)


plot_cells(
  cds3, 
  color_cells_by="cluster", 
  label_groups_by_cluster=F,
  group_label_size = 5,
  show_trajectory_graph  = F
)
plot_cells(
  cds3, 
  color_cells_by="partition", 
  label_groups_by_cluster=F,
  group_label_size = 5,
  show_trajectory_graph  = F
)

## manual set partition ---------
pmc_zone <- c('41', '9', '31', '43', '37', '44', '25', '47')

stopifnot(all(names(cds3@clusters$UMAP$clusters) == names(cds3@clusters$UMAP$partitions)))
levels(cds3@clusters$UMAP$partitions) <- 1:2
cds3@clusters$UMAP$partitions[  cds3@clusters$UMAP$clusters  %in% pmc_zone] <- 1
cds3@clusters$UMAP$partitions[! cds3@clusters$UMAP$clusters  %in% pmc_zone] <- 2

cds3@clusters$UMAP$partitions 

## split tumor  to OLFM4+ tumor and others tumor


# replace cluster
DimPlot(sc_epi2_subset, group.by = "RNA_snn_res.0.2", reduction = "harmonyumap")
DimPlot(sc_epi2_subset, group.by = "RNA_snn_res.0.4", reduction = "harmonyumap")
DimPlot(sc_epi2_subset, group.by = "RNA_snn_res.0.6", reduction = "harmonyumap")
DimPlot(sc_epi2_subset, group.by = "RNA_snn_res.0.8", reduction = "harmonyumap")



# partitions is supercluster, can be used to control trajectory count
# cluster control how trajectory is learnning
# cds3@clusters$UMAP$clusters = pData(cds3)$RNA_snn_res.0.8
# cds3@clusters$UMAP$partitions <- sc_epi2_subset$annotation

# learn graph ---------------

use_partition = T
cds_run_params <- list(
  minimal_branch_len = 30,
  nn.k = 100,
  euclidean_distance_ratio = 1
)

cds3 <- learn_graph(
  cds3,
  learn_graph_control = cds_run_params,
  close_loop = F,
  use_partition = use_partition,
  verbose = T
)


cds3_p1 <- plot_cells(
  cds3,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype3",
  group_cells_by      = "sub_celltype3",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

cds3_p1_2 <-
  cds3_p1 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # scale_color_manual(values = fibro_cols_sc) +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )

ggsave(
  cds3_p1_2,
  filename = file.path(outdir_monocle3, "monocle trajectory - annotation.pdf"),
  width  = 8,
  height = 6
)
# ggsave(
#   add_minicoord(cds3_p1_2),
#   filename = file.path(outdir_monocle3, "monocle trajectory - annotation - v2.pdf"),
#   width  = 8,
#   height = 6
# )

cds3_p2 <- plot_cells(
  cds3,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype",
  group_cells_by      = "sub_celltype",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

cds3_p2_2 <-
  cds3_p2 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # scale_color_manual(values = fibro_cols_sc) +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )

ggsave(
  cds3_p2_2,
  filename = file.path(outdir_monocle3, "monocle trajectory - annotation - v2.pdf"),
  width  = 8,
  height = 6
)

#get_earliest_principal_node <- function(cds3, time_bin="Pre_immature B"){
get_earliest_principal_node <- function(cds3, time_bin="0", group.by = "annotation"){
  stopifnot(group.by %in% colnames(colData(cds3)))
  stopifnot(time_bin %in% cds3[[group.by]])
  
  cell_ids <- which(colData(cds3)[, group.by] == time_bin)
  
  closest_vertex <- cds3@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds3), ])
  root_pr_nodes  <- igraph::V(principal_graph(cds3)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

cds3 <- 
  order_cells(
    cds3, 
    root_pr_nodes = get_earliest_principal_node(
      cds3,
      time_bin = "S1 fibroblasts",
      group.by = "annotation"
    )
  )

# saveRDS(cds3, file.path(outdir_monocle3, "monocle3.rds"))




# allorder = function(cds3,allnode = NULL,cell = NULL,outpath = NULL){
#   if(!is.null(allnode)){
#     cds3 <- order_cells(cds3,reduction_method = "UMAP",root_pr_nodes = allnode)
#   }else if(!is.null(cell)){
#     cds3 <- order_cells(cds3,reduction_method = "UMAP",root_cells = cell)
#   }
# }

cds3_p2 <- 
  plot_cells(
    cds3,
    reduction_method = "UMAP",
    color_cells_by = "pseudotime",
    label_cell_groups=FALSE,
    label_leaves=FALSE,
    label_branch_points=FALSE,
    trajectory_graph_segment_size = 1,
    graph_label_size=3,
    cell_size = 0.5
  ) 


cds3_p2_2 <-
  cds3_p2 +
  labs(
    title = "Monocle3-pseudotime",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # theme_classic() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )

ggsave(
  cds3_p2_2, 
  filename = file.path(outdir_monocle3, "monocle trajectory - pseudotime.pdf"),
  width = 8,
  height = 6
)

ggsave(
  add_minicoord(cds3_p2_2), 
  filename = file.path(outdir_monocle3, "monocle trajectory - pseudotime - v2.pdf"),
  width = 8,
  height = 6
)



# branch test
branch_genes <- graph_test(
  cds3,
  neighbor_graph = "principal_graph",  # 基于轨迹图
  cores = 1                            # 并行计算
)

# 筛选显著基因（按 q 值排序）
sig_genes <- 
  branch_genes %>% 
  filter(q_value < 0.05) %>% 
  arrange(desc(morans_I))

# 查看前10个基因
head(sig_genes, 10)

sig_genes



## plot PAGA
paga_pseudotime_df <- read.csv("monocle3/PAGA pseudotime.csv", row.names = 1)

sc_epi2_subset <- sc_epi2_subset %>% AddMetaData(paga_pseudotime_df)

p_PAGA <-
  FeaturePlot(
    sc_epi2_subset, 
    features = "dpt_pseudotime", 
    reduction = "harmonyumap"
  ) +
  viridis::scale_color_viridis(
    option = "C", 
    na.value = NA
  ) +
  ggtitle('PAGA-pseudotime')
p_PAGA$data <- p_PAGA$data %>% arrange(dpt_pseudotime)
p_PAGA



ggsave(
  p_PAGA, 
  filename = file.path(outdir_monocle3, "PAGA trajectory - pseudotime.pdf"),
  width = 8,
  height = 6
)

ggsave(
  add_minicoord(p_PAGA), 
  filename = file.path(outdir_monocle3, "PAGA trajectory - pseudotime - v2.pdf"),
  width = 8,
  height = 6
)



