library(monocle3)
library(Seurat)
library(tidyverse)

out_dir_260709 <- "20260709 monocle3 - annov7"
dir.create(out_dir_260709)

all_data_full <- readRDS('all_data_full_v7.rds')
meta_roi_full <- all_data_full[[]]


# 3 samples ---------------- 
cell_idx_3samples <- all_data_full$project == c('A246283', 'A255670', 'A234244')
epi_cells_in_subtype <- c(
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
)
sc_epi_full <- all_data_full[, cell_idx_3samples & (all_data_full$sub_celltype2_v7 %in% epi_cells_in_subtype )]
DefaultAssay(sc_epi_full) <- "Spatial.008um"

# QC
sc_epi_full <- sc_epi_full[,
    sc_epi_full$nCount_Spatial.008um > 100 &
    # sc_epi_full$site %in% c('body') &
    (! sc_epi_full$group %in% c("TLS", "other"))
]
sc_epi_full

# umap 
sc_epi_full <- FindVariableFeatures(sc_epi_full, nfeatures = 3000)
sc_epi_full <- ScaleData(
  sc_epi_full,
  vars.to.regress = c(
    'nCount_Spatial.008um'
  )
)

sc_epi_full <- RunPCA(sc_epi_full, features = VariableFeatures(sc_epi_full))

# pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(sc_epi_full, ndims = 30)
# dev.off()

dims <- 10

sc_epi_full <- RunUMAP(sc_epi_full, dims = 1:dims, return.model = TRUE, assay  = 'sketch')


# pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_epi_full, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi_full, group.by = "sub_celltype2", label = TRUE, reduction = "umap") + NoLegend()
# dev.off()


sc_epi_full <- sc_epi_full %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

sc_epi_full <- RunUMAP(
  sc_epi_full,
  dims = 1:dims,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(
  file = file.path(out_dir_260709, "umap plot of 3samples.pdf"),
  height = 6,
  width = 8
)
DimPlot(sc_epi_full, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi_full, group.by = "sub_celltype2", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



sc_epi_ls <- sc_epi_full %>% SplitObject('project')

# saveRDS(sc_epi_ls, file.path(out_dir_260709, 'sc_epi_ls.rds'))
# new monocle3 object ------------------
cds_ls <-
  sc_epi_ls %>%
  map(function(sc){
    expr_dat_monocle3 <- as(LayerData(JoinLayers(sc), "counts"), "sparseMatrix")
    
    p_data_monocle3 <- sc[[]] 
    f_data_monocle3 <- data.frame(
      gene_short_name = rownames(sc),
      row.names = rownames(sc)
    )
    
    cds3 <- new_cell_data_set(
      expr_dat_monocle3,
      cell_metadata = p_data_monocle3,
      gene_metadata = f_data_monocle3
    )
    
    cds3 <- estimate_size_factors(cds3)
    # add umap
    reducedDims(cds3)$UMAP <- Embeddings(sc[['harmonyumap']])
    
    
    
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
    
    
    cds3
  })




p_dr_by_monocle_ls <-
  cds_ls  %>%
  imap(function(cds, name){
    p1 <- plot_cells(
      cds,
      color_cells_by="cluster", 
      label_groups_by_cluster=F,
      group_label_size = 5,
      show_trajectory_graph  = F
    )
    p2 <- plot_cells(
      cds, 
      color_cells_by="partition", 
      label_groups_by_cluster=F,
      group_label_size = 5,
      show_trajectory_graph  = F
    )
    p1 + p2 + patchwork::plot_annotation(title = name)
  })
  
pdf(
  file = file.path(out_dir_260709, "umap plot of 3samples - monocle.pdf"),
  height = 6,
  width  = 15
)
p_dr_by_monocle_ls %>% walk(print)
dev.off()


# learn graph ---------------
cds_graph_ls <- 
  cds_ls %>%
  map(function(cds){
    use_partition = T
    cds_run_params <- list(
      minimal_branch_len = 30,
      nn.k = 100,
      euclidean_distance_ratio = 1
    )
    
    cds <- learn_graph(
      cds,
      learn_graph_control = cds_run_params,
      close_loop = F,
      use_partition = use_partition,
      verbose = T
    )
    cds
  })

# saveRDS(cds_graph_ls, file.path(out_dir_260709, 'cds_graph_ls.rds'))



p_graph_by_monocle_ls <- 
  cds_graph_ls %>%
  imap(function(cds, name){
    cds3_p1 <- plot_cells(
      cds,
      reduction_method    = "UMAP",
      color_cells_by      = "sub_celltype_v7",
      group_cells_by      = "sub_celltype_v7",
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
        title = name,
        x = "harmonyumap_1",
        y = "harmonyumap_2"
      ) +
      # scale_color_manual(values = fibro_cols_sc) +
      theme(
        plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black")
      )
  })


pdf(
  file = file.path(out_dir_260709, "monocle trajectory - annotation - sub_celltype.pdf"),
  width  = 8,
  height = 6
)
p_graph_by_monocle_ls %>% walk(print)
dev.off()


p_graph_by_monocle_ls2 <- 
  cds_graph_ls %>%
  imap(function(cds, name){
    cds3_p2 <- plot_cells(
      cds,
      reduction_method    = "UMAP",
      color_cells_by      = "sub_celltype2_v7",
      group_cells_by      = "sub_celltype2_v7",
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
        title = name,
        x = "harmonyumap_1",
        y = "harmonyumap_2"
      ) +
      # scale_color_manual(values = fibro_cols_sc) +
      theme(
        plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
        axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black")
      )
  })



pdf(
  file = file.path(out_dir_260709, "monocle trajectory - annotation - sub_celltype2.pdf"),
  width  = 8,
  height = 6
)
p_graph_by_monocle_ls2 %>% walk(print)
dev.off()



