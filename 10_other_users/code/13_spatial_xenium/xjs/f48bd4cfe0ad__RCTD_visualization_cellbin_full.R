library(tidyverse)
library(Seurat)

# set default ---------
rctd_out_dir <- 'RCTD_cellbin_full//'

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")




# data loading ----------
rctd_result_ls <- readRDS(str_glue("{rctd_out_dir}/rctd_result_ls.rds"))


sc_data <- readRDS('RCTD/SC2510275964.rds')
names(allcolour) <- sc_data$cell_type %>% unique()

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data@meta.data %>% head()

meta_rctd <- 
  rctd_result_ls %>% 
  map(function(rctd_result){
    rctd_result@results$weights
  }) %>% 
  map(as.matrix) %>%
  map(as.data.frame) %>%
  bind_rows() 

# add cell weights, and replace NA to 0
all_cells <- colnames(meta_rctd)
all_data <- all_data %>% AddMetaData(meta_rctd)
rctd_na_idx <- is.na(all_data@meta.data[all_cells])
all_data@meta.data[all_cells][rctd_na_idx] <- 0


# spatial 
cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

some_genes <- all_cells # c("Cd3d", "Col1a1", "Epcam", "Pecam1", "Cd68")
annotation_col <- NULL # 'first_type_annotation'
plot_dat <-
  FetchData(
    all_data, 
    vars = c(annotation_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')





# spatial plot of some gene ---------
p_sp_gene_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    some_genes %>%
      setNames(., .) %>% # .[1] %>%
      map(function(x){
        # sub
        dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
        
        # add scale to quantile
        q_cutoff <- 0.9
        n <- 0
        while( near(quantile(dat2[[x]], q_cutoff), 0) ){
          q_cutoff <- q_cutoff + 9/10^(n + 2)
          print(q_cutoff)
          n <- n + 1
          if(n > 3) break
        }
        
        if(n <= 3){
          dat2[[x]] <- pmin(dat2[[x]], quantile(dat2[[x]], q_cutoff))
        }
        
        
        # gene value range
        val_rng <- dat2[[x]] %>% range()
        height <- 5
        width <- 5
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(fill = "gray60", height = height, width = width) + 
          geom_tile(
            aes(fill = .data[[x]], alpha = .data[[x]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          scale_fill_gradient(low = "white", high = "red",limits = val_rng)+
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(str_glue("{sam} - {x}"))
        p
      })
  })


# export 
p_sp_gene_ls %>%
  iwalk(function(p_ls, sam){
    pdf(
      file.path(rctd_out_dir, str_glue("2. sptial plot of some gene in sample {sam}.pdf") ),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })



# test by dr and cluster ------------------
sc_rctd <- 
  CreateSeuratObject(
    CreateAssay5Object(
      data = t(all_data[[all_cells]])
    ),
    meta.data = all_data[[]]
  )

sc_rctd <- 
  sc_rctd %>%
  ScaleData() %>%
  RunPCA(
    npcs = length(all_cells) - 1, 
    features = all_cells
  )

ElbowPlot(sc_rctd)

sc_rctd <- sc_rctd %>% RunUMAP(dims = 1:5)
DimPlot(sc_rctd, group.by = "sample")

sc_rctd <-
  sc_rctd %>%
  FindNeighbors(dims = 1:5) %>%
  FindClusters(resolution = 0.3)
DimPlot(sc_rctd, group.by = "RNA_snn_res.0.3")

sc_rctd2 <- sc_rctd[, sc_rctd$RNA_snn_res.0.3 != "1"]
sc_rctd2 <-
  sc_rctd2  %>%
  FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8))



DimPlot(sc_rctd2, group.by = "RNA_snn_res.1")

# RNA_snn_res.0.2 RNA_snn_res.0.4 RNA_snn_res.0.6 RNA_snn_res.0.8
pdf(
  file.path(rctd_out_dir, "0. umap of recluster.pdf"),
  height = 6,
  width = 8
)
DimPlot(sc_rctd2, group.by = "RNA_snn_res.0.2")
DimPlot(sc_rctd2, group.by = "RNA_snn_res.0.4")
DimPlot(sc_rctd2, group.by = "RNA_snn_res.0.6")
DimPlot(sc_rctd2, group.by = "RNA_snn_res.0.8")
dev.off()


p_dot <- 
  DotPlot(
    sc_rctd2,
    group.by = "RNA_snn_res.0.2",
    features = all_cells,
  ) +
  scale_color_gradientn(colours = c("white", "gray99","lightblue", "blue")) +
  RotatedAxis()

ggsave(
  p_dot,
  filename = file.path(rctd_out_dir, "0. dotplot of recluster - res0.2.pdf"),
  height = 6,
  width = 8
)

# sc_rctd2$RNA_snn_res.0.2 %>% table

cell_match_ls <- c(
  '0' = "NK cells",
  '1' = "B cells",
  '2' = "AT2_uncertain",
  '3' = "Endothelial cells",
  '4' = "T cells",
  '5' = "Ciliated cells",
  '6' = "AT1",
  '7' = "Myeloid",
  '8' = "Fibroblast"
)

# Idents(sc_rctd2) <- "RNA_snn_res.0.2"
# sc_rctd2 <- RenameIdents(sc_rctd2, cell_match_ls)
# sc_rctd2$annotation <- Idents(sc_rctd2)


sc_rctd2$annotation <- 
  sc_rctd2$RNA_snn_res.0.2 %>%
  fct_recode(!!!structure(
    names = cell_match_ls,
    names(cell_match_ls)
  ))

DimPlot(
  sc_rctd2, 
  group.by = "annotation"
)


saveRDS(sc_rctd2, file = file.path(rctd_out_dir, "sc_rctd2.rds"))


# cell props 
p_col_props <-
  sc_rctd2[[]] %>%
  group_by(sample, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation)) +
  geom_col() +
  ggsci::scale_fill_igv()

ggsave(
  p_col_props,
  filename = file.path(rctd_out_dir, "0. cell props of recluster - res0.2.pdf"),
  height = 6,
  width = 5
)


# spatial plot - celltype ------------------
cells_position_ls2 <- 
  sc_rctd2[[c('x', 'y', 'sample', 'annotation')]] %>%
  dplyr::rename(SPATIAL_1 = x, SPATIAL_2 = y) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

annotation_col2 <- 'annotation'

p_sp_ct_ls2 <- 
  cells_position_ls2 %>% # .[1] %>%
  imap(function(coords, sam){
    # dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    dat2 <- coords
    height <- 5
    width <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[annotation_col2]]),
        height = height,
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p

  })


pdf(
  file.path(rctd_out_dir, "1. spatial plot of celltype in each sample - using re-cluster.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls2 %>% walk(print)
dev.off()








