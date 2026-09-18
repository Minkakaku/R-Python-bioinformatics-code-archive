library(Seurat)
library(tidyverse)


out_dir2 <- "density"
dir.create(out_dir2)

# cells_position_ls
# cell annotation
all_data <- readRDS("adata_filter.rds")
all_data[['RNA']] <- as(all_data[['RNA']], 'Assay5')

meta_data_anno <- read.csv('adata_celltype_obs.csv', row.names = 1)


all_data <- 
  all_data %>%
  AddMetaData(meta_data_anno['celltype'])

cells_position_ls <- 
  cbind(
    Embeddings(all_data[['spatial']]),
    all_data[['sample']]
  ) %>%
  rownames_to_column('cell_id') %>%
  dplyr::rename(
    x_centroid = SPATIAL_1,
    y_centroid = SPATIAL_2
  ) %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

annotation_col <- "celltype"
density_cell_to_plot <- 
  c(
    "Fibroblasts", 
    "Tumor",
    "T cell" 
  )

kde_dat_tidy_all <-
  cells_position_ls %>%
  map(function(x){
    sp_dat <- 
      x[c('cell_id', 'x_centroid', 'y_centroid')] %>% 
      inner_join(
        rownames_to_column(all_data[[annotation_col]], 'cell_id' ),
        by = "cell_id"
      ) %>%
      dplyr::filter(!is.na(.data[[annotation_col]]))
    
    # kde
    kde_dat_tidy_3cells <-
      density_cell_to_plot%>%
      map(function(x){
        sp_dat_sub <- sp_dat[sp_dat[[annotation_col]] == x, ]
        
        kde_dat <- MASS::kde2d(sp_dat_sub$x_centroid, sp_dat_sub$y_centroid)
        # contour(kde_dat)
        
        kde_dat_tidy <- kde_dat$z %>% reshape2::melt()
        kde_dat_tidy$x <- kde_dat$x[kde_dat_tidy$Var1]
        kde_dat_tidy$y <- kde_dat$y[kde_dat_tidy$Var2]
        contour(kde_dat)
        kde_dat_tidy$cell <- x
        kde_dat_tidy
      }) %>%
      bind_rows()
    kde_dat_tidy_3cells
  }) %>%
  enframe("project") %>%
  unnest(value)

kde_dat_tidy_all <-
  kde_dat_tidy_all %>%
  group_by(cell) %>%
  mutate(value_scaled = scales::rescale(value, to = c(0, 1)))


p_density_proj_ls <-
  cells_position_ls %>%
  imap(function(coords, name){
    coords <- coords %>% mutate(x = x_centroid, y = y_centroid)
    
    val_rng <-
      kde_dat_tidy_all %>% 
      group_by(cell) %>%
      nest() %>%
      deframe() %>%
      map(~ range(.x$value_scaled))
      
    
    plot_dat_celltype <-
      all_data[[annotation_col]] %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(.data[[annotation_col]] %in% density_cell_to_plot)
    plot_dat_celltype <- coords %>% inner_join(plot_dat_celltype, by = "cell_id")
    
    
    # bg
    p_bg <- 
      coords %>%
      ggplot() +
      geom_tile(aes(x = x, y = y), height = 10, width = 10, fill = "gray90") +
      theme_void() +
      coord_equal()
    
    # add all celltype
    p_add_celltype <- 
      p_bg +
      geom_tile(
        data = plot_dat_celltype,
        aes(x = x, y = y, fill = .data[[annotation_col]]),
        height = 15, 
        width  = 15,
        alpha = 0.35
      ) +
      ggsci::scale_fill_igv()
    
    # # add density
    # line_type <- c("dotted", "dashed", "twodash")
    # names(line_type) <- density_cell_to_plot
    
    kde_dat_used <- kde_dat_tidy_all[kde_dat_tidy_all$project %in% name, ]
    p_ls <- 
      density_cell_to_plot %>%
      setNames(., .) %>%
      map(function(cell){
        p_add_celltype +
          geom_contour(
            data = kde_dat_used %>% dplyr::filter(cell == .env[['cell']]),
            aes(x = x, y = y, z = value_scaled, colour = after_stat(level))
          ) + 
          scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu")), limits = val_rng[[cell]]) +
          ggtitle(cell)
      })

    # patchwork::wrap_plots(p_ls)
    p_ls
  })


# export

p_density_proj_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file.path(out_dir2, str_glue("density plot of 3 cells - {name}.pdf") ),
      height = 11,
      width  = 11
    )
    p_ls %>% walk(print)
    dev.off()
  })


