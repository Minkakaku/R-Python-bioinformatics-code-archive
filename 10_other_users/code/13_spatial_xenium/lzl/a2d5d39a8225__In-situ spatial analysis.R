library(tidyverse)
library(Seurat)

# set default ---------
out_dir <- 'OUT_DIR'
dir.create(out_dir)

# data loading ----------
all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data@meta.data %>% head()



# spatial dat
# when all_data is transformed from a h5ad, spatial coord can be acquired by `Embeddings(all_data[['spatial']])`
# or, can be acquired from all_data@images... when all_data is VisiumHD
# or, can be acquired from cells.parquet from Xenium orignial directory when all_data is Xenium
some_genes <- c("Cd3d", "Col1a1", "Epcam", "Pecam1", "Cd68")
annotation_col <- 'annotation'

cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    all_data, 
    vars = c(annotation_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[annotation_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "1. spatial plot of celltype in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()


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
      file.path(out_dir, str_glue("2. sptial plot of some gene in sample {sam}.pdf") ),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })
