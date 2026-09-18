# spatial plot of some gene ------
sp_outdir <- str_glue("spatial plot - {Date}")
dir.create(sp_outdir)

sp_genes <- c("Piezo1", "Mertk", "Hif1a", "Zbtb16", "Ece1", "Usp32", "Ube2r2", "Foxo1", "Atp2a2", "Cd38", "Fkbp5")
# sp_genes <- c("Piezo1")
sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    dat <-
      cells_boundary %>%
      right_join(dat_expr, by = "cell_id")
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = .data[[gene]], group = cell_id))  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })

# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_outdir, str_glue("spatial plot {grp}.pdf")),
      height = 20,
      width  = 20
    )
    print(p)
    dev.off()
  })










