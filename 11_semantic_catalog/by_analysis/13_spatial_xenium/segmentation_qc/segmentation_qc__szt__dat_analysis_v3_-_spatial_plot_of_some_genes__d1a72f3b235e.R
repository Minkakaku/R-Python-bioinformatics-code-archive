library(tidyverse)
library(Seurat)


all_data <- readRDS("all_data_projected.rds")
all_data <- all_data %>% JoinLayers()
DefaultAssay(all_data) <- "RNA"

sp_out_dir <- "20251026 spatial plot"
dir.create(sp_out_dir)


# 1 spatial plot of gene -----------
cells_boundary_ls <- readRDS('cells_boundary_ls.rds')

gene_list <- c(
  'Cxcl10',
  'Hspa1b',
  'Hspa1a',
  'Slfn4',
  'Dyrk3',
  'Oasl1',
  'Il1rn'
)

sp_p_gene_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
  
    all_meta <- all_data[[]] %>% rownames_to_column("cell_id") %>% dplyr::select(-x, -y)
    all_gene <- 
      all_data %>% 
      FetchData(vars = gene_list) %>%
      rownames_to_column("cell_id") 
    
    dat <- 
      cells_boundary %>%
      inner_join(all_meta, by = "cell_id") %>%
      inner_join(all_gene, by = 'cell_id')
    
    dat_sample <-
      dat %>% 
      mutate(roi = str_c(roi, group, sep = ":")) %>% 
      group_by(roi) %>% 
      summarise(x = min(x), y = min(y))
    
    gene_list %>%
      setNames(., .) %>%
      map(function(gene){
        # val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        # print(val_rng)
        # dat[[gene]] <- pmin(dat[[gene]], val_rng[2])
        # dat[[gene]] <- pmax(dat[[gene]], val_rng[1])
        
        p <- dat %>%
          ggplot(aes(x = x, y = y)) +
          geom_polygon(
            aes(fill = .data[[gene]], alpha = .data[[gene]], group = cell_id),
            show.legend = c(fill = T, alpha = F)
          )  +
          geom_text(data = dat_sample, aes(x = x, y = y, label = roi)) +
          scale_y_reverse() +
          # coord_fixed(ratio = 1/image_wh_r) +
          coord_equal() +
          scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds"), transform = "log1p") +
          theme_void() +
          ggtitle(gene)
        p
      })
  })


# export image
sp_p_gene_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })


# gene expression
# dat_gene <- 
#   all_data %>%
#   FetchData(vars = c('project', gene_list)) %>%
#   rownames_to_column("cell_id") 
# 
# 
# 
# dat_gene$Cxcl10 %>% density() %>% plot()
# 
# dat_gene 



p_vln_gene <-
  all_data %>%
  VlnPlot(
    features = gene_list,
    group.by = "project",
    # pt.size = 0
  )

pdf(
  file = file.path(sp_out_dir, "vlnplot of selected gene.pdf"),
  height = 12,
  width = 15
)
p_vln_gene %>% print
dev.off()


