library(Seurat)
library(tidyverse)

MultiSE_metadata_final <- read.table('/mnt/NexData_file/home/xjs/xutong_spatialEcoTyper/integrated_SE/MultiSE_metadata_final.tsv', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
MultiSE_metadata_final %>% head()
MultiSE_metadata_final %>% dim()
MultiSE_metadata_final$SE %>% is.na() %>% table()
MultiSE_metadata_final$SE %>% table()

all_data <- readRDS("~/R_OUTPUT/all_data_4slides_sample_filtered_reUMAP.rds")

all_data2 <- all_data[, all_data$annotation != 'Skeletal muscle cell']
all_data2$annotation %>% table()
all_data2@meta.data %>% head()
all_data2_meta <- 
  all_data2[[]] %>% 
  rownames_to_column('CID') %>% 
  left_join(MultiSE_metadata_final, by = 'CID') %>% 
  column_to_rownames('CID')
all_data2_meta %>% head()
all_data2 <- AddMetaData(all_data2, all_data2_meta['SE'])
all_data2@meta.data %>% head()
all_data2$SE %>% is.na() %>% table()




# spatial plot ---------
## all SE --------
all_data_meta <- 
  all_data2[[]] %>% 
  rownames_to_column('cell_id')
all_data_meta$SE <- factor(all_data_meta$SE, levels = c(str_c('SE', seq(1:10))))
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')


out_dir <- 'spatialEcoTyper/'
ct_sp_out_dir <- file.path(out_dir, 'SE spatial plot')
dir.create(ct_sp_out_dir)

ct_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = SE, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })

ct_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir, str_glue("SE spatial plot - {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


## split SE -------------
cn_sp_out_dir2 <- file.path(out_dir, 'SE spatial plot - split') 
dir.create(cn_sp_out_dir2)

cn_sp_p_ls2 <- 
  cells_boundary_ls %>% # .[1] %>%
  imap(function(cells_boundary, sample){ # browser()
    # cn_color <- structure(
    #   names = paste0("CN", 1:10),
    #   ggsci::pal_igv()(10)
    # ) %>%
    #   c(
    #     Others = "gray90"
    #   )
    # image_wh_r <- image_wh_r_ls[[sample]]
    
    anno_data   <- all_data_meta[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    cn_color <- structure(
      names = dat$SE %>% levels(),
      ggsci::pal_igv()(10)
    ) %>%
      c(
        Others = "gray90"
      )
    
    p_ls <-
      dat$SE %>% levels() %>% # .[1] %>%
      set_names(., .) %>%
      map(function(SE){
        dat2 <- dat
        
        dat2$SE <- as.character(dat2$SE)
        dat2$SE[dat2$SE != SE] <- "Others"
        
        p <- dat2 %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = SE, group = cell_id))  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cn_color) +
          theme_void() +
          labs(fill = NULL, title = paste(sample, group_info, SE, sep = "_") )
        p
      })
    
    p_ls
  })


# export image
cn_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(cn_sp_out_dir2, str_glue("SE spatial plot {name} - split by SE.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })




# correlation ------------
all_data_meta_cn <- readRDS("CN_result_10CN_4slides/all_data_meta_cn.rds")
all_data_meta_cn %>% dim()
all_data_meta %>% dim()
all_data_meta2 <- 
  all_data_meta[c('cell_id', 'SE')] %>% 
  left_join(all_data_meta_cn[c('cell_id', 'CN')] , by = 'cell_id')
all_data_meta2 %>% head

tab <- table(all_data_meta2$SE, all_data_meta2$CN)
print(tab)

p1 <- pheatmap::pheatmap(
  tab, 
  display_numbers = TRUE,          # 显示频数
  number_format = "%.0f",          # 数字格式
  main = "SE vs CN frequency heatmap",
  # scale = 'row', 
  color = colorRampPalette(c("white", "steelblue"))(100),
  fontsize_number = 12
)
ggsave(file.path(out_dir, 'SE vs CN frequency heatmap - raw.pdf'), p1, width = 8, height = 6)

p2 <- pheatmap::pheatmap(
  tab, 
  display_numbers = TRUE,          # 显示频数
  number_format = "%.0f",          # 数字格式
  main = "SE vs CN frequency heatmap",
  scale = 'row', 
  color = colorRampPalette(c("white", "steelblue"))(100),
  fontsize_number = 12
)
ggsave(file.path(out_dir, 'SE vs CN frequency heatmap - scale by SE.pdf'), p2, width = 8, height = 6)

p3 <- pheatmap::pheatmap(
  tab, 
  display_numbers = TRUE,          # 显示频数
  number_format = "%.0f",          # 数字格式
  main = "SE vs CN frequency heatmap",
  scale = 'column', 
  color = colorRampPalette(c("white", "steelblue"))(100),
  fontsize_number = 12
)
ggsave(file.path(out_dir, 'SE vs CN frequency heatmap - scale by CN.pdf'), p3, width = 8, height = 6)
