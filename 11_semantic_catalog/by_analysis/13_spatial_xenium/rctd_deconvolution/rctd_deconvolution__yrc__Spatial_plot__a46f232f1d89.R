library(Seurat)
library(tidyverse)
library(tibble)
library(data.table)

## spatial plot -----------------------------------------------------------------
spatial_theme <- 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(color = "black", size = 1, fill = NA), #边框
        #axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.title = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.title = element_text(hjust = 0.5, size = 20))


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")

cells_boundary_ls <- readRDS('./R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')
all_data <- readRDS('./R_OUTPUT/ALL.RDS')
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = allcolour) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })

grp <- names(cells_boundary_ls)
p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("/mnt/NexData_file/home/yrc/XE2511046228/26.6.3-邻域分析/spatial plot",  
                  str_glue("spatial dim plot {grp}_cellType.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  
  
  


## highlight -----------
### Tumor ----------
CAF <- readRDS('./RCTD_CAF/all_data_sp_rctd_anno.rds')
CAF@meta.data$first_type <- as.character(CAF@meta.data$first_type) %>% replace_na("unkonwn")
all_data@meta.data$annotation <- as.character(all_data@meta.data$annotation)
all_data@meta.data[rownames(all_data@meta.data)%in%rownames(CAF@meta.data),]$annotation <- CAF@meta.data$first_type

all_data$annotation %>% table()
all_data$annotation_highlight <- fct_collapse(
  all_data$annotation, 
  'Others' = unique(all_data$annotation)[!unique(all_data$annotation) %in% c('M-apCAF','Epithelial cells',"T cells")]
)
all_data$annotation_highlight %>% table()
all_data_meta <- all_data[[]] %>% rownames_to_column('cell_id')


p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    # if (!x %in% colnames(all_data_meta)) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation_highlight, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c("lightgrey","#40E0D0","#DC143C","#32CD32")) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  })


grp <- names(cells_boundary_ls)
p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("/mnt/NexData_file/home/yrc/XE2511046228/26.6.3-邻域分析/spatial plot",  
                  str_glue("spatial dim plot {grp}_highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  
