library(Seurat)
library(tidyverse)
library(tibble)
library(data.table)

out_dir <- 'DSS_spatial_plot'
dir.create(out_dir)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
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

seurat_ls <- readRDS('R_OUTPUT/seurat_DSS_HSgroup_ls.rds')
all_data <- merge(seurat_ls[[1]], seurat_ls[-1])
all_data$DSS %>% table()
all_data$annotation %>% table()
all_data$subtype %>% table()
all_data$subtype_re %>% table()
all_data$subtype_re2 %>% table()

all_data$annotation_full <- ifelse(all_data$group == 'HS_liver', as.character(all_data$subtype_re), as.character(all_data$subtype))
all_data$annotation_full %>% table()
saveRDS(all_data, 'R_OUTPUT/all_data_DSS.rds')

cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')





# target gene + celltype highlight ------------
target_gene_ls <- c(
  'Ripk3'
)
names(target_gene_ls) <- target_gene_ls
target_gene_ls %in% rownames(all_data)

unique(all_data$annotation_full[all_data$DSS == 'TRUE'])
unique(all_data$annotation_full)
# all_data$annotation7 <- fct_collapse(
#   all_data$annotation_full, 
#   'Others' = unique(all_data$annotation_full[all_data$DSS == 'FALSE'])
# )
# all_data$annotation7 %>% table()


all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
data_expr <- 
  all_data %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_meta, data_expr, by = 'cell_id')


p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) {
        if (!sample %in% data_expr_full$sample) return(NULL)
        
        # browser()
        dat <-
          cells_boundary %>%
          left_join(data_expr_full, by = 'cell_id') %>% 
          mutate(
            Annotation = 
              ifelse(DSS == 'TRUE', as.character(annotation_full), 'non-DSS') %>% 
              str_replace_all(c('Netrophil' = 'Neutrophil')), 
            Annotation = factor(Annotation, levels = c('non-DSS', unique(as.character(annotation_full))))
          )
       
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(fill = .data[[gen]], color = Annotation, group = cell_id), 
            size = 0.1,                # 边界粗细
            alpha = 0.8                # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 0.05,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = c("lightgrey", allcolour)) +
          guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(gen, " Expression"),
              barwidth = 0.5,
              barheight = 10
            ),
            color = guide_legend(
              title = "Cell Type",
              override.aes = list(
                size = 1.5,            # 图例中边界粗细
                fill = NA              # 图例中不显示填充
              )
            )
          ) + 
          coord_equal() +
          scale_y_reverse() + 
          spatial_theme
        
        p
      }) %>% 
      purrr::discard(function(x) length(x) < 1)
  })


p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(out_dir,  
                  str_glue("spatial dim plot {grp}_target gene - DSS_highlight.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 

