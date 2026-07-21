library(tidyverse)
library(Seurat)

out_dir <- 'target_gene_result'
dir.create(out_dir)

# MET expr ------------
all_data <- readRDS('all_data_annotation_filtered.rds')
all_data$annotation %>% table()
# all_data_pure <- subset(JoinLayers(all_data), annotation != 'low-nCounts')
# all_data_pure$annotation %>% table()
# saveRDS(all_data_pure, 'all_data_annotation_filtered.rds')


## boxplot ---------
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


target_gene_ls <- 'MET'
all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
data_expr <-  
  all_data %>% 
  FetchData(vars = unlist(target_gene_ls), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_meta, data_expr, by = 'cell_id')

boxplot <- 
  data_expr_full %>% 
  select(annotation, all_of(target_gene_ls)) %>% 
  pivot_longer(all_of(target_gene_ls), names_to = 'gene', values_to = 'expr') %>% 
  # filter(expr > 0) %>% 
  ggplot(aes(x = annotation, y = expr, fill = annotation)) + 
  geom_boxplot(outliers = F) + 
  geom_jitter(
    width = 0.2, 
    height = 0, 
    alpha = 0.5, 
    size = 1, 
    color = "steelblue"
  ) + 
  ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
  facet_wrap(~ gene, scale = "free_y", ncol = 5) + 
  scale_fill_manual(values = allcolour) + 
  labs(y = "Expression Level") +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12), 
    axis.text.x = element_text(angle = 90)
  )

ggsave(file.path(out_dir, 'Boxplot_MET_by_annotation.pdf'), width = 8, height = 8)


boxplot <- 
  data_expr_full %>% 
  select(annotation, all_of(target_gene_ls)) %>% 
  pivot_longer(all_of(target_gene_ls), names_to = 'gene', values_to = 'expr') %>% 
  filter(expr > 0) %>% 
  ggplot(aes(x = annotation, y = expr, fill = annotation)) + 
  geom_boxplot(outliers = F) + 
  geom_jitter(
    width = 0.2, 
    height = 0, 
    alpha = 0.5, 
    size = 1, 
    color = "steelblue"
  ) + 
  ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
  facet_wrap(~ gene, scale = "free_y", ncol = 5) + 
  scale_fill_manual(values = allcolour) + 
  labs(y = "Expression Level") +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12), 
    axis.text.x = element_text(angle = 90)
  )

ggsave(file.path(out_dir, 'Boxplot_MET_by_annotation_move0.pdf'), width = 8, height = 8)


## spatial plot -----------
cells_boundary_ls <- cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
names(target_gene_ls) <- target_gene_ls
target_gene_ls %in% rownames(all_data)

p_dimPlot_ls <- 
  target_gene_ls %>% 
  map(function(gen){
    cells_boundary_ls %>% 
      imap(function(cells_boundary, sample) { # browser()
        # if (!x %in% colnames(all_data_meta)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) + 
          geom_polygon(
            aes(
              fill = .data[[gen]],
              color = annotation,
              group = cell_id
              ),
            size = 0.1,                # 边界粗细
            alpha = 0.2               # 填充透明度 
          ) + 
          ggtitle(str_glue(sample, ' - ', gen)) + 
          # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
          # 基因表达颜色梯度
          scale_fill_gradient2(
            low = "blue",              # 低表达
            mid = "white",             # 中等表达
            high = "red",              # 高表达
            midpoint = 0.1,  # 中间值
            na.value = "lightgrey"     # NA值的颜色
          ) +
          scale_color_manual(values = allcolour) +
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
      })
  })

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(out_dir,  
                  str_glue("spatial dim plot {grp}_annotation.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 


marker_p_ls <- 
  target_gene_ls %>%
  map(function(x) {
    cells_boundary_ls %>%
      imap(function(cells_boundary, sample) { # browser()
        dat <-
          cells_boundary %>%
          right_join(data_expr_full, by = "cell_id")
        
        p <- dat %>%
          arrange(.data[[x]]) %>% 
          as.data.frame() %>%
          ggplot(aes(x = vertex_x, y = -vertex_y)) + 
          geom_polygon(aes(fill = .data[[x]], group = cell_id))  +
          # scattermore::geom_scattermore(aes(alpha = .data[[x]]), pointsize = 0.5, show.legend = TRUE) +
          guides(alpha = "none") +
          theme_void() +
          coord_equal() +
          ggtitle(sample)+ 
          scale_fill_gradientn(
            colours = c("blue", "white", "red"),
            # limits = c(0, ceiling(max(dat_expr[-1]))),  # 统一颜色范围
            limits = c(0, 0.1), 
            name = "Expression"
          ) + 
        guides(
          fill = guide_colorbar(
            title = str_glue(x, " Expression"),
            barwidth = 0.5,
            barheight = 10
          )
        )
        p
        
      })
  })


marker_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(out_dir,  
                  str_glue("spatial dim plot {grp}.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 
