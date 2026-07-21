#Epithelial subtype
#Epithelial subtype 在infercnv里的亚分型中选择resolution 0.2
Epithelial<-readRDS('epi_subcluster/epi_sub.rds')
Idents(Epithelial)<-"RNA_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
FeaturePlot(Epithelial,c("ESR1","PGR","ERBB2","EGFR","ERBB3"),reduction = 'harmonyumap',max.cutoff = 1)
DotPlot(Epithelial,c("ESR1","PGR","ERBB2","EGFR","ERBB3"),group.by = 'RNA_snn_res.0.2',col.min = 0) + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu")))
top10 <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DotPlot(Epithelial,features = unique(top10$gene) ,assay='RNA')+RotatedAxis()


pdf('/home/rstudio/epi_subcluster/RNA_snn_res.0.2_cluster.pdf', height = 6, width = 8)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2_cluster", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

scales::show_col(RColorBrewer::brewer.pal(Inf, "RdYlBu"))

# cell props barplot  --------

Epithelial_meta <-
  Epithelial[[]] %>%
  rownames_to_column("cell_id")

Epithelial_meta$project <- factor(Epithelial_meta$project,levels=c("GX_1","MXH_1","HXF_1","LS_1","GYF_1","GX_2","MXH_2","HXF_2","LS_2","GYF_2","GX_3","MXH_3","HXF_3","LS_3","GYF_3"))
Epithelial_meta$project <- factor(Epithelial_meta$project,levels=c("GX_1","GX_2","GX_3","MXH_1","MXH_2","MXH_3","HXF_1","HXF_2","HXF_3","LS_1","LS_2","LS_3","GYF_1","GYF_2","GYF_3"))

p_col_proj <-
  Epithelial_meta %>%
  group_by(project, RNA_snn_res.0.2_cluster) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = RNA_snn_res.0.2_cluster)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(dr_dir, "Epithelial project prop.pdf"),
  height = 5,
  width  = 8
)


##映射------
library(tibble)
library(data.table)
meta.data <-
  Epithelial[['RNA_snn_res.0.2']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation3 = ifelse((as.character(annotation)=='Epithelial'),paste0('cluster_',as.character(RNA_snn_res.0.2)), as.character(annotation)))%>%
  dplyr::select(cells, annotation3) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data['annotation3'])

all_data_meta<- all_data[[]]

##spatial
sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project) {
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray90") +
      scale_y_reverse() +
      coord_equal() +
      theme_void()
    
    # 修正后的代码
    p_ls <- p +
      geom_polygon(
        data = dat[which(dat$annotation == 'Epithelial'), ], 
        aes(x = vertex_x, y = vertex_y, fill = annotation3, group = cell_id)
      ) +
      scale_fill_manual(values = plan)
      # scale_fill_manual(values = c(
      #   "#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF", 
      #   "#BA6338FF", "#802268FF", "#6BD76BFF", "#D595A7FF", "#5DB1DDFF", 
      #   "#7A65A5FF", "#E4AF69FF"
      # ))
    
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path('modify/plot modify/epi_spatial', str_glue("spatial plot {grp} Epithelial.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })
