# 2025-12-09 spatial plot ---------------------
sp_output_filt <- "ROI_analysis/spatial plot"
dir.create(sp_output_filt, recursive = T)

all_data_meta <-
  all_data2[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan1 <- c("gray80", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan1)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - ROI.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })





### 添加一列新注释，区分ROI -------
all_data2@meta.data <- all_data2@meta.data %>% 
  mutate(annotation = as.character(annotation)) %>% 
  mutate(annotation_ROI = ifelse(ROI == "Other", "Other", annotation))

table(all_data2$annotation_ROI)


plan2 <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
           "#7F95D1", "#dee2d1", "gray80",  "#a5dff9", "#e8a0b8", "#ffeead", 
           "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
           "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
           "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
           "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
           "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
           "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
           "#005995", "#fa625f", "#f17d80")

all_data_meta <-
  all_data2[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_ROI, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan2)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })



# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_ROI.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })






### spatial plot of some gene ------
sp_genes <- c("Mapt")
sp_genes <- c("App")

sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data2, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    
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

# export
pdf(file.path("ROI_analysis/spatial plot", str_glue("spatial plot - Tau.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()




sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data2_filt, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    
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

# export
pdf(file.path("ROI_analysis/spatial plot", str_glue("spatial plot - Tau - ROI.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()



pdf(file.path("ROI_analysis/spatial plot", str_glue("Dot plot by group - Tau and App.pdf")), height = 6, width = 8)
DotPlot(all_data2, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()

pdf(file.path("ROI_analysis/spatial plot", str_glue("Dot plot by group - Tau and App - ROI.pdf")), height = 6, width = 8)
DotPlot(all_data2_filt, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()


pdf(file.path("ROI_analysis/spatial plot", str_glue("Dot plot by annotation - Tau and App.pdf")), height = 6, width = 8)
DotPlot(all_data2, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()

pdf(file.path("ROI_analysis/spatial plot", str_glue("Dot plot by annotation - Tau and App - ROI.pdf")), height = 6, width = 8)
DotPlot(all_data2_filt, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()






pdf(file.path("ROI_analysis/spatial plot", str_glue("Vlnplot by group - Tau and App.pdf")), height = 5, width = 10)
VlnPlot(all_data2, features = c("Mapt", "App"), group.by = "group", pt.size = 0)
dev.off()

pdf(file.path("ROI_analysis/spatial plot", str_glue("Vlnplot by group - Tau and App - ROI.pdf")), height = 5, width = 10)
VlnPlot(all_data2_filt, features = c("Mapt", "App"), group.by = "group", pt.size = 0)
dev.off()


pdf(file.path("ROI_analysis/spatial plot", str_glue("Vlnplot by annotation - Tau and App.pdf")), height = 5, width = 10)
VlnPlot(all_data2, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0)
dev.off()

pdf(file.path("ROI_analysis/spatial plot", str_glue("Vlnplot by annotation - Tau and App - ROI.pdf")), height = 5, width = 10)
VlnPlot(all_data2_filt, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0)
dev.off()












# 2025-12-31 spatial plot ---------------------
Date <- "2025-12-31"
sp_output_filt <- str_glue("ROI_analysis - {Date}/spatial plot")
dir.create(sp_output_filt, recursive = T)

all_data <- readRDS("data_for_analysis/all_data_h5ad_ROI_1230.rds")

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
cells_boundary_ls_filt <- cells_boundary_ls[!names(cells_boundary_ls) %in% c("B1", "B4")]
names(cells_boundary_ls_filt)

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan1 <- c("#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
           "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
           "#737495", "#bccf3d", "#de4307", "#f29c2b", "gray80", "#f6d04d", "#8bc24c", 
           "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
           "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
           "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
           "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
           "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
           "#005995", "#fa625f", "#f17d80")

sp_p_ls_CA1 <-
  cells_boundary_ls_filt %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI_CA1, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan1)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })

sp_p_ls_CA3 <-
  cells_boundary_ls_filt %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI_CA3, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan1)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })

sp_p_ls_DG <-
  cells_boundary_ls_filt %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI_DG, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan1)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })

sp_p_ls_Hippocampus <-
  cells_boundary_ls_filt %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI_Hippocampus, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan1)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })

# export image
sp_p_ls_CA1 %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - CA1.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })

sp_p_ls_CA3 %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - CA3.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })

sp_p_ls_DG %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - DG.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })

sp_p_ls_Hippocampus %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - Hippocampus.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })





### 添加一列新注释，区分ROI -------
all_data@meta.data <- all_data@meta.data %>% 
  mutate(annotation = as.character(annotation)) %>% 
  mutate(annotation_CA1 = ifelse(ROI_CA1 == "Other", "Other", annotation))
table(all_data$annotation_CA1)

all_data@meta.data <- all_data@meta.data %>% 
  mutate(annotation = as.character(annotation)) %>% 
  mutate(annotation_CA3 = ifelse(ROI_CA3 == "Other", "Other", annotation))
table(all_data$annotation_CA3)

all_data@meta.data <- all_data@meta.data %>% 
  mutate(annotation = as.character(annotation)) %>% 
  mutate(annotation_DG = ifelse(ROI_DG == "Other", "Other", annotation))
table(all_data$annotation_DG)

all_data@meta.data <- all_data@meta.data %>% 
  mutate(annotation = as.character(annotation)) %>% 
  mutate(annotation_Hippocampus = ifelse(ROI_Hippocampus == "Other", "Other", annotation))
table(all_data$annotation_Hippocampus)




all_data$ROI <- all_data$ROI_CA1

library(tibble)
library(dplyr)
meta.data.t <- 
  all_data[[c('ROI_CA3', 'ROI_DG')]] %>% 
  rownames_to_column('cells') %>%  
  right_join(
    all_data[['ROI_CA1']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    ROI_1 = ifelse((as.character(ROI_CA1) == 'Other'), as.character(ROI_CA3), as.character(ROI_CA1))) %>% 
  mutate(
    ROI = ifelse((as.character(ROI_1) == "Other"), as.character(ROI_DG), as.character(ROI_1))) %>% 
  dplyr::select(cells, ROI) %>%
  column_to_rownames('cells')

meta.data.t$ROI <- sub(".*?_", "", meta.data.t$ROI)

head(meta.data.t)
table(meta.data.t$ROI)

all_data <- AddMetaData(all_data, meta.data.t["ROI"])

pdf(file.path("mapping_annotation ROI umap.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "ROI", reduction = 'umap', pt.size = 2, repel = TRUE, cols = plan) 
dev.off()

saveRDS(all_data, file = "data_for_analysis/all_data_h5ad_ROI_0107.rds")



plan2 <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
           "#7F95D1", "gray80",  "#a5dff9", "#dee2d1", "#e8a0b8", "#ffeead", 
           "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
           "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
           "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
           "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
           "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
           "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
           "#005995", "#fa625f", "#f17d80")


Date <- "2025-12-31"
sp_output_filt <- str_glue("ROI_analysis - {Date}/spatial plot ROI annotation")
dir.create(sp_output_filt, recursive = T)

# 使用颜色名称而不是十六进制代码
cols_plan <- c(
  "Astrocytes" = "#ea7070",           # 红色
  "Endothelial" = "#fdc4b6",          # 浅橙色
  "GABAergic neurons" = "#e59572",    # 橙色
  "Glutamatergic neurons" = "#2694ab", # 青色
  "Microglial" = "#96ceb4",           # 浅青色
  "Oligo" = "#600473",                # 紫色
  "OPC" = "#7F95D1",                  # 蓝色
  "Other" = "gray80",                 # 灰色
  "Unknown" = "#a5dff9"               # 浅蓝色
)

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls_CA1 <-
  cells_boundary_ls_filt %>%
  imap(function(cells_boundary, project){ # browser()
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_CA1, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = cols_plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })

sp_p_ls_CA3 <- 
  cells_boundary_ls_filt %>% 
  imap(function(cells_boundary, project){ # browser()
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation_CA3, group = cell_id))  + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = cols_plan) + 
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })

sp_p_ls_DG <- 
  cells_boundary_ls_filt %>% 
  imap(function(cells_boundary, project){ # browser()
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation_DG, group = cell_id))  + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = cols_plan) + 
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) + 
      theme_void()
    p
  })

sp_p_ls_Hippocampus <- 
  cells_boundary_ls_filt %>% 
  imap(function(cells_boundary, project){ # browser()
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation_Hippocampus, group = cell_id)) + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = cols_plan) + 
      # coord_fixed(ratio = 1/image_wh_r) + 
      # ggsci::scale_fill_igv() + 
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) + 
      theme_void()
    p
  })



# export image
sp_p_ls_CA1 %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_CA1.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })

sp_p_ls_CA3 %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_CA3.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })

sp_p_ls_DG %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_DG.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })

sp_p_ls_Hippocampus %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_Hippocampus.pdf")),
      height = 12,
      width  = 12
    )
    print(p)
    dev.off()
  })






### spatial plot of some gene ------
Date <- "2025-12-31"
sp_output_filt <- str_glue("ROI_analysis - {Date}/spatial plot gene")
dir.create(sp_output_filt, recursive = T)

sp_genes <- c("Mapt")
sp_genes <- c("App")

sp_p_ls <- 
  cells_boundary_ls_filt %>%
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

# export
pdf(file.path(sp_output_filt, str_glue("spatial plot - App.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()




sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data2_filt, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    
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

# export
pdf(file.path("ROI_analysis/spatial plot", str_glue("spatial plot - Tau - ROI.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()



plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#e8a0b8", "#a5dff9", "#dee2d1", "#ffeead", "gray80",  
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


### annotation_CA1 ---------------
sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_CA1' # "annotation_CA3", "annotation_DG", "annotation_Hippocampus"
sp_genes <- c("Mapt")

sp_p_ls_Mapt <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")

    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })

sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_CA1' # "annotation_CA3", "annotation_DG", "annotation_Hippocampus"
sp_genes <- c("App")

sp_p_ls_App <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })


# export
pdf(file.path(sp_output_filt, str_glue("spatial plot CA1 - Mapt.pdf")), height = 25, width = 25)
print(sp_p_ls_Mapt)
dev.off()

pdf(file.path(sp_output_filt, str_glue("spatial plot CA1 - App.pdf")), height = 25, width = 25)
print(sp_p_ls_App)
dev.off()








### annotation_CA3 ---------------
sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_CA3' # "annotation_CA1", "annotation_DG", "annotation_Hippocampus"
sp_genes <- c("Mapt")

sp_p_ls_Mapt <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })

sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_CA3' # "annotation_CA1", "annotation_DG", "annotation_Hippocampus"
sp_genes <- c("App")

sp_p_ls_App <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })


# export
pdf(file.path(sp_output_filt, str_glue("spatial plot CA3 - Mapt.pdf")), height = 25, width = 25)
print(sp_p_ls_Mapt)
dev.off()

pdf(file.path(sp_output_filt, str_glue("spatial plot CA3 - App.pdf")), height = 25, width = 25)
print(sp_p_ls_App)
dev.off()










### annotation_DG ---------------
sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_DG' # "annotation_CA1", "annotation_CA3", "annotation_Hippocampus"
sp_genes <- c("Mapt")

sp_p_ls_Mapt <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })



sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_DG' # "annotation_CA1", "annotation_CA3","annotation_DG", "annotation_Hippocampus"
sp_genes <- c("App")

sp_p_ls_App <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })


# export
pdf(file.path(sp_output_filt, str_glue("spatial plot DG - Mapt.pdf")), height = 25, width = 25)
print(sp_p_ls_Mapt)
dev.off()

pdf(file.path(sp_output_filt, str_glue("spatial plot DG - App.pdf")), height = 25, width = 25)
print(sp_p_ls_App)
dev.off()











### annotation_Hippocampus ---------------
sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_Hippocampus' # "annotation_CA1", "annotation_CA3","annotation_DG", "annotation_Hippocampus"
sp_genes <- c("Mapt")

sp_p_ls_Mapt <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })



sp_cells <- c("Glutamatergic neurons", "Astrocytes", "GABAergic neurons", "Endothelial", "OPC", "Oligo", "Microglial", "Unknown")
annotation_col <- 'annotation_Hippocampus' # "annotation_CA1", "annotation_CA3","annotation_DG", "annotation_Hippocampus"
sp_genes <- c("App")

sp_p_ls_App <-
  cells_boundary_ls_filt %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id") %>% 
      right_join(dat_expr, by = "cell_id")
    
    all_cells <- dat[[annotation_col]] %>% unique()
    others_cells <- all_cells %>% setdiff(sp_cells)
    # color_for_cells <- c(
    #   # structure(names = sp_cells, plan[seq_along(sp_cells)]), # sp_cells
    #   structure(names = sp_cells, plan[1:8]), # sp_cells
    #   structure(names = others_cells, rep('gray80', times = length(others_cells))), # others cells
    #   Others = 'gray80'
    # )
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        # change expr of other cells to NA/0
        dat[[annotation_col]][ dat[[annotation_col]] %in% others_cells ] <- 'Others'
        dat[[gene]][ dat[[annotation_col]] == 'Others'] <- NA
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(
              fill = .data[[gene]], 
              # color = .data[[annotation_col]],
              group = cell_id), 
            linewidth = 0.2, 
            alpha = 0.8
          )  +
          scale_y_reverse() +
          coord_equal() +
          # scale_color_manual(values = color_for_cells, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3))) +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })


# export
pdf(file.path(sp_output_filt, str_glue("spatial plot Hippocampus - Mapt.pdf")), height = 25, width = 25)
print(sp_p_ls_Mapt)
dev.off()

pdf(file.path(sp_output_filt, str_glue("spatial plot Hippocampus - App.pdf")), height = 25, width = 25)
print(sp_p_ls_App)
dev.off()









all_data_CA1 <- all_data[, !all_data$ROI_CA1 %in% c("Other")]
all_data_CA3 <- all_data[, !all_data$ROI_CA3 %in% c("Other")]
all_data_DG <- all_data[, !all_data$ROI_DG %in% c("Other")]
all_data_Hippocampus <- all_data[, !all_data$ROI_Hippocampus %in% c("Other")]
dim(all_data_CA1)
dim(all_data_CA3)
dim(all_data_DG)
dim(all_data_Hippocampus)


### 气泡图 ------------
pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by group - Mapt and App.pdf")), height = 6, width = 8)
DotPlot(all_data, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by group - Mapt and App - CA1.pdf")), height = 6, width = 8)
DotPlot(all_data_CA1, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by group - Mapt and App - CA3.pdf")), height = 6, width = 8)
DotPlot(all_data_CA3, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by group - Mapt and App - DG.pdf")), height = 6, width = 8)
DotPlot(all_data_DG, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by group - Mapt and App - Hippocampus.pdf")), height = 6, width = 8)
DotPlot(all_data_Hippocampus, features = c("Mapt", "App"), dot.scale = 10, group.by = "group")
dev.off()




pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by annotation - Tau and App.pdf")), height = 6, width = 8)
DotPlot(all_data, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by annotation - Tau and App - CA1.pdf")), height = 6, width = 8)
DotPlot(all_data_CA1, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by annotation - Tau and App - CA3.pdf")), height = 6, width = 8)
DotPlot(all_data_CA3, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by annotation - Tau and App - DG.pdf")), height = 6, width = 8)
DotPlot(all_data_DG, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Dot plot by annotation - Tau and App - Hippocampus.pdf")), height = 6, width = 8)
DotPlot(all_data_Hippocampus, features = c("Mapt", "App"), dot.scale = 10, group.by = "annotation")
dev.off()





### 小提琴图 ------------
pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by group - Tau and App.pdf")), height = 5, width = 10)
VlnPlot(all_data, features = c("Mapt", "App"), group.by = "group", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by group - Tau and App - CA1.pdf")), height = 5, width = 10)
VlnPlot(all_data_CA1, features = c("Mapt", "App"), group.by = "group", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by group - Tau and App - CA3.pdf")), height = 5, width = 10)
VlnPlot(all_data_CA3, features = c("Mapt", "App"), group.by = "group", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by group - Tau and App - DG.pdf")), height = 5, width = 10)
VlnPlot(all_data_DG, features = c("Mapt", "App"), group.by = "group", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by group - Tau and App - Hippocampus.pdf")), height = 5, width = 10)
VlnPlot(all_data_Hippocampus, features = c("Mapt", "App"), group.by = "group", pt.size = 0, cols = plan)
dev.off()



pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by annotation - Tau and App.pdf")), height = 5, width = 12)
VlnPlot(all_data, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by annotation - Tau and App - CA1.pdf")), height = 5, width = 12)
VlnPlot(all_data_CA1, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by annotation - Tau and App - CA3.pdf")), height = 5, width = 12)
VlnPlot(all_data_CA3, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by annotation - Tau and App - DG.pdf")), height = 5, width = 12)
VlnPlot(all_data_DG, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0, cols = plan)
dev.off()

pdf(file.path("ROI_analysis - 2025-12-31/spatial plot gene/", str_glue("Vlnplot by annotation - Tau and App - Hippocampus.pdf")), height = 5, width = 12)
VlnPlot(all_data_Hippocampus, features = c("Mapt", "App"), group.by = "annotation", pt.size = 0, cols = plan)
dev.off()



















