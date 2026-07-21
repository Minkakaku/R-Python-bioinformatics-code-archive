#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2_group.rds --species human -g Group -c annotation2 --ident1 Mild --ident2 Normal -o Mild_vs_Normal_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2_group.rds --species human -g Group -c annotation2 --ident1 Severe --ident2 Normal -o Severe_vs_Normal_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2_group.rds --species human -g Group -c annotation2 --ident1 Severe --ident2 Mild -o Severe_vs_Mild_DEG &

#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2_group.rds --species human -g Group -c annotation2 --ident1 Mild --ident2 Normal -o Mild_vs_Normal_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2_group.rds --species human -g Group -c annotation2 --ident1 Severe --ident2 Normal -o Severe_vs_Normal_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2_group.rds --species human -g Group -c annotation2 --ident1 Severe --ident2 Mild -o Severe_vs_Mild_DEG &




# 映射到空间切片
spatial_out_dir <- "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/spatial" 
dir.create(spatial_out_dir)
spatial_out_dir_1 <- "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/spatial/cell_spatial_1" 
dir.create(spatial_out_dir_1)


all_data <- readRDS("/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2.rds")

group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
group_levels <- c("Normal", "Mild", "Severe")
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


all_data@meta.data <- all_data@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map),
    Group = factor(Group, levels = group_levels)
  )

table(all_data$Group)

all_data$core_id <- all_data$cell_id2
head(all_data)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
cells_boundary_ls <- map(cells_boundary_ls, ~mutate(.x, core_id = str_remove(cell_id, "^TMA_")))
#cells_boundary_df <- cells_boundary_ls %>% bind_rows()
#cells_boundary_df <- cells_boundary_df %>%
#  mutate(core_id = str_remove(cell_id, "^TMA_")) # 移除开头 TMA_
head(cells_boundary_ls)

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') 
head(all_data_meta)

all_data_meta$sample %>% table()


# 1. 定义需要展示的细胞类型（对应table里非Others）
keep <- c("CXCL12+Fibroblast", "Podocyte", "PTS1", "PTS2", "PTS3")

# 2. 提取IGV调色板前10个颜色，匹配10种细胞
igv_col <- ggsci::pal_igv()(length(keep))
# 命名映射：细胞名 = IGV颜色
cell_color_map <- setNames(igv_col, keep)
# 追加Others灰色
cell_color_map[["Others"]] <- "#dcdcdc"

cell_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "Group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "core_id")
    
    dat <- dat %>%
      mutate(plot_group = ifelse(annotation2 %in% keep, annotation2, "Others"))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation2, group = core_id))  +
      scale_y_reverse() +
      coord_equal() + 
      # 自定义颜色：为keep列表中的细胞分配颜色，"Others"设为灰色
      scale_fill_manual(
        values = cell_color_map,
        limits = c(keep, "Others")
      ) +
      theme_void() +
      labs(fill = NULL, title = sample)
    p
  })

cell_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(spatial_out_dir_1, str_glue("Spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })



spatial_out_dir_2 <- "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/spatial/cell_spatial_2" 
dir.create(spatial_out_dir_2)


cell_sp_p_ls2 <- 
  cells_boundary_ls %>% # .[1] %>%
  imap(function(cells_boundary, sample){ # browser()
    
    anno_data   <- all_data_meta[c("sample", "Group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "core_id")
    
    p_ls <-
      keep %>% # .[1] %>%
      set_names(., .) %>%
      map(function(cn){
        dat2 <- dat
        
        dat2$annotation2 <- as.character(dat2$annotation2)
        dat2$annotation2[dat2$annotation2 != cn] <- "Others"
        
        p <- dat2 %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = annotation2, group = core_id))  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cell_color_map) +
          theme_void() +
          labs(fill = NULL, title = paste(sample, group_info, cn, sep = "_") )
        p
      })
    
    p_ls
  })


# export image
cell_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(spatial_out_dir_2, str_glue("Spatial plot {name} - split by celltype.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })




# 基因空间图输出目录
spatial_out_dir_3 <- "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/spatial/gene_spatial" 
dir.create(spatial_out_dir_3, recursive = TRUE, showWarnings = FALSE)

some_genes <- c("HOXB7", "ECRG4", "RGN")
annotation_col <- 'annotation2'

# ---------------------- 修复FetchData报错：分开提取表达+meta ----------------------
# 1. 提取基因表达
plot_dat <- FetchData(all_data, vars = c(annotation_col, some_genes)) %>%
  rownames_to_column("core_id") %>%
  # 删掉 "-" 及后面所有字符，还原成 aaaaagpe-1
  mutate(core_id = str_remove(core_id, "-[^-]+$"))

# ---------------------- 基于cells_boundary_ls绘制细胞多边形基因图 ----------------------
p_sp_gene_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ 
    # 边界数据 + 基因表达合并
    dat <- cells_boundary %>%
      right_join(plot_dat, by = "core_id")
    # 统计无边界坐标的细胞数量
    #na_coord_count <- sum(is.na(dat$vertex_x))
    #message(str_glue("样本 {sample} 仅存在表达、无细胞边界的细胞数量：{na_coord_count}"))
    # 遍历每个目标基因生成单张图
    gene_p_ls <-
      some_genes %>%
      setNames(., .) %>%
      map(function(gene_x){
        dat2 <- dat
        
        # 分位数截断极值（保留你原有防极端高表达逻辑）
        q_cutoff <- 0.9
        n <- 0
        while( near(quantile(dat2[[gene_x]], q_cutoff, na.rm = TRUE), 0) ){
          q_cutoff <- q_cutoff + 9/10^(n + 2)
          print(q_cutoff)
          n <- n + 1
          if(n > 3) break
        }
        if(n <= 3){
          dat2[[gene_x]] <- pmin(dat2[[gene_x]], quantile(dat2[[gene_x]], q_cutoff, na.rm = TRUE))
        }
        
        # 表达量范围
        val_rng <- range(dat2[[gene_x]], na.rm = TRUE)
        
        # 细胞多边形绘图，替换原来geom_tile
        p <- dat2 %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(fill = .data[[gene_x]], alpha = .data[[gene_x]], group = core_id),
            linewidth = 0.03 # 细胞轮廓细线，可选调0隐藏边框
          ) +
          scale_y_reverse() +
          coord_equal() +
          # 红白连续渐变，图例仅保留fill
          scale_fill_gradient(
            low = "blue", high = "red",
            limits = val_rng,
            name = "Expression"
          ) +
          scale_alpha(range = c(0.2, 1), guide = "none") + # alpha不显示图例
          theme_void() +
          theme(panel.background = element_rect(fill = "white", color = NA)) +
          ggtitle(str_glue("{sample} - {gene_x}"))
        return(p)
      })
    return(gene_p_ls)
  })

# ---------------------- 导出PDF（和细胞分型图逻辑完全一致） ----------------------
p_sp_gene_ls %>% 
  iwalk(function(p_ls, sam_name){
    # 空坐标跳过
    if(all(is.na(p_ls[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p_ls[[1]]$data$vertex_x))) {
      warning(str_glue("样本 {sam_name} vertex_x存在NA值"), immediate. = TRUE)
    }
    # 输出PDF，尺寸和细胞边界图统一20×20
    pdf(
      file.path(spatial_out_dir_3, str_glue("Sptial plot of gene in sample {sam_name}.pdf") ),
      width  = 20,
      height = 20
    )
    p_ls %>% walk(print)
    dev.off()
  })



spatial_out_dir_4 <- "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/spatial/gene_spatial_2" 
dir.create(spatial_out_dir_4, recursive = TRUE, showWarnings = FALSE)
# 指定目标细胞类型
target_cell_type <- "CXCL12+Fibroblast"

# ---------------------- 修复FetchData报错：分开提取表达+meta ----------------------
# 1. 提取基因表达
plot_dat <- FetchData(all_data, vars = c(annotation_col, some_genes)) %>%
  rownames_to_column("core_id") %>%
  # 删掉 "-" 及后面所有字符，还原成 aaaaagpe-1
  mutate(core_id = str_remove(core_id, "-[^-]+$"))

# ---------------------- 基于cells_boundary_ls绘制细胞多边形基因图 ----------------------
p_sp_gene_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ 
    # 边界数据 + 基因表达合并
    dat <- cells_boundary %>%
      right_join(plot_dat, by = "core_id")
    
    dat <- dat %>%
      mutate(across(all_of(some_genes), ~ifelse(.data[[annotation_col]] == target_cell_type, .x, 0)))
    # 统计无边界坐标的细胞数量
    #na_coord_count <- sum(is.na(dat$vertex_x))
    #message(str_glue("样本 {sample} 仅存在表达、无细胞边界的细胞数量：{na_coord_count}"))
    # 遍历每个目标基因生成单张图
    gene_p_ls <-
      some_genes %>%
      setNames(., .) %>%
      map(function(gene_x){
        dat2 <- dat
        
        # 分位数截断极值（保留你原有防极端高表达逻辑）
        q_cutoff <- 0.9
        n <- 0
        while( near(quantile(dat2[[gene_x]], q_cutoff, na.rm = TRUE), 0) ){
          q_cutoff <- q_cutoff + 9/10^(n + 2)
          print(q_cutoff)
          n <- n + 1
          if(n > 3) break
        }
        if(n <= 3){
          dat2[[gene_x]] <- pmin(dat2[[gene_x]], quantile(dat2[[gene_x]], q_cutoff, na.rm = TRUE))
        }
        
        # 表达量范围
        val_rng <- range(dat2[[gene_x]], na.rm = TRUE)
        
        # 细胞多边形绘图，替换原来geom_tile
        p <- dat2 %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(
            aes(fill = .data[[gene_x]], alpha = .data[[gene_x]], group = core_id),
            linewidth = 0.03 # 细胞轮廓细线，可选调0隐藏边框
          ) +
          scale_y_reverse() +
          coord_equal() +
          # 红白连续渐变，图例仅保留fill
          scale_fill_gradient(
            low = "blue", high = "red",
            limits = val_rng,
            name = "Expression"
          ) +
          scale_alpha(range = c(0.2, 1), guide = "none") + # alpha不显示图例
          theme_void() +
          theme(panel.background = element_rect(fill = "white", color = NA)) +
          ggtitle(str_glue("{sample} - {gene_x} in {target_cell_type}"))
        return(p)
      })
    return(gene_p_ls)
  })

# ---------------------- 导出PDF（和细胞分型图逻辑完全一致） ----------------------
p_sp_gene_ls %>% 
  iwalk(function(p_ls, sam_name){
    # 空坐标跳过
    if(all(is.na(p_ls[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p_ls[[1]]$data$vertex_x))) {
      warning(str_glue("样本 {sam_name} vertex_x存在NA值"), immediate. = TRUE)
    }
    # 输出PDF，尺寸和细胞边界图统一20×20
    pdf(
      file.path(spatial_out_dir_4, str_glue("Sptial plot of gene in CXCL12+Fibroblast {sam_name}.pdf") ),
      width  = 20,
      height = 20
    )
    p_ls %>% walk(print)
    dev.off()
  })


# 加载包
library(Seurat)
library(tidyverse)
library(dplyr)
library(tibble)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

all_data <- readRDS("/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/alldata/alldata_filter_annotation2.rds")

all_data@meta.data <- all_data@meta.data %>%
  rownames_to_column("cell_id")

group_map <- c(
  "N"       = "Normal",
  "M_igAN"  = "Mild",
  "S_igAN"  = "Severe"
)
group_levels <- c("Normal", "Mild", "Severe")
group_col <- c(
  "Normal"  = "#2E86AB",
  "Mild"    = "#A23B72",
  "Severe"  = "#F18F01"
)


all_data@meta.data <- all_data@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map),
    Group = factor(Group, levels = group_levels)
  )

table(all_data$Group)

# 1. 定义6个基因集列表
gene_sets <- list(
  "ECM_synthesis_assembly" = c("CCN1", "CCN2", "COL4A4", "FBLN5", "TIMP2", "CCN4", "COL4A2"),
  "Anti_fibrosis_protective" = c("ECRG4", "NR4A1", "SOCS2", "SOCS3", "GATA6", "FOXO3"),
  "Pro_fibrosis_pro_inflammation" = c("LTBP1", "TNC", "MYLK", "EDNRB", "PDGFRA", "PDGFRB", "CTNNB1", "ZEB2"),
  "Proliferation_survival_translation" = c("EEF1G", "EIF3C", "HNRNPM", "HNRNPUL1", "NONO", "HSPD1", "EIF4A1"),
  "Migration_contractility" = c("MYLK", "TUBB", "ACTN4", "PALLD", "SNAI2"),
  "Metabolism_stress_response" = c("FMO3", "CYP1B1", "AOC3", "PRNP", "PHB2", "SOD1")
)

# 2. 分组配色（对应附图蓝色=正常，橙色系=病变）
group_color <- c(
  Normal = "#4477CC",
  Mild = "#EE9955",
  Severe = "#DD5533"
)

# 3. 输出目录
heatmap_out <- "/mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/spatial/gene_set_heatmap"
dir.create(heatmap_out, recursive = T, showWarnings = F)

# 4. 提取表达矩阵+样本分组注释
## 提取标准化表达矩阵（data槽，log归一化表达，和附图逻辑一致）
expr_mat  <- all_data %>% JoinLayers() %>% LayerData(layer = "data")


## 样本元数据：样本名+分组
sample_meta <- all_data@meta.data %>%
  distinct(sample, Group) %>%
  arrange(Group, sample) # 固定样本顺序

# 顶部样本分组注释条（复刻图上方彩色横条）
top_anno <- HeatmapAnnotation(
  Group = sample_meta$Group,
  col = list(Group = group_color),
  show_legend = T,
  annotation_name_side = "left",
  height = unit(3, "mm")
)

# 5. 循环每个基因集，批量出图
imap(gene_sets, function(gene_vec, set_name) {
  # 过滤存在于表达矩阵的基因（去除缺失基因）
  gene_valid <- intersect(gene_vec, rownames(expr_mat))
  if (length(gene_valid) == 0) {
    message(str_glue("【{set_name}】无有效基因，跳过"))
    return(NULL)
  }
  
  # 提取当前基因集表达矩阵，按样本取均值（关键：附图是样本水平汇总，不是单细胞！）
  expr_long <- expr_mat[gene_valid, , drop=F] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell_id", values_to = "expr") %>%
    left_join(all_data@meta.data %>% select(cell_id, sample), by = "cell_id")
  
  # 按样本+基因取平均表达（和原图一致：每个格子是该样本内该基因平均表达）
  sample_gene_mean <- expr_long %>%
    group_by(sample, gene) %>%
    summarise(mean_expr = mean(expr, na.rm = T), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = mean_expr) %>%
    column_to_rownames("gene")
  
  # 匹配样本顺序，保证注释和矩阵列对齐
  sample_gene_mean <- sample_gene_mean[, sample_meta$sample, drop=F]
  
  # --- 关键修改1：转换为矩阵 ---
  sample_gene_mean_mat <- as.matrix(sample_gene_mean)
  
  # 检查是否有表达值全为0的基因，如果有，加一个小偏移避免全0行
  # 但实际上热图可以处理，只是颜色可能不显示
  
  # --- 关键修改2：动态计算热图高度，设置最小值 ---
  n_genes <- nrow(sample_gene_mean_mat)
  # 每个基因至少0.3cm，但总高度至少3cm
  row_height <- max(1, 10 / n_genes)
  total_height <- max(5, n_genes * row_height + 1)  # +1 为标题留空间
  
  # 颜色标尺
  col_fun <- colorRamp2(
    c(min(sample_gene_mean_mat), max(sample_gene_mean_mat)), 
    c("white", "red")
  )
  
  # 绘制热图
  ht <- Heatmap(
    matrix = sample_gene_mean_mat,  # 使用矩阵
    name = "Expression",
    col = col_fun,
    top_annotation = top_anno,
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    use_raster = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_width = unit(16, "cm"),
    heatmap_height = unit(total_height, "cm"),  # 使用计算的高度
    column_title = str_glue("{set_name} gene expression heatmap"),
    column_title_gp = gpar(fontsize = 14)
  )
  
  # 导出PDF
  pdf_path <- file.path(heatmap_out, str_glue("Heatmap_{set_name}.pdf"))
  pdf(pdf_path, width = 10, height = max(5, total_height + 1))
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  message(str_glue("完成：{pdf_path}"))
})
