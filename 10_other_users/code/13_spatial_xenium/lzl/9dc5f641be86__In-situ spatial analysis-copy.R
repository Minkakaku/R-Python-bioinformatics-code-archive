library(tidyverse)
library(Seurat)

# set default ---------
out_dir <- 'Celltype_Spatial_1'
dir.create(out_dir)

# data loading ----------
#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
#all_data@meta.data %>% head()

all_data <- readRDS('/mnt/NexData_file/home/lzl/case/XE2602040886-02/output/data/adata_filtered_add_celltype.rds')

all_data$celltype <- as.character(all_data$celltype)
all_data$celltype[all_data$celltype == "Column cytotrophoblast"] <- "Column Cytotrophoblast"
all_data$celltype <- as.factor(all_data$celltype)


all_data$annotation <- as.character(all_data$annotation)
all_data$annotation[all_data$annotation == "Column cytotrophoblast"] <- "Column Cytotrophoblast"
all_data$annotation <- as.factor(all_data$annotation)


#all_data_1 <- readRDS('/mnt/NexData_file/home/lzl/case/XE2602040886-02/output/data/adata_celltype_all.rds')
# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")

# normalize
all_data <- all_data %>% NormalizeData()

#all_data@assays$RNA@layers$
 


# spatial dat
# when all_data is transformed from a h5ad, spatial coord can be acquired by `Embeddings(all_data[['spatial']])`
# or, can be acquired from all_data@images... when all_data is VisiumHD
# or, can be acquired from cells.parquet from Xenium orignial directory when all_data is Xenium
some_genes <- c("PRG2","TNFSF10","HCAR2","HPGD","HCAR3","FSTL1","LTBP1","MXRA7","MCAM","FLT4","LAIR2","LOX","MMP12")
rownames(expr_matrix)
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
    vars = c(annotation_col,some_genes)
  ) %>%
  rownames_to_column('cell_id')



# plot celltype ------------------
cell_colors <- c(
  "Column Cytotrophoblast" = "#1f77b4",
  "Erythroblast" = "#ff7f0e",
  "Extravillous Trophoblast" = "#279e68",
  "Hofbauer Cell" = "#d62728",
  "Macrophage" = "#aa40fc",
  "Myofibroblast-like Fibroblast" = "#8c564b",
  "Renin-Angiotensin System Fibroblast" = "#e377c2",
  "Syncytiotrophoblast" = "#b5bd61",
  "Syncytiotrophoblast Precursor" = "#17becf",
  "Vascular Endothelial Cell(VEC)" = "#aec7e8",
  "Villous Cytotrophoblast" = "#ffbb78",
  "Villous Mesenchymal Stem Cell(CV-MSCs)" = "#98df8a"
)


p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 8
    width <- 8
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
      scale_fill_manual(values = cell_colors) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })

#p_sp_ct_ls[4]
pdf(
  file.path(out_dir, "1. spatial plot of celltype in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()

#####alltype循环
all_cell_types <- names(cell_colors)

# 遍历每种细胞类型
for(target_cell_type in all_cell_types){
  
  message("正在绘制：", target_cell_type)
  
  # 强制锁定当前细胞类型和颜色（解决循环+imap冲突！）
  current_cell <- target_cell_type
  current_color <- cell_colors[[current_cell]]
  
  p_sp_ct_ls <- cells_position_ls %>%
    imap(function(coords, sam){
      
      # 强制内部使用锁定的变量，不受循环影响
      ct <- current_cell
      #cc <- current_color
      
      dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
      
      # 高亮分组
      dat2$highlight <- ifelse(dat2[[annotation_col]] == ct, ct, "Other")
      
      dat2$highlight <- factor(
        dat2$highlight,
        levels = c(ct, "Other")
      )
      
      # 绘图
      p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
        geom_tile(aes(fill = highlight,alpha = highlight), height = 8, width = 8) +
        theme_void() +
        coord_equal() +
        scale_y_reverse() +
        ##透明度
        scale_alpha_manual(
          values = setNames(c(1, 0.3), c(ct, "Other"))  # 目标不透明，其他15%透明
        ) +
        # 核心：强制命名颜色（解决ggplot不认变量问题）
        scale_fill_manual(
          values = setNames(c(current_color, "gray60"), c(ct, "Other")),
          guide = guide_legend(override.aes = list(alpha = c(1, 1)))
        ) +
        
        theme(
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle(sam)
      
      return(p)
    })
  
  # 输出PDF
  pdf(
    file.path(out_dir, paste0("spatial plot of ", target_cell_type, ".pdf")),
    height = 8, width = 8
  )
  walk(p_sp_ct_ls, print)
  dev.off()
}

message("全部完成！已输出每个细胞类型独立PDF！")
out_dir <- "EVT Gene list DEG"
file.remove(list.files(path = out_dir, pattern = "sptial.*\\.pdf$", full.names = TRUE))
# spatial plot of some gene ---------
some_genes <- c("DIO2" ,"LAIR2" ,"MFAP5" , "HPGD"  , "HSD3B1")
p_sp_gene_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    some_genes %>%
      setNames(., .) %>% # .[1] %>%
      map(function(x){
        # sub
        dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
        
        stopifnot(x %in% colnames(dat2))

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
        height <- 8
        width <- 8
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(fill = "gray60", height = height, width = width) + 
          geom_tile(
            aes(fill = .data[[x]], alpha = .data[[x]]), #
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

p_sp_gene_ls[[1]]
# export 
p_sp_gene_ls %>%
  iwalk(function(p_ls, sam){
    pdf(
      file.path(out_dir, str_glue("sptial plot of EVT marker in sample {sam}.pdf") ),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })




gene_colors <- c(
  "PRG2"    = "#00E5FF",  # 亮青 (Cyan-ish) - 类似 DAPI/Alexa488 混合感
  "TNFSF10" = "#39FF14",  # 荧光绿 (Neon Green) - 经典 FITC 通道
  "HCAR2"   = "#FFEA00",  # 纯黄 (Pure Yellow) - 高亮，不偏橙
  "HPGD"    = "#FF0055",  # 荧光粉 (Hot Pink) - 类似 TRITC/mCherry
  "HCAR3"   = "#FF9100",  # 亮橙 (Vivid Orange) - 区分于黄色和粉色
  "FSTL1"   = "#D500F9",  # 电紫 (Electric Purple) - 类似 Cy5/Alexa647
  "LTBP1"   = "#00B0FF",  # 天蓝 (Sky Blue) - 区分于青色，更偏蓝
  "MXRA7"   = "#FF3D00",  # 朱红 (Vermilion) - 高饱和，不偏黄
  "MCAM"    = "#1DE9B6",  # 薄荷青 (Mint) - 区分于纯绿，带一点蓝调
  "FLT4"    = "#FFFF00",  # 正黄 (Pure Yellow) - 如果和HCAR2太近，可改用 #B2FF59 (黄绿)
  "LAIR2"   = "#7C4DFF",  # 靛紫 (Indigo) - 区分于电紫，更偏蓝紫
  "LOX"     = "#F50057",  # 玫红 (Rose) - 区分于荧光粉，稍微偏冷
  "MMP12"   = "#FFAB00"   # 琥珀 (Amber) - 区分于纯黄和橙色，更有质感
)
# 修改绘图代码，为每个基因使用指定颜色
p_sp_gene_ls <-
  cells_position_ls %>%
  imap(function(coords, sam){
    some_genes %>%
      setNames(., .) %>%
      map(function(x){
        dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
        stopifnot(x %in% colnames(dat2))
        
        # quantile cutoff
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
        
        val_rng <- dat2[[x]] %>% range()
        height <- 8
        width <- 8
        
        # 获取当前基因的颜色
        gene_color <- gene_colors[x]
        
        p <- dat2 %>%
          ggplot(aes(x = SPATIAL_1, y = SPATIAL_2)) +
          geom_tile(fill = "gray60", height = height, width = width, alpha = 0.3) + 
          geom_tile(
            aes(fill = .data[[x]], alpha = .data[[x]]),
            height = height, 
            width = width
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # 使用基因指定颜色：白→基因色
          scale_fill_gradient(low = "white", high = gene_color, limits = val_rng) +
          scale_alpha_continuous(range = c(0.25, 1), guide = "none") +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(str_glue("{sam} - {x}"))
        
        p
      })
  })

# 导出
p_sp_gene_ls %>%
  iwalk(function(p_ls, sam){
    pdf(
      file.path(out_dir, str_glue("spatial plot of EVT deg gene in sample {sam}.pdf")),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })

##
# C组findallmarker差异基因和monocle分化基因交集
c("MCAM","FLT4","LAIR2","LOX","MMP12")
# T组findallmarker差异基因和monocle分化基因交集
c("PRG2","TNFSF10","HCAR2","HPGD","HCAR3","FSTL1","LTBP1","MXRA7")

#
Layers(all_data)
all_data[["RNA"]]@layers$some_genes <- NULL
all_data <- JoinLayers(
  all_data
)
meta <- all_data@meta.data
expr_matrix <- GetAssayData(all_data, assay = "RNA", slot = "data")

meta_sub <- subset(meta, annotation == "Extravillous Trophoblast")
cells_sub <- rownames(meta_sub)

# 表达矩阵也只取这些细胞
expr_sub <- expr_matrix[, cells_sub]

# ==============================================================================
# 1. 批量做两组差异检验（Wilcoxon 秩和检验）
# ==============================================================================
diff_results <- lapply(some_genes, function(g) {
  if (!g %in% rownames(expr_sub)) {
    warning(paste(g, "不在表达矩阵中！"))
    return(NULL)
  }
  
  expr <- expr_sub[g, ] %>% as.vector()
  df <- data.frame(
    gene = g,
    group = meta_sub$group,
    expr = expr,
    check.names = F
  )
  
  # 两组Wilcox检验
  groups <- unique(meta_sub$group) # 这里我帮你修正为 meta_sub（更严谨）
  if (length(groups) != 2) stop("必须只有两个分组！")
  
  g1 <- df$expr[df$group == groups[1]]
  g2 <- df$expr[df$group == groups[2]]
  
  wt <- wilcox.test(g1, g2)
  
  # ===================== 【新增】计算两组中位数 =====================
  median1 <- round(median(g1), 4)
  median2 <- round(median(g2), 4)
  
  data.frame(
    gene = g,
    group1 = groups[1],
    group2 = groups[2],
    median1 = median1,    # 新增：第一组中位数
    median2 = median2,    # 新增：第二组中位数
    pval = wt$p.value,
    padj = p.adjust(wt$p.value, method = "BH"),
    stringsAsFactors = F
  )
}) %>% bind_rows()

head(diff_results)


expression_ratio_stats <- lapply(some_genes, function(g) {
  
  if (!g %in% rownames(expr_sub)) return(NULL)
  
  # 该基因在该细胞中的表达
  exp <- as.numeric(expr_sub[g, ])
  group <- meta_sub$group
  group1 <- unique(group)[1]
  group2 <- unique(group)[2]
  
  # 阳性细胞比例（>0 就算表达）
  ratio1 <- sum(exp[group == group1] > 0) / sum(group == group1)
  ratio2 <- sum(exp[group == group2] > 0) / sum(group == group2)
  
  # 百分比
  pct1 <- round(ratio1 * 100, 2)
  pct2 <- round(ratio2 * 100, 2)
  
  # 差异
  diff_pct <- round(pct2 - pct1, 2)
  
  data.frame(
    gene = g,
    cell_type = "Extravillous Trophoblast",
    group1 = group1,
    group2 = group2,
    pct_group1 = pct1,
    pct_group2 = pct2,
    diff_pct = diff_pct  # 正数=group2更高，负数=group1更高
  )
}) %>% bind_rows()

final_stats <- merge(
  diff_results,          # 差异检验表
  expression_ratio_stats,# 阳性细胞比例表
  by = c("gene", "group1", "group2"),  # 三列共同匹配
  all = TRUE             # 保留所有基因
)

out_dir <- 'EVT Gene list DEG'
dir.create(out_dir)
write.csv(final_stats, file.path(out_dir, "Gene_list_deg.csv"), row.names = F)


expression_ratio_stats_sample <- lapply(some_genes, function(g) {
  
  # 基因不在矩阵中则跳过
  if (!g %in% rownames(expr_sub)) return(NULL)
  
  # 提取该基因的表达值
  exp <- as.numeric(expr_sub[g, ])
  sample <- meta_sub$sample  # 按样本分组
  all_samples <- unique(sample)
  
  # 对每个样本分别计算阳性比例
  sample_stats_list <- lapply(all_samples, function(s) {
    n_total <- sum(sample == s)
    n_positive <- sum(exp[sample == s] > 0)
    pct <- round(n_positive / n_total * 100, 2)
    
    data.frame(
      gene = g,
      cell_type = "Extravillous Trophoblast",
      sample = s,          # 样本名
      total_cells = n_total,   # 该样本总细胞数
      positive_cells = n_positive, # 阳性细胞数
      pct_positive = pct   # 阳性比例（%）
    )
  })
  
  # 合并当前基因的所有样本结果
  do.call(rbind, sample_stats_list)
}) %>% bind_rows()
expression_ratio_stats_sample <- expression_ratio_stats_sample %>%
  arrange(gene, sample)
######等高线
target_genes <- c("MMP12", "LOX", "TNFSF10","PRG2")
rownames(expr_matrix)
annotation_col <- 'annotation'
HRAS_position_ls <- map(cells_position_ls, function(sample_df) {
  
  # 3.1 从 Seurat 对象中提取当前样本细胞的表达量
  # 注意：这里假设 sample_df$cell_id 中的 ID 与 all_data 的 colnames 完全匹配
  # FetchData 会自动处理 cell_id 的匹配
  
  # 提取表达矩阵子集
  gene_expr_matrix <- FetchData(all_data, 
                                vars = c(annotation_col,target_genes), 
                                cells = sample_df$cell_id)
  
  # 3.2 将提取的矩阵转换为 data.frame 并加上 cell_id 列以便合并
  # rownames(gene_expr_matrix) 就是 cell_id
  expr_df <- gene_expr_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "cell_id")
  
  # 3.3 将表达量数据 (expr_df) 左连接合并到原始坐标数据 (sample_df) 中
  # 使用 left_join 保证原始坐标数据的顺序和行数不变
  final_df <- left_join(sample_df, expr_df, by = "cell_id")
  
  return(final_df)
})

p_hras_contour_ls <- 
  HRAS_position_ls %>% 
  imap(function(coords, sam){
    
    # --- 1. 准备当前样本的绘图数据 ---
    # coords 就是当前样本的 HRAS_position_ls[[i]] (即之前的 dat)
    all_cells <- coords %>% remove_missing()
    hras_positive_cells <- all_cells %>% filter(TNFSF10 > 0)
    
    # 如果该样本没有阳性细胞，直接返回空，防止后续报错
    if (nrow(hras_positive_cells) == 0) return(NULL)
    
    # --- 2. 计算当前样本的二维核密度估计 (KDE) ---
    kde_result <- MASS::kde2d(
      x = hras_positive_cells$SPATIAL_1, 
      y = -hras_positive_cells$SPATIAL_2, 
      n = 100            
    )
    
    # 转换为长格式数据框
    kde_df <- data.frame(
      x = rep(kde_result$x, times = length(kde_result$y)),
      y = rep(kde_result$y, each = length(kde_result$x)),
      z = as.vector(kde_result$z)
    )
    
    # --- 3. 绘制当前样本的叠加图 ---
    p <- ggplot() +
      # 【第一层：背景】细胞类型注释底图
      geom_tile(
        data = all_cells,
        aes(x = SPATIAL_1, y = -SPATIAL_2, fill = annotation), # Y轴同步翻转
        height = 8,
        width = 8
      ) +
      scale_fill_manual(values = cell_colors) +
      
      # 【第二层：前景】HRAS 阳性细胞的密度等高线
      geom_contour(data = kde_df, 
                   aes(x = x, y = y, z = z, color = after_stat(level)), 
                   size = 0.8) +
      
      # 设置配色、主题与标签
      scale_color_viridis_c(option = "C") +
      theme_void() +
      coord_equal() +
      labs(fill = "Celltype", color = "Density") +
      theme(
        panel.background = element_rect(fill = "black", color = NA), # 保持黑色背景
        legend.key.size = unit(0.5, 'cm'),       
        legend.text = element_text(size = 10),   
        legend.title = element_text(size = 12)   
      ) +
      ggtitle(sam)+ # 加上样本名称作为标题
      guides(
        fill = guide_legend(order = 1),       # 离散图例用 guide_legend
        color = guide_colorbar(order = 2)     # 连续图例用 guide_colorbar
      )
    
    return(p)
  })
pdf(
  file.path(out_dir, "contour_plot_TNFSF10 in each sample.pdf"),
  height = 6,
  width = 10
)
p_hras_contour_ls %>% walk(print)
dev.off()
out_dir <- "EVT Gene list DEG"

####