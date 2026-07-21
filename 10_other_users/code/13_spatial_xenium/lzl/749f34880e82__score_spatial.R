library(Seurat)
library(tidyverse)
library(future)
# set default ---------
out_dir <- 'Score_spatial'
dir.create(out_dir)

all_data <- readRDS('/mnt/NexData_file/home/lzl/case/XE2602120662-01/output/data/adata_filtered_add_celltype.rds')

# change assay
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")
all_data$annotation <- gsub("Netrophil", "Neutrophil", all_data$annotation)

# normalize
plan(sequential) 
all_data <- all_data %>% NormalizeData()

Score1 <- c(
  "EPCAM", "TACSTD2", "CLDN4", "MUC1",
  "KRT7", "KRT8", "KRT18", "KRT19", "KRT20",
  "UPK1A", "UPK1B", "UPK2", "UPK3A",
  "GATA3", "FOXA1", "PPARG"
)
Score1 <- stringr::str_to_title(tolower(Score1))
Score1 <- Score1[Score1%in%rownames(all_data)]
Score2 <- c(
  "EPCAM", "TACSTD2", "CLDN3", "CLDN4", "CLDN7",
  "KRT7", "KRT8", "KRT18", "KRT19",
  "MUC1", "AGR2", "S100A14", "SFN",
  "GATA3", "FOXA1"
)
Score2 <- stringr::str_to_title(tolower(Score2))
Score2 <- Score2[Score2%in%rownames(all_data)]

#Urothelial / epithelial score
Score3 <- c('Epcam','Krt8','Krt18','Krt19','Krt5','Krt14','Krt20','Upk1a','Upk1b','Upk2','Upk3a','Gata3','Foxa1','Pparg')
Score3 <- Score3[Score3%in%rownames(all_data)]
#MB49-like / mesenchymal tumor-like score
Score4 <- c('Vim','Fn1','Cdh2','Zeb1','Zeb2','Snai1','Snai2','Twist1','Mmp2','Itga5','Itgb1','Tnc','Serpina3n','Lgals3','Axl','Postn','S100a4')
Score4 <- Score4[Score4%in%rownames(all_data)]
#Proliferation score
Score5 <- c('Mki67','Top2a','Pcna','Mcm2','Mcm5','Mcm6','Cdk1','Ccnb1','Ccna2')
Score5 <- Score5[Score5%in%rownames(all_data)]
#focus_gene
focus_gene <- c("EPCAM","KRT8","KRT18","KRT19",
                "KRT5","KRT14",
                "COL1A1","DCN","LUM",
                "VIM","FN1",
                "PTPRC","CD3D","NKG7",
                "PECAM1","VWF")
focus_gene <- stringr::str_to_title(tolower(focus_gene))
focus_gene <- focus_gene[focus_gene%in%rownames(all_data)]
####
all_data <- JoinLayers(all_data)
all_data <- AddModuleScore(
  object = all_data,
  features = list(Score1, Score2,Score3,Score4,Score5),
  name = c("Score1", "Score2","Score3","Score4","Score5"),  # 列名会变成 Score11, Score21（Seurat自动加数字）
  ctrl = 100  # 控制基因集大小，默认100
)
colnames(all_data@meta.data)[36:40] <- c("Score1","Score2","epithelial score","mesenchymal tumor-like score","Proliferation score")

####
cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()
####
plot_dat <-
  FetchData(
    all_data, 
    vars = c("Score1","Score2","epithelial score","mesenchymal tumor-like score","Proliferation score")  # 两个得分列名
  ) %>%
  rownames_to_column('cell_id')

# 定义颜色方案
score_colors <- list(
  `epithelial score` = c("grey90", "gold", "orange", "red", "darkred"),
  `mesenchymal tumor-like score` = c("grey90",  "skyblue", "dodgerblue", "blue", "navy"),
  `Proliferation score` = c("grey90", "#FDE725", "#73D055", "#238A8D", "#440154")
)

# 绘制每个样本的两个得分分布
# 遍历每个样本，每个样本生成一个PDF（内含两页：两个得分）
iwalk(cells_position_ls, function(coords, sam){
  
  # 准备该样本的数据
  dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
  
  # 打开PDF，文件名包含样本名
  pdf(
    file.path(out_dir, paste0("1. spatial plot of new scores - ", sam, ".pdf")),
    height = 8,
    width = 10
  )
  
  # 遍历两个得分，每页一个
  for (score_name in names(score_colors)) {
    
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, y = SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[score_name]]),
        height = 8, 
        width = 8
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      scale_fill_gradientn(
        colors = score_colors[[score_name]],
        na.value = "grey90",
        name = score_name
      ) +
      theme(
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right"
      ) +
      ggtitle(paste0(sam, " - ", score_name))
    
    print(p)
  }
  
  dev.off()
})


# 设置要绘制的基因
focus_gene <- c("Epcam", "Col1a1", "Lum", "Fn1", "Ptprc", "Cd3d", "Pecam1", "Vwf")
#target_gene <- "Epcam"  # 替换为你的目标基因名
spatial_colors <- c("#00008B", "#0000CD", "#4169E1", "#00BFFF", 
                    "#00FF7F", "#ADFF2F", "#FFFF00", "#FFD700", "#FF4500")

# 遍历每个基因，每个基因生成一个PDF
for (target_gene in focus_gene) {
  
  # 提取基因表达数据
  plot_dat <-
    FetchData(
      all_data, 
      vars = c(target_gene)
    ) %>%
    rownames_to_column('cell_id')
  
  # 绘制空间表达分布
  p_sp_gene_ls <-
    cells_position_ls %>%
    imap(function(coords, sam){
      dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
      
      p <-
        dat2 %>%
        ggplot(aes(x = SPATIAL_1, y = SPATIAL_2)) +
        geom_tile(
          aes(fill = .data[[target_gene]]),
          height = 8, 
          width = 8
        ) +
        theme_void() +
        coord_equal() +
        scale_y_reverse() +
        # SpatialFeaturePlot 经典渐变
        scale_fill_gradientn(
          colors = spatial_colors,
          na.value = "grey90",
          name = target_gene
        ) +
        theme(
          panel.background = element_rect(fill = NA, color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right"
        ) +
        ggtitle(paste0(sam, " - ", target_gene))
      
      p
    })
  
  # 输出PDF
  pdf(
    file.path(out_dir, paste0("spatial plot of ", target_gene, " expression.pdf")),
    height = 8,
    width = 10
  )
  p_sp_gene_ls %>% walk(print)
  dev.off()
  
  message("完成: ", target_gene)
}
###celltype
out_dir <- "Celltype_spatial"
dir.create(out_dir)
annotation_col <- 'annotation'
all_data$annotation[all_data$annotation == "Mesothelial"] <- "Tumor"
table(all_data$annotation)
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
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')

unique(all_data$annotation)
table(all_data$annotation)
cell_colors <- c(
  "Endothelial" = "#1f77b4",
  "Epithelial" = "#ff7f0e",
  "Fibroblast" = "#279e68",
  "Fibroblastic reticular cell(FRC)" = "#d62728",
  "Lipid-Processing Fibroblast" = "#aa40fc",
  "Macrophage" = "#8c564b",
  "Tumor" = "#e377c2",
  "Monocyte" = "#b5bd61",
  "Neutrophil" = "#17becf",
  "Other" = "#aec7e8",
  "Smooth muscle cell" = "#ffbb78",
  "T cell" = "#98df8a"
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
      
      dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
      
      # 【核心新增】：如果当前样本中没有目标细胞类型，直接跳过该样本
      if (!ct %in% dat2[[annotation_col]]) {
        message("⏭️ 样本 [", sam, "] 不含 ", ct, "，已跳过...")
        return(NULL) # 在 map 函数中，返回 NULL 代表跳过此元素
      }
      
      # 高亮分组（此时可以安全地使用 ifelse，因为目标细胞必然存在）
      dat2$highlight <- ifelse(dat2[[annotation_col]] == ct, ct, "Others")
      
      dat2$highlight <- factor(
        dat2$highlight,
        levels = c(ct, "Others")
      )
      
      # 绘图
      p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
        geom_tile(aes(fill = highlight,alpha = highlight), height = 8, width = 8) +
        theme_void() +
        coord_equal() +
        scale_y_reverse() +
        
        ##透明度
        scale_alpha_manual(
          values = setNames(c(1, 0.3), c(ct, "Others")) 
        ) +
        
        # 核心：强制命名颜色（解决ggplot不认变量问题）
        scale_fill_manual(
          values = setNames(c(current_color, "gray60"), c(ct, "Others")),
          guide = guide_legend(override.aes = list(alpha = c(1, 1)))
        ) +
        
        theme(
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle(sam)
      
      return(p)
    }) %>% 
    compact() # 【核心新增】：自动剔除列表中的 NULL 元素，防止后续 walk 报错
  
  # 输出PDF
  pdf(
    file.path(out_dir, paste0("spatial plot of ", target_cell_type, ".pdf")),
    height = 8, width = 8
  )
  walk(p_sp_ct_ls, print)
  dev.off()
}

message("全部完成！已输出每个细胞类型独立PDF！")

######
immnue_celltype <- c("T cell", "Neutrophil", "Macrophage", "Monocyte")

plot_dat <-
  FetchData(
    all_data, 
    vars = c("annotation")  # 获取原始注释，不是二分的 immune
  ) %>%
  rownames_to_column('cell_id')

p_sp_immnue_ls <- cells_position_ls %>%
  imap(function(coords, sam){
    
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    # 创建颜色映射：免疫细胞用原色，其他用 gray60
    dat2$color_group <- ifelse(
      dat2$annotation %in% immnue_celltype,
      dat2$annotation,           # 保留原始细胞类型名
      "Others"                   # 其他统一
    )
    
    # 构建颜色向量
    plot_colors <- cell_colors[immnue_celltype]  # 提取免疫细胞颜色
    plot_colors["Others"] <- "gray60"             # 其他颜色
    
    # 只保留 dat2 中实际存在的类型
    present_types <- unique(dat2$color_group)
    plot_colors <- plot_colors[names(plot_colors) %in% present_types]
    
    # 透明度：免疫细胞不透明，其他半透明
    alpha_values <- setNames(
      ifelse(names(plot_colors) == "Others", 0.3, 1),
      names(plot_colors)
    )
    
    # 确保因子顺序：免疫细胞在前，Others 最后
    type_levels <- c(setdiff(present_types, "Others"), "Others")
    dat2$color_group <- factor(dat2$color_group, levels = type_levels)
    
    p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
      geom_tile(aes(fill = color_group, alpha = color_group), height = 8, width = 8) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      
      scale_alpha_manual(values = alpha_values) +
      
      scale_fill_manual(
        values = plot_colors,
        guide = guide_legend(override.aes = list(alpha = 1))  # 图例统一不透明
      ) +
      
      theme(
        panel.background = element_rect(fill = NA, color = NA),
        plot.title = element_text(hjust = 0.5, size = 16)
      ) +
      ggtitle(sam)
    
    return(p)
  })

pdf(
  file.path(out_dir, "spatial plot of immune.pdf"),
  height = 8, width = 8
)
walk(p_sp_immnue_ls, print)
dev.off()

dir.create("rds")
saveRDS(all_data,"./rds/alldata_update_cell.rds")



#######
df = all_data@meta.data
#df = read.csv("/lvdata/wzb/scRNA/FWSC20240734/ReAll/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('annotation','group')])

# 计算所有样本 × 所有亚型的频数
summary_t <- as.data.frame(summary)

# 按分组，计算每个样本内的比例
summary_prop <- summary_t %>%
  group_by(group) %>%
  mutate(total = sum(Freq), prop = Freq/total) %>%
  ungroup()

pdf(
  file.path(out_dir, "celltype_prop - all group.pdf"),
  width = 10,   # 加宽以容纳多个样本
  height = 8
)

p <- ggplot(summary_prop, aes(x = group, y = prop, fill = annotation)) +
  geom_col(position = "stack", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors) +
  labs(
    x = "",
    y = "Cell proportion",
    fill = "Cell type",
    title = "Celltype proportion across group"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),  # x轴标签倾斜防重叠
    axis.text.y = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(p)
dev.off()


df = all_data@meta.data
#df = read.csv("/lvdata/wzb/scRNA/FWSC20240734/ReAll/celltype/meta.xls",sep = "\t",row.names=1)
summary <- table(df[,c('annotation','project')])

# 计算所有样本 × 所有亚型的频数
summary_t <- as.data.frame(summary)

# 按分组，计算每个样本内的比例
summary_prop <- summary_t %>%
  group_by(project) %>%
  mutate(total = sum(Freq), prop = Freq/total) %>%
  ungroup()

pdf(
  file.path(out_dir, "celltype_prop - all sample.pdf"),
  width = 10,   # 加宽以容纳多个样本
  height = 8
)

p <- ggplot(summary_prop, aes(x = project, y = prop, fill = annotation)) +
  geom_col(position = "stack", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors) +
  labs(
    x = "",
    y = "Cell proportion",
    fill = "Cell type",
    title = "Celltype proportion across sample"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),  # x轴标签倾斜防重叠
    axis.text.y = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(p)
dev.off()
