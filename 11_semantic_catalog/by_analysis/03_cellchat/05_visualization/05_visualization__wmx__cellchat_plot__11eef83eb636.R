library(CellChat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

## ALl 大群 ================================
## 1. 参数

rdata_file <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/All/tmp_object.RData"   
outdir <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/All/Microglia_subtype_compare"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


target_cell <- "Oligodendrocyte"

source_cells <- c("Microglia_Chka"  , "Microglia_Cmss1" , "Microglia_Ccl4"  , "Microglia_Csmd3" , "Microglia_Hspa1b", "Microglia_Tle4")

## 如果你的细胞名不是这个，先运行：
## levels(cellchat@idents)
## 然后把上面三个名字改成完全一致


## 2. 读取 CellChat 对象

load(rdata_file)

## 自动寻找 CellChat 对象
obj_names <- ls()
cellchat <- object.list$all

message("读取到 CellChat 对象：", 'All')

cell_levels <- levels(cellchat@idents)

message("当前 CellChat 细胞类型如下：")
print(cell_levels)

if (!all(c(source_cells, target_cell) %in% cell_levels)) {
  stop(
    "source 或 target 名称不在 cellchat@idents 中，请检查：\n",
    "source = ", source_cells, "\n",
    "target_cell = ", target_cell, "\n\n",
    "当前可用细胞类型：\n",
    paste(cell_levels, collapse = "\n")
  )
}


## 3. 提取 microglia → oligodendrocyte 通讯结果

df.net <- subsetCommunication(
  cellchat,
  sources.use = source_cells,
  targets.use = target_cell
)

if (nrow(df.net) == 0) {
  stop("没有提取到指定 source → target 的通讯结果，请检查细胞名或 CellChat 是否已完成 computeCommunProb/filterCommunication")
}

df.net <- df.net %>%
  mutate(
    source_group = case_when(
      TRUE ~ as.character(source)
    ),
    LR_pair = paste0(ligand, " - ", receptor)
  )


## 4. 总通讯数量 / 总通讯强度比较
df.summary <- df.net %>%
  group_by(source_group) %>%
  summarise(
    interaction_number = n(),
    interaction_strength = sum(prob, na.rm = TRUE),
    mean_strength = mean(prob, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  df.summary,
  file = file.path(outdir, "0.microglia_subtype_to_oligodendrocyte_summary.csv"),
  row.names = FALSE
)

pC1 <- ggplot(df.summary, aes(x = source_group, y = interaction_number, fill = source_group)) +
  geom_col(width = 0.65) +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Number of interactions",
    title = "Microglia → Oligodendrocyte interaction number"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

pC2 <- ggplot(df.summary, aes(x = source_group, y = interaction_strength, fill = source_group)) +
  geom_col(width = 0.65) +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Total interaction strength",
    title = "Microglia → Oligodendrocyte interaction strength"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(outdir, "0.total_interaction_number_strength.pdf"), 
       plot = pC1 + pC2, 
       device = cairo_pdf, 
       width = 10, height = 5)

ggsave(
  filename = file.path(outdir, "0.total_interaction_number_strength.png"),
  plot = pC1 + pC2,
  width = 10,
  height = 5,
  dpi = 300
)


## 5. 差异 signaling pathway dotplot
df.pathway <- df.net %>%
  group_by(source_group, pathway_name) %>%
  summarise(
    pathway_strength = sum(prob, na.rm = TRUE),
    interaction_number = n(),
    .groups = "drop"
  )

df.pathway.wide <- df.pathway %>%
  select(source_group, pathway_name, pathway_strength) %>%
  pivot_wider(
    names_from = source_group,
    values_from = pathway_strength,
    values_fill = 0
  )

source_cols <- setdiff(
  colnames(df.pathway.wide),
  "pathway_name"
)
df.pathway.wide <- df.pathway.wide %>%
  rowwise() %>%
  mutate(
    max_strength = max(c_across(all_of(source_cols))),
    min_strength = min(c_across(all_of(source_cols))),
    strength_difference = max_strength - min_strength,
    mean_strength = mean(
      c_across(all_of(source_cols))
    ),
    # 最大来源
    max_source =source_cols[which.max(c_across(all_of(source_cols)))],
    # 最小来源
    min_source =source_cols[which.min(c_across(all_of(source_cols)))],
    
  ) %>%
  ungroup()

write.csv(
  df.pathway.wide,
  file = file.path(outdir, "01.pathway_strength_subtype_compare.csv"),
  row.names = FALSE
)

top_pathways <- df.pathway.wide %>%
  filter(strength_difference > 0) %>%
  slice_head(n = 30) %>%
  pull(pathway_name)

df.pathway.plot <- df.pathway %>%
  filter(pathway_name %in% top_pathways) %>%
  left_join(
    df.pathway.wide %>% select(pathway_name, strength_difference),
    by = "pathway_name"
  ) %>%
  mutate(
    pathway_name = factor(pathway_name, levels = rev(top_pathways)),
    source_group = factor(source_group, levels = source_cells)
  )

library(RColorBrewer)

pD <-ggplot(
  df.pathway.plot,
  aes(
    x = source_group,
    y = pathway_name,
    size = pathway_strength,
    color = pathway_strength
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(
    colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
    name = "Pathway strength"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "Pathway strength"
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = NULL,
    title = "Pathway-level signaling: microglia → oligodendrocyte"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(outdir, "01.pathway_dotplot_microglia_to_oligodendrocyte.pdf"), 
       plot = pD, 
       device = cairo_pdf, 
       width = 7, height = 8)

ggsave(
  filename = file.path(outdir, "01.pathway_dotplot_microglia_to_oligodendrocyte.png"),
  plot = pD,
  width = 7,
  height = 8,
  dpi = 300
)


## 6. ligand-receptor pair bubble plot

df.lr.wide <- df.net %>%
  group_by(
    source_group,
    LR_pair,
    pathway_name,
    ligand,
    receptor
  ) %>%
  summarise(
    LR_strength = sum(prob, na.rm = TRUE),
    pval_min = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = source_group,
    values_from = LR_strength,
    values_fill = 0
  )
source_cols <- setdiff(
  colnames(df.lr.wide),
  c(
    "LR_pair",
    "pathway_name",
    "ligand",
    "receptor",
    "pval_min"
  )
)

df.lr.wide <- df.lr.wide %>%
  rowwise() %>%
  mutate(
    # 最大/最小强度
    max_strength = max(c_across(all_of(source_cols)), na.rm = TRUE),
    min_strength = min(c_across(all_of(source_cols)),na.rm = TRUE),
    # 差异范围
    strength_difference =max_strength - min_strength,
    # 平均强度
    mean_strength =mean(c_across(all_of(source_cols)), na.rm = TRUE),
    # 最大来源
    max_source =source_cols[which.max(c_across(all_of(source_cols)))],
    # 最小来源
    min_source =source_cols[which.min(c_across(all_of(source_cols)))],
    # 平均strength
    mean_strength = mean(
      c_across(all_of(source_cols))
    ),
  ) %>%
  ungroup() %>%
  arrange(
    desc(abs(strength_difference))
  )

write.csv(
  df.lr.wide,
  file = file.path(outdir, "02.LR_pair_microglia_subtype_compare.csv"),
  row.names = FALSE
)

top_lr <- df.lr.wide %>%
  slice_head(n = 40) %>%
  pull(LR_pair)

df.lr.plot <- df.net %>%
  group_by(source_group, LR_pair, pathway_name) %>%
  summarise(
    LR_strength = sum(prob, na.rm = TRUE),
    pval_min = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(LR_pair %in% top_lr) %>%
  mutate(
    LR_pair = factor(LR_pair, levels = rev(top_lr)),
    source_group = factor(source_group, levels = source_cells),
    neglog10_pval = -log10(pval_min + 1e-300)
  )

pE <- ggplot(
  df.lr.plot,
  aes(
    x = source_group,
    y = LR_pair,
    size = LR_strength ,
    color = LR_strength
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(
    colors =c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
    name = "LR_strength"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "LR_strength"
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = NULL,
    title = "Ligand-receptor pairs: microglia → oligodendrocyte"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(outdir, "02.LR_pair_bubble_microglia_to_oligodendrocyte.pdf"), 
      plot = pE, 
      device = cairo_pdf, 
      width = 7, height = 8)

ggsave(
  filename = file.path(outdir, "02.LR_pair_bubble_microglia_to_oligodendrocyte.png"),
  plot = pE,
  width = 8,
  height = 10,
  dpi = 300
)


## 7. 重点 pathway 的 circle plot

focus_pathways <- df.pathway.wide %>% pull(pathway_name)

## 只保留三个节点：high microglia、low microglia、oligodendrocyte
vertex.receiver <- match(target_cell, levels(cellchat@idents))

for (pw in focus_pathways) {
  
  message("绘制 pathway circle plot: ", pw)
  
  pdf(
    file = file.path(outdir, paste0("03.circle_", pw, "_microglia_to_oligodendrocyte.pdf")),
    width = 6,
    height = 6
  )
  
  tryCatch({
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      sources.use = source_cells,
      targets.use = target_cell,
      layout = "circle",
      vertex.receiver = vertex.receiver,
      vertex.label.cex = 0.8,
      edge.width.max = 8,
      signaling.name = paste0(pw, ": microglia → oligodendrocyte")
    )
  }, error = function(e) {
    message("pathway 绘制失败: ", pw)
    message(e$message)
  })
  
  dev.off()
  
  png(
    file = file.path(outdir, paste0("03.circle_", pw, "_microglia_to_oligodendrocyte.png")),
    width = 600,
    height = 600
  )
  
  tryCatch({
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      sources.use = source_cells,
      targets.use = target_cell,
      layout = "circle",
      vertex.receiver = vertex.receiver,
      vertex.label.cex = 0.8,
      edge.width.max = 8,
      signaling.name = paste0(pw, ": microglia → oligodendrocyte")
    )
  }, error = function(e) {
    message("pathway 绘制失败: ", pw)
    message(e$message)
  })
  
  dev.off()
}

## All 亚群 -----------------------------------------

## 1. 参数

rdata_file <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/All_sub_oli/tmp_object.RData"

outdir <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/All_sub_oli/Microglia_subtype_compare_Oligo_subtypes"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

source_cells <- c("Microglia_Chka"  , "Microglia_Cmss1" , "Microglia_Ccl4"  , "Microglia_Csmd3" , "Microglia_Hspa1b", "Microglia_Tle4")

## 2. 读取 CellChat 对象

load(rdata_file)

cellchat <- object.list$all

cell_levels <- levels(cellchat@idents)

message("当前 CellChat 细胞类型如下：")
print(cell_levels)

## 3. 设置 Oligodendrocyte 亚群

## 方法 A：手动指定，推荐最稳
target_cells <- c(
  "Oligodendrocyte_Kirrel3",
  "Oligodendrocyte_Grm7",
  "Oligodendrocyte_Cdkn1a",
  "Oligodendrocyte_Selenop",
  "Oligodendrocyte_Il33",
  "Oligodendrocyte_S100b",
  "Oligodendrocyte_Hebp1"
)


message("本次分析 target cells 为：")
print(target_cells)

if (!all(source_cells %in% cell_levels)) {
  stop(
    "source 名称不在 cellchat@idents 中，请检查：\n",
    paste(source_cells, collapse = ", "), "\n\n",
    "当前可用细胞类型：\n",
    paste(cell_levels, collapse = "\n")
  )
}

if (!all(target_cells %in% cell_levels)) {
  stop(
    "部分 target 名称不在 cellchat@idents 中，请检查：\n",
    paste(setdiff(target_cells, cell_levels), collapse = ", "), "\n\n",
    "当前可用细胞类型：\n",
    paste(cell_levels, collapse = "\n")
  )
}

## 文件名安全处理函数，避免 T/NK 这种名字导致路径错误
safe_name <- function(x) {
  x %>%
    gsub("[/\\\\:*?\"<>| ]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}

## 4. 定义单个 target 的分析函数

run_one_target <- function(cellchat, target_cell, outdir_target) {
  
  dir.create(outdir_target, showWarnings = FALSE, recursive = TRUE)
  
  message("====================================")
  message("开始分析 target: ", target_cell)
  message("输出目录: ", outdir_target)
  
  ## ----------------
  ## 提取通讯结果
  ## ----------------
  
  df.net <- subsetCommunication(
    cellchat,
    sources.use = source_cells,
    targets.use = target_cell
  )
  
  if (nrow(df.net) == 0) {
    message("跳过：没有提取到 ", paste(source_cells, collapse = "/"), " → ", target_cell, " 的通讯结果")
    return(NULL)
  }
  
  df.net <- df.net %>%
    mutate(
      source_group = case_when(
        TRUE ~ as.character(source)
      ),
      target_group = as.character(target),
      LR_pair = paste0(ligand, " - ", receptor)
    )
  
  
  ## ================================
  ## C. 总通讯数量 / 总通讯强度
  ## ================================
  
  df.summary <- df.net %>%
    group_by(source_group) %>%
    summarise(
      interaction_number = n(),
      interaction_strength = sum(prob, na.rm = TRUE),
      mean_strength = mean(prob, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(
    df.summary,
    file = file.path(outdir_target, paste0("0.", safe_name(target_cell), "_summary.csv")),
    row.names = FALSE
  )
  
  pC1 <- ggplot(df.summary, aes(x = source_group, y = interaction_number, fill = source_group)) +
    geom_col(width = 0.65) +
    theme_classic(base_size = 10) +
    labs(
      x = NULL,
      y = "Number of interactions",
      title = paste0("Microglia → ", target_cell, "\ninteraction number")
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  
  pC2 <- ggplot(df.summary, aes(x = source_group, y = interaction_strength, fill = source_group)) +
    geom_col(width = 0.65) +
    theme_classic(base_size = 10) +
    labs(
      x = NULL,
      y = "Total interaction strength",
      title = paste0("Microglia → ", target_cell, "\ninteraction strength")
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  
  ggsave(
    file.path(outdir_target, paste0("0.", safe_name(target_cell), "_total_interaction_number_strength.pdf")),
    plot = pC1 + pC2,
    device = cairo_pdf,
    width = 10,
    height = 5
  )
  
  ggsave(
    file.path(outdir_target, paste0("0.", safe_name(target_cell), "_total_interaction_number_strength.png")),
    plot = pC1 + pC2,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  ## ================================
  ## D. pathway-level dotplot
  ## ================================
  
  df.pathway <- df.net %>%
    group_by(source_group, pathway_name) %>%
    summarise(
      pathway_strength = sum(prob, na.rm = TRUE),
      interaction_number = n(),
      .groups = "drop"
    )
  
  df.pathway.wide <- df.pathway %>%
    select(source_group, pathway_name, pathway_strength) %>%
    pivot_wider(
      names_from = source_group,
      values_from = pathway_strength,
      values_fill = 0
    )
  
  source_cols <- setdiff(
    colnames(df.pathway.wide),
    "pathway_name"
  )
  df.pathway.wide <- df.pathway.wide %>%
    rowwise() %>%
    mutate(
      max_strength = max(c_across(all_of(source_cols))),
      min_strength = min(c_across(all_of(source_cols))),
      strength_difference = max_strength - min_strength,
      mean_strength = mean(
        c_across(all_of(source_cols))
      ),
      # 最大来源
      max_source =source_cols[which.max(c_across(all_of(source_cols)))],
      # 最小来源
      min_source =source_cols[which.min(c_across(all_of(source_cols)))],
      
    ) %>%
    ungroup()
  
  write.csv(
    df.pathway.wide,
    file = file.path(outdir_target, paste0("01.", safe_name(target_cell), "_pathway_strength_microglia_subtype_compare.csv")),
    row.names = FALSE
  )
  
  top_pathways <- df.pathway.wide %>%
    slice_head(n = 30) %>%
    pull(pathway_name)
  
  df.pathway.plot <- df.pathway %>%
    filter(pathway_name %in% top_pathways) %>%
    left_join(
      df.pathway.wide %>%  select(pathway_name, strength_difference),
      by = "pathway_name"
    ) %>%
    mutate(
      pathway_name = factor(pathway_name, levels = rev(top_pathways)),
      source_group = factor(source_group, levels = source_cells)
    )
  
  pD <- ggplot(
    df.pathway.plot,
    aes(
      x = source_group,
      y = pathway_name,
      size = pathway_strength,
      color = pathway_strength
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
      name = "Pathway strength"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "Pathway strength"
    ) +
    theme_classic(base_size = 13) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Pathway signaling: microglia → ", target_cell)
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(outdir_target, paste0("01.", safe_name(target_cell), "_pathway_dotplot.pdf")),
    plot = pD,
    device = cairo_pdf,
    width = 8,
    height = 8
  )
  
  ggsave(
    file.path(outdir_target, paste0("01.", safe_name(target_cell), "_pathway_dotplot.png")),
    plot = pD,
    width = 8,
    height = 8,
    dpi = 300
  )
  
  ## ================================
  ## E. LR pair bubble plot
  ## ================================
  
  df.lr.wide <- df.net %>%
    group_by(
      source_group,
      LR_pair,
      pathway_name,
      ligand,
      receptor
    ) %>%
    summarise(
      LR_strength = sum(prob, na.rm = TRUE),
      pval_min = min(pval, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = source_group,
      values_from = LR_strength,
      values_fill = 0
    )
  source_cols <- setdiff(
    colnames(df.lr.wide),
    c(
      "LR_pair",
      "pathway_name",
      "ligand",
      "receptor",
      "pval_min"
    )
  )
  
  df.lr.wide <- df.lr.wide %>%
    rowwise() %>%
    mutate(
      # 最大/最小强度
      max_strength = max(c_across(all_of(source_cols)), na.rm = TRUE),
      min_strength = min(c_across(all_of(source_cols)),na.rm = TRUE),
      # 差异范围
      strength_difference =max_strength - min_strength,
      # 平均强度
      mean_strength =mean(c_across(all_of(source_cols)), na.rm = TRUE),
      # 最大来源
      max_source =source_cols[which.max(c_across(all_of(source_cols)))],
      # 最小来源
      min_source =source_cols[which.min(c_across(all_of(source_cols)))],
      # 平均strength
      mean_strength = mean(
        c_across(all_of(source_cols))
      ),
    ) %>%
    ungroup() %>%
    arrange(
      desc(abs(strength_difference))
    )
  
  write.csv(
    df.lr.wide,
    file = file.path(outdir_target, paste0("02.", safe_name(target_cell), "_LR_pair_microglia_subtype_compare.csv")),
    row.names = FALSE
  )
  
  top_lr <- df.lr.wide %>%
    slice_head(n = 40) %>%
    pull(LR_pair)
  
  df.lr.plot <- df.net %>%
    group_by(source_group, LR_pair, pathway_name) %>%
    summarise(
      LR_strength = sum(prob, na.rm = TRUE),
      pval_min = min(pval, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(LR_pair %in% top_lr) %>%
    mutate(
      LR_pair = factor(LR_pair, levels = rev(top_lr)),
      source_group = factor(source_group, levels = source_cells),
      neglog10_pval = -log10(pval_min + 1e-300)
    )
  
  pE <- ggplot(
    df.lr.plot,
    aes(
      x = source_group,
      y = LR_pair,
      size = LR_strength,
      color = LR_strength
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
      name = "LR strength"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "LR strength"
    ) +
    theme_classic(base_size = 13) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Ligand-receptor pairs: microglia → ", target_cell)
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(outdir_target, paste0("02.", safe_name(target_cell), "_LR_pair_bubble.pdf")),
    plot = pE,
    device = cairo_pdf,
    width = 8,
    height = 10
  )
  
  ggsave(
    file.path(outdir_target, paste0("02.", safe_name(target_cell), "_LR_pair_bubble.png")),
    plot = pE,
    width = 8,
    height = 10,
    dpi = 300
  )
  
  ## ================================
  ## F. 重点 pathway circle plot
  ## ================================
  
  focus_pathways <- df.pathway.wide %>%
    pull(pathway_name)
  
  vertex.receiver <- match(target_cell, levels(cellchat@idents))
  
  for (pw in focus_pathways) {
    
    message("绘制 pathway circle plot: ", target_cell, " - ", pw)
    
    pdf(
      file = file.path(outdir_target, paste0("03.circle_", safe_name(target_cell), "_", safe_name(pw), ".pdf")),
      width = 6,
      height = 6
    )
    
    tryCatch({
      netVisual_aggregate(
        cellchat,
        signaling = pw,
        sources.use = source_cells,
        targets.use = target_cell,
        layout = "circle",
        vertex.receiver = vertex.receiver,
        vertex.label.cex = 0.8,
        edge.width.max = 8,
        signaling.name = paste0(pw, ": microglia → ", target_cell)
      )
    }, error = function(e) {
      message("pathway 绘制失败: ", target_cell, " - ", pw)
      message(e$message)
    })
    
    dev.off()
    
    png(
      file = file.path(outdir_target, paste0("03.circle_", safe_name(target_cell), "_", safe_name(pw), ".png")),
      width = 800,
      height = 800,
      res = 120
    )
    
    tryCatch({
      netVisual_aggregate(
        cellchat,
        signaling = pw,
        sources.use = source_cells,
        targets.use = target_cell,
        layout = "circle",
        vertex.receiver = vertex.receiver,
        vertex.label.cex = 0.8,
        edge.width.max = 8,
        signaling.name = paste0(pw, ": microglia → ", target_cell)
      )
    }, error = function(e) {
      message("pathway 绘制失败: ", target_cell, " - ", pw)
      message(e$message)
    })
    
    dev.off()
  }
  
  message("完成 target: ", target_cell)
  
  return(
    list(
      df.net = df.net,
      df.summary = df.summary,
      df.pathway.wide = df.pathway.wide,
      df.lr.wide = df.lr.wide
    )
  )
}

## 5. 循环所有 Oligodendrocyte 亚群

all_results <- list()

for (tc in target_cells) {
  
  outdir_target <- file.path(outdir, safe_name(tc))
  
  res <- run_one_target(
    cellchat = cellchat,
    target_cell = tc,
    outdir_target = outdir_target
  )
  
  all_results[[tc]] <- res
}


message("全部 Oligodendrocyte 亚群分析完成，结果保存在：", outdir)


## 

# 组间比较筛选 ，大群 ---------------------------
## 1. 参数设置
rdata_g1_g2 <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/Compare/M_15d_M_3d/tmp_object.RData"
rdata_g1_g3 <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/Compare/M_30d_M_3d/tmp_object.RData"

outdir <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/Compare/Microglia_to_Olig_group_compare/"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

source_cells <- c("Microglia_Chka"  , "Microglia_Cmss1" , "Microglia_Ccl4"  , "Microglia_Csmd3" , "Microglia_Hspa1b", "Microglia_Tle4")

## 如果是大类 Olig
target_cells <- c("Oligodendrocyte")

top_pathway_n <- 30
top_lr_n <- 40

safe_name <- function(x) {
  x %>%
    gsub("[/\\\\:*?\"<>| ]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}

## 2. 读取两个 RData 中的 object.list

load_cellchat_list <- function(rdata_file) {
  e <- new.env()
  load(rdata_file, envir = e)
  
  if (!exists("object.list", envir = e)) {
    stop("这个 RData 中没有 object.list: ", rdata_file)
  }
  
  obj.list <- get("object.list", envir = e)
  return(obj.list)
}

object.list.g12 <- load_cellchat_list(rdata_g1_g2)
object.list.g13 <- load_cellchat_list(rdata_g1_g3)

message("RData1 object.list names:")
print(names(object.list.g12))

message("RData2 object.list names:")
print(names(object.list.g13))
cellchat_g1 <- object.list.g12$M_3d
cellchat_g2 <- object.list.g12$M_15d
cellchat_g3 <- object.list.g13$M_30d

## 3. 提取 source → target 通讯结果

extract_comm_one <- function(cellchat, group_name, source_cells, target_cells) {
  
  cell_levels <- levels(cellchat@idents)
  
  if (!all(source_cells %in% cell_levels)) {
    stop(
      group_name, " 中 source_cells 不存在：",
      paste(setdiff(source_cells, cell_levels), collapse = ", "),
      "\n当前可用细胞类型：\n",
      paste(cell_levels, collapse = "\n")
    )
  }
  
  if (!all(target_cells %in% cell_levels)) {
    stop(
      group_name, " 中 target_cells 不存在：",
      paste(setdiff(target_cells, cell_levels), collapse = ", "),
      "\n当前可用细胞类型：\n",
      paste(cell_levels, collapse = "\n")
    )
  }
  
  df <- subsetCommunication(
    cellchat,
    sources.use = source_cells,
    targets.use = target_cells
  )
  
  if (nrow(df) == 0) {
    warning(group_name, " 没有提取到通讯结果")
    return(NULL)
  }
  
  df <- df %>%
    mutate(
      group = group_name,
      source_group = as.character(source),
      target_group = as.character(target),
      LR_pair = paste0(ligand, " - ", receptor)
    )
  
  return(df)
}

df.g1 <- extract_comm_one(cellchat_g1, "M_3d", source_cells, target_cells)
df.g2 <- extract_comm_one(cellchat_g2, "M_15d", source_cells, target_cells)
df.g3 <- extract_comm_one(cellchat_g3, "M_30d", source_cells, target_cells)

df.all <- bind_rows(df.g1, df.g2, df.g3)

write.csv(
  df.all,
  file = file.path(outdir, "00.all_group_source_to_target_communication.csv"),
  row.names = FALSE
)

## 4. 总通讯强度 summary

df.summary.all <- df.all %>%
  group_by(group, source_group, target_group) %>%
  summarise(
    interaction_number = n(),
    interaction_strength = sum(prob, na.rm = TRUE),
    mean_strength = mean(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = c("M_3d", "M_15d", "M_30d")),
    source_group = factor(source_group, levels = source_cells)
  )

write.csv(
  df.summary.all,
  file = file.path(outdir, "00.summary_interaction_number_strength_all_groups.csv"),
  row.names = FALSE
)

p_summary_strength <- ggplot(
  df.summary.all,
  aes(
    x = group,
    y = interaction_strength,
    fill = source_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  # facet_wrap(~ group, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Total interaction strength",
    fill = NULL,
    title = "Microglia → Oligodendrocyte interaction strength"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(outdir, "00.summary_total_interaction_strength.pdf"),
  plot = p_summary_strength,
  device = cairo_pdf,
  width = max(8, length(target_cells) * 3),
  height = 5
)

ggsave(
  file.path(outdir, "00.summary_total_interaction_strength.png"),
  plot = p_summary_strength,
  width = max(8, length(target_cells) * 3),
  height = 5,
  dpi = 300
)

p_summary_number <- ggplot(
  df.summary.all,
  aes(
    x = group,
    y = interaction_number,
    fill = source_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  # facet_wrap(~ target_group, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Number of interactions",
    fill = NULL,
    title = "Microglia → Oligodendrocyte interaction number"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(outdir, "00.summary_interaction_number.pdf"),
  plot = p_summary_number,
  device = cairo_pdf,
  width = max(8, length(target_cells) * 3),
  height = 5
)

ggsave(
  file.path(outdir, "00.summary_interaction_number.png"),
  plot = p_summary_number,
  width = max(8, length(target_cells) * 3),
  height = 5,
  dpi = 300
)

## 6. pathway 层面统计

df.pathway.all <- df.all %>%
  group_by(group, source_group, target_group, pathway_name) %>%
  summarise(
    pathway_strength = sum(prob, na.rm = TRUE),
    interaction_number = n(),
    .groups = "drop"
  )

df.pathway.compare <- df.pathway.all %>%
  select(group, source_group, target_group, pathway_name, pathway_strength) %>%
  pivot_wider(
    names_from = group,
    values_from = pathway_strength,
    values_fill = 0
  )

df.pathway.compare <- df.pathway.compare %>%
  mutate(
    diff_M15d_vs_M3d = M_15d - M_3d,
    diff_M30d_vs_M3d = M_30d - M_3d,
    
    increased_in_M15d = diff_M15d_vs_M3d > 0,
    increased_in_M30d =  diff_M30d_vs_M3d > 0,
    increased_in_M15d_and_M30d = diff_M15d_vs_M3d > 0 &  diff_M30d_vs_M3d > 0,
    
    decreased_in_M15d = diff_M15d_vs_M3d < 0,
    decreased_in_M30d = diff_M30d_vs_M3d < 0,
    decreased_in_M15d_and_M30d = diff_M15d_vs_M3d < 0 & diff_M30d_vs_M3d < 0,
    
  ) %>%
  arrange(source_group, target_group)

write.csv(
  df.pathway.compare,
  file = file.path(outdir, "01.pathway_M15d_M30d_vs_M3d_compare_with_flags.csv"),
  row.names = FALSE
)


## 7. pathway bubble plot

plot_pathway_enriched_type <- function(
    source_use,
    enriched_type = c("M15d", "M30d", "both"),
    df.pathway.compare,
    df.pathway.all,
    outdir,
    top_n = 30
) {
  
  enriched_type <- match.arg(enriched_type)
  
  if (enriched_type == "M15d") {
    df.enrich <- df.pathway.compare %>%
      filter(source_group == source_use, increased_in_M15d) %>%
      arrange(desc(increased_in_M15d))
    
    title_suffix <- "increased M_15d"
    file_suffix <- "increased_M_15d"
    
  } else if (enriched_type == "M30d") {
    df.enrich <- df.pathway.compare %>%
      filter(source_group == source_use, increased_in_M30d) %>%
      arrange(desc(increased_in_M15d))
    
    title_suffix <- "increased M_30d"
    file_suffix <- "increased_M_30d"
    
  } else {
    df.enrich <- df.pathway.compare %>%
      filter(source_group == source_use, increased_in_M15d_and_M30d )
    
    title_suffix <- "increased M_15d_and_M_30d "
    file_suffix <- "increased_M_15d_and_M_30d "
  }
  
  enriched_pathways <- df.enrich %>%
    slice_head(n = top_n) %>%
    pull(pathway_name) %>%
    unique()
  
  if (length(enriched_pathways) == 0) {
    message(source_use, " 没有找到 ", enriched_type, " 增强的 pathway")
    return(NULL)
  }
  
  df.plot <- df.pathway.all %>%
    filter(
      source_group == source_use,
      pathway_name %in% enriched_pathways
    ) %>%
    mutate(
      group = factor(group, levels = c("M_3d", "M_15d", "M_30d")),
      pathway_name = factor(pathway_name, levels = rev(enriched_pathways))
    )
  
  p <- ggplot(
    df.plot,
    aes(
      x = group,
      y = pathway_name,
      size = pathway_strength,
      color = pathway_strength
    )
  ) +
    geom_point(alpha = 0.9) +
    # facet_wrap(~ target_group, nrow = 1) +
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
      name = "Pathway strength"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "Pathway strength"
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(source_use, " → Oligodendrocyte: pathways ", title_suffix)
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(outdir, paste0("01.", safe_name(source_use), "_to_Olig_pathway_", file_suffix, ".pdf")),
    plot = p,
    device = cairo_pdf,
    width = 10,
    height = 8
  )
  
  ggsave(
    file.path(outdir, paste0("01.", safe_name(source_use), "_to_Olige_pathway_", file_suffix, ".png")),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  return(p)
}

for(soucesss in source_cells){
  p.pathway.high.g2 <- plot_pathway_enriched_type(
    soucesss, "M15d", df.pathway.compare, df.pathway.all, outdir, top_pathway_n
  )

  p.pathway.high.g3 <- plot_pathway_enriched_type(
    soucesss, "M30d", df.pathway.compare, df.pathway.all, outdir, top_pathway_n
  )

  p.pathway.high.both <- plot_pathway_enriched_type(
    soucesss, "both", df.pathway.compare, df.pathway.all, outdir, top_pathway_n
  )
}
## 8. LR pair 层面统计

df.lr.all <- df.all %>%
  group_by(group, source_group, target_group, pathway_name, ligand, receptor, LR_pair) %>%
  summarise(
    LR_strength = sum(prob, na.rm = TRUE),
    pval_min = min(pval, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  df.lr.all,
  file = file.path(outdir, "02.LR_strength_all_groups.csv"),
  row.names = FALSE
)

df.lr.compare <- df.lr.all %>%
  select(
    group,
    source_group,
    target_group,
    pathway_name,
    ligand,
    receptor,
    LR_pair,
    LR_strength
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = LR_strength,
    values_fill = 0
  )


df.lr.compare <- df.lr.compare %>%
  mutate(
    diff_M15d_vs_M3d = M_15d - M_3d,
    diff_M30d_vs_M3d = M_30d - M_3d,
    
    increased_in_M15d = diff_M15d_vs_M3d > 0,
    increased_in_M30d =  diff_M30d_vs_M3d > 0,
    increased_in_M15d_and_M30d = diff_M15d_vs_M3d > 0 &  diff_M30d_vs_M3d > 0,
    
    decreased_in_M15d = diff_M15d_vs_M3d < 0,
    decreased_in_M30d = diff_M30d_vs_M3d < 0,
    decreased_in_M15d_and_M30d = diff_M15d_vs_M3d < 0 & diff_M30d_vs_M3d < 0,
    
  ) %>%
  arrange(source_group, target_group)

write.csv(
  df.lr.compare,
  file = file.path(outdir, "02.LR_group2_group3_vs_group1_compare_with_flags.csv"),
  row.names = FALSE
)


## 9. LR bubble plot

plot_lr_enriched_type <- function(
    source_use,
    enriched_type = c("M15d", "M30d", "both"),
    df.lr.compare,
    df.lr.all,
    outdir,
    top_n = 40
) {
  
  enriched_type <- match.arg(enriched_type)
  
  if (enriched_type == "M15d") {
    df.enrich <- df.lr.compare %>%
      filter(source_group == source_use, increased_in_M15d) %>%
      arrange(desc(increased_in_M15d))
    
    title_suffix <- "increased M_15d"
    file_suffix <- "increased_M_15d"
    
  } else if (enriched_type == "M30d") {
    df.enrich <-df.lr.compare %>%
      filter(source_group == source_use, increased_in_M30d) %>%
      arrange(desc(increased_in_M15d))
    
    title_suffix <- "increased M_30d"
    file_suffix <- "increased_M_30d"
    
  } else {
    df.enrich <-df.lr.compare %>%
      filter(source_group == source_use, increased_in_M15d_and_M30d )
    
    title_suffix <- "increased M_15d_and_M_30d "
    file_suffix <- "increased_M_15d_and_M_30d "
  }
  enriched_lr <- df.enrich$LR_pair
  
  if (length(enriched_lr) == 0) {
    message(source_use, " 没有找到 ", enriched_type, " 增强的 LR pair")
    return(NULL)
  }
  
  df.plot <- df.lr.all %>%
    filter(
      source_group == source_use,
      LR_pair %in% enriched_lr
    ) %>%
    mutate(
      group = factor(group, levels = c("M_3d", "M_15d", "M_30d")),
      LR_pair = factor(LR_pair, levels = rev(enriched_lr))
    )
  
  p <- ggplot(
    df.plot,
    aes(
      x = group,
      y = LR_pair,
      size = LR_strength,
      color = LR_strength
    )
  ) +
    geom_point(alpha = 0.9) +
    # facet_wrap(~ target_group, nrow = 1) +
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
      name = "LR strength"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "LR strength"
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(source_use, " → Oligodendrocyte: LR pairs ", title_suffix)
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(outdir, paste0("02.", safe_name(source_use), "_to_Olig_LR_", file_suffix, ".pdf")),
    plot = p,
    device = cairo_pdf,
    width = max(8, length(target_cells) * 3),
    height = 10
  )
  
  ggsave(
    file.path(outdir, paste0("02.", safe_name(source_use), "_to_Olig_LR_", file_suffix, ".png")),
    plot = p,
    width = max(8, length(target_cells) * 3),
    height = 10,
    dpi = 300
  )
  
  return(p)
}
for(soucesss in source_cells){
  p.lr.high.g2 <- plot_lr_enriched_type(
    soucesss, "M15d", df.lr.compare, df.lr.all, outdir, top_lr_n
  )

  p.lr.high.g3 <- plot_lr_enriched_type(
    soucesss, "M30d", df.lr.compare, df.lr.all, outdir, top_lr_n
  )

  p.lr.high.both <- plot_lr_enriched_type(
    soucesss, "both", df.lr.compare, df.lr.all, outdir, top_lr_n
  )
}

## 9. 指定 pathway：按 group 分别绘制 circle plot 并拼图

library(grid)
library(gridExtra)

pathway_use <-  union(df.pathway.compare$pathway_name[df.pathway.compare$increased_in_M15d ==TRUE],df.pathway.compare$pathway_name[df.pathway.compare$increased_in_M30d ==TRUE])


plot_one_pathway_three_groups_base <- function(
    pathway_use,
    source_cells,
    target_cell,
    cellchat_g1,
    cellchat_g2,
    cellchat_g3,
    outdir,
    group_labels = c("M_3d", "M_15d", "M_30d"),
    edge.width.max = 2,
    vertex.label.cex = 1.0,
    width = 18,
    height = 6
) {
  
  outdir_circle <- file.path(
    outdir,
    paste0(
      "03.circle_compare_Microglia_subtypes",
      "_to_",
      safe_name(target_cell)
    )
  )
  dir.create(outdir_circle, showWarnings = FALSE, recursive = TRUE)
  
  pdf_file <- file.path(
    outdir_circle,
    paste0(
      "circle_compare_",
      safe_name(pathway_use), "_",
      "Microglia_subtypes",
      "_to_",
      safe_name(target_cell),
      "_M3d_M15d_M30d.pdf"
    )
  )
  
  png_file <- file.path(
    outdir_circle,
    paste0(
      "circle_compare_",
      safe_name(pathway_use), "_",
      "Microglia_subtypes",
      "_to_",
      safe_name(target_cell),
      "_M3d_M15d_M30d.png"
    )
  )
  
  plot_one <- function(cellchat, group_label) {
    
    cell_levels <- levels(cellchat@idents)
    
    if (!all(source_cells %in% cell_levels)) {
      message(group_label, " 缺少 source: ", paste(setdiff(source_cells, cell_levels), collapse = ", "))
      plot.new()
      title(main = paste0(group_label, "\nmissing source"))
      return(NULL)
    }
    
    if (!target_cell %in% cell_levels) {
      message(group_label, " 缺少 target: ", target_cell)
      plot.new()
      title(main = paste0(group_label, "\nmissing target"))
      return(NULL)
    }
    
    if (!pathway_use %in% cellchat@netP$pathways) {
      message(group_label, " 缺少 pathway: ", pathway_use)
      plot.new()
      title(main = paste0(group_label, "\nNo ", pathway_use))
      return(NULL)
    }
    
    vertex.receiver <- match(target_cell, cell_levels)
    
    tryCatch({
      netVisual_aggregate(
        cellchat,
        signaling = pathway_use,
        sources.use = source_cells,
        targets.use = target_cell,
        layout = "circle",
        vertex.receiver = vertex.receiver,
        vertex.label.cex = vertex.label.cex,
        edge.width.max = edge.width.max,
        signaling.name = paste0(
          group_label, "\n",
          pathway_use
        )
      )
    }, error = function(e) {
      message("绘图失败: ", group_label, " | ", pathway_use)
      message(e$message)
      plot.new()
      title(main = paste0(group_label, "\n", pathway_use, " failed"))
    })
  }
  
  ## -------------------------
  ## 保存 PDF
  ## -------------------------
  
  pdf(pdf_file, width = width, height = height)
  par(mfrow = c(1, 3), xpd = TRUE)
  
  plot_one(cellchat_g1, group_labels[1])
  plot_one(cellchat_g2, group_labels[2])
  plot_one(cellchat_g3, group_labels[3])
  
  dev.off()
  
  ## -------------------------
  ## 保存 PNG
  ## -------------------------
  
  png(png_file, width = width * 300, height = height * 300, res = 300)
  par(mfrow = c(1, 3), xpd = TRUE)
  
  plot_one(cellchat_g1, group_labels[1])
  plot_one(cellchat_g2, group_labels[2])
  plot_one(cellchat_g3, group_labels[3])
  
  dev.off()
  
  message("保存完成：", pdf_file)
  message("保存完成：", png_file)
}

## 运行：指定 pathway
src = source_cells 
tg = target_cell


for (pw in pathway_use) {
  plot_one_pathway_three_groups_base(
    pathway_use = pw,
    source_cells = source_cells,
    target_cell = target_cell,
    cellchat_g1 = cellchat_g1,
    cellchat_g2 = cellchat_g2,
    cellchat_g3 = cellchat_g3,
    outdir = outdir,
    group_labels = c("M_3d", "M_15d", "M_30d"),
    edge.width.max = 2,
    vertex.label.cex = 1.0,
    width = 18,
    height = 6
  )
}



# 组间比较筛选，亚群 ----------

rdata_g1_g2 <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/Compare_sub_oli/M_15d_M_3d/tmp_object.RData"
rdata_g1_g3 <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/Compare_sub_oli/M_30d_M_3d/tmp_object.RData"

outdir <- "/mnt/NexData_file/home/wmx/Project/SC2603312591-01/shouhou9/cellchat/Compare_sub_oli/Microglia_subtypes_to_Olig_group_compare"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

source_cells <- c("Microglia_Chka"  , "Microglia_Cmss1" , "Microglia_Ccl4"  , "Microglia_Csmd3" , "Microglia_Hspa1b", "Microglia_Tle4")


target_cells <-  c(
  "Oligodendrocyte_Kirrel3",
  "Oligodendrocyte_Grm7",
  "Oligodendrocyte_Cdkn1a",
  "Oligodendrocyte_Selenop",
  "Oligodendrocyte_Il33",
  "Oligodendrocyte_S100b",
  "Oligodendrocyte_Hebp1"
)

top_pathway_n <- 30
top_lr_n <- 40

safe_name <- function(x) {
  x %>%
    gsub("[/\\\\:*?\"<>| ]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}

## 2. 读取两个 RData 中的 object.list

load_cellchat_list <- function(rdata_file) {
  e <- new.env()
  load(rdata_file, envir = e)
  
  if (!exists("object.list", envir = e)) {
    stop("这个 RData 中没有 object.list: ", rdata_file)
  }
  
  obj.list <- get("object.list", envir = e)
  return(obj.list)
}

object.list.g12 <- load_cellchat_list(rdata_g1_g2)
object.list.g13 <- load_cellchat_list(rdata_g1_g3)

cellchat_g1 <- object.list.g12$M_3d
cellchat_g2 <- object.list.g12$M_15d
cellchat_g3 <- object.list.g13$M_30d

## 3. 提取 source → target 通讯结果

extract_comm_one <- function(cellchat, group_name, source_cells, target_cells) {
  
  cell_levels <- levels(cellchat@idents)
  
  if (!all(source_cells %in% cell_levels)) {
    stop(
      group_name, " 中 source_cells 不存在：",
      paste(setdiff(source_cells, cell_levels), collapse = ", "),
      "\n当前可用细胞类型：\n",
      paste(cell_levels, collapse = "\n")
    )
  }
  
  if (!all(target_cells %in% cell_levels)) {
    stop(
      group_name, " 中 target_cells 不存在：",
      paste(setdiff(target_cells, cell_levels), collapse = ", "),
      "\n当前可用细胞类型：\n",
      paste(cell_levels, collapse = "\n")
    )
  }
  
  df <- subsetCommunication(
    cellchat,
    sources.use = source_cells,
    targets.use = target_cells
  )
  
  if (nrow(df) == 0) {
    warning(group_name, " 没有提取到通讯结果")
    return(NULL)
  }
  
  df <- df %>%
    mutate(
      group = group_name,
      source_group = as.character(source),
      target_group = as.character(target),
      LR_pair = paste0(ligand, " - ", receptor)
    )
  
  return(df)
}

df.g1 <- extract_comm_one(cellchat_g1, "M_3d",  source_cells, target_cells)
df.g2 <- extract_comm_one(cellchat_g2, "M_15d", source_cells, target_cells)
df.g3 <- extract_comm_one(cellchat_g3, "M_30d", source_cells, target_cells)

df.all <- bind_rows(df.g1, df.g2, df.g3)

write.csv(
  df.all,
  file = file.path(outdir, "00.all_group_source_to_target_communication.csv"),
  row.names = FALSE
)


## 4. 总通讯强度 summary

df.summary.all <- df.all %>%
  group_by(group, source_group, target_group) %>%
  summarise(
    interaction_number = n(),
    interaction_strength = sum(prob, na.rm = TRUE),
    mean_strength = mean(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = c("M_3d", "M_15d", "M_30d")),
    source_group = factor(source_group, levels = source_cells)
  )

write.csv(
  df.summary.all,
  file = file.path(outdir, "00.summary_interaction_number_strength_all_groups.csv"),
  row.names = FALSE
)

p_summary_strength <- ggplot(
  df.summary.all,
  aes(
    x = group,
    y = interaction_strength,
    fill = source_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  facet_wrap(~ target_group, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Total interaction strength",
    fill = NULL,
    title = "Microglia → Oligodendrocyte subtypes interaction strength"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(outdir, "00.summary_total_interaction_strength.pdf"),
  plot = p_summary_strength,
  device = cairo_pdf,
  width = max(8, length(target_cells) * 4),
  height = 5
)

ggsave(
  file.path(outdir, "00.summary_total_interaction_strength.png"),
  plot = p_summary_strength,
  width = max(8, length(target_cells) * 4),
  height = 5,
  dpi = 300
)

p_summary_number <- ggplot(
  df.summary.all,
  aes(
    x = group,
    y = interaction_number,
    fill = source_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  facet_wrap(~ target_group, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Number of interactions",
    fill = NULL,
    title = "Microglia → Oligodendrocyte subtypes interaction number"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(outdir, "00.summary_interaction_number.pdf"),
  plot = p_summary_number,
  device = cairo_pdf,
  width = max(8, length(target_cells) * 4),
  height = 5
)

ggsave(
  file.path(outdir, "00.summary_interaction_number.png"),
  plot = p_summary_number,
  width = max(8, length(target_cells) * 4),
  height = 5,
  dpi = 300
)

## 5. pathway 层面统计

df.pathway.all <- df.all %>%
  group_by(group, source_group, target_group, pathway_name) %>%
  summarise(
    pathway_strength = sum(prob, na.rm = TRUE),
    interaction_number = n(),
    .groups = "drop"
  )

write.csv(
  df.pathway.all,
  file = file.path(outdir, "01.pathway_strength_all_groups.csv"),
  row.names = FALSE
)

df.pathway.compare <- df.pathway.all %>%
  select(group, source_group, target_group, pathway_name, pathway_strength) %>%
  pivot_wider(
    names_from = group,
    values_from = pathway_strength,
    values_fill = 0
  )

if (!"M_3d" %in% colnames(df.pathway.compare)) df.pathway.compare$M_3d <- 0
if (!"M_15d" %in% colnames(df.pathway.compare)) df.pathway.compare$M_15d <- 0
if (!"M_30d" %in% colnames(df.pathway.compare)) df.pathway.compare$M_30d <- 0

df.pathway.compare <- df.pathway.compare %>%
  mutate(
    diff_M15d_vs_M3d = M_15d - M_3d,
    diff_M30d_vs_M3d = M_30d - M_3d,
    
    increased_in_M15d = diff_M15d_vs_M3d > 0,
    increased_in_M30d = diff_M30d_vs_M3d > 0,
    increased_in_M15d_and_M30d = diff_M15d_vs_M3d > 0 & diff_M30d_vs_M3d > 0,
    
    decreased_in_M15d = diff_M15d_vs_M3d < 0,
    decreased_in_M30d = diff_M30d_vs_M3d < 0,
    decreased_in_M15d_and_M30d = diff_M15d_vs_M3d < 0 & diff_M30d_vs_M3d < 0,
    
  ) %>%
  arrange(source_group, target_group)

write.csv(
  df.pathway.compare,
  file = file.path(outdir, "01.pathway_M15d_M30d_vs_M3d_compare_with_flags.csv"),
  row.names = FALSE
)


## 7. pathway bubble plot

plot_pathway_enriched_type <- function(
    source_use,
    enriched_type = c("M15d", "M30d", "both"),
    df.pathway.compare,
    df.pathway.all,
    outdir,
    top_n = 30
) {
  
  enriched_type <- match.arg(enriched_type)
  
  if (enriched_type == "M15d") {
    df.enrich <- df.pathway.compare %>%
      filter(source_group == source_use, increased_in_M15d) %>%
      arrange(desc(diff_M15d_vs_M3d))
    
    title_suffix <- "increased in M_15d"
    file_suffix <- "increased_M_15d"
    
  } else if (enriched_type == "M30d") {
    df.enrich <- df.pathway.compare %>%
      filter(source_group == source_use, increased_in_M30d) %>%
      arrange(desc(diff_M30d_vs_M3d))
    
    title_suffix <- "increased in M_30d"
    file_suffix <- "increased_M_30d"
    
  } else {
    df.enrich <- df.pathway.compare %>%
      filter(source_group == source_use, increased_in_M15d_and_M30d) 
    
    title_suffix <- "increased in M_15d and M_30d"
    file_suffix <- "increased_M_15d_and_M_30d"
  }
  
  enriched_key <- df.enrich %>%
    group_by(target_group) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(key = paste(target_group, pathway_name, sep = "___")) %>%
    pull(key) %>%
    unique()
  
  if (length(enriched_key) == 0) {
    message(source_use, " 没有找到 ", enriched_type, " 增强的 pathway")
    return(NULL)
  }
  
  df.plot <- df.pathway.all %>%
    mutate(key = paste(target_group, pathway_name, sep = "___")) %>%
    filter(
      source_group == source_use,
      key %in% enriched_key
    ) %>%
    group_by(target_group) %>%
    mutate(
      pathway_name = factor(
        pathway_name,
        levels = rev(unique(pathway_name[order(pathway_strength, decreasing = TRUE)]))
      )
    ) %>%
    ungroup() %>%
    mutate(
      group = factor(group, levels = c("M_3d", "M_15d", "M_30d"))
    )
  
  p <- ggplot(
    df.plot,
    aes(
      x = group,
      y = pathway_name,
      size = pathway_strength,
      color = pathway_strength
    )
  ) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ target_group, nrow = 1, scales = "free_y") +
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
      name = "Pathway strength"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "Pathway strength"
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(source_use, " → Oligodendrocyte subtypes: pathways ", title_suffix)
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(outdir, paste0("01.", safe_name(source_use), "_to_Olig_pathway_", file_suffix, ".pdf")),
    plot = p,
    device = cairo_pdf,
    width = max(10, length(target_cells) * 4),
    height = 8
  )
  
  ggsave(
    file.path(outdir, paste0("01.", safe_name(source_use), "_to_Olig_pathway_", file_suffix, ".png")),
    plot = p,
    width = max(10, length(target_cells) * 4),
    height = 8,
    dpi = 300
  )
  
  return(p)
}

for(soucesss in source_cells){
  p.pathway.high.M15d <- plot_pathway_enriched_type(soucesss, "M15d", df.pathway.compare, df.pathway.all, outdir, top_pathway_n)
  p.pathway.high.M30d <- plot_pathway_enriched_type(soucesss, "M30d", df.pathway.compare, df.pathway.all, outdir, top_pathway_n)
  p.pathway.high.both <- plot_pathway_enriched_type(soucesss, "both", df.pathway.compare, df.pathway.all, outdir, top_pathway_n)
}
## 8. LR pair 层面统计

df.lr.all <- df.all %>%
  group_by(group, source_group, target_group, pathway_name, ligand, receptor, LR_pair) %>%
  summarise(
    LR_strength = sum(prob, na.rm = TRUE),
    pval_min = min(pval, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  df.lr.all,
  file = file.path(outdir, "02.LR_strength_all_groups.csv"),
  row.names = FALSE
)

df.lr.compare <- df.lr.all %>%
  select(
    group,
    source_group,
    target_group,
    pathway_name,
    ligand,
    receptor,
    LR_pair,
    LR_strength
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = LR_strength,
    values_fill = 0
  )

if (!"M_3d" %in% colnames(df.lr.compare)) df.lr.compare$M_3d <- 0
if (!"M_15d" %in% colnames(df.lr.compare)) df.lr.compare$M_15d <- 0
if (!"M_30d" %in% colnames(df.lr.compare)) df.lr.compare$M_30d <- 0

df.lr.compare <- df.lr.compare %>%
  mutate(
    diff_M15d_vs_M3d = M_15d - M_3d,
    diff_M30d_vs_M3d = M_30d - M_3d,

    increased_in_M15d = diff_M15d_vs_M3d > 0,
    increased_in_M30d = diff_M30d_vs_M3d > 0,
    increased_in_M15d_and_M30d = diff_M15d_vs_M3d > 0 & diff_M30d_vs_M3d > 0,
    
    decreased_in_M15d = diff_M15d_vs_M3d < 0,
    decreased_in_M30d = diff_M30d_vs_M3d < 0,
    decreased_in_M15d_and_M30d = diff_M15d_vs_M3d < 0 & diff_M30d_vs_M3d < 0,
    
  ) %>%
  arrange(source_group, target_group)

write.csv(
  df.lr.compare,
  file = file.path(outdir, "02.LR_M15d_M30d_vs_M3d_compare_with_flags.csv"),
  row.names = FALSE
)


## 10. LR bubble plot

plot_lr_enriched_type <- function(
    source_use,
    enriched_type = c("M15d", "M30d", "both"),
    df.lr.compare,
    df.lr.all,
    outdir,
    top_n = 40
) {
  
  enriched_type <- match.arg(enriched_type)
  
  if (enriched_type == "M15d") {
    df.enrich <- df.lr.compare %>%
      filter(source_group == source_use, increased_in_M15d) %>%
      arrange(desc(diff_M15d_vs_M3d))
    
    title_suffix <- "increased in M_15d"
    file_suffix <- "increased_M_15d"
    
  } else if (enriched_type == "M30d") {
    df.enrich <- df.lr.compare %>%
      filter(source_group == source_use, increased_in_M30d) %>%
      arrange(desc(diff_M30d_vs_M3d))
    
    title_suffix <- "increased in M_30d"
    file_suffix <- "increased_M_30d"
    
  } else {
    df.enrich <- df.lr.compare %>%
      filter(source_group == source_use, increased_in_M15d_and_M30d) 
    
    title_suffix <- "increased in M_15d and M_30d"
    file_suffix <- "increased_M_15d_and_M_30d"
  }
  
  enriched_key <- df.enrich %>%
    group_by(target_group) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(key = paste(target_group, LR_pair, sep = "___")) %>%
    pull(key) %>%
    unique()
  
  if (length(enriched_key) == 0) {
    message(source_use, " 没有找到 ", enriched_type, " 增强的 LR pair")
    return(NULL)
  }
  
  df.plot <- df.lr.all %>%
    mutate(key = paste(target_group, LR_pair, sep = "___")) %>%
    filter(
      source_group == source_use,
      key %in% enriched_key
    ) %>%
    group_by(target_group) %>%
    mutate(
      LR_pair = factor(
        LR_pair,
        levels = rev(unique(LR_pair[order(LR_strength, decreasing = TRUE)]))
      )
    ) %>%
    ungroup() %>%
    mutate(
      group = factor(group, levels = c("M_3d", "M_15d", "M_30d"))
    )
  
  p <- ggplot(
    df.plot,
    aes(
      x = group,
      y = LR_pair,
      size = LR_strength,
      color = LR_strength
    )
  ) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ target_group, nrow = 1, scales = "free_y") +
    scale_color_gradientn(
      colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#f4a582", "#b2182b"),
      name = "LR strength"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      name = "LR strength"
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(source_use, " → Oligodendrocyte subtypes: LR pairs ", title_suffix)
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(outdir, paste0("02.", safe_name(source_use), "_to_Olig_LR_", file_suffix, ".pdf")),
    plot = p,
    device = cairo_pdf,
    width = max(10, length(target_cells) * 4),
    height = 10
  )
  
  ggsave(
    file.path(outdir, paste0("02.", safe_name(source_use), "_to_Olig_LR_", file_suffix, ".png")),
    plot = p,
    width = max(10, length(target_cells) * 4),
    height = 10,
    dpi = 300
  )
  
  return(p)
}

for(soucesss in source_cells){
  p.lr.high.M15d <- plot_lr_enriched_type(soucesss, "M15d", df.lr.compare, df.lr.all, outdir, top_lr_n)
  p.lr.high.M30d <- plot_lr_enriched_type(soucesss, "M30d", df.lr.compare, df.lr.all, outdir, top_lr_n)
  p.lr.high.both <- plot_lr_enriched_type(soucesss, "both", df.lr.compare, df.lr.all, outdir, top_lr_n)
}

## 11. circle plot：多个 target 同时展示

pathway_use <- union(
  df.pathway.compare$pathway_name[df.pathway.compare$increased_in_M15d == TRUE],
  df.pathway.compare$pathway_name[df.pathway.compare$increased_in_M30d == TRUE]
)


plot_one_pathway_three_groups_base_multi_target <- function(
    pathway_use,
    source_cells,
    target_cells,
    cellchat_g1,
    cellchat_g2,
    cellchat_g3,
    outdir,
    group_labels = c("M_3d", "M_15d", "M_30d"),
    edge.width.max = 2,
    vertex.label.cex = 0.8,
    width = 18,
    height = 6
) {
  
  source_label <- "Microglia_subtypes"
  target_label <- paste(target_cells, collapse = "_")
  
  outdir_circle <- file.path(
    outdir,
    paste0(
      "03.circle_compare_Microglia_subtypes",
      "_to_",
      'Oligodendrocyte_subtype'
    )
  )
  dir.create(outdir_circle, showWarnings = FALSE, recursive = TRUE)
  
  pdf_file <- file.path(
    outdir_circle,
    paste0(
      "circle_compare_",
      safe_name(pathway_use), "_",
      safe_name(source_label),
      "_to_",
      'Oligodendrocyte_subtype',
      "_M3d_M15d_M30d.pdf"
    )
  )
  
  png_file <- file.path(
    outdir_circle,
    paste0(
      "circle_compare_",
      safe_name(pathway_use), "_",
      safe_name(source_label),
      "_to_",
      'Oligodendrocyte_subtype',
      "_M3d_M15d_M30d.png"
    )
  )
  
  plot_one <- function(cellchat, group_label) {
    
    cell_levels <- levels(cellchat@idents)
    
    missing_sources <- setdiff(source_cells, cell_levels)
    missing_targets <- setdiff(target_cells, cell_levels)
    
    if (length(missing_sources) > 0) {
      message(group_label, " 缺少 source: ", paste(missing_sources, collapse = ", "))
      plot.new()
      title(main = paste0(group_label, "\nmissing source"))
      return(NULL)
    }
    
    if (length(missing_targets) > 0) {
      message(group_label, " 缺少 target: ", paste(missing_targets, collapse = ", "))
      plot.new()
      title(main = paste0(group_label, "\nmissing target"))
      return(NULL)
    }
    
    if (!pathway_use %in% cellchat@netP$pathways) {
      message(group_label, " 缺少 pathway: ", pathway_use)
      plot.new()
      title(main = paste0(group_label, "\nNo ", pathway_use))
      return(NULL)
    }
    
    vertex.receiver <- match(target_cells, cell_levels)
    
    tryCatch({
      netVisual_aggregate(
        cellchat,
        signaling = pathway_use,
        sources.use = source_cells,
        targets.use = target_cells,
        layout = "circle",
        vertex.receiver = vertex.receiver,
        vertex.label.cex = vertex.label.cex,
        edge.width.max = edge.width.max,
        signaling.name = paste0(
          group_label, "\n",
          pathway_use
        )
      )
    }, error = function(e) {
      message("绘图失败: ", group_label, " | ", pathway_use)
      message(e$message)
      plot.new()
      title(main = paste0(group_label, "\n", pathway_use, " failed"))
    })
  }
  
  pdf(pdf_file, width = width, height = height)
  par(mfrow = c(1, 3), xpd = TRUE)
  plot_one(cellchat_g1, group_labels[1])
  plot_one(cellchat_g2, group_labels[2])
  plot_one(cellchat_g3, group_labels[3])
  dev.off()
  
  png(png_file, width = width * 300, height = height * 300, res = 300)
  par(mfrow = c(1, 3), xpd = TRUE)
  plot_one(cellchat_g1, group_labels[1])
  plot_one(cellchat_g2, group_labels[2])
  plot_one(cellchat_g3, group_labels[3])
  dev.off()
  
  message("保存完成：", pdf_file)
  message("保存完成：", png_file)
}

for (pw in pathway_use) {
  plot_one_pathway_three_groups_base_multi_target(
    pathway_use = pw,
    source_cells = source_cells,
    target_cells = target_cells,
    cellchat_g1 = cellchat_g1,
    cellchat_g2 = cellchat_g2,
    cellchat_g3 = cellchat_g3,
    outdir = outdir,
    group_labels = c("M_3d", "M_15d", "M_30d"),
    edge.width.max = 2,
    vertex.label.cex = 0.8,
    width = 18,
    height = 6
  )
}
