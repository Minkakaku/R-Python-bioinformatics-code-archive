# 这一部分是对cellchat的结果图进行调整
# 因为细胞类型太多了，且名字很长，先对报告里的图进行调整
# 针对T --> 上皮，上皮 --> B 绘制一些图

.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/lib/R/library/")
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(pdftools,lib.loc = "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/lib/R/library/")
library(stringr)
library(jpeg)
library(rmarkdown)

# ============================================================
# 1. 参数设置
# ============================================================

#  CellChat 输出的总目录
root_dir <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou19/cellchat/"

# 是否把 PDF 转成 PNG
convert_pdf_to_png <- TRUE

# source 细胞类型：例如展示 T 细胞到上皮
# 可以用正则匹配
source_pattern <- c("CD4_Tfh","CD8_Tm")

# target 细胞类型：上皮相关
target_pattern  <- c("CEACAM7+ Tumor","IGKC+ Tumor","KRT16+ Tumor","MUC5B+ Epi","MMP20+ Tumor","GZMM+ Tumor","GDPD2+ Tumor")

measure_list <- c("count", "weight")
heatmap_color <- c("white", "red2")


# ============================================================
# 小函数：安全保存 PDF
# ============================================================

safe_pdf <- function(file, width = 8, height = 7, expr) {
  pdf(file, width = width, height = height)
  tryCatch(
    {
      force(expr)
    },
    error = function(e) {
      message("[ERROR] ", file)
      message(e$message)
    },
    finally = {
      dev.off()
    }
  )
}

# ============================================================
# pathway整体气泡图：一个方向只画一张完整pathway bubble
# ============================================================
plot_pathway_bubble_all <- function(cellchat_obj, group_name, outdir,
                                    senders, receivers,
                                    direction_name,
                                    title_name) {
  
  df.netp <- subsetCommunication(cellchat_obj, slot.name = "netP")
  
  message("    netP columns: ", paste(colnames(df.netp), collapse = ", "))
  
  # 自动识别 pathway 列
  pathway_col <- intersect(
    c("pathway", "pathway_name", "signaling", "pathway.name"),
    colnames(df.netp)
  )
  
  if (length(pathway_col) == 0) {
    message("[SKIP] netP 中没有找到 pathway/pathway_name/signaling 列")
    write.table(
      df.netp,
      paste0(outdir, "/", group_name, "_", direction_name, "_netP_raw.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    return(NULL)
  }
  
  pathway_col <- pathway_col[1]
  
  # 自动识别 prob 列
  prob_col <- intersect(
    c("prob", "probability", "communication_probability"),
    colnames(df.netp)
  )
  
  if (length(prob_col) == 0) {
    message("[SKIP] netP 中没有找到 prob 列")
    write.table(
      df.netp,
      paste0(outdir, "/", group_name, "_", direction_name, "_netP_raw.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    return(NULL)
  }
  
  prob_col <- prob_col[1]
  
  # 自动识别 pval 列
  pval_col <- intersect(
    c("pval", "p.value", "pvalue", "p"),
    colnames(df.netp)
  )
  
  if (length(pval_col) == 0) {
    df.netp$pval_use <- NA
    pval_col <- "pval_use"
  } else {
    pval_col <- pval_col[1]
  }
  
  # 统一字段
  df.netp$pathway_use <- as.character(df.netp[[pathway_col]])
  df.netp$prob_use <- as.numeric(df.netp[[prob_col]])
  df.netp$pval_use <- as.numeric(df.netp[[pval_col]])
  
  df.pathway <- df.netp[
    df.netp$source %in% senders & df.netp$target %in% receivers,
  ]
  
  if (nrow(df.pathway) == 0) {
    message("[SKIP] no pathway communication for ", group_name, " ", direction_name)
    return(NULL)
  }
  
  df.pathway$interaction <- paste0(df.pathway$source, " -> ", df.pathway$target)
  
  # 保存 pathway 表
  write.table(
    df.pathway,
    paste0(outdir, "/", group_name, "_", direction_name, "_all_pathway.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # pathway 排序：按照总通讯强度
  pathway_order <- df.pathway %>%
    dplyr::group_by(pathway_use) %>%
    dplyr::summarise(
      total_prob = sum(prob_use, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(total_prob) %>%
    dplyr::pull(pathway_use)
  
  # x轴排序：按照总通讯强度
  interaction_order <- df.pathway %>%
    dplyr::group_by(interaction) %>%
    dplyr::summarise(
      total_prob = sum(prob_use, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(total_prob) %>%
    dplyr::pull(interaction)
  
  df.pathway$pathway_use <- factor(df.pathway$pathway_use, levels = pathway_order)
  df.pathway$interaction <- factor(df.pathway$interaction, levels = interaction_order)
  
  df.pathway$pval_plot <- df.pathway$pval_use
  df.pathway$pval_plot <- df.pathway$pval_use
  df.pathway$pval_plot[is.na(df.pathway$pval_plot)] <- 1
  
  df.pathway$p_group <- dplyr::case_when(
    df.pathway$pval_plot < 0.01 ~ "p < 0.01",
    df.pathway$pval_plot < 0.05 ~ "p < 0.05",
    TRUE ~ "p > 0.05"
  )
  
  df.pathway$p_group <- factor(
    df.pathway$p_group,
    levels = c("p < 0.01", "p < 0.05", "p > 0.05")
  )
  
  p <- ggplot(
    df.pathway,
    aes(
      x = interaction,
      y = pathway_use
    )
  ) +
    geom_point(
      aes(
        color = prob_use,
        size = p_group
      )
    ) +
    scale_color_gradientn(
      colours = c(
        "#3B4CC0",
        "#5AA9A6",
        "#A6D96A",
        "#FEE08B",
        "#F46D43",
        "#B2182B"
      ),
      name = "Commun. Prob.",
      na.value = "grey90"
    ) +
    scale_size_manual(
      name = "p-value",
      values = c(
        "p < 0.01" = 4.2,
        "p < 0.05" = 2.8,
        "p > 0.05" = 1.4
      ),
      breaks = c("p < 0.01", "p < 0.05", "p > 0.05")
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 8,
        color = "black"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.35, "cm")
    ) +
    labs(
      title = paste0(group_name, " ", title_name, " pathway-level bubble")
    )
  h_use <- max(5, 0.28 * length(unique(df.pathway$pathway_use)) + 2)
  w_use <- max(8, 0.55 * length(unique(df.pathway$interaction)) + 4)
  safe_pdf(
    file = paste0(outdir, "/", group_name, "_", direction_name, "_all_pathway_bubble.pdf"),
    width = w_use,
    height = h_use,
    expr = {
      print(p)
    }
  )
  
  ggsave(
    paste0(outdir, "/", group_name, "_", direction_name, "_all_pathway_bubble.png"),
    p,
    width = w_use,
    height = h_use,
    dpi = 300
  )
  
  message("    pathway number: ", length(unique(df.pathway$pathway_use)))
}
# ============================================================
# 新增函数：source -> target 和 target -> source 双向绘图
# 1. Number circle
# 2. Weight circle
# 3. LR pair bubble
# 4. pathway整体bubble
# 5. LR表和pathway表
# ============================================================

plot_source_to_target <- function(cellchat_obj, group_name, outdir,
                                  source_use = NULL,
                                  target_use = NULL) {
  
  mat_count <- cellchat_obj@net$count
  mat_weight <- cellchat_obj@net$weight
  
  all_celltypes <- rownames(mat_count)
  
  source_use <- intersect(source_use, all_celltypes)
  target_use <- intersect(target_use, all_celltypes)
  
  message("[", group_name, "] source_use: ", paste(source_use, collapse = ", "))
  message("[", group_name, "] target_use: ", paste(target_use, collapse = ", "))
  
  if (length(source_use) == 0 || length(target_use) == 0) {
    message("[SKIP] source 或 target 没有匹配到细胞类型")
    return(NULL)
  }
  
  plot_one_direction <- function(senders, receivers, direction_name, title_name) {
    
    message("  [PLOT] ", direction_name)
    message("    senders: ", paste(senders, collapse = ", "))
    message("    receivers: ", paste(receivers, collapse = ", "))
    
    # ----------------------------
    # count matrix
    # ----------------------------
    
    mat_count2 <- matrix(
      0,
      nrow = nrow(mat_count),
      ncol = ncol(mat_count),
      dimnames = dimnames(mat_count)
    )
    
    mat_count2[senders, receivers] <- mat_count[senders, receivers, drop = FALSE]
    
    # ----------------------------
    # weight matrix
    # ----------------------------
    
    mat_weight2 <- matrix(
      0,
      nrow = nrow(mat_weight),
      ncol = ncol(mat_weight),
      dimnames = dimnames(mat_weight)
    )
    
    mat_weight2[senders, receivers] <- mat_weight[senders, receivers, drop = FALSE]
    
    # ----------------------------
    # Number circle
    # ----------------------------
    
    safe_pdf(
      file = paste0(outdir, "/", group_name, "_", direction_name, "_Number_netVisual_circle.pdf"),
      width = 9,
      height = 8,
      expr = {
        netVisual_circle(
          mat_count2,
          vertex.weight = 6,
          weight.scale = TRUE,
          label.edge = FALSE,
          edge.weight.max = max(mat_count, na.rm = TRUE),
          edge.width.max = 10,
          vertex.label.cex = 0.6,
          title.name = paste0(group_name, " ", title_name, " Number")
        )
      }
    )
    
    # ----------------------------
    # Weight circle
    # ----------------------------
    
    safe_pdf(
      file = paste0(outdir, "/", group_name, "_", direction_name, "_Weight_netVisual_circle.pdf"),
      width = 12,
      height = 10,
      expr = {
        netVisual_circle(
          mat_weight2,
          vertex.weight = 6,
          weight.scale = TRUE,
          label.edge = FALSE,
          edge.weight.max = max(mat_weight, na.rm = TRUE),
          edge.width.max = 10,
          vertex.label.cex = 0.6,
          title.name = paste0(group_name, " ", title_name, " Weight")
        )
      }
    )
    
    # ----------------------------
    # ligand-receptor table
    # ----------------------------
    
    df.net <- subsetCommunication(cellchat_obj)
    
    df.st <- df.net[
      df.net$source %in% senders & df.net$target %in% receivers,
    ]
    
    write.table(
      df.st,
      paste0(outdir, "/", group_name, "_", direction_name, "_lr.txt"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    # ----------------------------
    # LR pair bubble plot
    # ----------------------------
    
    if (nrow(df.st) > 0) {
      
      safe_pdf(
        file = paste0(outdir, "/", group_name, "_", direction_name, "_bubble.pdf"),
        width = 10,
        height = 12,
        expr = {
          p <- netVisual_bubble(
            cellchat_obj,
            sources.use = senders,
            targets.use = receivers,
            remove.isolate = FALSE,
            angle.x = 45
          )
          print(p)
        }
      )
      
    } else {
      message("    [SKIP] no ligand-receptor pairs for ", direction_name)
    }
    
    # ----------------------------
    # pathway-level table + 整体pathway bubble
    # 不再单独画每个pathway
    # ----------------------------
    
    plot_pathway_bubble_all(
      cellchat_obj = cellchat_obj,
      group_name = group_name,
      outdir = outdir,
      senders = senders,
      receivers = receivers,
      direction_name = direction_name,
      title_name = title_name
    )
  }
  
  # ============================================================
  # 1. source -> target
  # ============================================================
  
  plot_one_direction(
    senders = source_use,
    receivers = target_use,
    direction_name = "source_to_target",
    title_name = "source -> target"
  )
  
  # ============================================================
  # 2. target -> source
  # ============================================================
  
  plot_one_direction(
    senders = target_use,
    receivers = source_use,
    direction_name = "target_to_source",
    title_name = "target -> source"
  )
}
# ============================================================
# 主循环：读取所有 tmp_object.RData
# ============================================================

tmp_files <- list.files(
  root_dir,
  pattern = "^tmp_object\\.RData$",
  recursive = TRUE,
  full.names = TRUE
)

message("[INFO] 找到 tmp_object.RData 数量：", length(tmp_files))

for (tmp_file in tmp_files) {
  
  message("====================================")
  message("[LOAD] ", tmp_file)
  
  load(tmp_file)  # 读入 object.list
  
  old_outdir <- dirname(tmp_file)
  outdir <- file.path(old_outdir, "new_plot")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  message("[OUTDIR] ", outdir)
  
  group_names <- names(object.list)
  if (is.null(group_names)) {
    group_names <- paste0("group", seq_along(object.list))
  }
  
  common_cell_types <- Reduce(
    intersect,
    lapply(object.list, function(x) levels(droplevels(x@idents)))
  )
  
  num_celltypes <- length(common_cell_types)
  if (num_celltypes <= 4) {
    n <- c(1, num_celltypes)
  } else {
    n <- c(ceiling(num_celltypes / 2), 3)
  }
  
  if (length(common_cell_types) <= 6) {
    labeledge <- TRUE
  } else {
    labeledge <- FALSE
  }
  
  if (length(object.list) == 1) {
    title_Number <- paste0("Number of interactions - ", names(object.list)[1])
    title_Weight <- paste0("Weight of interactions - ", names(object.list)[1])
    title_name <- names(object.list)[1]
  } else {
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    title_Number <- c(
      paste0("Number of interactions - ", names(object.list)[1]),
      paste0("Number of interactions - ", names(object.list)[2])
    )
    title_Weight <- c(
      paste0("Weight of interactions - ", names(object.list)[1]),
      paste0("Weight of interactions - ", names(object.list)[2])
    )
    title_name <- c(names(object.list)[1], names(object.list)[2])
  }
  
  title_list <- list(title_Number, title_Weight)
  
  # ========================================================
  # 1. 原来的 group heatmap，文件名不变
  # ========================================================
  
  for (i in seq_along(measure_list)) {
    
    gg1 <- netVisual_heatmap(
      object.list[[1]],
      measure = measure_list[i],
      color.heatmap = heatmap_color,
      title.name = title_list[[i]][1]
    )
    
    if (length(object.list) == 1) {
      gg <- gg1
    } else {
      gg2 <- netVisual_heatmap(
        object.list[[2]],
        measure = measure_list[i],
        color.heatmap = heatmap_color,
        title.name = title_list[[i]][2]
      )
      gg <- gg1 + gg2
    }
    
    if (i == 1) {
      label <- "Number"
    } else {
      label <- "Weight"
    }
    
    safe_pdf(
      file = paste0(outdir, "/group_", label, "_netVisual_heatmap.pdf"),
      width = ifelse(length(object.list) == 1, 10, 20),
      height = 10,
      expr = {
        print(gg)
      }
    )
  }
  

  # ========================================================
  # 2. 新增：source -> target 图
  # 文件名是新增的，不影响原来文件
  # ========================================================
  
  for (j in seq_along(object.list)) {
    
    plot_source_to_target(
      cellchat_obj = object.list[[j]],
      group_name = names(object.list)[j],
      outdir = outdir,
      source_use = source_pattern,
      target_use =  target_pattern
    )
  }
  
  outdir1 <- paste0(outdir,'/B_Epith/')
  dir.create(outdir1,recursive = T)
  source_pattern1 <- c("Naive_B" ,"Activated_B","Memory_B")
  for (j in seq_along(object.list)) {
    
    plot_source_to_target(
      cellchat_obj = object.list[[j]],
      group_name = names(object.list)[j],
      outdir = outdir1,
      source_use = source_pattern1,
      target_use =  target_pattern
    )
  }
  
  # ========================================================
  # 3. 补充剩余图片
  # ========================================================
  if (length(object.list) == 1) {
    
    pathway.union_intersection <- object.list[[1]]@netP$pathways
    height_h <- ceiling(length(pathway.union_intersection) / 10) * 2.5
    safe_pdf(
      file = paste0(outdir, "/incoming.pdf"),
      width = 15,
      height = 0.5*height_h+1,
      expr = {
        ht <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = 'incoming', signaling = pathway.union_intersection,title = 'all',width = 15, height = height_h)
        print(ht)
      }
    )
    safe_pdf(
      file = paste0(outdir, "/outgoing.pdf"),
      width = 15,
      height = 0.5*height_h+1,
      expr = {
        ht <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = 'outgoing', signaling = pathway.union_intersection,title = 'all',width = 15, height = height_h)
        print(ht)
      }
    )
    
  }
  
  if (length(object.list) >= 2) {
    
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    safe_pdf(
      file = paste0(outdir, "/Interactions_diff_heatmap.pdf"),
      width = 17,
      height = 12,
      expr = {
        gg1 <- netVisual_heatmap(cellchat, comparison = c(2, 1))
        gg2 <- netVisual_heatmap(cellchat, comparison = c(2, 1), measure = "weight")
        print(gg1 + gg2)
      }
    )
    pathway.union_intersection <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
    height_h <- ceiling(length(pathway.union_intersection) / 10) * 2.5
    safe_pdf(
      file = paste0(outdir, "/incoming.pdf"),
      width = 20,
      height = 0.5*height_h+1,
      expr = {
        ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = 'incoming', signaling = pathway.union_intersection,title = names(object.list)[1],width = 15, height = height_h)
        ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = 'incoming',signaling = pathway.union_intersection,title = names(object.list)[2],width = 15, height = height_h)
        ht <- ht1 + ht2
        print(ht)
      }
    )
    safe_pdf(
      file = paste0(outdir, "/outgoing.pdf"),
      width = 20,
      height = 0.5*height_h+1,
      expr = {
        ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = 'outgoing', signaling = pathway.union_intersection,title = names(object.list)[1],width = 15, height = height_h)
        ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = 'outgoing',signaling = pathway.union_intersection,title = names(object.list)[2],width = 15, height = height_h)
        ht <- ht1 + ht2
        print(ht)
      }
    )
    
    safe_pdf(
      file = paste0(outdir, "/source_target_LR_bubble_compare.pdf"),
      width = 12,
      height = 12,
      expr = {
        
        p <- netVisual_bubble(
          cellchat,
          sources.use = source_pattern,
          targets.use = target_pattern,
          comparison = c(1,2),
          remove.isolate = FALSE,
          angle.x = 45
        )
        
        print(p)
      }
    )
    safe_pdf(
      file = paste0(outdir, "/target_source_LR_bubble_compare.pdf"),
      width = 12,
      height = 12,
      expr = {
        
        p <- netVisual_bubble(
          cellchat,
          sources.use = target_pattern,
          targets.use = source_pattern,
          comparison = c(1,2),
          remove.isolate = FALSE,
          angle.x = 45
        )
        
        print(p)
      }
    )
    safe_pdf(
      file = paste0(outdir1, "/source_target_LR_bubble_compare.pdf"),
      width = 12,
      height = 12,
      expr = {
        
        p <- netVisual_bubble(
          cellchat,
          sources.use = source_pattern1,
          targets.use = target_pattern,
          comparison = c(1,2),
          remove.isolate = FALSE,
          angle.x = 45
        )
        
        print(p)
      }
    )
    safe_pdf(
      file = paste0(outdir1, "/target_source_LR_bubble_compare.pdf"),
      width = 12,
      height = 12,
      expr = {
        
        p <- netVisual_bubble(
          cellchat,
          sources.use = target_pattern,
          targets.use = source_pattern1,
          comparison = c(1,2),
          remove.isolate = FALSE,
          angle.x = 45
        )
        
        print(p)
      }
    )
    
  }
  
  # ========================================================
  # 4. PDF 转 PNG，文件名仍然对应原 PDF
  # ========================================================
  
  pdf_files <- list.files(
    path = outdir,
    pattern = "\\.pdf$",
    full.names = TRUE
  )
  pdf_files1 <- list.files(
    path = outdir1,
    pattern = "\\.pdf$",
    full.names = TRUE
  )
  pdf_files <- c(pdf_files1,pdf_files)
  
  for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    tryCatch(
      {
        pdf_convert(
          file_name,
          dpi = 300,
          filenames = sub(
            "\\.pdf$",
            ".png",
            file_name
          )
        )
      },
      error = function(e) {
        message("[WARN] PDF 转 PNG 失败：", file_name)
        message(e$message)
      }
    )
  }
  
}



