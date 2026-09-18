# ============================================================
# Re-root Monocle2 pseudotime
# 只修改 root_state / reverse，不重新跑 reduceDimension
# ============================================================

library('igraph', lib.loc = '/mnt/NexData_file/home/wmx/envs/scRNA/lib/R/library/')
.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle2/lib/R/library/")

suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(reshape2)
  library(getopt)
})

options(stringsAsFactors = FALSE, warn = -1)

# ============================================================
# 参数
# ============================================================

command <- matrix(c(
  'rds',       'r', 1, "character",
  'outpath',   'o', 1, "character",
  'RootState', 'S', 1, "numeric",
  'Reverse',   'v', 2, "logical",
  'plotcols',  'P', 2, "character"
), byrow = TRUE, ncol = 4)

args <- getopt(command)

if (is.null(args$rds)) {
  stop("Please provide monocle rds file by -r")
}

if (is.null(args$outpath)) {
  stop("Please provide output path by -o")
}

if (is.null(args$RootState)) {
  stop("Please provide RootState by -S")
}

if (is.null(args$Reverse)) {
  args$Reverse <- FALSE
}

plotcols <- args$plotcols
if (is.null(plotcols)) {
  plotcols <- c("State", "Pseudotime", "Cluster", "sub_cell_type")
} else {
  plotcols <- unlist(strsplit(plotcols, ","))
  plotcols <- trimws(plotcols)
}

cds_file <- args$rds
outdir <- args$outpath
RootState <- args$RootState
Reverse <- args$Reverse

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

message("========================================")
message("[Input monocle rds] ", cds_file)
message("[Out path] ", outdir)
message("[RootState] ", RootState)
message("[Reverse] ", Reverse)
message("[Plot columns] ", paste(plotcols, collapse = ","))
message("========================================")

# ============================================================
# 通用函数
# ============================================================

safe_name <- function(x) {
  x <- gsub("[/\\\\:;,*?\"<>| ]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

plot_theme <- function() {
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )
}

save_plot_pdf_png <- function(gg, prefix, width = 8, height = 6) {
  ggsave(paste0(prefix, ".png"), gg, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(prefix, ".pdf"), gg, width = width, height = height, limitsize = FALSE)
}

get_real_colname <- function(cds, col) {
  meta_cols <- colnames(pData(cds))
  hit <- meta_cols[tolower(meta_cols) == tolower(col)]
  
  if (length(hit) > 0) {
    return(hit[1])
  } else {
    return(NA)
  }
}

get_nrow_facet <- function(x) {
  n <- length(unique(na.omit(x)))
  if (n <= 1) return(1)
  out <- floor(sqrt(n))
  if (out < 1) out <- 1
  return(out)
}

# ============================================================
# 绘图函数
# ============================================================

PsePlotByCols <- function(cds,
                          plotcols = c("State", "Pseudotime"),
                          do_facet = TRUE,
                          do_tree = TRUE,
                          do_density = TRUE) {
  
  meta <- pData(cds)
  
  message("========================================")
  message("[Plot] Available columns:")
  message(paste(colnames(meta), collapse = ", "))
  message("========================================")
  
  for (col in plotcols) {
    
    real_col <- get_real_colname(cds, col)
    
    if (is.na(real_col)) {
      message("[Skip] Column not found: ", col)
      next
    }
    
    safe_col <- safe_name(real_col)
    
    message("----------------------------------------")
    message("[Plot] color_by = ", real_col)
    
    gg <- plot_cell_trajectory(cds, color_by = real_col) +
      plot_theme()
    
    save_plot_pdf_png(
      gg = gg,
      prefix = paste0("Pseudotime.", safe_col),
      width = 8,
      height = 6
    )
    
    x <- meta[[real_col]]
    n_unique <- length(unique(na.omit(x)))
    
    is_discrete <- is.factor(x) ||
      is.character(x) ||
      real_col %in% c("State", "Cluster") ||
      n_unique <= 30
    
    if (do_facet && is_discrete && n_unique > 1) {
      
      nrow_use <- get_nrow_facet(x)
      ncol_use <- ceiling(n_unique / nrow_use)
      
      gg_facet <- plot_cell_trajectory(cds, color_by = real_col) +
        facet_wrap(as.formula(paste0("~", real_col)), nrow = nrow_use) +
        plot_theme()
      
      save_plot_pdf_png(
        gg = gg_facet,
        prefix = paste0("Pseudotime.", safe_col, ".Each"),
        width = max(8, 4 * ncol_use),
        height = max(5, 3 * nrow_use)
      )
    }
    
    if (do_tree) {
      
      gg_tree <- plot_complex_cell_trajectory(
        cds,
        color_by = real_col,
        show_branch_points = TRUE,
        cell_size = 0.5,
        cell_link_size = 0.3
      ) +
        plot_theme()
      
      save_plot_pdf_png(
        gg = gg_tree,
        prefix = paste0("Pseudotime.", safe_col, ".tree"),
        width = 8,
        height = 6
      )
    }
    
    if (do_density &&
        "Pseudotime" %in% colnames(meta) &&
        is_discrete &&
        real_col != "Pseudotime" &&
        n_unique > 1) {
      
      den_df <- data.frame(
        Pseudotime = meta$Pseudotime,
        group = meta[[real_col]]
      )
      
      den_df <- den_df[
        !is.na(den_df$Pseudotime) & !is.na(den_df$group),
      ]
      
      if (nrow(den_df) > 0 && length(unique(den_df$group)) > 1) {
        
        gg_den <- ggplot(den_df, aes(x = Pseudotime, fill = group)) +
          geom_density(alpha = 0.5) +
          theme_cowplot() +
          theme(legend.position = "right")
        
        save_plot_pdf_png(
          gg = gg_den,
          prefix = paste0("PseudotimeDensity.", safe_col),
          width = 8,
          height = 5
        )
        
        n_group <- length(unique(den_df$group))
        nrow_use <- floor(sqrt(n_group))
        if (nrow_use < 1) nrow_use <- 1
        ncol_use <- ceiling(n_group / nrow_use)
        
        gg_den_facet <- ggplot(den_df, aes(x = Pseudotime, fill = group)) +
          geom_density(alpha = 0.5) +
          facet_wrap(~group, nrow = nrow_use) +
          theme_cowplot() +
          theme(legend.position = "right")
        
        save_plot_pdf_png(
          gg = gg_den_facet,
          prefix = paste0("PseudotimeDensity.", safe_col, ".Each"),
          width = max(8, 4 * ncol_use),
          height = max(5, 2.5 * nrow_use)
        )
      }
    }
  }
  
  if (do_density && "Pseudotime" %in% colnames(meta)) {
    
    den_df <- data.frame(Pseudotime = meta$Pseudotime)
    den_df <- den_df[!is.na(den_df$Pseudotime), , drop = FALSE]
    
    if (nrow(den_df) > 0) {
      gg_den_all <- ggplot(den_df, aes(x = Pseudotime)) +
        geom_density(alpha = 0.5, fill = "blue") +
        theme_cowplot()
      
      save_plot_pdf_png(
        gg = gg_den_all,
        prefix = "PseudotimeDensity",
        width = 8,
        height = 5
      )
    }
  }
  
  message("[Plot] All selected columns finished.")
}

# ============================================================
# 主流程
# ============================================================

message("[Read] ", cds_file)
cds <- readRDS(cds_file)

message("========================================")
message("[Before reorder] State table:")
print(table(pData(cds)$State))

if ("sub_cell_type" %in% colnames(pData(cds))) {
  message("[Before reorder] State vs sub_cell_type:")
  print(table(pData(cds)$State, pData(cds)$sub_cell_type))
}

if ("cell_type" %in% colnames(pData(cds))) {
  message("[Before reorder] State vs cell_type:")
  print(table(pData(cds)$State, pData(cds)$cell_type))
}

message("========================================")
message("[Reorder cells]")
cds <- orderCells(
  cds,
  root_state = RootState,
  reverse = Reverse
)

message("[After reorder] Pseudotime summary:")
print(summary(pData(cds)$Pseudotime))

# ============================================================
# 保存统计表
# ============================================================

if ("Cluster" %in% colnames(pData(cds)) &&
    "State" %in% colnames(pData(cds))) {
  
  stat_table <- table(pData(cds)$Cluster, pData(cds)$State) %>%
    as.data.frame() %>%
    dcast(Var1 ~ Var2)
  
  colnames(stat_table)[1] <- "Cluster"
  
  write.table(
    stat_table,
    file = "Pseudotime.Statistics.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# 保存新的 pseudotime 信息
pData(cds) %>%
  data.frame(CellName = row.names(pData(cds))) %>%
  write.table(
    file = "Pseudotime.Summary.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

# 保存新的 monocle 对象
saveRDS(cds, file = "Pseudotime.monocle.reRoot.rds")

# 重新画图
PsePlotByCols(
  cds = cds,
  plotcols = plotcols,
  do_facet = TRUE,
  do_tree = TRUE,
  do_density = TRUE
)

message("========================================")
message("[Done]")
message("[Output] ", outdir)
message("========================================")