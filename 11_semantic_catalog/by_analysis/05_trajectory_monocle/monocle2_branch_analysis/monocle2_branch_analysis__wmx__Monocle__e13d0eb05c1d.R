# Seurat4Pseudotime
# ------------------------------------------------------------
# Usage example:
# Rscript monocle_pseudotime.R \
#   -r input.rds \
#   -l 0.001 \
#   -o ./output/ \
#   -P cell_type,state,pseudotime,group,sample
# ------------------------------------------------------------

library('igraph', lib.loc = '/mnt/NexData_file/home/wmx/envs/scRNA/lib/R/library/') # 加载2.0.3的igraph不然ordercells会报错
.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle2/lib/R/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(ggplot2)
  library(getopt)
  library(reshape2)
  library(cowplot)
  library(Matrix)
})

options(stringsAsFactors = FALSE, warn = -1)

# ============================================================
# 参数
# ============================================================

command <- matrix(c(
  'rds',                'r', 1, "character",
  'gene',               'g', 2, "character",
  'lambda',             'l', 1, "numeric",
  'outpath',            'o', 1, "character",
  'RootState',          'S', 2, "numeric",
  'Reverse',            'v', 2, "logical",
  'celllist',           'C', 2, "character",
  'clusters',           'c', 2, "character",
  'RandomCells',        'R', 2, "logical",
  'CellsPercent',       'p', 2, "numeric",
  'MinCellsPerCluster', 'M', 2, "numeric",
  'plotcols',           'P', 2, "character"
), byrow = TRUE, ncol = 4)

args <- getopt(command)

if (is.null(args$rds)) {
  stop("Please provide input rds file by -r")
}

if (is.null(args$outpath)) {
  stop("Please provide output path by -o")
}

if (is.null(args$lambda)) {
  args$lambda <- 0.001
}

if (is.null(args$gene)) {
  args$gene <- NULL
}

if (is.null(args$RootState)) {
  args$RootState <- NULL
}

if (is.null(args$Reverse)) {
  args$Reverse <- NULL
}

if (is.null(args$CellsPercent)) {
  args$CellsPercent <- 0.1
}

if (is.null(args$MinCellsPerCluster)) {
  args$MinCellsPerCluster <- 1000
}

if (is.null(args$RandomCells)) {
  args$RandomCells <- FALSE
}

rdsFile <- args$rds
geneFile <- args$gene
outpath <- args$outpath
lambda <- args$lambda
RootState <- args$RootState
Reverse <- args$Reverse
celllist <- args$celllist

cluster <- args$clusters
if (!is.null(cluster)) {
  cluster <- unlist(strsplit(cluster, ","))
  cluster <- trimws(cluster)
}

RandomCells <- args$RandomCells
cellsPercent <- args$CellsPercent
minCellsPerCluster <- args$MinCellsPerCluster

plotcols <- args$plotcols
if (is.null(plotcols)) {
  plotcols <- c("State", "Cluster", "Pseudotime")
} else {
  plotcols <- unlist(strsplit(plotcols, ","))
  plotcols <- trimws(plotcols)
}

message("========================================")
message("[Input rds] ", rdsFile)
message("[Gene file] ", ifelse(is.null(geneFile), "NULL", geneFile))
message("[Out path] ", outpath)
message("[Lambda] ", lambda)
message("[RootState] ", ifelse(is.null(RootState), "NULL", RootState))
message("[Reverse] ", ifelse(is.null(Reverse), "NULL", Reverse))
message("[Cluster use] ", ifelse(is.null(cluster), "NULL", paste(cluster, collapse = ",")))
message("[RandomCells] ", RandomCells)
message("[CellsPercent] ", cellsPercent)
message("[MinCellsPerCluster] ", minCellsPerCluster)
message("[Plot columns] ", paste(plotcols, collapse = ","))
message("========================================")

outdir <- file.path(outpath, "Pseudotime")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

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

get_real_colname <- function(seuCDS, col) {
  meta_cols <- colnames(seuCDS@phenoData@data)
  hit <- meta_cols[tolower(meta_cols) == tolower(col)]
  
  if (length(hit) > 0) {
    return(hit[1])
  } else {
    return(NA)
  }
}

get_nrow_facet <- function(x) {
  n <- length(unique(na.omit(x)))
  if (n <= 1) {
    return(1)
  }
  out <- floor(sqrt(n))
  if (out < 1) {
    out <- 1
  }
  return(out)
}

# ============================================================
# 随机抽样函数
# ============================================================

setNum <- function(number, cellsPercent, minCellsPerCluster) {
  selectNum <- round(number * cellsPercent)
  
  if (selectNum >= minCellsPerCluster) {
    finalNum <- selectNum
  } else {
    finalNum <- min(number, minCellsPerCluster)
  }
  
  return(finalNum)
}

randomSample <- function(cluster, size = NULL, cellsPercent, minCellsPerCluster) {
  samples <- c()
  ids <- unique(cluster)
  
  for (i in ids) {
    clusterInfo <- cluster[cluster == i]
    
    if (is.null(size)) {
      finalNum <- setNum(
        number = length(clusterInfo),
        cellsPercent = cellsPercent,
        minCellsPerCluster = minCellsPerCluster
      )
    } else {
      finalNum <- size
    }
    
    finalNum <- min(finalNum, length(clusterInfo))
    insample <- sample(names(clusterInfo), finalNum, replace = FALSE)
    samples <- c(samples, insample)
  }
  
  return(samples)
}

# ============================================================
# Seurat 预处理
# ============================================================

Psepre <- function(seuset,
                   celllist = NULL,
                   cluster = NULL,
                   RandomCells = FALSE,
                   cellsPercent = 0.1,
                   minCellsPerCluster = 1000) {
  
  message("[Psepre] Original cells: ", ncol(seuset))
  
  if (!is.null(celllist)) {
    message("[Psepre] Use celllist: ", celllist)
    cells_use <- read.table(celllist, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1]
    cells_use <- intersect(cells_use, colnames(seuset))
    message("[Psepre] Cells from celllist found in object: ", length(cells_use))
    seuset <- subset(seuset, cells = cells_use)
  }
  
  if (!is.null(cluster)) {
    message("[Psepre] Subset clusters: ", paste(cluster, collapse = ","))
    keep_cells <- names(seuset@active.ident)[seuset@active.ident %in% cluster]
    message("[Psepre] Cells in selected clusters: ", length(keep_cells))
    seuset <- subset(seuset, cells = keep_cells)
  }
  
  if (isTRUE(RandomCells)) {
    message("[Psepre] Random sampling cells")
    
    newcelllist <- randomSample(
      cluster = seuset@active.ident,
      cellsPercent = cellsPercent,
      minCellsPerCluster = minCellsPerCluster
    )
    
    seuset <- subset(seuset, cells = as.character(newcelllist))
    
    write.table(
      data.frame(Cell = colnames(seuset), Cluster = seuset@active.ident),
      file = "RandomCells.txt",
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    message("[Psepre] Cells after random sampling: ", ncol(seuset))
  }
  
  seuset <- AddMetaData(
    object = seuset,
    metadata = seuset@active.ident,
    col.name = "Cluster"
  )
  
  message("[Psepre] Final cells: ", ncol(seuset))
  
  return(seuset)
}

# ============================================================
# Seurat 转 Monocle2 CDS
# ============================================================

ImportSeurat4CDS <- function(otherCDS) {
  
  message("[ImportSeurat4CDS] Get count matrix")
  
  data <- GetAssayData(otherCDS, slot = "counts", assay = "RNA")
  
  if (class(data) == "data.frame") {
    data <- as(as.matrix(data), "sparseMatrix")
  }
  
  pd <- tryCatch(
    {
      new("AnnotatedDataFrame", data = otherCDS@meta.data)
    },
    error = function(e) {
      pData <- data.frame(
        cell_id = colnames(data),
        row.names = colnames(data)
      )
      message("[ImportSeurat4CDS] This Seurat object doesn't provide metadata")
      new("AnnotatedDataFrame", data = pData)
    }
  )
  
  if (length(setdiff(colnames(data), rownames(pd))) > 0) {
    data <- data[, rownames(pd)]
  }
  
  fd <- data.frame(
    gene_short_name = row.names(data),
    row.names = row.names(data)
  )
  fd <- new("AnnotatedDataFrame", data = fd)
  
  expressionFamily <- negbinomial.size()
  
  monocle_cds <- newCellDataSet(
    data,
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0,
    expressionFamily = expressionFamily
  )
  
  if ("Monocle" %in% names(otherCDS@misc)) {
    otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
    otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
    monocle_cds <- otherCDS@misc$Monocle
    mist_list <- otherCDS
  } else {
    mist_list <- otherCDS
  }
  
  monocle_cds@auxClusteringData$seurat <- mist_list
  
  return(monocle_cds)
}

# ============================================================
# 按指定列画 pseudotime 图
# ============================================================

PsePlotByCols <- function(seuCDS,
                          plotcols = c("State", "Cluster", "Pseudotime"),
                          do_facet = TRUE,
                          do_tree = TRUE,
                          do_density = TRUE) {
  
  meta <- seuCDS@phenoData@data
  
  message("========================================")
  message("[Plot] Available columns:")
  message(paste(colnames(meta), collapse = ", "))
  message("========================================")
  
  for (col in plotcols) {
    
    real_col <- get_real_colname(seuCDS, col)
    
    if (is.na(real_col)) {
      message("[Skip] Column not found: ", col)
      next
    }
    
    safe_col <- safe_name(real_col)
    
    message("----------------------------------------")
    message("[Plot] color_by = ", real_col)
    
    # 普通 trajectory
    gg <- plot_cell_trajectory(seuCDS, color_by = real_col) +
      plot_theme()
    
    save_plot_pdf_png(
      gg = gg,
      prefix = paste0("Pseudotime.", safe_col),
      width = 8,
      height = 6
    )
    
    # 判断是否为离散变量
    x <- meta[[real_col]]
    n_unique <- length(unique(na.omit(x)))
    
    is_discrete <- is.factor(x) ||
      is.character(x) ||
      real_col %in% c("State", "Cluster") ||
      n_unique <= 30
    
    # 分面图
    if (do_facet && is_discrete && n_unique > 1) {
      
      nrow_use <- get_nrow_facet(x)
      ncol_use <- ceiling(n_unique / nrow_use)
      
      gg_facet <- plot_cell_trajectory(seuCDS, color_by = real_col) +
        facet_wrap(as.formula(paste0("~", real_col)), nrow = nrow_use) +
        plot_theme()
      
      save_plot_pdf_png(
        gg = gg_facet,
        prefix = paste0("Pseudotime.", safe_col, ".Each"),
        width = max(8, 4 * ncol_use),
        height = max(5, 3 * nrow_use)
      )
    }
    
    # complex tree
    if (do_tree) {
      
      gg_tree <- plot_complex_cell_trajectory(
        seuCDS,
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
    
    # pseudotime density
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
        if (nrow_use < 1) {
          nrow_use <- 1
        }
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
  
  # 总体 pseudotime density
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
# 主函数
# ============================================================

mainfun <- function(seuset,
                    geneFile,
                    lambda,
                    RootState,
                    Reverse,
                    celllist,
                    cluster,
                    RandomCells,
                    cellsPercent,
                    minCellsPerCluster,
                    plotcols) {
  
  message("========================================")
  message("[Main] Read Seurat RDS")
  seuset <- readRDS(seuset)
  message("[Main] Seurat cells: ", ncol(seuset))
  message("[Main] Seurat genes: ", nrow(seuset))
  
  # 这里默认用 leiden 作为 active.ident
  # 如果你想改成 cell_type，就改成：seuset@active.ident <- seuset$cell_type
  if ("leiden" %in% colnames(seuset@meta.data)) {
    seuset@active.ident <- as.factor(seuset$leiden)
    message("[Main] Use active.ident = leiden")
  } else if ("cell_type" %in% colnames(seuset@meta.data)) {
    seuset@active.ident <- as.factor(seuset$cell_type)
    message("[Main] leiden not found, use active.ident = cell_type")
  } else {
    message("[Main] Use original active.ident")
  }
  
  seuset <- Psepre(
    seuset = seuset,
    celllist = celllist,
    cluster = cluster,
    RandomCells = RandomCells,
    cellsPercent = cellsPercent,
    minCellsPerCluster = minCellsPerCluster
  )
  
  message("[Main] Psepre done")
  
  # ordering genes
  if (!is.null(geneFile)) {
    message("[Main] Use gene file: ", geneFile)
    genelist <- unique(read.table(geneFile, sep = "\t", header = TRUE)[, 1])
  } else {
    message("[Main] Find variable features")
    seuset <- FindVariableFeatures(
      seuset,
      selection.method = "vst",
      nfeatures = 2000
    )
    genelist <- VariableFeatures(seuset)
  }
  
  genelist <- intersect(genelist, rownames(seuset))
  message("[Main] Ordering genes used: ", length(genelist))
  
  if (length(genelist) < 10) {
    stop("Too few ordering genes found in Seurat object")
  }
  
  lambda_use <- lambda * ncol(seuset)
  message("[Main] lambda after scaling: ", lambda_use)
  
  message("[Main] Build monocle CDS")
  
  seuCDS <- ImportSeurat4CDS(seuset) %>%
    setOrderingFilter(genelist) %>%
    estimateSizeFactors() %>%
    reduceDimension(
      reduction_method = "DDRTree",
      lambda = lambda_use
      # 如果细胞数量较少报错，可尝试打开：
      # auto_param_selection = FALSE,
      # ncenter = 2
    ) %>%
    orderCells()
  
  message("[Main] orderCells done")
  
  # 修改起点或反向
  if (!is.null(RootState) || !is.null(Reverse)) {
    message("[Main] Reorder cells by RootState / Reverse")
    seuCDS <- orderCells(
      seuCDS,
      root_state = RootState,
      reverse = Reverse
    )
  }
  
  # 保存统计表
  if ("Cluster" %in% colnames(seuCDS@phenoData@data) &&
      "State" %in% colnames(seuCDS@phenoData@data)) {
    
    stat_table <- table(seuCDS$Cluster, seuCDS$State) %>%
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
  
  # 保存细胞 pseudotime 信息
  seuCDS@phenoData@data %>%
    data.frame(CellName = row.names(seuCDS@phenoData@data)) %>%
    write.table(
      file = "Pseudotime.Summary.txt",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  
  # 按指定列画图
  PsePlotByCols(
    seuCDS = seuCDS,
    plotcols = plotcols,
    do_facet = TRUE,
    do_tree = TRUE,
    do_density = TRUE
  )
  
  # 保存 monocle 对象
  saveRDS(seuCDS, file = "Pseudotime.monocle.rds")
  
  message("========================================")
  message("[Main] All done")
  message("[Main] Results saved in: ", getwd())
  message("========================================")
}

# ============================================================
# 运行
# ============================================================

mainfun(
  seuset = rdsFile,
  geneFile = geneFile,
  lambda = lambda,
  RootState = RootState,
  Reverse = Reverse,
  celllist = celllist,
  cluster = cluster,
  RandomCells = RandomCells,
  cellsPercent = cellsPercent,
  minCellsPerCluster = minCellsPerCluster,
  plotcols = plotcols
)