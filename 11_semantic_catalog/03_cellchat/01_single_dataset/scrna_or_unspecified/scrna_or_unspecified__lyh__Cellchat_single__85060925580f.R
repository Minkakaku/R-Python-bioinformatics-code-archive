#!/usr/bin/env Rscript

# ==============================================================================
# CellChat Analysis Pipeline (Refactored for Production)
# Engineer: Liu Yuhua
# Logic Source: Merged Cellchat_single.R & cellchat_analysis.R
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(NMF)
  library(ggalluvial)
  library(ComplexHeatmap)
  library(qs)
  library(future)
  library(reticulate)
})

# ==============================================================================
# 1. Command Line Arguments Setup
# ==============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input Seurat object path (.rds, .qs, or .RData)", metavar = "file"),
  make_option(c("-o", "--outdir"), type = "character", default = "./CellChat_Output", 
              help = "Output directory path [default= %default]", metavar = "path"),
  make_option(c("-s", "--species"), type = "character", default = "human", 
              help = "Species: 'human', 'mouse' or 'zebrafish' [default= %default]", metavar = "character"),
  make_option(c("-c", "--celltype_col"), type = "character", default = "celltype", 
              help = "Metadata column name for cell type annotation [default= %default]", metavar = "character"),
  make_option(c("--subset_col"), type = "character", default = NULL, 
              help = "Metadata column to subset data (e.g., 'Group') [optional]", metavar = "character"),
  make_option(c("--subset_group"), type = "character", default = NULL, 
              help = "Value in subset_col to keep (e.g., 'HR') [optional]", metavar = "character"),
  make_option(c("--python_path"), type = "character", default = "/usr/local/bin/python3", 
              help = "Path to python interpreter [default= %default]", metavar = "path"),
  make_option(c("--cores"), type = "integer", default = 4, 
              help = "Number of cores for parallel processing [default= %default]", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be supplied.", call. = FALSE)
}

# ==============================================================================
# 2. Configuration & Custom Functions
# ==============================================================================

# Configure Python
message(paste0(">>> Setting Python interpreter: ", opt$python_path))
use_python(opt$python_path, required = TRUE)

# Helper to convert PDF to PNG (ImageMagick)
run_convert <- function(pdf_file) {
  if (!file.exists(pdf_file)) return()
  png_file <- sub(".pdf$", ".png", pdf_file)
  cmd <- paste("convert -density 300", pdf_file, png_file)
  try(system(cmd), silent = TRUE)
}

# Custom NMF Function (Preserved from original script for visualization consistency)
myidentifyCommunicationPatterns <- function (object, slot.name = "netP", pattern = c("outgoing", 
    "incoming"), k = NULL, heatmap.show = TRUE, 
    color.use = NULL, color.heatmap = "Spectral", title.legend = "Contributions", 
    width = 4, height = 6, font.size = 8) 
{
    pattern <- match.arg(pattern)
    prob <- methods::slot(object, slot.name)$prob
    
    if (pattern == "outgoing") {
        data_sender <- apply(prob, c(1, 3), sum)
        data_sender = sweep(data_sender, 2L, apply(data_sender, 
            2, function(x) max(x, na.rm = TRUE)), "/", check.margin = FALSE)
        data0 = as.matrix(data_sender)
    } else if (pattern == "incoming") {
        data_receiver <- apply(prob, c(2, 3), sum)
        data_receiver = sweep(data_receiver, 2L, apply(data_receiver, 
            2, function(x) max(x, na.rm = TRUE)), "/", check.margin = FALSE)
        data0 = as.matrix(data_receiver)
    }
    
    options(warn = -1)
    data <- data0
    data <- data[rowSums(data) != 0, ]
    
    # Filter columns with low variance/counts
    if(nrow(data) >= 4){
        cols_to_remove <- apply(data, 2, function(x) length(unique(x)) < 4)
    } else {
        cols_to_remove <- apply(data, 2, function(x) all(x == 0 | x == 1))
    }
    
    cols_to_keep <- !cols_to_remove
    data_filt <- data[, cols_to_keep, drop = FALSE]
    data <- data_filt
    data <- na.omit(data)
    if(nrow(data) > 0) data <- data[rowSums(data) > 0, , drop=FALSE]
    
    # Safety checks
    if(ncol(data) == 0 || nrow(data) == 0) return(object)
    if (is.null(k)) stop("Please provide k!")
    
    # Critical: Adjust K if dimensions are too small (from analysis logic)
    if(ncol(data) < k || nrow(data) < k) {
        message(paste0("   [Warning] Data dim (",nrow(data),"x",ncol(data),") smaller than k=",k,". Adjusting k."))
        k <- min(ncol(data), nrow(data))
        if(k < 2) return(object) 
    }

    outs_NMF <- NMF::nmf(data, rank = k, method = "lee", seed = "nndsvd")
    W <- scaleMat(outs_NMF@fit@W, "r1")
    H <- scaleMat(outs_NMF@fit@H, "c1")
    colnames(W) <- paste0("Pattern ", seq(1, ncol(W)))
    rownames(H) <- paste0("Pattern ", seq(1, nrow(H)))
    
    if (heatmap.show) {
        net <- W
        if (is.null(color.use)) {
            color.use <- scPalette(length(rownames(net)))
        }
        color.heatmap = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, 
            name = color.heatmap))))(1000)
        df <- data.frame(group = rownames(net))
        rownames(df) <- rownames(net)
        cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
        row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned), 
            which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
            simple_anno_size = grid::unit(0.2, "cm"))
        
        ht1 = Heatmap(net, col = color.heatmap, na_col = "white", 
            name = "Contribution", left_annotation = row_annotation, 
            cluster_rows = T, cluster_columns = F, clustering_method_rows = "average", 
            row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
            column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                "cm"), height = unit(height, "cm"), show_heatmap_legend = F, 
            column_title = "Cell patterns", column_title_gp = gpar(fontsize = 10))
        
        net <- t(H)
        ht2 = Heatmap(net, col = color.heatmap, na_col = "white", 
            name = "Contribution", cluster_rows = T, cluster_columns = F, 
            clustering_method_rows = "average", row_names_side = "left", 
            row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
            column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                "cm"), height = unit(height, "cm"), column_title = "Communication patterns", 
            column_title_gp = gpar(fontsize = 10), heatmap_legend_param = list(title = title.legend, 
                title_gp = gpar(fontsize = 8, fontface = "plain"), 
                title_position = "leftcenter-rot", border = NA, 
                at = c(round(min(net, na.rm = T), digits = 1), 
                  round(max(net, na.rm = T), digits = 1)), legend_height = unit(20, 
                  "mm"), labels_gp = gpar(fontsize = 6), grid_width = unit(2, 
                  "mm")))
        
        draw(ht1 + ht2) 
    }
    
    data_W <- as.data.frame(as.table(W))
    colnames(data_W) <- c("CellGroup", "Pattern", "Contribution")
    data_H <- as.data.frame(as.table(H))
    colnames(data_H) <- c("Pattern", "Signaling", "Contribution")
    res.pattern = list(cell = data_W, signaling = data_H)
    methods::slot(object, slot.name)$pattern[[pattern]] <- list(data = data0, 
        pattern = res.pattern)
    return(object)
}

my_netAnalysis_signalingRole_heatmap <- function (object, signaling = NULL, pattern = c("outgoing", "incoming", 
    "all"), slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
    title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
    cluster.rows = FALSE, cluster.cols = FALSE) 
{
    pattern <- match.arg(pattern)
    if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
    }
    else if (pattern == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
    }
    else if (pattern == "all") {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
    }
    if (is.null(title)) {
        title <- paste0(legend.name, " signaling patterns")
    }
    else {
        title <- paste0(paste0(legend.name, " signaling patterns"), 
            " - ", title)
    }
    if (!is.null(signaling)) {
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
        idx <- match(rownames(mat1), signaling)
        mat[idx[!is.na(idx)], ] <- mat1
        dimnames(mat) <- list(signaling, colnames(mat1))
    }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
    if (is.null(color.use)) {
        color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
        name = color.heatmap))))(1000)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    names(color.use) <- colnames(mat)
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
        simple_anno_size = grid::unit(0.2, "cm"))
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
        border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
        show_annotation_name = FALSE)
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
        values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
            length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
        show_annotation_name = FALSE)
    if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
        legend.break <- max(mat, na.rm = T)
    }
    else {
        legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
            round(max(mat, na.rm = T), digits = 1))
    }
    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
        name = "Relative strength", bottom_annotation = col_annotation, 
        top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
        cluster_columns = cluster.rows, row_names_side = "left", 
        row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
        column_names_gp = gpar(fontsize = font.size), width = unit(width, 
            "cm"), height = unit(height, "cm"), column_title = title, 
        column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
            fontface = "plain"), title_position = "leftcenter-rot", 
            border = NA, at = legend.break, legend_height = unit(20, 
                "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                "mm")))
    return(ht1)
}


# ==============================================================================
# 3. Main Workflow
# ==============================================================================

# Create output dir
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

message(">>> Loading input data: ", opt$input)

# Load Data based on extension
if (grepl("\\.qs$", opt$input)) {
  seurat_obj <- qread(opt$input)
} else if (grepl("\\.rds$", opt$input)) {
  seurat_obj <- readRDS(opt$input)
} else {
  env <- new.env()
  load(opt$input, envir = env)
  obj_names <- ls(env)
  target_obj <- NULL
  for(nm in obj_names){
    if(inherits(get(nm, envir = env), "Seurat")){
      target_obj <- nm
      break
    }
  }
  if(is.null(target_obj)) stop("No Seurat object found in RData file.")
  seurat_obj <- get(target_obj, envir = env)
  rm(env)
}

#seurat_obj <- UpdateSeuratObject(seurat_obj)

# Set Identity
if (!(opt$celltype_col %in% colnames(seurat_obj@meta.data))) {
  stop(paste0("Column '", opt$celltype_col, "' not found in metadata."))
}
seurat_obj <- SetIdent(seurat_obj, value = opt$celltype_col)

# Subset Data
if (!is.null(opt$subset_col) && !is.null(opt$subset_group)) {
  message(paste0(">>> Subsetting data: ", opt$subset_col, " == ", opt$subset_group))
  if (!(opt$subset_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste0("Subset column '", opt$subset_col, "' not found."))
  }
  seurat_obj <- subset(seurat_obj, subset = !!sym(opt$subset_col) == opt$subset_group)
}

message(paste0(">>> Cells for analysis: ", ncol(seurat_obj)))

# Initialize CellChat
message(">>> Creating CellChat Object...")
cellchat <- createCellChat(object = seurat_obj, group.by = opt$celltype_col)

# Set Database
if (tolower(opt$species) == "human") {
  CellChatDB <- CellChatDB.human
  message(">>> Using Human Database")
} else if (tolower(opt$species) == "zebrafish") {
  CellChatDB <- CellChatDB.zebrafish
  message(">>> Using Zebrafish Database")
} else {
  CellChatDB <- CellChatDB.mouse
  message(">>> Using Mouse Database")
}

select_class <- unlist(strsplit("Secreted Signaling,ECM-Receptor,Cell-Cell Contact", ","))
CellChatDB.use <- subsetDB(CellChatDB, search = select_class, key = "annotation", non_protein = TRUE)
cellchat@DB <- CellChatDB.use

# Pre-processing
message(">>> Pre-processing and Inference...")
options(future.globals.maxSize = 100000 * 1024^2)
future::plan("multisession", workers = opt$cores)

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Inference (Updated: Added population.size=TRUE from analysis script)
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ==============================================================================
# 4. Basic Visualizations (Overview)
# ==============================================================================
message(">>> Generating Basic Plots...")
viz_dir <- file.path(opt$outdir, "Overview")
if(!dir.exists(viz_dir)) dir.create(viz_dir, recursive = TRUE)

groupSize <- as.numeric(table(cellchat@idents))

# 4.1 Global Interactions
f1 <- file.path(viz_dir, "network_interactions_count.pdf")
pdf(f1, width=6, height=7)
par(mar = c(1, 1, 1, 1), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
run_convert(f1)

f2 <- file.path(viz_dir, "network_interactions_weight.pdf")
pdf(f2, width=6, height=7)
par(mar = c(1, 1, 1, 1), xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
run_convert(f2)

# 4.2 Split Views (Interaction per Cell Type) 
# [Logic from cellchat_analysis.R]
mat <- cellchat@net$weight
cells_num <- length(levels(cellchat@idents))

# Precise dynamic layout logic from analysis script
if(cells_num <= 3){
    nrow=1; ncol=cells_num; width=10; height=5
} else if(cells_num > 3 && cells_num <= 4){
    nrow=2; ncol=2; width=8; height=6
} else if(cells_num > 4 && cells_num <= 6){
    nrow=2; ncol=3; width=12; height=6
} else if(cells_num > 6 && cells_num <= 9){
    nrow=3; ncol=3; width=12; height=8
} else if(cells_num > 9 && cells_num <= 10){
    nrow=2; ncol=5; width=22; height=6
} else if(cells_num > 10 && cells_num <= 12){
    nrow=3; ncol=4; width=16; height=10
} else if(cells_num > 12 && cells_num <= 15){
    nrow=3; ncol=5; width=22; height=10
} else if(cells_num > 15 && cells_num <= 16){
    nrow=4; ncol=4; width=16; height=12
} else if(cells_num > 16 && cells_num <= 20){
    nrow=4; ncol=5; width=22; height=12
} else {
    nrow=5; ncol=6; width=26; height=15
}

f3 <- file.path(viz_dir, "network_interactions_cells.pdf")
pdf(f3, width=width, height=height)
par(mfrow = c(nrow,ncol), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
run_convert(f3)

# ==============================================================================
# 5. Pathway Level Analysis & Heatmaps
# ==============================================================================
message(">>> Generating Pathway Plots...")
path_dir <- file.path(opt$outdir, "Pathways")
if(!dir.exists(path_dir)) dir.create(path_dir, recursive = TRUE)

pathways <- cellchat@netP$pathways
seq_num <- ceiling(cells_num/2) # for hierarchy plot splitting

for (i in 1:length(pathways)) {
  path <- pathways[i]
  
  # Hierarchy (Added vertex.receiver logic from analysis)
  f_hier <- file.path(path_dir, paste0(path, ".hierarchy.pdf"))
  pdf(f_hier, width=12, height=8) 
  plot_res <- tryCatch({
    vertex.receiver = seq(1, seq_num)
    netVisual_aggregate(cellchat, signaling = path, vertex.receiver = vertex.receiver, layout="hierarchy")
    TRUE
  }, error = function(e) FALSE)
  dev.off()
  if(plot_res) run_convert(f_hier) else if(file.exists(f_hier)) file.remove(f_hier)
  
  # Circle
  f_circ <- file.path(path_dir, paste0(path, ".circle.pdf"))
  pdf(f_circ, width=8, height=8)
  plot_res <- tryCatch({
    netVisual_aggregate(cellchat, signaling = path, layout="circle")
    TRUE
  }, error = function(e) FALSE)
  dev.off()
  if(plot_res) run_convert(f_circ) else if(file.exists(f_circ)) file.remove(f_circ)
  
  # Heatmap
  f_heat <- file.path(path_dir, paste0(path, ".heatmap.pdf"))
  pdf(f_heat, width=8, height=7)
  plot_res <- tryCatch({
    print(netVisual_heatmap(cellchat, signaling = path, color.heatmap = "Reds"))
    TRUE
  }, error = function(e) FALSE)
  dev.off()
  if(plot_res) run_convert(f_heat) else if(file.exists(f_heat)) file.remove(f_heat)
  
  # Contribution
  f_cont <- file.path(path_dir, paste0(path, ".contribution.pdf"))
  pdf(f_cont, width=8, height=6)
  plot_res <- tryCatch({
    print(netAnalysis_contribution(cellchat, signaling = path))
    TRUE
  }, error = function(e) FALSE)
  dev.off()
  if(plot_res) run_convert(f_cont) else if(file.exists(f_cont)) file.remove(f_cont)
}

# Bubble Plots (Optimized loop with interaction limit)
bubble_dir <- file.path(opt$outdir, "BubblePlots")
if(!dir.exists(bubble_dir)) dir.create(bubble_dir, recursive = TRUE)

for (i in 1:cells_num) {
  source_group <- levels(cellchat@idents)[i]
  f_bub <- file.path(bubble_dir, paste0(source_group, ".bubble.pdf"))
  
  tryCatch({
    # Logic to limit pairs if too many, preventing unreadable plots
    df.net_tmp <- subsetCommunication(cellchat, sources.use = source_group, targets.use = c(1:cells_num))
    
    # 增加非空判断，防止某种细胞完全没有互作导致报错
    if(nrow(df.net_tmp) == 0) {
        message(paste("Skipping", source_group, "- No interactions found."))
        next
    }
    
    len_inter <- length(unique(df.net_tmp$interaction_name))
    
    # Calc Dynamic Height: Base 4 + small increment per interaction
    bub_height <- 4 + (len_inter * 0.15)
    bub_height <- min(max(bub_height, 6), 25) # Clamp constraints
    
    if(len_inter > 50){
       # --- FIX START ---
       # 取前50个并将向量转换为数据框，指定列名为 interaction_name
       interaction_vec <- unique(df.net_tmp$interaction_name)[1:50]
       pairLR.use.df <- data.frame(interaction_name = interaction_vec)
       
       p <- netVisual_bubble(cellchat, sources.use = source_group, targets.use = c(1:cells_num), 
                             remove.isolate = FALSE, pairLR.use = pairLR.use.df)
       # --- FIX END ---
    } else {
       p <- netVisual_bubble(cellchat, sources.use = source_group, targets.use = c(1:cells_num), 
                             remove.isolate = FALSE)
    }
    
    ggsave(f_bub, p, height=bub_height, width=8)
    run_convert(f_bub)
  }, error = function(e) { message(paste("Bubble plot skipped for", source_group, ":", e$message)) })
}

# ==============================================================================
# 6. NMF Pattern Analysis (Improved K Selection with 2nd Derivative)
# ==============================================================================
message(">>> Running NMF Pattern Recognition...")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

nmf_dir <- file.path(opt$outdir, "Patterns")
if(!dir.exists(nmf_dir)) dir.create(nmf_dir, recursive = TRUE)

patterns <- c("outgoing", "incoming")

for (pat in patterns) {
  message(paste0("    Processing ", pat, " patterns..."))
  out_sub <- file.path(nmf_dir, pat)
  if(!dir.exists(out_sub)) dir.create(out_sub, recursive = TRUE)
  
  # Select K
  res_selectK <- selectK(cellchat, pattern = pat)
  f_sel <- file.path(out_sub, "selectK_plot.pdf")
  ggsave(f_sel, res_selectK$plot, width=6, height=4)
  run_convert(f_sel)
  
  # --- Improved Automated K selection (using Change Rate/Acceleration from analysis) ---
  df_k <- res_selectK$data
  
  # Cophenetic Logic
  coph_scores <- df_k$score[df_k$Measure == "Cophenetic"]
  k_vals <- df_k$k[df_k$Measure == "Cophenetic"]
  
  if(length(k_vals) > 2) {
      changes_co <- diff(coph_scores)
      change_rate_co <- c(NA, diff(changes_co))
      # +1 because diff reduces length
      k_co <- k_vals[which.min(change_rate_co) + 1] 
  } else {
      k_co <- 3
  }
  if(length(k_co) == 0) k_co <- 3
  
  # Silhouette Logic
  sil_scores <- df_k$score[df_k$Measure == "Silhouette"]
  if(length(k_vals) > 2) {
      changes_si <- diff(sil_scores)
      change_rate_si <- c(NA, diff(changes_si))
      k_si <- k_vals[which.min(change_rate_si) + 1]
  } else {
      k_si <- 3
  }
  if(length(k_si) == 0) k_si <- 3
  
  # Logic to pick best K
  selected_k <- min(k_co, k_si, na.rm = TRUE)
  
  # Limit K to dimensions of data (Prevent Rank > Dimension error)
  prob_mtx <- methods::slot(cellchat, "netP")$prob
  if(pat == "outgoing") data_check <- apply(prob_mtx, c(1, 3), sum)
  else data_check <- apply(prob_mtx, c(2, 3), sum)
  
  # Remove zero rows/cols roughly to estimate rank validity
  valid_dims <- min(dim(data_check))
  selected_k <- min(selected_k, valid_dims - 1)
  if(selected_k < 2) selected_k <- 2
  
  message(paste0("    Selected K for ", pat, ": ", selected_k, " (Co:", k_co, ", Si:", k_si, ")"))
  
  # Dynamic Height for NMF Pattern Plot
  nmf_plot_height <- max(6, 3 + (length(pathways) * 0.3))
  
  # Run Pattern Heatmap
  f_pat <- file.path(out_sub, paste0("CommunicationPatterns_", pat, ".pdf"))
  pdf(f_pat, width=10, height=nmf_plot_height)
  plot_res <- tryCatch({
      # Pass calculated height to custom function
      cellchat <- myidentifyCommunicationPatterns(cellchat, pattern = pat, k = selected_k, 
                                                  height=nmf_plot_height, heatmap.show = TRUE)
      TRUE
  }, error = function(e) { message(paste0("NMF Heatmap Error: ", e$message)); FALSE })
  dev.off()
  if(plot_res) run_convert(f_pat) else if(file.exists(f_pat)) file.remove(f_pat)
  
  # River Plot (converted logic from analysis: width=6, height=pdf_height)
  # Analysis: pdf_height <- (6+(length(pathways) * 0.3))/2.54
  river_height <- (6 + (length(pathways) * 0.3)) / 2.54
  river_height <- max(river_height, 4) # Minimum height check
  
  f_riv <- file.path(out_sub, paste0("River_", pat, ".pdf"))
  pdf(f_riv, width=6, height=river_height)
  plot_res <- tryCatch({
    netAnalysis_river(cellchat, pattern = pat)
    TRUE
  }, error = function(e) FALSE)
  dev.off()
  if(plot_res) run_convert(f_riv) else if(file.exists(f_riv)) file.remove(f_riv)
  
  # Dot Plot
  # Analysis: pdf_height <- (6+(cells_num * 0.5))/2.54
  # Analysis: pdf_width <- (10 +(length(pathways) * 0.3))/2.54
  dot_height <- (6 + (cells_num * 0.5)) / 2.54
  dot_width <- (10 + (length(pathways) * 0.3)) / 2.54
  
  f_dot <- file.path(out_sub, paste0("Dot_", pat, ".pdf"))
  pdf(f_dot, width=dot_width, height=dot_height)
  plot_res <- tryCatch({
    netAnalysis_dot(cellchat, pattern = pat)
    TRUE
  }, error = function(e) FALSE)
  dev.off()
  if(plot_res) run_convert(f_dot) else if(file.exists(f_dot)) file.remove(f_dot)
}

# ==============================================================================
# 7. Signaling Roles (2D Space)
# ==============================================================================
role_dir <- file.path(opt$outdir, "SignalingRoles")
if(!dir.exists(role_dir)) dir.create(role_dir, recursive = TRUE)

# 7.1 Scatter Plot
f_scat <- file.path(role_dir, "SignalingRole_Scatter.pdf")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(f_scat, gg1, width=6.5, height=5.5)
run_convert(f_scat)

# 7.2 Heatmaps with Dynamic Sizing (Logic from cellchat_analysis.R)
n_pathways_centr <- length(cellchat@netP$centr)
if(n_pathways_centr == 0) n_pathways_centr <- 10

# Formula directly from analysis script logic:
# pdf_width <- 10 + (cells_num * 0.5)
# pdf_height <- 8 + (length(pathways) * 0.3)
# The analysis script opens pdf with width=(pdf_width/2.54+2), height=(pdf_height/2.54+2)
raw_width <- 10 + (cells_num * 0.5)
raw_height <- 8 + (n_pathways_centr * 0.3)

calc_width  <- (raw_width / 2.54) + 2
calc_height <- (raw_height / 2.54) + 2

message(paste0(">>> Dynamic Role Heatmap Size: ", round(calc_width, 2), " x ", round(calc_height, 2), " inches"))

for (pat in c("outgoing", "incoming", "all")) {
  f_role <- file.path(role_dir, paste0("SignalingRole_", pat, "_heatmap.pdf"))
  
  pdf(f_role, width = calc_width, height = calc_height)
  
  plot_res <- tryCatch({
    # The custom function uses 'unit' internally. 
    # Analysis script passes raw_width/raw_height to the function's width/height params
    # while setting the PDF to calc_width/calc_height.
    
    draw(my_netAnalysis_signalingRole_heatmap(
      cellchat, 
      pattern = pat, 
      width = raw_width,   # passing the raw "cm-like" value
      height = raw_height, # passing the raw "cm-like" value
      font.size = 10 
    ))
    TRUE
  }, error = function(e) {
    message(paste("   [Error] Plotting", pat, "failed:", e$message))
    FALSE
  })
  dev.off()
  
  if(plot_res) run_convert(f_role) else if(file.exists(f_role)) file.remove(f_role)
}

# ==============================================================================
# 8. Manifold and Classification Learning Analysis (REQUIRES PYTHON)
# ==============================================================================
message(">>> Running Manifold and Classification Learning (Python)...")
group_dir <- file.path(opt$outdir, "SignalingGroups")
if(!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)

for(type in c("functional", "structural")) {
  message(paste0("    Processing ", type, " similarity..."))
  
  tryCatch({
    # Compute similarity and embedding
    cellchat <- computeNetSimilarity(cellchat, type = type)
    cellchat <- netEmbedding(cellchat, type = type, umap.method = 'uwot') # Added umap.method='uwot'
    cellchat <- netClustering(cellchat, type = type, do.plot = FALSE)
    
    # Visualizations
    p1 <- netVisual_embedding(cellchat, type = type, label.size = 3.5)
    f_embed <- file.path(group_dir, paste0("SignalingGroups_", type, ".pdf"))
    ggsave(f_embed, plot=p1, height=4, width=5)
    run_convert(f_embed)
    
    p2 <- netVisual_embeddingZoomIn(cellchat, type = type, nCol = 2)
    f_zoom <- file.path(group_dir, paste0("SignalingGroups_", type, "_split.pdf"))
    ggsave(f_zoom, plot=p2, height=5, width=6)
    run_convert(f_zoom)
    
  }, error = function(e) {
    message(paste0("Error in ", type, " analysis (check Python env): ", e$message))
  })
}

# ==============================================================================
# 9. Save Final Object & Tables
# ==============================================================================
save_path <- file.path(opt$outdir, "cellchat_object.rds")
message(">>> Saving CellChat object to: ", save_path)
saveRDS(cellchat, save_path)

message(">>> Exporting communication network table...")
df.net <- subsetCommunication(cellchat)

# 1. Save as TXT
txt_path <- file.path(opt$outdir, "Communication_info.txt")
write.table(df.net, txt_path, sep="\t", row.names=F, quote=F)

# 2. Save as CSV (Added requirement)
csv_path <- file.path(opt$outdir, "Communication_info.csv")
write.csv(df.net, csv_path, row.names = FALSE)

message(">>> Analysis Complete!")