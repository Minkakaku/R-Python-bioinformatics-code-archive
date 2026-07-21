#!/usr/bin/env Rscript

# ==============================================================================
# CellChat Comparative Analysis Pipeline (Final Fix: PNG Conversion & TXT Export)
# Engineer: Liu Yuhua
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(qs)
  library(future)
  library(reticulate)
})

# ==============================================================================
# 1. Arguments & Setup
# ==============================================================================

option_list <- list(
  make_option(c("--input_ctrl"), type = "character", default = NULL, help = "Control RDS"),
  make_option(c("--input_treat"), type = "character", default = NULL, help = "Treat RDS"),
  make_option(c("--name_ctrl"), type = "character", default = "Control", help = "Control Name"),
  make_option(c("--name_treat"), type = "character", default = "Treat", help = "Treat Name"),
  make_option(c("-o", "--outdir"), type = "character", default = "./CellChat_Compare_Output", help = "Output Dir"),
  make_option(c("--python_path"), type = "character", default = "/mnt/NexData_file/home/lyh/Software/conda_env/lyh/bin/python", help = "Python Path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_ctrl) || is.null(opt$input_treat)) {
  stop("Input files missing.", call. = FALSE)
}

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# Python setup
use_python(opt$python_path, required = TRUE)

# Helper: Convert PDF to PNG
run_convert <- function(pdf_file) {
  if (!file.exists(pdf_file)) return()
  png_file <- sub(".pdf$", ".png", pdf_file)
  try(system(paste("convert -density 300", pdf_file, png_file)), silent = TRUE)
}

write_status_note <- function(note_file, lines) {
  writeLines(enc2utf8(lines), con = note_file)
}

# Helper: Custom Heatmap
my_netAnalysis_signalingRole_heatmap <- function (object, signaling = NULL, pattern = c("outgoing", "incoming", "all"), slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE) {
    pattern <- match.arg(pattern)
    if (length(slot(object, slot.name)$centr) == 0) stop("Please run `netAnalysis_computeCentrality` first!")
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") { mat <- t(outgoing); legend.name <- "Outgoing" }
    else if (pattern == "incoming") { mat <- t(incoming); legend.name <- "Incoming" }
    else { mat <- t(outgoing + incoming); legend.name <- "Overall" }
    
    if (is.null(title)) title <- paste0(legend.name, " signaling patterns")
    else title <- paste0(paste0(legend.name, " signaling patterns"), " - ", title)
    
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
    
    if (is.null(color.use)) color.use <- scPalette(length(colnames(mat)))
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(1000)
    
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = setNames(color.use, colnames(mat))), which = "column", show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = grid::unit(0.2, "cm"))
    
    pSum <- rowSums(mat.ori); pSum[is.na(pSum)] <- 0
    pSum_log <- -1/log(pSum); pSum_log[is.infinite(pSum_log) | pSum_log < 0] <- 0
    ha1 = rowAnnotation(Strength = anno_barplot(pSum_log, border = FALSE), show_annotation_name = FALSE)
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE, gp = gpar(fill = color.use, col = color.use)), show_annotation_name = FALSE)

    Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength", bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, cluster_columns = cluster.rows, row_names_side = "left", row_names_gp = gpar(fontsize = font.size), column_names_gp = gpar(fontsize = font.size), width = unit(width, "cm"), height = unit(height, "cm"), column_title = title, column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90)
}

# ==============================================================================
# 2. Load & Merge
# ==============================================================================
message(">>> Loading datasets...")
cellchat_ctrl <- readRDS(opt$input_ctrl)
cellchat_treat <- readRDS(opt$input_treat)

if (exists("updateCellChat")) {
  message(">>> Updating CellChat objects...")
  tryCatch({
    cellchat_ctrl <- updateCellChat(cellchat_ctrl)
    cellchat_treat <- updateCellChat(cellchat_treat)
  }, error = function(e) {
    message(">>> Update failed or not needed, proceeding...")
  })
}

# Harmonize
common_idents <- union(levels(cellchat_ctrl@idents), levels(cellchat_treat@idents))
cellchat_ctrl <- liftCellChat(cellchat_ctrl, common_idents)
cellchat_treat <- liftCellChat(cellchat_treat, common_idents)

object.list <- list(cellchat_ctrl, cellchat_treat)
names(object.list) <- c(opt$name_ctrl, opt$name_treat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat, file.path(opt$outdir, "cellchat_merged.rds"))

# ==============================================================================
# 3. Overview Plots (Fixed PNG conversions)
# ==============================================================================
message(">>> Generating Overview Plots...")

# Barplots
f_bar <- file.path(opt$outdir, "compareInteractions.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
ggsave(f_bar, gg1 + gg2, height=3, width=7)
run_convert(f_bar) # ADDED

# Diff Circle
f_diff_count <- file.path(opt$outdir, "diffInteraction_count.pdf")
pdf(f_diff_count, height=7.5, width=6)
par(xpd=TRUE); netVisual_diffInteraction(cellchat, weight.scale = T); dev.off()
run_convert(f_diff_count) # ADDED

f_diff_weight <- file.path(opt$outdir, "diffInteraction_weight.pdf")
pdf(f_diff_weight, height=7.5, width=6)
par(xpd=TRUE); netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight"); dev.off()
run_convert(f_diff_weight) # ADDED

# Diff Heatmap
f_diff_hm <- file.path(opt$outdir, "diffInteraction_count_heatmap.pdf")
pdf(f_diff_hm, height=6, width=6)
netVisual_heatmap(cellchat); dev.off()
run_convert(f_diff_hm) # ADDED

# Group Circle
f_group <- file.path(opt$outdir, "interactions_group.pdf")
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(f_group, width=12, height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interactions - ", names(object.list)[i]))
}
dev.off()
run_convert(f_group) # ADDED

# Scatter
f_scatter <- file.path(opt$outdir, "signalingRole_scatter.pdf")
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
ggsave(f_scatter, patchwork::wrap_plots(plots = gg), height=3, width=7)
run_convert(f_scatter) # ADDED

# ==============================================================================
# 4. Signaling Role Heatmaps
# ==============================================================================
message(">>> Generating Comparative Heatmaps...")
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
cells_num <- length(common_idents)

heatmap_w_cm <- 10 + (cells_num * 0.5)
heatmap_h_cm <- 6 + (length(pathway.union) * 0.3)

for (pat in c("outgoing", "incoming", "all")){
  f_hm <- file.path(opt$outdir, paste0("signalingRole_heatmap_", pat, ".pdf"))
  pdf(f_hm, width = (heatmap_w_cm/2.54 + 2), height = (heatmap_h_cm/2.54 + 2))
  tryCatch({
    ht1 = my_netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = pat, signaling = pathway.union, title = names(object.list)[1], width = 5, height = heatmap_h_cm)
    ht2 = my_netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = pat, signaling = pathway.union, title = names(object.list)[2], width = 5, height = heatmap_h_cm)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  }, error = function(e) message(paste("Error heatmap:", pat, e$message)))
  dev.off()
  run_convert(f_hm)
}

# ==============================================================================
# 5. Differential Bubble Plots (Fixed TXT export)
# ==============================================================================
message(">>> Generating Differential Bubble Plots...")

dataset.name <- names(object.list)
pos.dataset <- dataset.name[2]
features.name <- paste0(pos.dataset, ".merged")

deg_dir <- file.path(opt$outdir, "DEG_Network")
if(!dir.exists(deg_dir)) dir.create(deg_dir)
deg_status_file <- file.path(deg_dir, "deg_status.txt")

deg_ok <- TRUE
deg_error <- NULL
net.up <- NULL
net.down <- NULL

tryCatch({
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05, min.cells = 10)
  net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
  net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset, ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(cellchat, net = net, datasets = dataset.name[1], ligand.logFC = -0.05, receptor.logFC = NULL)
}, error = function(e) {
  deg_ok <<- FALSE
  deg_error <<- e$message
})

if (!deg_ok) {
  message(paste0(">>> DEG/Bubble module skipped: ", deg_error))
  write_status_note(deg_status_file, c(
    "Status: skipped",
    paste0("Reason: ", deg_error)
  ))
} else {
  if (is.null(net.up)) net.up <- data.frame()
  if (is.null(net.down)) net.down <- data.frame()

  # Save ALL tables as both CSV and TXT
  write.csv(net.up, file.path(deg_dir, "netMappingDEG_UP.csv"), row.names=F)
  write.table(net.up, file.path(deg_dir, "netMappingDEG_UP.txt"), sep="\t", row.names=F, quote=F)

  write.csv(net.down, file.path(deg_dir, "netMappingDEG_DOWN.csv"), row.names=F)
  write.table(net.down, file.path(deg_dir, "netMappingDEG_DOWN.txt"), sep="\t", row.names=F, quote=F)

  write_status_note(deg_status_file, c(
    "Status: completed",
    paste0("UP rows: ", nrow(net.up)),
    paste0("DOWN rows: ", nrow(net.down))
  ))
}

cell_type <- common_idents

dir_bub <- file.path(opt$outdir, "Diff_BubblePlots")
if(!dir.exists(dir_bub)) dir.create(dir_bub)

if (deg_ok && (nrow(net.up) > 0 || nrow(net.down) > 0) && length(cell_type) >= 5) {
  message(">>> Cell types >= 5, generating individual bubble plots...")
  
  for (i in 1:length(cell_type)) {
    ct <- cell_type[i]
    p_Inc <- tryCatch({
      netVisual_bubble(cellchat, sources.use = ct, targets.use = cell_type, comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased in ", dataset.name[2]), angle.x = 45, remove.isolate = TRUE)
    }, error = function(e) NULL)
    
    p_Dec <- tryCatch({
      netVisual_bubble(cellchat, sources.use = ct, targets.use = cell_type, comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased in ", dataset.name[2]), angle.x = 45, remove.isolate = TRUE)
    }, error = function(e) NULL)
    
    p_final <- NULL
    if (!is.null(p_Inc) || !is.null(p_Dec)) {
       if(!is.null(p_Inc) && !is.null(p_Dec)) {
          p_final <- p_Inc + p_Dec
          w_calc <- (length(unique(p_Inc$data$source.target)) * 0.2 + 4) * 2
          h_calc <- max(length(unique(p_Inc$data$interaction_name_2)), length(unique(p_Dec$data$interaction_name_2))) * 0.08 + 3
       } else if (!is.null(p_Inc)) {
          p_final <- p_Inc
          w_calc <- length(unique(p_Inc$data$source.target)) * 0.2 + 4
          h_calc <- length(unique(p_Inc$data$interaction_name_2)) * 0.08 + 3
       } else {
          p_final <- p_Dec
          w_calc <- length(unique(p_Dec$data$source.target)) * 0.2 + 4
          h_calc <- length(unique(p_Dec$data$interaction_name_2)) * 0.08 + 3
       }
       
       f_name <- file.path(dir_bub, paste0(ct, "_Up-Down_Regulation.pdf"))
       ggsave(f_name, p_final, width=w_calc, height=h_calc, limitsize=FALSE)
       run_convert(f_name)
    }
  }
} else if (deg_ok && (nrow(net.up) > 0 || nrow(net.down) > 0)) {
  message(">>> Cell types < 5, generating combined bubble plot...")
  p_Inc <- tryCatch({
      netVisual_bubble(cellchat, sources.use = cell_type, targets.use = cell_type, comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased in ", dataset.name[2]), angle.x = 45, remove.isolate = TRUE)
    }, error = function(e) NULL)
  p_Dec <- tryCatch({
      netVisual_bubble(cellchat, sources.use = cell_type, targets.use = cell_type, comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased in ", dataset.name[2]), angle.x = 45, remove.isolate = TRUE)
    }, error = function(e) NULL)
  
  if (!is.null(p_Inc) || !is.null(p_Dec)) {
       if(!is.null(p_Inc) && !is.null(p_Dec)) {
          p_final <- p_Inc + p_Dec
          w_calc <- (length(unique(p_Inc$data$source.target)) * 0.2 + 4) * 2
          h_calc <- max(length(unique(p_Inc$data$interaction_name_2)), length(unique(p_Dec$data$interaction_name_2))) * 0.08 + 3
       } else if (!is.null(p_Inc)) {
          p_final <- p_Inc
          w_calc <- length(unique(p_Inc$data$source.target)) * 0.2 + 4
          h_calc <- length(unique(p_Inc$data$interaction_name_2)) * 0.08 + 3
       } else {
          p_final <- p_Dec
          w_calc <- length(unique(p_Dec$data$source.target)) * 0.2 + 4
          h_calc <- length(unique(p_Dec$data$interaction_name_2)) * 0.08 + 3
       }
       f_name <- file.path(dir_bub, "AllCells_Up-Down_Regulation.pdf")
       ggsave(f_name, p_final, width=w_calc, height=h_calc, limitsize=FALSE)
       run_convert(f_name)
  }
} else {
  message(">>> Differential bubble plots skipped: no significant differential signaling interactions found.")
}

# ==============================================================================
# 6. Shared Pathway Details
# ==============================================================================
message(">>> Generating Shared Pathway Plots...")
inter_path <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
dir_path <- file.path(opt$outdir, "Shared_Pathways")
if(!dir.exists(dir_path)) dir.create(dir_path)

if(length(inter_path) > 0){
  for (path in inter_path){
    # Circle
    f_circ <- file.path(dir_path, paste0(path, "_circle.pdf"))
    tryCatch({
        weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path)
        pdf(f_circ, width=10, height=5)
        par(mfrow = c(1,2), xpd=TRUE)
        for (i in 1:length(object.list)) netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]))
        dev.off()
        run_convert(f_circ)
    }, error = function(e) NULL)
    
    # Heatmap
    f_heat <- file.path(dir_path, paste0(path, "_heatmap.pdf"))
    tryCatch({
      pdf(f_heat, width=10, height=5)
      ht <- list()
      for (i in 1:length(object.list)) ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = path, color.heatmap = "Reds", title.name = paste(path, names(object.list)[i]))
      ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
      dev.off()
      run_convert(f_heat)
    }, error = function(e) NULL)
  }
}

message(">>> Comparison Analysis Complete!")
