library(CellChat)
library(ggplot2)
library(tidyverse)

rm(list = ls())
setwd("/mnt/NexData_file/home/wmx/Project/AE2408071595-01/shouhou7/")
outdir = "/mnt/NexData_file/home/wmx/Project/AE2408071595-01/shouhou7/TNSF_plot"
dir.create(outdir,showWarnings = F)

# 加载cellchat结果
load("/mnt/NexData_file/home/wmx/Project/AE2408071595-01/shouhou7/cellchat/HLH_iMCD-TAFRO/tmp_object.RData")

data <- read.delim2( "/mnt/NexData_file/home/wmx/Project/AE2408071595-01/shouhou5/final_cellchat_interaction.txt",header = T)
genes <- c("TNFSF13B", "TNFSF12", "TNFSF10", "TNFSF13", "TNFSF14", "TNFSF8")
lr_network <- data  |> filter(receptor.symbol %in% genes | ligand.symbol %in% genes )

# =========================
# bubble group comparison 
# =========================
cellchat_merge <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

pathway_use <- Reduce(
  intersect,
  list(
    cellchat_merge@netP$HLH$pathways,
    cellchat_merge@netP$`iMCD-TAFRO`$pathways,
    unique(lr_network$pathway_name)
  )
)

pdf(
  file.path(outdir, "group_comparison_bubble.pdf"),
  width = 15,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_merge,
    signaling = pathway_use,
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()

png(
  file.path(outdir, "group_comparison_bubble.png"),
  width = 15, height = 8, units = "in", res = 300
)

print(
  netVisual_bubble(
    cellchat_merge,
    signaling = pathway_use,
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()

# 挑选nonclassical 发出，其它为target的图
pdf(
  file.path(outdir, "NCM-Other_group_comparison_bubble.pdf"),
  width = 8,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_merge,
    signaling = pathway_use,
    sources.use = "Nonclassical Monocyte",
    targets.use = setdiff(rownames((object.list$HLH@net$weight)),"Nonclassical Monocyte"),
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()

png(
  file.path(outdir, "NCM-Other_group_comparison_bubble.png"),
  width = 8, height = 8, units = "in", res = 300
)

print(
  netVisual_bubble(
    cellchat_merge,
    signaling = pathway_use,
    sources.use = "Nonclassical Monocyte",
    targets.use = setdiff(rownames((object.list$HLH@net$weight)),"Nonclassical Monocyte"),
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()

# 挑选发出为其它，接受为non classical
pdf(
  file.path(outdir, "Other-NCM_group_comparison_bubble.pdf"),
  width = 8,
  height = 8
)

print(
  netVisual_bubble(
    cellchat_merge,
    signaling = pathway_use,
    targets.use = "Nonclassical Monocyte",
    sources.use  = setdiff(rownames((object.list$HLH@net$weight)),"Nonclassical Monocyte"),
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()

png(
  file.path(outdir, "Other-NCM_group_comparison_bubble.png"),
  width = 8, height = 8, units = "in", res = 300
)

print(
  netVisual_bubble(
    cellchat_merge,
    signaling = pathway_use,
    targets.use = "Nonclassical Monocyte",
    sources.use  = setdiff(rownames((object.list$HLH@net$weight)),"Nonclassical Monocyte"),
    comparison = c(1, 2),
    remove.isolate = FALSE,
    angle.x = 90
  )
)

dev.off()

# =========================
# CCL ,IFN-II 通路的层次图绘制
# =========================

Reduce(
  intersect,
  list(
    cellchat_merge@netP$HLH$pathways,
    cellchat_merge@netP$`iMCD-TAFRO`$pathways,
    c("CCL","IFN-II")
  )
)

outdir = "/mnt/NexData_file/home/wmx/Project/AE2408071595-01/shouhou7/CCL_IFN_plot"
dir.create(outdir,showWarnings = F)

reorder_cellchat <- function(cellchat, new_order) {
  
  old_order <- levels(cellchat@idents)
  
  # 检查细胞类型是否一致
  if (!setequal(old_order, new_order)) {
    cat("old_order:\n")
    print(old_order)
    cat("new_order:\n")
    print(new_order)
    stop("new_order 和 cellchat@idents 的细胞类型不完全一致")
  }
  
  # 修改 idents 顺序
  cellchat@idents <- factor(cellchat@idents, levels = new_order)
  
  # 修改 @net 中的矩阵或三维数组
  for (slot_name in names(cellchat@net)) {
    x <- cellchat@net[[slot_name]]
    
    if (is.matrix(x) || is.data.frame(x)) {
      if (!is.null(rownames(x)) && !is.null(colnames(x))) {
        if (all(new_order %in% rownames(x)) && all(new_order %in% colnames(x))) {
          cellchat@net[[slot_name]] <- x[new_order, new_order, drop = FALSE]
        }
      }
    }
    
    if (length(dim(x)) == 3) {
      dn <- dimnames(x)
      if (!is.null(dn[[1]]) && !is.null(dn[[2]])) {
        if (all(new_order %in% dn[[1]]) && all(new_order %in% dn[[2]])) {
          cellchat@net[[slot_name]] <- x[new_order, new_order, , drop = FALSE]
        }
      }
    }
  }
  
  # 修改 @netP 中的三维数组
  for (slot_name in names(cellchat@netP)) {
    x <- cellchat@netP[[slot_name]]
    
    if (length(dim(x)) == 3) {
      dn <- dimnames(x)
      if (!is.null(dn[[1]]) && !is.null(dn[[2]])) {
        if (all(new_order %in% dn[[1]]) && all(new_order %in% dn[[2]])) {
          cellchat@netP[[slot_name]] <- x[new_order, new_order, , drop = FALSE]
        }
      }
    }
  }
  
  return(cellchat)
}
object.list <- lapply(object.list, function(x) {
  
  celltypes <- rownames(x@net$weight)
  
  new_order <- c(
    "Nonclassical Monocyte",
    setdiff(celltypes, "Nonclassical Monocyte")
  )
  
  x <- reorder_cellchat(x, new_order = new_order)
  
  return(x)
})


for (nm in names(object.list)) {
  cellchat <- object.list[[nm]]
  for(sig in c("CCL","IFN-II")){
    
    # 层次图
    png(
      file.path(outdir, paste0(nm,"_",sig,"_hierarchy.png")),
      width = 14, height = 9, units = "in", res = 300
    )
    
    print(netVisual_aggregate(object =cellchat, 
                              signaling = sig,  
                              vertex.receiver = c(1,2,3,4,5),
                              layout = "hierarchy"))
    dev.off()
    
    pdf(
      file.path(outdir, paste0(nm,"_",sig,"_hierarchy.pdf")),
      width = 14,
      height = 9
    )
    print(netVisual_aggregate(object =cellchat, 
                              signaling = sig,  
                              vertex.receiver = c(1,2,3,4,5),
                              layout = "hierarchy"))
    dev.off()
    
    # heatmap图
    cellchat <- netAnalysis_computeCentrality(object = cellchat, slot.name = "netP") 
    pdf(
      file.path(outdir, paste0(nm,"_",sig,"_heatmap.pdf")),
      width = 10,
      height = 7
    )
    netAnalysis_signalingRole_network(cellchat,
                                      signaling = sig,
                                      width = 10, height = 3.5,
                                      font.size = 10)
    dev.off()
    png(
      file.path(outdir, paste0(nm,"_",sig,"_heatmap.png")),
      width = 10, height = 7, units = "in", res = 300
    )
    netAnalysis_signalingRole_network(cellchat,
                                      signaling = sig,
                                      width = 10, height = 3.5,
                                      font.size = 10)
    dev.off()
  }
}


