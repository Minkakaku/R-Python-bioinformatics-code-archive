setwd('/mnt/NexData_file/home/zwl/Project/AE2404160288/require2/cellchat/')
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(argparser)
library(stringr)
library(jpeg)
rm(list=ls())
gc()


vertex_cal = function (object, signaling, slot.name = "netP"){
  centr0 <- slot(object, slot.name)$centr[signaling]
  centr0 <- centr0[[1]]
  mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
  mat <- t(mat)
  rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
  mat <- mat[c("outdeg","indeg"),]
  mat <- as.data.frame(t(mat))
  top3_indices <- head(order(mat$outdeg, decreasing = TRUE), 3)
}


top5_pathway =function(RData,outdir_A){
  
  load(RData)
  cellchat <- mergeCellChat(object.list,add.names = names(object.list))
  # 1. 识别所有的巨噬细胞亚群 (根据你的 metadata 或 cellchat 对象)
  # 假设 cellchat@idents 存储了细胞类型
  all_clusters <- levels(cellchat@idents[["joint"]])
  mac_clusters <- all_clusters[str_detect(all_clusters, "Macro")]
  

  # 2. 开始循环处理每个亚群
  for (mac_cell in mac_clusters) {
    outdir = paste0(outdir_A,mac_cell,'/')
    print(outdir)
    if(!dir.exists(outdir)) dir.create(outdir,recursive = T)
    # 定义两种角色：作为 Sender 和 作为 Target
    roles <- list(sender = "source", target = "target")
    
    for (role_name in names(roles)) {
      col_name <- roles[[role_name]]
      
      # 提取并过滤数据
      plot_df <- object.list %>% 
        map(subsetCommunication) %>%
        map_dfr(function(x) {
          # 动态过滤：当前巨噬细胞亚群在该角色位置
          x_filtered <- x %>% dplyr::filter(!!sym(col_name) == mac_cell)
          return(x_filtered)
        }, .id = "group")
      
      if (nrow(plot_df) == 0) next # 如果没数据则跳过
      
      # 筛选 Top 5 通路 (基于 prob 均值或最大值)
      top_pathways <- plot_df %>%
        group_by(pathway_name) %>%
        summarise(max_prob = max(prob)) %>%
        slice_max(max_prob, n = 5) %>%
        pull(pathway_name)
      
      # 绘制top5通路的热图等等
      # 组1
      for(k in str_c(top_pathways)){
        pairLR_top <- tryCatch(extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE), error=function(e) NULL)
        if (!is.null(pairLR_top) && nrow(pairLR_top) > 0) {
          pairLR_top <- pairLR_top[1,]
          vertex_receiver <- vertex_cal(object.list[[1]],k)
          pdf(paste0(outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
          gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 16, height = 2.5, font.size = 10)
          dev.off()
          
          pdf(paste0(outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 8)
          print(netAnalysis_contribution(object.list[[1]], signaling = k))
          dev.off()
          
          # pdf(paste0(outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
          # netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
          # dev.off()
        }else{
          print(paste0(k,' not express in ',names(object.list)[1]))
        }
      }
      # 组2
      for(k in str_c(top_pathways)){
        pairLR_top <- tryCatch(extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE), error=function(e) NULL)
        if (!is.null(pairLR_top) && nrow(pairLR_top) > 0) {
          pairLR_top <- pairLR_top[1,]
          vertex_receiver <- vertex_cal(object.list[[2]],k)
          pdf(paste0(outdir,k,"-",names(object.list)[2],"-centrality_heatmap.pdf"),width = 8,height= 3)
          gg <- netAnalysis_signalingRole_network(object.list[[2]], signaling = k, width = 16, height = 2.5, font.size = 10)
          dev.off()
          
          pdf(paste0(outdir,k,"-",names(object.list)[2],"-contribution_barplot.pdf"),width = 8,height= 8)
          print(netAnalysis_contribution(object.list[[2]], signaling = k))
          dev.off()
          
          # pdf(paste0(outdir,k,"-",names(object.list)[2],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
          # netVisual_individual(object.list[[2]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
          # dev.off()
        }else{
          print(paste0(k,' not express in ',names(object.list)[2]))
        }
      }
      
      
      # 最终绘图数据：仅保留 Top 5 通路相关的受配体对
      final_df <- plot_df %>%
        filter(pathway_name %in% top_pathways) %>%
        mutate(cell_pair = str_c(source, target, sep = " -> \n"))
      
      # 3. 绘图
      p <- ggplot(final_df, aes(x = group, y = interaction_name_2)) +
        geom_point(aes(color = prob, size = pval)) +
        facet_wrap(~ cell_pair, scales = "free_x", nrow = 1) +
        theme_classic() +
        scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) +
        scale_size(transform = "reverse", range = c(1, 4)) +
        labs(title = paste(mac_cell, "as", role_name),
             x = "Group", y = "Interaction (Top 5 Pathways)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8, face = "bold"))
      
      # 4. 动态保存
      # 文件名示例: Macrophages_as_sender_dotplot.pdf
      safe_name <- str_replace_all(mac_cell, " ", "_")
      file_base <- paste0(outdir, safe_name, "_as_", role_name, "_top5_dotplot")
      
      ggsave(paste0(file_base, ".png"), p, width = 17, height = 8, dpi = 150)
      ggsave(paste0(file_base, ".pdf"), p, width = 17, height = 8)
      
      print(paste("Finished:", mac_cell, "as", role_name))
    }
  }
}















files1 <- list.files('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/cellchat/myeloid_Fibro',
                     pattern = "tmp_object.RData", recursive = TRUE, full.names = TRUE)
files2 <- list.files('/mnt/NexData_file/home/zwl/Project/AE2404160288/require1/cellchat/Treg_myeloid',
                     pattern = "tmp_object.RData", recursive = TRUE, full.names = TRUE)
files3 <- list.files('/mnt/NexData_file/home/zwl/Project/AE2404160288/require2/cellchat/Macro',
                     pattern = "tmp_object.RData", recursive = TRUE, full.names = TRUE)
files4 <- list.files('/mnt/NexData_file/home/zwl/Project/AE2404160288/require2/cellchat/myeloid_Epi',
                     pattern = "tmp_object.RData", recursive = TRUE, full.names = TRUE)
files = c(files1,files2,files3,files4)
for (fi in files[7:12]) {
  outdir = paste0(dirname(fi),'/top_Path/')
  if(!dir.exists(outdir)) dir.create(outdir,recursive = T)
  top5_pathway(fi,outdir)
}



