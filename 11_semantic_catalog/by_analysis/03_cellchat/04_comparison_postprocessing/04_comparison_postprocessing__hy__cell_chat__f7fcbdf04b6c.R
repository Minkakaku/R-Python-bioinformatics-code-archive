rm(list = ls())

library(cowplot)

library(Seurat)

library(RColorBrewer)

library(CellChat)

library(ggplot2)

library(ggalluvial)

library(svglite)

library(patchwork)


library(CellChat)

options(stringsAsFactors = FALSE)

outdir = "/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/RCTD/compare/MPC/comp_1/CN7_CN8/cell_chat/CN2"

inrds = "/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/RCTD/compare/MPC/comp_1/CN7_CN8/cell_chat/Cn2_celltype_anno.rds"

species = 'human'
library(ComplexHeatmap)

setwd(outdir)

load("/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/RCTD/compare/MPC/comp_1/CN7_CN8/cell_chat/CN2/cellchat_object_list.RData")

load("/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/RCTD/compare/MPC/comp_1/CN7_CN8/cell_chat/CN2/cellchat_object.RData")

group = c("MPC","TNBC")


plot_save = function(filename,plt,dpi=300,w=7,h=7){

  ggsave(paste0(filename,'.pdf'),plt,height = h,width = w)

  

}

fig_save = function(filename,plt,dpi=300,w=7,h=7){

  pdf(paste0(filename,'.pdf'),height = h,width = w)

  draw(plt)

  dev.off()
}

recor_save = function(filename,plt,dpi=300,w=7,h=7){

  pdf(paste0(filename,'.pdf'),height = h,width = w)

  replayPlot(plt)

  dev.off()

}


for (i in 1:length(group)) {

  for (j in (1:length(group))[(1:length(group)) !=i]) {

    compare.list = c(i,j)

    

    # cellchat@netP$

    gg3 = netVisual_diffInteraction(cellchat,  comparison =compare.list,

                                    weight.scale = T,measure = "count")

    dev.off()

    recor_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.net"),gg3)

    gg4 = netVisual_diffInteraction(cellchat, comparison = compare.list,

                                    weight.scale = T, measure = "weight")

    recor_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.net"),gg4)

    gg5 = netVisual_heatmap(cellchat, comparison = compare.list)

    fig_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.heatmap"),gg5)

    gg6 = netVisual_heatmap(cellchat, measure = "weight", comparison = compare.list)

    fig_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.heatmap"),gg6)

    

    gg7 = rankNet(cellchat,

                  mode = "comparison",

                  measure = "weight",

                  sources.use = NULL,

                  targets.use = NULL, comparison =compare.list ,

                  stacked = T,

                  do.stat = TRUE)

    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.weight"),gg7)

    gg8 = rankNet(cellchat,

                  mode = "comparison",

                  measure = "weight",

                  sources.use = NULL,

                  targets.use = NULL, comparison = compare.list,

                  stacked = F,

                  do.stat = TRUE)

    plot_save(paste0(group[i],"VS",group[j],"_Information flow.weight"),gg8)

    

    gg9 = rankNet(cellchat,

                  mode = "comparison",

                  measure = "count",

                  sources.use = NULL,

                  targets.use = NULL,

                  comparison =compare.list,

                  stacked = T,

                  do.stat = TRUE)

    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.count"),gg9)

    gg10 = rankNet(cellchat,

                   mode = "comparison",

                   measure = "count",

                   sources.use = NULL,

                   targets.use = NULL, comparison = compare.list,

                   stacked = F,

                   do.stat = TRUE)

    plot_save(paste0(group[i],"VS",group[j],"_Information flow.count"),gg10)

    

    

    

    dif_pathway = as.character(gg8$data[c(1:3,(nrow(gg8$data)-3):nrow(gg8$data)),1])

    

    # 比较群体间的传入传出区别

    library(ComplexHeatmap)

    pathway.list = lapply(cellchat@netP, function(x){

      unique(x$pathways)

    }) %>% unlist() %>% unique()

    # group.max = 

    ht1 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "outgoing", 

                                            signaling = pathway.list,

                                            title = names(cellchat.objet.list.new)[i], 

                                            width = 10, height = length(pathway.list)/3)

    ht2 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "outgoing",

                                            signaling = pathway.list,

                                            title = names(cellchat.objet.list.new)[j],

                                            width = 10, height = length(pathway.list)/3)

    fig_save(paste0(group[i],"VS",group[j],'_Outgoing compare.heatmap'),plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),

             dpi = 300,w=10,h=length(pathway.list)/3)

    

    ht3 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "incoming", 

                                            signaling = pathway.list,

                                            title = names(cellchat.objet.list.new)[i], 

                                            width = 10, height = length(pathway.list)/3)

    ht4 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "incoming",

                                            signaling = pathway.list,

                                            title = names(cellchat.objet.list.new)[j], 

                                            width = 10, height = length(pathway.list)/3)

    fig_save(paste0(group[i],"VS",group[j],'_Incoming compare.heatmap'),plt = draw(ht3 + ht4, ht_gap = unit(0.5, "cm")),

             dpi = 300,w=10,h=length(pathway.list)/3)

    

    ht5 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "all", 

                                            signaling = pathway.list,

                                            title = names(cellchat.objet.list.new)[i], 

                                            width = 10, height = length(pathway.list)/3)

    ht6 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "all",

                                            signaling = pathway.list,

                                            title = names(cellchat.objet.list.new)[j], 

                                            width = 10, height = length(pathway.list)/3)

    fig_save(paste0(group[i],"VS",group[j],'_All compare.heatmap'),plt = draw(ht5 + ht6, ht_gap = unit(0.5, "cm")),

             dpi = 300,w=10,h=length(pathway.list)/3)

    

    

  }

}
