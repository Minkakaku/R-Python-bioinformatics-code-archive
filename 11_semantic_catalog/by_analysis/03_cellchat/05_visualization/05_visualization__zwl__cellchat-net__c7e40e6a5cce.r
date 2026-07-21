setwd('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require5')

library(CellChat)
rm(list = ls())
gc()

for (mi in c('IR+NGF_IR','NGF_Control')) {
  load(paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require3/ALL/cellchat/Early/',mi,'/tmp_object.RData'))
  outdir ='/mnt/NexData_file/home/zwl/Project/FWSC20252451/require5/Early/'
  for (i in names(object.list)) {
    cellchat <- aggregateNet(object.list[[i]],
                             sources.use = c("Neurons", "OPC", "Microglia"),
                             targets.use = c("Neurons", "OPC", "Microglia"),)
    groupSize <- sapply(c("Neurons", "OPC", "Microglia"), function(ct) {
      sum(cellchat@meta$cell_type == ct)
    })
    
    pdf(paste0(outdir,i,"_subcelltype_interactions_netplot.pdf"),width = 12,height = 7)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
  }
  
  load(paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require3/ALL/cellchat/Late/',mi,'/tmp_object.RData'))
  outdir ='/mnt/NexData_file/home/zwl/Project/FWSC20252451/require5/Late/'
  for (i in names(object.list)) {
    cellchat <- aggregateNet(object.list[[i]],
                             sources.use = c("Neurons", "OPC", "Microglia"),
                             targets.use = c("Neurons", "OPC", "Microglia"),)
    groupSize <- sapply(c("Neurons", "OPC", "Microglia"), function(ct) {
      sum(cellchat@meta$cell_type == ct)
    })
    
    pdf(paste0(outdir,i,"_subcelltype_interactions_netplot.pdf"),width = 12,height = 7)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
  }
  
  
}


setwd('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require5')

library(CellChat)
rm(list = ls())
gc()

for (mi in c('IR+NGF_IR','NGF_Control')) {
  load(paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require3/Neuron_Mic/Early/',mi,'/tmp_object.RData'))
  outdir ='/mnt/NexData_file/home/zwl/Project/FWSC20252451/require5/Neuron_Mic/Early/'

  for (i in names(object.list)) {
    cellchat <- aggregateNet(object.list[[i]])
    groupSize <- table(cellchat@meta$cell_type)
    
    pdf(paste0(outdir,i,"_interactions_netplot.pdf"),width = 12,height = 7)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
  }
  
  load(paste0('/mnt/NexData_file/home/zwl/Project/FWSC20252451/require3/Neuron_Mic/Late/',mi,'/tmp_object.RData'))
  outdir ='/mnt/NexData_file/home/zwl/Project/FWSC20252451/require5/Neuron_Mic/Late/'
  for (i in names(object.list)) {
    cellchat <- aggregateNet(object.list[[i]])
    groupSize <- table(cellchat@meta$cell_type)
    
    pdf(paste0(outdir,i,"_interactions_netplot.pdf"),width = 12,height = 7)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
  }
  
  
}
