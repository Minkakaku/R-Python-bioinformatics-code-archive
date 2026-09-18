.libPaths( "/home/szt/miniforge3/envs/rpy/lib/R/library")
options(stringsAsFactors = FALSE)
library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
setwd("/home/wyt/test/visiumhd/cellchat/")

cols <- c(brewer.pal(9,"Set1")[c(1:5,7,6,8,9)],brewer.pal(9,"Pastel1")[c(7,2,3,5,1)],brewer.pal(12,"Set3")[c(-2,-4,-6,-9,-11)],brewer.pal(8,"Set2")[c(3,4,7)],brewer.pal(8,"Pastel2")[8],brewer.pal(12,"Paired")[c(5,6,9,10,12)])

# gro <- "muti_group"
# infle = ":q!"
# path = "/home/wyt/test/protocol/cellchat_Spatial/visiumHD/"
# assay = "Spatial"
# seg_type = "visiumHD"
# spot.size = 8
# species = "mouse"
# cell_type = "celltype"
# sample_label = c("KD1W","KD2W","KD3W","KD4W")
# fov_label = c("KD1WL","KD2WL","KD3WL","KD4WL")



parser = argparse::ArgumentParser(description="celllchat_Spatial analysis")
parser$add_argument('--gro')
parser$add_argument('--infle')
parser$add_argument('--path')
parser$add_argument('--assay')
parser$add_argument('--seg_type')
parser$add_argument('--spot.size',type="integer")
parser$add_argument('--species')
parser$add_argument('--cell_type')
parser$add_argument('--sample_label',type="character",nargs="+")
parser$add_argument('--fov_label',type="character",nargs="+")

args = parser$parse_args()
gro = args$gro
fov = args$fov
infle  = args$infle
assay  = args$assay
path = args$path
seg_type = args$seg_type
spot.size = args$spot.size
species = args$species
cell_type = args$cell_type
sample_label =  args$sample_label
fov_label =  args$fov_label

dir.create(paste0(path,"/",gro))
setwd(paste0(path,"/",gro))
spatial <- readRDS(infle)
spatial <- subset(spatial,sample %in% sample_label)
for (i in fov_label){
  DefaultBoundary(spatial@images[[i]]) <- "centroid" 
  print(i)
}
spatial@meta.data[[cell_type]] <- as.factor(spatial@meta.data[[cell_type]])
Idents(spatial) <- cell_type
data.input = Seurat::GetAssayData(spatial, layer  = "data", assay = assay) 
meta = data.frame(labels = Idents(spatial), 
                  slices = spatial@meta.data$sample,
                  samples = as.factor(spatial@meta.data$sample)) 
meta$labels <- factor(meta$labels, levels = levels(Idents(spatial)))
meta$slices <- factor(meta$slices, levels = sample_label)

spatial_total.locs <- NULL
for( i in fov_label){
  spatial.locs = Seurat::GetTissueCoordinates(spatial[[i]])
  spatial.locs$cell <- NULL
  spatial_total.locs = rbind(spatial_total.locs,spatial.locs)
}


spatial.factors_total <- NULL
if (seg_type == "visiumHD") {
  for( i in fov_label){
    scale_factors = spatial@images[[i]]@scale.factors
    spot.size = spot.size 
    conversion.factor = spot.size/scale_factors$spot
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    spatial.factors_total <- rbind(spatial.factors_total, spatial.factors)
  }
} else if (seg_type == "xenium") {
  for( i in fov_label){
    spot.size = 10
    conversion.factor = 1
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    spatial.factors_total <- rbind(spatial.factors_total, spatial.factors)
  }
}

rownames(spatial.factors_total) <- sample_label

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial_total.locs, spatial.factors  = spatial.factors_total)

if (species == "mouse"){
  CellChatDB <- CellChatDB.mouse
} else if(species == "human"){
  CellChatDB <-  CellChatDB.human
}
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
print("identifyOverExpressedInteractions")
options(future.globals.maxSize = 16 * 1024^2 * 1024) 
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,contact.knn.k = 6)
##########################################################################
cat('computeCommunProb have done')
cellchat <- filterCommunication(cellchat, min.cells = 10)
#提取推断出的通信网络
cellchat <- computeCommunProbPathway(cellchat)
#在信号通路水平上推断细胞间通信
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat, paste0(gro,"_cellchat_netP.rds"))

pathways <- cellchat@netP$pathways
write.table(pathways,paste0(gro,"_pathways.txt"))

df.net <- subsetCommunication(cellchat)
write.csv(df.net,paste0(gro,"_df.net.csv"))

pdf(paste0(gro,"_interactions_Number_weights_strength.pdf"),12,6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

hei <- ceiling(length(unique(meta$label))/4)* 5.5

pdf(paste0(gro,"celltype_interactions_Number_weights_strength.pdf"),24,hei)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


library(NMF)
library(ggalluvial)

pdf(paste0(gro,"_out_income.pdf"),12,6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()

pdf(paste0(gro,"_outgoing.pdf"),12,6)
out_Patterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = out_Patterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

pdf(paste0(gro,"_income.pdf"),12,6)
input_Patterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = input_Patterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()


dir.create(paste0(path,"/",gro,"/aggregate"))

for (i in pathways) {
  for( j in sample_label){
    tryCatch({
      p <- netVisual_aggregate(cellchat, signaling = i, layout = "spatial", 
                               sample.use = j,edge.width.max = 2, vertex.size.max = 1, 
                               alpha.image = 0.1, vertex.label.cex = 3.5)
      p1 <- p + theme_bw() + xlab("") + ylab("") + 
        theme(axis.text = element_blank(), axis.ticks = element_blank())
      ggsave(p1, filename = paste0(path, "/", gro, "/aggregate/", 'cellchat_aggregate', i,"_",j, ".png"),
             width = 15, height = 15, units = "in", dpi = 300)
    },error = function(e) {
      message(paste("Error processing pathway", i,"fov",j, ":", e$message))
    }
    )
  }
}


options(future.globals.maxSize = 10 * 1024^3) 
pairs =as.character(unique(df.net$interaction_name))
dir.create(paste0(path,"/",gro,"/pair"))
for (j in pairs){
  int <- length(strsplit(j, "_")[[1]])
  for (z in sample_label){
    p2 <- spatialFeaturePlot(cellchat, pairLR.use = j, sample.use = z,point.size = 0.3, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
    for (i in 1:length(p2)){
      p2[[i]] <- p2[[i]] + theme_bw()+ xlab("")+ylab("")+ theme(axis.text = element_blank(),axis.ticks=element_blank()) + coord_fixed()  
    }
    ggsave(p2,filename = paste0(path,"/",gro,"/pair/",'Pairs_',j,"sample_",z,".png"),width = 16,height = 8,units = "in",dpi = 300)
    if (int == 2){
      p3 <- spatialFeaturePlot(cellchat, pairLR.use = j,sample.use = z, point.size = 0.3, do.binary = T, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
      p4 <- p3 +  theme_bw()+ xlab("")+ylab("")+ theme(axis.text = element_blank(),axis.ticks=element_blank()) + coord_fixed()
      ggsave(p4,filename = paste0(path,"/",gro,"/pair/",'Pairs_binary',j,"sample_",z,".png"),width = 16,height = 8,units = "in",dpi = 300)
    }
  }
}

cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "functional") #????uwotʵ??umap, ?????Դ???ֻ??ʹ?ò???umap.method = 'uwot'
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space

pdf(paste0(gro,"_functional.pdf"),12,6)
p1 <- netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
p2 <- netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
p1 + p2
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space

pdf(paste0(gro,"_structural.pdf"),12,6)
p1 <- netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
p2 <- netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
p1 + p2
dev.off()

saveRDS(cellchat, paste0(gro,"_cellchat_fun_str.rds"))


