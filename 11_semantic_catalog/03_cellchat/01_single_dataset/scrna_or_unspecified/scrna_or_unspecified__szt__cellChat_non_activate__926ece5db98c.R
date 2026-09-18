
#devtools::install_github("sqjin/CellChat")
library(patchwork)
library(CellChat)
options(stringsAsFactors = FALSE)
#############################################################################
harmony <- readRDS("non_activate.rds")

data <- harmony@assays$RNA@data
meta <- harmony@meta.data

cell.Normal <- rownames(meta)[meta$AUC_TLS_Group == "non_activate"]
#############################################################################

#第一部分：CellChat对象的数据输入、处理及初始化

#############################################################################
# Prepare input data for CelChat analysis
data.input <- data[,cell.Normal]
meta.input <- meta[cell.Normal,]

unique(meta.input$subType)
#-------------------------------------------
#创建CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "subType")

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#设置配体受体交互数据库
# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
table(CellChatDB$interaction$annotation)

# use Secreted Signaling
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

# set the used database in the object
#cellchat@DB <- CellChatDB.use
cellchat@DB <- CellChatDB
#-------------------------------------------
#预处理用于细胞通信分析的表达数据
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
#############################################################################

#第二部分：细胞通信网络的推断

#############################################################################
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
groupSize
cellchat <- filterCommunication(cellchat, min.cells = 10)

#，以轻松访问推断感兴趣的细胞-细胞通信
subsetCommunication(cellchat)
#subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

#在信号通路级别推断细胞-细胞通信
cellchat <- computeCommunProbPathway(cellchat)

#计算整合的细胞通信网络
cellchat <- aggregateNet(cellchat)

#可视化整合的细胞通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#由于细胞通信网络复杂，我们可以检查每个细胞组发送的信号。
mat <- cellchat@net$weight

pdf("non_activate.pdf", width = 16, height = 14)
par(mfrow = c(4,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
#############################################################################


#第三部分：细胞通信网络系统分析

#############################################################################
#计算和可视化网络中心分数
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pathways.show <- c("CXCL")
netAnalysis_signalingRole_scatter(cellchat)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

saveRDS(cellchat, file = "cellchat_non_activate.rds")

