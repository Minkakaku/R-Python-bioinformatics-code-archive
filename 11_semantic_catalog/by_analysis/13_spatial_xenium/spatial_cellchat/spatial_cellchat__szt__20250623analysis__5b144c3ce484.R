all_data@meta.data %>% head()
###大类
all_data@meta.data$annotation %>% unique()
#B cells  pDCs  T cells Endothelial  Macrophage Fibroblasts        
#Plasma  MergeDC Epithelial Smooth muscle cells
###亚类
all_data@meta.data$annotation2 %>% unique()
#"B cells" "pDCs" "CD4_Tfh"  "CD8_Tem" 
#"Naive_T"  "Endothelial"   "T_like_Macrophage"   "Fibroblasts" 
#"Tfh"  "B_like_Macrophage"   "cycling macrophage"  "Treg"
#"CD8_naive""proliferating T" "M2"   "Plasma"  
#"M1"   "MergeDC"  "Epithelial"   "CD4_naive" 
#"DC_like_Macrophage"  "Epithelial_like_Macrophage" "CD8_MAIT" "Smooth muscle cells"  
#"Fibroblast_like_Macrophage" "Mast_like_Macrophage"  "Th2"
table(all_data@meta.data$annotation ,all_data@meta.data$celltype)
###分组
all_data@meta.data$group.y %>% unique()
# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","Ni_LN","Nmi_LN","N1_LN"
###样本
all_data@meta.data$sample %>% unique()

na_count <- sum(is.na(all_data@meta.data$annotation2))
cat("annotation2 列中 NA 的数量：", na_count, "\n")

p <- DimPlot(all_data,group.by = "annotation2")
ggsave(
  p,
  filename = file.path("cell.png"),
  height = 6,
  width  = 12
)

# 
# # 转换为 data.table
# library(data.table)
# dt <- as.data.table(all_data@meta.data)
# 
# # 创建 annotation3 并填充
# dt[, annotation3 := annotation2]
# dt[is.na(annotation2), annotation3 := annotation]
# 
# # 转换回原始数据结构
# all_data@meta.data <- as.data.frame(dt)
# 
# # 验证
# na_before <- sum(is.na(all_data@meta.data$annotation2))
# na_after <- sum(is.na(all_data@meta.data$annotation3))
# cat("填充前后 NA 数量变化:", na_before, "→", na_after, "\n")
# 
# all_data@meta.data$annotation3



all_data <- readRDS('/home/rstudio/R_OUTPUT/all_data.rds')
Mac <- readRDS('/home/rstudio/Macrophage_subtype/macrophage_sub.rds')
T_cell <- readRDS('/home/rstudio/T_subtype/All_T_subtype.rds')
Epi <- readRDS('/home/rstudio/Epi_subcluster/Epi.rds')
#映射——add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(T_cell[['subtype']] %>% rownames_to_column('cells'), 
                 Mac[['subtype']] %>% rownames_to_column('cells'))

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2= ifelse((as.character(annotation) %in% c('T cells','Macrophage')),
                          as.character(subtype), 
                          as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
all_data$annotation %>% table()
all_data$annotation2 %>% table()

all_data <-  JoinLayers(all_data,layers = 'data')
# 增加内存限制至15 GiB（根据实际需求调整）
options(future.globals.maxSize = 15 * 1024^3)  # 15 GiB = 15 × 1024^3 字节

# 重新执行标准化（确保已设置并行计划）
library(future)
plan("multisession", workers = 4)  # 根据CPU核心数调整workers
all_data <- NormalizeData(all_data)


# 1. 转换为字符型（确保类型一致）
all_data@meta.data$annotation2 <- as.charater(all_data@meta.data$annotation2)
# 2. 替换NA为"un"
all_data@meta.data$annotation2[is.na(all_data@meta.data$annotation2)] <- "un"
# # 3. 若需要恢复为因子类型（需包含"un"水平）
# if (is.factor(all_data@meta.data$annotation2)) {
#   all_data@meta.data$annotation2 <- factor(
#     all_data@meta.data$annotation2,
#     levels = unique(c(levels(all_data@meta.data$annotation2), "un"))
#   )
# }
# 4. 验证
na_count_after <- sum(is.na(all_data@meta.data$annotation2))
cat("替换后NA数量：", na_count_after, "\n")

all_data@meta.data$annotation2 <- as.factor(all_data@meta.data$annotation2)

saveRDS(all_data, file = file.path("R_OUTPUT", "all_data_normalize.rds"))

all_data<-readRDS('/home/rstudio/R_OUTPUT/all_data_normalize.rds')

#4.2Cellchat----
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
library(patchwork)

###Part I: Data input & processing and initialization of CellChat object
# load("/Output/seurat_result/cluster_result/Filter_Basicinfo.rData")
# if(!("seurat_clusters" %in% colnames(pbmc_filt@meta.data))){
#   pbmc_filt<-UpdateSeuratObject(pbmc_filt)
#   pbmc_filt[["seurat_clusters"]]=pbmc_filt@meta.data$res.0.8
# }else{
#   print(head(pbmc_filt@meta.data))
# }

#细胞注释
# cell_annot<-read.table("../cell_annot.txt",header=T,sep="\t",quote="",check.names=F)
# cell.names <- rownames(pbmc_filt@meta.data)
# cell.types <- cell_annot$celltype[match(cell.names, cell_annot$cells)]  
# pbmc_filt[["celltype"]] <- cell.types

# #筛选分析细胞
# cells=read.table(file="select_cells.txt",header=F,row.names=1,sep="\t",quote="",check.names=F)
# pbmc_filt<-subset(pbmc_filt,cells=rownames(cells))

all_data@meta.data$annotation2 %>% unique()

###肿瘤-基质----
alldataCC1 <- subset(all_data, 
                     subset = annotation %in% c("Epithelial", "Endothelial", "Fibroblasts"))

alldataCC1 <- alldataCC1 <- JoinLayers(alldataCC1,layers = 'data')


cellchat <- createCellChat(object = alldataCC1@assays$RNA@layers$data,meta=alldataCC1@meta.data,group.by='annotation')
cellchat <- setIdent(cellchat, ident.use = "annotation",levels=unique(alldataCC1@meta.data$annotation))
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human

p<-showDatabaseCategory(CellChatDB)
ggsave("Cellchat_result/CellChatDB.mmu.png",width=10,height=5,dpi=300)
ggsave("Cellchat_result/CellChatDB.mmu.pdf",width=10,height=5)

print(colnames(CellChatDB$interaction))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

###Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
###信息提取
df.net <- subsetCommunication(cellchat)
write.table(df.net,"Cellchat_result/Communication_info.txt",sep="\t",quote=F,row.names=F)

##Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
head(cellchat@net$count)
head(cellchat@net$weight)

pdf(file="Cellchat_result/network_interactions_count.pdf",height=8,width=7)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf(file="Cellchat_result/network_interactions_weight.pdf",height=8,width=7)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
cells_num<-length(levels(cellchat@idents))
if(cells_num<=3){
  nrow=1
  ncol=cells_num
  width=10
  height=5
}else if(cells_num>3 && cells_num<=4){
  nrow=2
  ncol=2
  width=8
  height=6
}else if(cells_num>4 && cells_num<=6){
  nrow=2
  ncol=3
  width=12
  height=6
}else if(cells_num>6 && cells_num<=9){
  nrow=3
  ncol=3
  width=12
  height=8
}else if(cells_num>9 && cells_num<=10){
  nrow=2
  ncol=5
  width=22
  height=6
}else if(cells_num>10 && cells_num<=12){
  nrow=3
  ncol=4
  width=16
  height=10
}else if(cells_num>12 && cells_num<=15){
  nrow=3
  ncol=5
  width=22
  height=10
}else if(cells_num>15 && cells_num<=16){
  nrow=4
  ncol=4
  width=16
  height=12
}else if(cells_num>16 && cells_num<=20){
  nrow=4
  ncol=5
  width=22
  height=12
}else{
  nrow=5
  ncol=6
  width=26
  height=15
}
pdf("Cellchat_result/network_interactions_cells.pdf",width=width,height=height)
par(mfrow = c(nrow,ncol), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

####Part III: Visualization of cell-cell communication network
###可视化
print(cellchat@netP$pathways)
for (i in 1:length(cellchat@netP$pathways)){
  pathways.show <- cellchat@netP$pathways[i]
  file=paste("Cellchat_result/pathway/",cellchat@netP$pathways[i],".pdf",sep="")
  vertex.receiver = seq(1,4)
  pdf(file=file,width=18,height=8)
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  dev.off()
  file2=paste("Cellchat_result/pathway/",cellchat@netP$pathways[i],".circle.pdf",sep="")
  pdf(file=file2,width=10,height=10)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
  file3=paste("Cellchat_result/pathway/",cellchat@netP$pathways[i],".chord.pdf",sep="")
  pdf(file=file3,width=10,height=10)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
  ####heatmap循环运行pdf得到的图片是空的
  #file4=paste("Cellchat_result/pathway/",cellchat@netP$pathways[i],".heatmap.pdf",sep="")
  #pdf(file=file4,width=12,height=10)
  #netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  #dev.off()
  ####contribution循环运行得到的图片是空的
  #file5=paste("Cellchat_result/pathway/",cellchat@netP$pathways[i],".contribution.pdf",sep="")
  #pdf(file=file5,width=12,height=10)
  #netAnalysis_contribution(cellchat, signaling = pathways.show)
  #dev.off()
}

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
for (i in 1:length(cellchat@netP$pathways)){
  pathways.show <- cellchat@netP$pathways[i]
  file=paste("Cellchat_result/pathway/",cellchat@netP$pathways[i],".signalingRole_network.pdf",sep="")
  pdf(file=file,width=8,height=6)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  dev.off()
}

gg1 <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(file="Cellchat_result/signalingRole_scatter.pdf",plot=gg1,width=10,height=10)

pdf(file="Cellchat_result/signalingRole_heatmap_outgoing.pdf",width=8,height=6)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()
pdf(file="Cellchat_result/signalingRole_heatmap_incoming.pdf",width=8,height=6)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

library(NMF)
library(ggalluvial)
pdf(file="Cellchat_result/selectK_outgoing.pdf",width=8,height=6)
selectK(cellchat, pattern = "outgoing")
dev.off()

#设置聚类的个数
nPatterns = 3
pdf(file="Cellchat_result/CommunicationPatterns_outgoing.pdf",width=12,height=6)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)  ##必须赋值到cellchat，否则后续无法绘图
dev.off()

pdf(file="Cellchat_result/netAnalysis_river_outgoing.pdf",width=12,height=6)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf(file="Cellchat_result/netAnalysis_dot_outgoing.pdf",width=10,height=6)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

pdf("Cellchat_result/selectK_incoming.pdf",width=8,height=6)
selectK(cellchat, pattern = "incoming")
dev.off()
#设置聚类的个数
nPatterns = 4
pdf("Cellchat_result/CommunicationPatterns_incoming.pdf",width=12,height=6)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()

pdf("Cellchat_result/netAnalysis_river_incoming.pdf",width=12,height=6)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

pdf("Cellchat_result/netAnalysis_dot_incoming.pdf",width=10,height=6)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

save(cellchat,file="Cellchat_result/cellchat.rData")


###肿瘤-免疫----
alldataCC2
###免疫-免疫----
alldataCC3



#4.3拟时序----
library(Seurat)
library(tidyverse)
library(monocle)
source("monocle fix bugs.R")

###4.3.1T细胞----
monocle_out_dir1 <- "/home/rstudio/monocle_result/T_cell" 
dir.create(monocle_out_dir1) 
sc_tcell <- T_cell
sc_tcell <- FindVariableFeatures(sc_tcell, nfeatures = 3000)
Idents(sc_tcell)<-"subtype"

D <- DimPlot(sc_tcell, reduction = 'harmonyumap') 
ggsave("monocle_result/T_cell/dimplot.png", width = 8, height = 5)

VariableFeatures(sc_tcell)

# subset
set.seed(1234)
meta_sc_tcell_subset <-
  sc_tcell[['subtype']] %>%
  rownames_to_column("cell_id") %>%
  group_by(subtype) %>%
  slice_sample(n = 5000)
# meta_sc_epi2_subset <- sc_epi2[['annotation_full']] %>% rownames_to_column("cell_id") 

# cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(sc_tcell), "counts"), "sparseMatrix")
p_data_monocle <- sc_tcell[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sc_tcell),
  row.names       = rownames(sc_tcell)
)

cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_sc_tcell_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_sc_tcell_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# cds <- clusterCells(cds)

# ##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

#cds <- reduceDimension(cds)
cds<-reduceDimension(cds,max_components=2,method ='DDRTree')
#undebug(orderCells)
cds <- orderCells(cds)

# plot ----------
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds,color_by="subtype", size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="group.y", size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)

cds <- orderCells(cds, root_state = 4)

saveRDS(cds, file = file.path(monocle_out_dir1, "cds.rds"))


# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds, color_by="Pseudotime", size=0.5,show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by="subtype",cell_size = 0.5, size=0.5,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by="group.y", size=0.5,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=0.5,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "group.y", cell_size = 0.5)#+ scale_color_manual(values = col_plan_2)



pdf(file.path(monocle_out_dir1, "cell trajectory.pdf"), width = 8, height = 6)

plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="group.y", cell_size =0.1,show_backbone=TRUE) +
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="subtype", cell_size =0.1, show_backbone=TRUE) + 
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir1, "cell trajectory - annotation - split by group.pdf"), width = 15, height = 12)
local({
  Biobase::pData(cds)$group.y <- all_data[['group.y']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$group)]
  
  plot_cell_trajectory(cds, color_by = "subtype",cell_size = 0.1,show_backbone=TRUE) +
    #scale_color_manual(values = col_plan_1) +
    facet_wrap(~ group.y) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()


pdf(file.path(monocle_out_dir1, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "subtype") +
  #scale_color_manual(values = col_plan_1) +
  facet_wrap(~ subtype)
dev.off()


pdf(file.path(monocle_out_dir1, "barplot tof cell props group by State - group.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, group.y) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = group.y))
dev.off()

pdf(file.path(monocle_out_dir1, "barplot of cell props group by State - annotation.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, subtype ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = subtype ))
dev.off()


pdf(file.path(monocle_out_dir1, "Dimplot of T_cell - group and cell annotation.pdf"), width = 8, height = 5)
DimPlot(
  sc_tcell,
  group.by = "subtype"
)

DimPlot(
  sc_tcell,
  group.by = "group.y"
)
dev.off()



# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
  dplyr::select(subtype, State) %>%
  group_by(subtype, State) %>%
  mutate(subtype = as.character(subtype)) %>%
  # mutate(annotation_full = factor(annotation_full, levels = c(
  #   'S1 fibroblasts', 
  #   'Peri-cryptal S2 fibroblasts',
  #   'Immediate early-response (IER) fibroblasts',
  #   'S3 fibroblasts',
  #   'Myofibroblasts'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype, y = freq)) +
  geom_col(aes(fill = State)) +
  #scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_out_dir1, "cell stat barplot.pdf"), height = 4, width = 8)




# de analysis -----------
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 10, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir1, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(monocle_out_dir1, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

#undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir1, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)  
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)


# BEAM analysis ------------------
#undebug(BEAM)
gene_for_beam <- VariableFeatures(sc_tcell)
BEAM_res <- BEAM(cds[gene_for_beam,], branch_point = 1, cores = 10) 

saveRDS(BEAM_res, file.path(monocle_out_dir, "BEAM_res.rds"))

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res, file.path(monocle_out_dir, "BEAM_res.csv"), row.names = F)

# plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
#                             branch_point = 1, #绘制的是哪个分支
#                             num_clusters = 4, #分成几个cluster，根据需要调整
#                             cores = 5,
#                             use_gene_short_name = T,
#                             show_rownames = T) #有632个gene，太多了

#选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()




library(Seurat)
library(tidyverse)
library(monocle)
source("monocle fix bugs.R")
###4.3.2Macrophage----
monocle_out_dir2 <- "/home/rstudio/monocle_result/Mac" 
dir.create(monocle_out_dir2) 
sc_Mac <- Mac
sc_Mac <- FindVariableFeatures(sc_Mac, nfeatures = 3000)
Idents(sc_Mac)<-"subtype"

D <- DimPlot(sc_Mac, reduction = 'harmonyumap') 
ggsave("monocle_result/Mac/dimplot.png", width = 8, height = 5)

VariableFeatures(sc_Mac)

# subset
set.seed(1234)
meta_sc_Mac_subset <-
  sc_Mac[['subtype']] %>%
  filter(subtype %in% c("M1","M2")) %>% 
  rownames_to_column("cell_id") %>%
  group_by(subtype) %>%
  slice_sample(n = 5000)
# meta_sc_epi2_subset <- sc_epi2[['annotation_full']] %>% rownames_to_column("cell_id") 

# cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(sc_Mac), "counts"), "sparseMatrix")
p_data_monocle <- sc_Mac[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sc_Mac),
  row.names       = rownames(sc_Mac)
)

cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_sc_Mac_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_sc_Mac_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# cds <- clusterCells(cds)

# ##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

#cds <- reduceDimension(cds)
cds<-reduceDimension(cds,max_components=2,method ='DDRTree')
#undebug(orderCells)
cds <- orderCells(cds)

# plot ----------
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds,color_by="subtype", size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="group.y", size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)

cds <- orderCells(cds, root_state = 4)

saveRDS(cds, file = file.path(monocle_out_dir2, "cds.rds"))

#所有巨噬:cds_allMac.rds;M1、M2:cds.rds
cds <- readRDS(file = file.path(monocle_out_dir2, "cds.rds"))

# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds, color_by="Pseudotime", size=0.5,show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by="subtype",cell_size = 0.5, size=0.5,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by="group.y", size=0.5,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=0.5,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "group.y", cell_size = 0.5)#+ scale_color_manual(values = col_plan_2)



pdf(file.path(monocle_out_dir2, "cell trajectory.pdf"), width = 8, height = 6)

plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="group.y", cell_size =0.1,show_backbone=TRUE) +
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="subtype", cell_size =0.1, show_backbone=TRUE) + 
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir2, "cell trajectory - annotation - split by group.pdf"), width = 15, height = 12)
local({
  Biobase::pData(cds)$group.y <- all_data[['group.y']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$group)]
  
  plot_cell_trajectory(cds, color_by = "subtype",cell_size = 0.1,show_backbone=TRUE) +
    #scale_color_manual(values = col_plan_1) +
    facet_wrap(~ group.y) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()


pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "subtype") +
  #scale_color_manual(values = col_plan_1) +
  facet_wrap(~ subtype)
dev.off()


pdf(file.path(monocle_out_dir2, "barplot tof cell props group by State - group.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, group.y) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = group.y))
dev.off()

pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, subtype ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = subtype ))
dev.off()


pdf(file.path(monocle_out_dir2, "Dimplot of T_cell - group and cell annotation.pdf"), width = 8, height = 5)
DimPlot(
  sc_Mac,
  group.by = "subtype"
)

DimPlot(
  sc_Mac,
  group.by = "group.y"
)
dev.off()



# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
  dplyr::select(subtype, State) %>%
  group_by(subtype, State) %>%
  mutate(subtype = as.character(subtype)) %>%
  # mutate(annotation_full = factor(annotation_full, levels = c(
  #   'S1 fibroblasts', 
  #   'Peri-cryptal S2 fibroblasts',
  #   'Immediate early-response (IER) fibroblasts',
  #   'S3 fibroblasts',
  #   'Myofibroblasts'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = subtype, y = freq)) +
  geom_col(aes(fill = State)) +
  #scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_out_dir2, "cell stat barplot.pdf"), height = 4, width = 8)




# de analysis -----------
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 10, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir2, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(monocle_out_dir2, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

#undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir2, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)  
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)


  
  
  
###4.3.3Epithelial----
monocle_out_dir3 <- "/home/rstudio/monocle_result/Epi" 
dir.create(monocle_out_dir3) 
sc_Epi <-Epi
sc_Epi <- FindVariableFeatures(sc_Epi, nfeatures = 3000)
Idents(sc_Epi)<-"RNA_snn_res.0.1"

D <- DimPlot(sc_Epi, reduction = 'harmonyumap') 
ggsave("monocle_result/Epi/dimplot.png", width = 8, height = 5)

VariableFeatures(sc_Epi)

# subset
set.seed(1234)
meta_sc_Epi_subset <-
  sc_Epi[['RNA_snn_res.0.1']] %>%
  rownames_to_column("cell_id") %>%
  group_by(RNA_snn_res.0.1) %>%
  slice_sample(n = 5000)
# meta_sc_epi2_subset <- sc_epi2[['annotation_full']] %>% rownames_to_column("cell_id") 

# cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(sc_Epi), "counts"), "sparseMatrix")
p_data_monocle <- sc_Epi[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sc_Epi),
  row.names       = rownames(sc_Epi)
)

cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_sc_Epi_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_sc_Epi_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# cds <- clusterCells(cds)

# ##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

#cds <- reduceDimension(cds)
cds<-reduceDimension(cds,max_components=2,method ='DDRTree')
#undebug(orderCells)
cds <- orderCells(cds)

# plot ----------
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds,color_by="RNA_snn_res.0.1", size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="group.y", size=1,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)

#cds <- orderCells(cds, root_state = 4)

saveRDS(cds, file = file.path(monocle_out_dir3, "cds.rds"))


# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds, color_by="Pseudotime", size=0.5,show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by="RNA_snn_res.0.1",cell_size = 0.5, size=0.5,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by="group.y", size=0.5,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=0.5,show_backbone=TRUE) #+ scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "RNA_snn_res.0.1")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "group.y", cell_size = 0.5)#+ scale_color_manual(values = col_plan_2)



pdf(file.path(monocle_out_dir3, "cell trajectory.pdf"), width = 8, height = 6)

plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="group.y", cell_size =0.1,show_backbone=TRUE) +
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="RNA_snn_res.0.1", cell_size =0.1, show_backbone=TRUE) + 
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  #scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir3, "cell trajectory - annotation - split by group.pdf"), width = 15, height = 12)
local({
  Biobase::pData(cds)$group.y <- all_data[['group.y']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$group)]
  
  plot_cell_trajectory(cds, color_by = "RNA_snn_res.0.1",cell_size = 0.1,show_backbone=TRUE) +
    #scale_color_manual(values = col_plan_1) +
    facet_wrap(~ group.y) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()


pdf(file.path(monocle_out_dir3, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "RNA_snn_res.0.1") +
  #scale_color_manual(values = col_plan_1) +
  facet_wrap(~ RNA_snn_res.0.1)
dev.off()


pdf(file.path(monocle_out_dir3, "barplot tof cell props group by State - group.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, group.y) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = group.y))
dev.off()

pdf(file.path(monocle_out_dir3, "barplot of cell props group by State - annotation.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, RNA_snn_res.0.1 ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = RNA_snn_res.0.1 ))
dev.off()


pdf(file.path(monocle_out_dir3, "Dimplot of - group and cell annotation.pdf"), width = 8, height = 5)
DimPlot(
  sc_Epi,
  group.by = "RNA_snn_res.0.1"
)

DimPlot(
  sc_Epi,
  group.by = "group.y"
)
dev.off()



# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
  dplyr::select(RNA_snn_res.0.1, State) %>%
  group_by(RNA_snn_res.0.1, State) %>%
  mutate(RNA_snn_res.0.1 = as.character(RNA_snn_res.0.1)) %>%
  # mutate(annotation_full = factor(annotation_full, levels = c(
  #   'S1 fibroblasts', 
  #   'Peri-cryptal S2 fibroblasts',
  #   'Immediate early-response (IER) fibroblasts',
  #   'S3 fibroblasts',
  #   'Myofibroblasts'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.1, y = freq)) +
  geom_col(aes(fill = State)) +
  #scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_out_dir3, "cell stat barplot.pdf"), height = 4, width = 8)




# de analysis -----------
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 10, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir3, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(monocle_out_dir3, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

#undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir3, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)  
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)













#4.4inferCNV----
# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


# 1.所有上皮细胞和down sample部分免疫细胞，跑inferCNV分析，
# 识别出正常肝细胞（肝母细胞/胆管上皮细胞？），恶性上皮细胞。
# 恶性上皮细胞跑cNMF分析，识别出上皮细胞程序。

library(infercnv)

###  epi subcluster --------------
Epi_out_dir <- "Epi_subcluster"
dir.create(Epi_out_dir)


Epi <- all_data %>% JoinLayers() %>% .[, .$annotation %in% "Epithelial"]


pdf(file.path(Epi_out_dir, 'nCount_Spatial_008um plot.pdf'), height = 6, width = 8)
VlnPlot(Epi, features = "nCount_RNA", group.by = 'project', pt.size = 0)
dev.off()


Epi[['RNA']] <- split(Epi[['RNA']], Epi$project, layers = "counts")


# cellcycle
Epi <- 
  Epi %>%
  CellCycleScoring(
    s.features = cc.genes.updated.2019$s.genes %>% str_to_title,
    g2m.features = cc.genes.updated.2019$g2m.genes %>% str_to_title
  )

# preprocess
Epi <- 
  Epi %>% 
  NormalizeData(scale.factor = floor(mean(all_data$nCount_RNA)))%>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(
    vars.to.regress = c(
      'nCount_RNA'
    )
  )



# sc_hep <- sc_hep %>% ScaleData(vars.to.regress = c('nCount_RNA', 'S.Score', 'G2M.Score'))

VariableFeaturePlot(Epi)

# VariableFeatures(sc_hep) <- VariableFeatures(sc_hep) %>% setdiff(c(ig_genes, tcr_genes, mt_genes, ribo_genes))



## pca -----
Epi <- Epi %>% RunPCA()

pdf(file.path(Epi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epi, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_Epi <- 
  Epi[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(Epi) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(Epi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 1:5)
dev.off()
pdf(file.path(Epi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 6:10)
dev.off()
pdf(file.path(Epi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 11:15)
dev.off()
pdf(file.path(Epi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Epi[, suset_cells_Epi], dims = 16:20)
dev.off()


## dimension reduction --------
epi_dims <-  15


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Epi <- RunUMAP(Epi, dims = 1:epi_dims)


pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{epi_dims} project.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "group.y", label = TRUE) + NoLegend()
dev.off()


Epi <- Epi %>% 
  harmony::RunHarmony(
    dims.use = 1:epi_dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Epi <- RunUMAP(
  Epi, 
  dims = 1:epi_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{epi_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "group.y", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# plot QC figure
pdf(file.path(Epi_out_dir, str_glue("qc - Dimplot dim{epi_dims}.pdf")), height = 6, width = 8)
FeaturePlot(Epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI')) + NoLegend()
FeaturePlot(Epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("qc - Vlnplot dim{epi_dims}.pdf")), height = 6, width = 8)
VlnPlot(Epi, features = c('nCount_RNA', 'nFeature_RNA', 'log10GenesPerUMI'), group.by = 'sample', pt.size = 0) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("FeaturePlot - dim{epi_dims} - CellCycle.pdf")), height = 6, width = 12)
FeaturePlot(Epi, features = c('S.Score', 'G2M.Score'), reduction = "harmonyumap")
dev.off()

pdf(file.path(Epi_out_dir, str_glue("DimPlot - dim{epi_dims} - Phase.pdf")), height = 5, width = 6)
DimPlot(Epi, group.by = c('Phase'), reduction = "harmonyumap")
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
Epi <-
  FindNeighbors(
    Epi,
    reduction = reduction.use,
    dims = 1:epi_dims,
    k.param = 10
  )

Epi <-
  FindClusters(
    Epi,
    resolution = c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 5, width = 6)
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_hep, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Epi, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# Dotplot -------
# epi cluster
marker_epi <- list(
  "AT1"     = c("AGER","EMP2","CAV1","RTKN2","MYL9","HOPX","PDPN","MSLN","CLIC5"),
  # "AT1-AT2" = c("EGFR","SHH","MUC1","PDPN"),
  "AT2"     = c("MUC1","SFTPC","NPC2","SFTPA1","NAPSA","SFTPA2","CTSH","PGC","SFTPB","SFTPD","ETV5","ABCA3"),
  basal     = c("KRT5","KRT17","MIR205HG","DLK2","BCL11A","KRT14","TP63","DAPL1"),
  'ciliated cells' = c("CCDC113","FOXJ1","DNAH3","DNAIC2","DNAAF1","LRRC6"),
  'club cells'     = c("SCGB1A1","SCGB3A2","CYP2F2","CCKAR")
)

# Hep cluster
marker_Hep <- list(
  "PC_hep"    = c("Oat","Glul","Cyp2c37"),
  "MZ_hep"    = "Hamp",
  "PP_hep"    = c("Hal","Sds","Gls2"),
  "Cyc_hep"   = c("Mki67","Top2a","Cenpa")
)


table(Epi$RNA_snn_res.0.1)
pdf(
  file.path(Epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - major marker - RNA_snn_res.0.1.pdf")), 
  height = 5,
  width = 20
)
# major marker
DotPlot(
  Epi,
  features = c(
    "Immune cell" = "Ptprc",
    "T cell" = "Il2rb",
    "T cell" = "Ms4a4b",
    "T cell" = "Il7r",
    "T cell" = "Lck",
    "T cell" = "Cd3d",
    "T cell" = "Cd3e",
    "B cell" = "Cd19",
    "B cell" = "Cd79a",
    "B cell" = "Ms4a1",
    "PC_hep" = "Oat",
    "PC_hep" = "Glul",
    "PC_hep" = "Cyp2c37",
    "MZ_hep" = "Hamp",
    "PP_hep" = "Hal",
    "PP_hep" = "Sds",
    "PP_hep" = "Gls2",
    "Cyc_hep" = "Mki67",
    "Cyc_hep" = "Top2a",
    "Cyc_hep" = "Cenpa",
    "HSCs" = "Reln",
    "HSCs" = "Dcn",
    "HSCs" = "Nr1h5",
    "Portal_fibroblasts" = "Aebp1",
    "Portal_fibroblasts" = "Mgp", # "IGJ", "JCHAIN", "XBP1"
    "Portal_fibroblasts" = "Igfbp5",
    "Fibroblasts" = "Pdgfra",
    "Fibroblasts" = "Col1a1", # "COL1A1","DCN","ACTA2"
    "Fibroblasts" = "Col1a2",
    "LSECs" = "Kdr",
    "LSECs" = "Flt1",
    "LVECs" = "Vwf",
    "LVECs" = "Edn1",
    "LVECs" = "Jam2",
    "Cholangiocytes" = "Spp1",
    "Cholangiocytes" = "Epcam",
    "Cholangiocytes" = "Glis3",
    "Kupffer" = "Clec4f", #  "KRT5"
    "Kupffer" = "Cd5l",
    "Kupffer" = "Trpm2",
    "Dendritic" = "Siglech", # "COL1A1","DCN","ACTA2"
    "Dendritic" = "Irf8",
    "Dendritic" = "Runx2",
    "Monocytes/non-Kupffer_cell_macrophages" = "Ccr2",
    "Monocytes/non-Kupffer_cell_macrophages" = "Myo1f", # "VWF","AQP1"
    "Monocytes/non-Kupffer_cell_macrophages" = "Gpnmb",
    "Pericytes" = "Rgs5",
    "Pericytes" = "Kcnj8",
    "Pericytes" = "Higd1b",
    "Pericytes" = "Abcc9",
    "MDSCs" = "Cd33",
    "MDSCs" = "Itgam"
  ),
  group.by = "RNA_snn_res.0.1"
) + 
  RotatedAxis()
dev.off()

pdf(
  file.path(Epi_out_dir, str_glue("Dotplot dim{epi_dims} cluster - {reduction.use} - marker_epi - RNA_snn_res.0.1.pdf")), 
  height = 5,
  width = 10
)
# liver hep marker
DotPlot(
  Epi,
  features = marker_epi,
  group.by = "RNA_snn_res.0.1"
) +
  RotatedAxis()
dev.off()

### dotplot ------
major_marker_ls<-c(
  "epithelial state" = "CDH1",
  "epithelial state" = "EPCAM",
  "epithelial state" = "CLDN3",
  "epithelial state" = "OCLN",
  "mesenchymal state" = "VIM",
  "mesenchymal state" = "CDH2",
  "mesenchymal state" = "FN1",
  "EMT" = "SNAI1",
  "EMT" = "SNAI2",
  "EMT" = "TWIST1",
  "EMT" = "TWIST2",
  "EMT" = "ZEB1",
  "EMT" = "ZEB2",
  "EMT" = "TGFB1",
  "EMT" = "MMPS"
)

celltype_marker<-c(
  "B cell"= "CD19", 
  "B cell"= "CD79A", 
  "B cell"= "CD79B", 
  "B cell"= "MS4A1",
  "Plasma" = "TNFRSF17", 
  "Plasma" = "MZB1",
  # "Plasma" = "XBP1",
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "CD4",
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  "Endothelial" = "PECAM1",
  "Endothelial" = "CLDN5", 
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT18",
  "Epithelial" = "CDH1",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Monocytes" = "ITGAM",
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  "Monocytes" = "CD68",
  "Macrophage" = "CSF1R",
  "Macrophage" = "MARCO", 
  "DC"= "FCER1A", 
  "DC"= "CD1C", 
  "DC"= "CD83", 
  "DC"= "LILRA4", 
  "DC"= "CLEC9A", 
  "DC"= "IDO1",
  # 'pDC' = 'CLEC4C',
  # 'pDC' = "LILRA4",
  # 'pDC' = 'IL3RA',
  # 'pDC' = 'IRF8',
  # 'mregDC' = 'CCR7',
  # 'mregDC' ='CCL19',
  # 'mregDC' ='LAMP3',
  # 'mregDC' ='CD83',
  "NK" = "NCAM1",
  "Mast cell" = "CTSG", 
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = 'HDC', 
  'Neutrophil'= "CSF3R", 
  'Neutrophil'= "FCGR3B", # "S100A8", "S100A9" 
  "Melanocy" = "MITF",
  "Melanocy" = "MLANA"
)

# dotplot
pdf(file.path(Epi_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Epi, group.by = "RNA_snn_res.0.1", features = celltype_marker)+RotatedAxis()
DotPlot(Epi, group.by = "RNA_snn_res.0.1", features = major_marker_ls)+RotatedAxis()
DotPlot(all_data, group.by = "res_0.8", features = major_marker_ls)+RotatedAxis()
dev.off()

# 保存中间文件
saveRDS(Epi, file.path(Epi_out_dir, "Epi.rds"))

#加入上皮再分群cluster到all_data----
Epi_anno_meta <-
  Epi[['RNA_snn_res.0.1']] %>% 
  rownames_to_column("cell_id") 

#B cells  pDCs  T cells Endothelial  Macrophage Fibroblasts        
#Plasma  MergeDC Epithelial Smooth muscle cells

all_data_anno_meta <- 
  all_data[[]] %>% dplyr::filter(nCount_RNA > 100) %>% .['annotation'] %>%
  # all_data[['annotation']] %>%
  rownames_to_column('cell_id') %>%
  # dplyr::filter(annotation %in% c("Club cells","TRB-SC","Basal","AT2","Plasma","Macrophage","T cells")) %>%
  dplyr::filter(annotation %in% c("Epithelial","T cells")) %>%
  left_join(Epi_anno_meta, by = "cell_id") %>% 
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(
    annotation %in% "Epithelial",
    RNA_snn_res.0.1,
    annotation
  )) %>% 
  dplyr::select(cell_id, annotation = annotation2)


# all_data <- all_data %>% AddMetaData(all_data_anno_meta)  
meta_data_infer_subset <-
  all_data_anno_meta %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.1)

### inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹

infercnv_out_dir = "infercnv_hepsub_variable_all500"
dir.create(infercnv_out_dir) 


# meta_data_infer <-
#   all_data[[]] %>%
#   rownames_to_column("cell_id") %>%
#   dplyr::filter(
#     annotation %in% c("Club cells","TRB-SC","Basal","AT2",'Endothelial',"Plasma","Macrophage","T cells")
#   ) %>%
#   dplyr::select(cell_id, annotation, Spatial.008um_snn_res.0.8) %>%
#   map_dfc(as.character) 

# unique(meta_data_infer$annotation)
# "Acinar cells"     "Ductal cells"     "Endocrine cells"  "Epithelial cells" "NKT cells"


# cell is too much
# less than 50000
# 取1/20细胞： 1179497 --> 47179.88
# 1662182 --> 83000
# meta_data_infer_subset <-
#   meta_data_infer %>%
#   group_by(annotation) %>%
#   slice_sample(prop = 0.5)

# 合并annotation
# meta_data_infer_subset <-
#   meta_data_infer_subset %>%
#   mutate(Spatial.008um_snn_res.0.8 = str_c('Cluster', Spatial.008um_snn_res.0.8, sep = "_")) %>%
#   mutate(annotation2 = ifelse(
#     annotation %in% c("Plasma", "Macrophage","T cells","Endothelial"),
#     annotation,
#     Spatial.008um_snn_res.0.8
#   )) %>%
#   dplyr::select(
#     cell_id,
#     annotation = annotation2
#   )

# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )


all_data <- FindVariableFeatures(all_data)

#undebug(CreateInfercnvObject)
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(all_data)[1:500], meta_data_infer_subset$cell_id] ,
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path('infercnv2', "gene_pos47.txt"),
  
  min_max_counts_per_cell = c(30, Inf),
  
  delim             = "\t",
  ref_group_names   = "T cells"
)




# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.001) } # 304 <--- 
infercnv_obj@count.data %>% rowMeans %>% { sum(. < 0.005) } # 218


infercnv_obj2  = infercnv::run(
  infercnv_obj,
  cutoff              =    0.001,
  out_dir    =    infercnv_out_dir,
  cluster_by_groups   =    TRUE,
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  png_res             =    300,
  # useRaster   = FALSE,
  output_format = 'pdf',
  no_plot = T #  figure margin too large ----
)
infercnv::plot_cnv(
  infercnv_obj2,
  out_dir           =  infercnv_out_dir, 
  output_filename   =  'infercnv',
  output_format     =  "pdf",
  png_res = 300, 
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)

saveRDS(infercnv_obj2, file.path(infercnv_out_dir, 'infercnv_obj2.rds'))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, 'infercnv_obj.rds'))

# # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
# infercnv_obj_default = infercnv::run(
#   infercnv_obj,
#   cutoff              =    0.01, 
#   out_dir    =    infercnv_out_dir,
#   cluster_by_groups   =    TRUE, 
#   plot_steps          =    FALSE,
#   denoise             =    TRUE,
#   HMM                 =    FALSE,
#   no_prelim_plot      =    TRUE,
#   png_res             =    300,
#   # useRaster   = FALSE,
#   output_format = 'pdf',
#   no_plot = T #  figure margin too large ----
# )
# 
# # C stack error
# infercnv::plot_cnv(
#   infercnv_obj_default,
#   out_dir           =  infercnv_out_dir,
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 60,
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )


# download inferCNV object from linux server
# plot inferCNV locally
# saveRDS(infercnv_obj_default, file.path(infercnv_out_dir, 'infercnv_obj_default.rds'))
# 
# infercnv_obj_default <- readRDS("inferCNV/infercnv_obj_default.rds")
# 
# infercnv::plot_cnv(
#   infercnv_obj_default,
#   out_dir           =  'inferCNV', 
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 300, 
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )
# 
infercnv_obj2_default_median_filtering <-
  infercnv_obj2 %>%
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           =  infercnv_out_dir,
    output_filename   =  'infercnv_median_filtering',
    # output_format     =  "pdf",
    output_format="pdf",
    cluster_by_groups =  TRUE,
    x.range           =  "auto",
    x.center          =  1,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering .rds'))

# expFile='expFile.txt'
# groupFiles='groupinfo.txt'
# geneFile='gene_pos.txt'
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
#                                     annotations_file=groupFiles,
#                                     delim="\t",
#                                     gene_order_file= geneFile,
#                                     ref_group_names=c('B_cells','plasma_cells','T_cells','microglial')) 
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics



# infercnv::plot_cnv(
#   infercnv_obj2,
#   out_dir           =  'infercnv_output_2',
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 300,
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.01, 
#                               out_dir= "infercnv_output_2",  
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               hclust_method="ward.D2", plot_steps=F,no_plot = T)
# infercnv::plot_cnv(infercnv_obj2, 
#                    out_dir='infercnv_output_2',
#                    output_filename='infercnv', 
#                    x.range="auto",
#                    x.center=1,
#                    title = "infercnv", 
#                    output_format="pdf",
#                    color_safe_pal = FALSE)
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.1, 
#                               out_dir= "infercnv_output_sd",  
#                               write_expr_matrix=T,
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               sd_amplifier=3,
#                               denoise=T,
#                               noise_logistic=TRUE,
#                               hclust_method="ward.D2", plot_steps=F)
# saveRDS(infercnv_obj2, file.path(infercnv_output, 'infercnv_obj2.rds'))
# saveRDS(infercnv_obj, file.path(infercnv_output, 'infercnv_obj.rds'))


# cnv score ------------
library(scales)
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}
# cnvScore <- function(data){
#   sd_num <- apply(data, 1, sd)
#   data <- data[sd_num > median(sd_num), ]
#   
#   data <- data %>% as.matrix() %>%
#     t() %>% 
#     scale() %>% 
#     t()
#   
#   cnv_score <- as.data.frame(colSums(data * data))
#   return(cnv_score)
# }
cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <-'CNV_score'
temp <-
  cnv_score %>% 
  rownames_to_column('cell_id') %>% 
  left_join(.,as.data.frame(meta_data_infer_subset),'cell_id')  

p_infer_cnv_score <-
  temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+geom_violin(aes(fill= annotation),color="NA",drop = FALSE)+
  ggsci::scale_fill_igv()+
  # scale_fill_manual(values = color_v)+
  theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height = 6,
  width  = 12
)

cols_brewer = c("#f94144","#f25c27","#eb770a","#f89f12","#f9c74f","#90be6d",
                "#43aa8b","#4d908e","#577590","#277da1","#967aa1",
                "#F2698D","#F391AE","#F3B9CF","#AF8EB8","#7B84B7","#3F6FAF",
                "#035AA6","#0473B4","#8491B7","#737AAB","#5E5E97")

p_infer_cnv_score <-
  temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+
  geom_boxplot(aes(fill= annotation),color="NA",outlier.colour = "black")+
  # ggsci::scale_fill_igv()+
  scale_fill_manual(values = cols_brewer)+
  theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_boxplot.pdf"),
  height = 6,
  width  = 12
)



# spatial plot --------------
# Idents(all_data) <- 'Spatial.008um'
# Idents(all_data) <- 'Spatial.qc.pass'

sc_hep_infer_anno <-
  sc_hep[['Spatial.008um_snn_res.0.4']] %>% 
  rownames_to_column("cell_id") %>%
  mutate(annotation_epi = ifelse(
    Spatial.008um_snn_res.0.4 %in% c('3', '6', '12' ,'13'),
    "Tumor cells",
    "Epithelial cells"
    # Spatial.008um_snn_res.0.4
  )) 

all_data_infer_anno <-
  all_data[['annotation']] %>%
  rownames_to_column("cell_id") %>%
  left_join(sc_hep_infer_anno, by = "cell_id") %>%
  mutate(annotation2 = ifelse(
    is.na(annotation_epi),
    as.character(annotation),
    as.character(annotation_epi)
  )) %>%
  column_to_rownames("cell_id")


sp_p_ls <-
  Images(all_data) %>%
  set_names(., .) %>%
  map(function(img){
    all_data %>%
      AddMetaData(all_data_infer_anno['annotation2']) %>%
      SpatialDimPlot(
        images = img, 
        group.by = "annotation2",
        pt.size.factor	= 2,
        image.scale = "lowres",
        crop = FALSE
      ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img)
  })


# export
pdf(
  file.path(hep_out_dir, "spatial plot of infercnv analysis.pdf"),
  height = 8,
  width  = 8
)
sp_p_ls %>% walk(print)
dev.off()





# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","Ni_LN","Nmi_LN","N1_LN"





