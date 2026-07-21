library(Seurat)

sc <- readRDS("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou1/h5_trans.rds")

exclude_types <- c(
  "Endothelial cells", 
  "Enteroendocrine", 
  "Macrophages", 
  "Mast cells", 
  "T cells",
  "Epithelium",
  "Tuft"
)

# 1. 提取原始数据（表达矩阵、元数据、降维结果等）
expr_matrix <- GetAssayData(sc, assay = "RNA", slot = "counts")  # 替换为你的实际assay（如"SCT"）
meta_data <- sc@meta.data
#reductions <- sc@reductions  # 保留降维结果（如UMAP、PCA）

# 2. 筛选需要保留的细胞
# 筛选表达矩阵和元数据
keep_cells <- colnames(sc)[! (sc$cell_type %in% exclude_types)]
expr_filtered <- expr_matrix[, keep_cells, drop = FALSE]
meta_filtered <- meta_data[keep_cells, , drop = FALSE]

# 3. 重新构建Seurat对象
sc_filtered <- CreateSeuratObject(
  counts = expr_filtered,
  meta.data = meta_filtered,
  project = "filtered"
)
# 替换特定标签
# 确保 cell_type 是字符而不是 factor
sc_filtered$cell_type <- as.character(sc_filtered$cell_type)

# 替换
sc_filtered$cell_type[sc_filtered$cell_type == "Myofibroblast 0"] <- "Schwann cells"

sc_filtered$cell_type[sc_filtered$cell_type == "Smooth muscle cells 0"] <- "SMC 0"
sc_filtered$cell_type[sc_filtered$cell_type == "Smooth muscle cells 1"] <- "SMC 1"
sc_filtered$cell_type[sc_filtered$cell_type == "Smooth muscle cells 2"] <- "SMC 1"
sc_filtered$cell_type[sc_filtered$cell_type == "Smooth muscle cells 3"] <- "SMC 2"
sc_filtered$cell_type[sc_filtered$cell_type == "Smooth muscle cells 4"] <- "SMC 0"
sc_filtered$cell_type[sc_filtered$cell_type == "Smooth muscle cells 5"] <- "SMC 3"

sc_filtered$cell_type[sc_filtered$cell_type == "SMC 3"] <- "cajal cells"
# 再次检查
table(sc_filtered$cell_type)

### 继续分不同子集进行评分可视化
Idents(sc_filtered) <- sc_filtered$cell_type
sc_myoFB <- subset(sc_filtered,idents=c("Schwann cells","Myofibroblast 1","Myofibroblast 2","Myofibroblast 3","Myofibroblast 4"))
sc_SMC <- subset(sc_filtered,idents=c("SMC 0","SMC 1","SMC 2","cajal cells"))
sc_final <- merge(sc_myoFB,sc_SMC)
table(sc_final$cell_type)

cell_type_colors <- c(
  "Myofibroblast 1" = "#377EB8",
  "Myofibroblast 2" = "#4DAF4A",
  "Myofibroblast 3" = "#984EA3",
  "Myofibroblast 4" = "#FF7F00",
  "Schwann cells" = "#FFFF33",
  "cajal cells" = "#E41A1C",
  "SMC 0" = "#A65628",
  "SMC 1" = "#F781BF",
  "SMC 2" = "#F44684"
)

cell_type_order <- c(
  "Schwann cells",
  "Myofibroblast 1",
  "Myofibroblast 2",
  "Myofibroblast 3",
  "Myofibroblast 4",
  "cajal cells",
  "SMC 0",
  "SMC 1",
  "SMC 2"
)


### ====常用函数：批量进行目标基因集addmodulescore ==== ###
genesets_addmodulescore <- function(gene_sets,gene_name,scRNA){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 100,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    cell_type = scRNA_fin$cell_type,  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(gene_name," Module Score"), x = "Cell Type", y = "Score")+
    scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou1/scscore/",gene_name,"_addmodulescore.png"),
         width = 16,height = 9,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0("/mnt/NexData_file/home/xzx/project/SC2507113630/shouhou1/scscore/",gene_name,"_addmodulescore.pdf"),
         width = 16,height = 9,dpi = 300,bg = "white")
  
  return(scRNA_fin)
}



Schwann_genes <- c("PLP1","S100B","UCHL1","SOX10","CDH19","CRYAB","PMP22","SCN7A","NGFR")
#sc_filtered <- genesets_addmodulescore(gene_sets = Schwann_genes,gene_name = "SchwannGenes",scRNA = sc_filtered)
sc_final <- genesets_addmodulescore(gene_sets = Schwann_genes,gene_name = "SchwannGenes",scRNA = sc_final)


cajal_genes <- list(c("KIT","ETV1","ANO1","VIM","PDGFRA","NOS1","HCN1","HCN2"))
cajal_genes <- c("KIT", "ANO1", "CACNA1H", "HCN4", "GJA1","ETV1","VIM","PDGFRA","NOS1","HCN1","HCN2")
#sc_filtered <- genesets_addmodulescore(gene_sets = cajal_genes,gene_name = "cajalGenes",scRNA = sc_filtered)
sc_final <- genesets_addmodulescore(gene_sets = cajal_genes,gene_name = "cajalGenes",scRNA = sc_final)



sc_myoFB <- genesets_addmodulescore(gene_sets = Schwann_genes,gene_name = "myoFB_SchwannGenes",scRNA = sc_myoFB)

sc_SMC <- genesets_addmodulescore(gene_sets = cajal_genes,gene_name = "SMC_cajalGenes",scRNA = sc_SMC)





####   特定基因集评分  ####
library(Seurat)
library(dplyr)
#library(tidyverse)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
#library(gt)
library(KEGGREST)
library(clusterProfiler)
#library(msigdbr)
library(data.table)
library(enrichplot)

####
gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/all_IonChannel_hsa.gmt") 
ion_channel <- gene_sets$gene
sc_final <- genesets_addmodulescore(gene_sets = ion_channel,gene_name = "ion_channel",scRNA = sc_final)
sc_SMC <- genesets_addmodulescore(gene_sets = ion_channel,gene_name = "SMC_ion_channel",scRNA = sc_SMC)
sc_myoFB <- genesets_addmodulescore(gene_sets = ion_channel,gene_name = "myoFB_ion_channel",scRNA = sc_myoFB)


####
# gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/Neuroactive_ligand_receptor.gmt") 
# Neuroactive_ligand_receptor <- gene_sets$gene
# sc_final <- genesets_addmodulescore(gene_sets = Neuroactive_ligand_receptor,gene_name = "Neuroactive_ligand_receptor",scRNA = sc_final)


####Schwann细胞相关通路
gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/GOBP_REGULATION_OF_MYELINATION.gmt") 
REGULATION_OF_MYELINATION <- gene_sets$gene
#sc_final <- genesets_addmodulescore(gene_sets = REGULATION_OF_MYELINATION,gene_name = "REGULATION_OF_MYELINATION",scRNA = sc_final)
sc_myoFB <- genesets_addmodulescore(gene_sets = REGULATION_OF_MYELINATION,gene_name = "myoFB_REGULATION_OF_MYELINATION",scRNA = sc_myoFB)

gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/GOBP_PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT.gmt") 
PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT <- gene_sets$gene
#sc_final <- genesets_addmodulescore(gene_sets = PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT,gene_name = "REGULATION_OF_MYELINATION",scRNA = sc_final)
sc_myoFB <- genesets_addmodulescore(gene_sets = PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT,gene_name = "myoFB_PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT",scRNA = sc_myoFB)

gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/GOBP_REGULATION_OF_AXON_GUIDANCE.gmt") 
REGULATION_OF_AXON_GUIDANCE <- gene_sets$gene
#sc_final <- genesets_addmodulescore(gene_sets = PERIPHERAL_NERVOUS_SYSTEM_AXON_ENSHEATHMENT,gene_name = "REGULATION_OF_MYELINATION",scRNA = sc_final)
sc_myoFB <- genesets_addmodulescore(gene_sets = REGULATION_OF_AXON_GUIDANCE,gene_name = "myoFB_REGULATION_OF_AXON_GUIDANCE",scRNA = sc_myoFB)


####cajal细胞相关通路
gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/selectedpathway.gmt") 
Calcium_signaling <- gene_sets[gene_sets$term == "Calcium signaling - Homo sapiens (human)",]
Calcium_signaling <- Calcium_signaling$gene
#sc_final <- genesets_addmodulescore(gene_sets = Calcium_signaling,gene_name = "Calcium_signaling",scRNA = sc_final)
sc_SMC <- genesets_addmodulescore(gene_sets = Calcium_signaling,gene_name = "SMC_Calcium_signaling",scRNA = sc_SMC)

Gap_junction <- gene_sets[gene_sets$term == "Gap junction - Homo sapiens (human)",]
Gap_junction <- Gap_junction$gene
sc_SMC <- genesets_addmodulescore(gene_sets = Gap_junction,gene_name = "SMC_Gap_junction",scRNA = sc_SMC)




gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/GOMF_VOLTAGE_GATED_CALCIUM_CHANNEL_ACTIVITY.gmt")
VOLTAGE_GATED_Ca <- gene_sets$gene
gene_sets <- read.gmt("./SC2507113630/shouhou1/scscore/GOMF_VOLTAGE_GATED_POTASSIUM_CHANNEL_ACTIVITY.gmt")
VOLTAGE_GATED_Ka <- gene_sets$gene
VOLTAGE_GATED <- c(VOLTAGE_GATED_Ca,VOLTAGE_GATED_Ka)

#sc_final <- genesets_addmodulescore(gene_sets = VOLTAGE_GATED,gene_name = "VOLTAGE_GATED_ION_CHANNEL_ACTIVITY",scRNA = sc_final)
sc_SMC <- genesets_addmodulescore(gene_sets = VOLTAGE_GATED,gene_name = "SMC_VOLTAGE_GATED_ION_CHANNEL_ACTIVITY",scRNA = sc_SMC)






#####################
table(sc_final$cell_type)
# 查看当前所有的 assay 名称
names(sc_final@assays)
saveRDS(sc_final, "./SC2507113630/shouhou1/results/cellchat_final.rds")
