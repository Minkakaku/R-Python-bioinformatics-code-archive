setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/")

library(Seurat)
library(ggplot2)

genesets_addmodulescore <- function(gene_sets,gene_name,scRNA,group_name,out_dir){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 200,
    #nbin = 5,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    cell_type = sc@meta.data[[group_name]],  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  #plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0(gene_name," Module Score"), x = "Cell Type", y = "Score")
  #+
  #  scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.png"),
         width = 10,height = 6,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.pdf"),
         width = 10,height = 6,dpi = 300,bg = "white")
  
  return(scRNA_fin)
}


##T
sc <- readRDS('./shouhou2/recluster_fin/T_cells/h5_trans.rds')

Immune_Checkpoint_Activation_Receptor_Gene <- c("CD28","ICOS","TNFRSF9","TNFRSF14","TNFSF14","TNFRSF4","CD27","TNFRSF18","CD226","CD244","ADORA2A","CD209","TNFRSF25")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "group",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/T_cells/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/T_cells/")

Immune_Checkpoint_Inhibition_Receptor<- c("LAG3","PDCD1","CTLA4","TMIGD2","BTLA","CD160","TIGIT","CD96","HAVCR2","AHR","CEACAM1","SIRPA","CD200R1","CD44")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "group",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/T_cells/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "cell_type",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/T_cells/")


##myeloid
sc <- readRDS('./shouhou2/recluster_fin/myeloid/h5_trans.rds')

Immune_Checkpoint_Activation_Receptor_Gene <- c("CD28","ICOS","TNFRSF9","TNFRSF14","TNFSF14","TNFRSF4","CD27","TNFRSF18","CD226","CD244","ADORA2A","CD209","TNFRSF25")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "group",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/myeloid/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/myeloid/")

Immune_Checkpoint_Inhibition_Receptor<- c("LAG3","PDCD1","CTLA4","TMIGD2","BTLA","CD160","TIGIT","CD96","HAVCR2","AHR","CEACAM1","SIRPA","CD200R1","CD44")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "group",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/myeloid/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou2/recluster_fin/myeloid/")


##B
sc <- readRDS('./shouhou3/Bcells/celltype/h5_trans.rds')

Immune_Checkpoint_Activation_Receptor_Gene <- c("CD28","ICOS","TNFRSF9","TNFRSF14","TNFSF14","TNFRSF4","CD27","TNFRSF18","CD226","CD244","ADORA2A","CD209","TNFRSF25")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "group",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Bcells/celltype/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "cell_type",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Bcells/celltype/")

Immune_Checkpoint_Inhibition_Receptor<- c("LAG3","PDCD1","CTLA4","TMIGD2","BTLA","CD160","TIGIT","CD96","HAVCR2","AHR","CEACAM1","SIRPA","CD200R1","CD44")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "group",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Bcells/celltype/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "cell_type",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Bcells/celltype/")

##FB
sc <- readRDS('./shouhou3/Fibroblasts/celltype/h5_trans.rds')

Immune_Checkpoint_Activation_Receptor_Gene <- c("CD28","ICOS","TNFRSF9","TNFRSF14","TNFSF14","TNFRSF4","CD27","TNFRSF18","CD226","CD244","ADORA2A","CD209","TNFRSF25")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "group",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Fibroblasts/celltype/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "cell_type",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Fibroblasts/celltype/")

Immune_Checkpoint_Inhibition_Receptor<- c("LAG3","PDCD1","CTLA4","TMIGD2","BTLA","CD160","TIGIT","CD96","HAVCR2","AHR","CEACAM1","SIRPA","CD200R1","CD44")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "group",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Fibroblasts/celltype/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "cell_type",
                              out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou3/Fibroblasts/celltype/")

