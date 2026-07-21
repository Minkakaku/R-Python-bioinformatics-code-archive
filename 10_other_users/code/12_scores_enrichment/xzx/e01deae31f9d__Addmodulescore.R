setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/")
library(Seurat)
library(ggplot2)


sc <- readRDS('./h5_trans.rds')

genesets_addmodulescore <- function(gene_sets,gene_name,scRNA,group_name,out_dir){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 500,
    nbin = 5,
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

Immune_Checkpoint_Activation_Receptor_Gene <- c("CD28","ICOS","TNFRSF9","TNFRSF14","TNFSF14","TNFRSF4","CD27","TNFRSF18","CD226","CD244","ADORA2A","CD209","TNFRSF25")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "group",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Activation_Receptor_Gene,gene_name = "Immune_Checkpoint_Activation_Receptor_Gene",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/")

Immune_Checkpoint_Inhibition_Receptor<- c("LAG3","PDCD1","CTLA4","TMIGD2","BTLA","CD160","TIGIT","CD96","HAVCR2","AHR","CEACAM1","SIRPA","CD200R1","CD44")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "group",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/")
sc <- genesets_addmodulescore(gene_sets = Immune_Checkpoint_Inhibition_Receptor,gene_name = "Immune_Checkpoint_Inhibition_Receptor",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/SC2507294017/")

MSC_like <- c("ENG", "MCAM", "THY1", "NT5E", "VIM", "PDGFRB")
sc <- genesets_addmodulescore(gene_sets = MSC_like,gene_name = "MSC_like",scRNA = sc,group_name = "cluster",out_dir = "/mnt/NexData_file/home/xzx/project/FWSC20253283/")
sc <- genesets_addmodulescore(gene_sets = MSC_like,gene_name = "MSC_like",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/FWSC20253283/")

Pluripotent_like <- c("POU5F1", "KLF4","NR5A2","MYC")
sc <- genesets_addmodulescore(gene_sets = Pluripotent_like,gene_name = "Pluripotent_like",scRNA = sc,group_name = "cluster",out_dir = "/mnt/NexData_file/home/xzx/project/FWSC20253283/")
sc <- genesets_addmodulescore(gene_sets = Pluripotent_like,gene_name = "Pluripotent_like",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/FWSC20253283/")

Senescent_like <- c("CDKN2A", "CDKN1A", "TGFB1", "TGM2", "GLB1")
sc <- genesets_addmodulescore(gene_sets = Senescent_like,gene_name = "Senescent_like",scRNA = sc,group_name = "cluster",out_dir = "/mnt/NexData_file/home/xzx/project/FWSC20253283/")
sc <- genesets_addmodulescore(gene_sets = Senescent_like,gene_name = "Senescent_like",scRNA = sc,group_name = "cell_type",out_dir = "/mnt/NexData_file/home/xzx/project/FWSC20253283/")
