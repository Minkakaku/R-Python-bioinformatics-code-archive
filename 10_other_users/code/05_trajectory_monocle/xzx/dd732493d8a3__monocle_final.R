setwd("/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/")

rm(list = ls())

library(monocle)
library(Seurat)
library(ggplot2)
#library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)
library(getopt)



#读取参数和文件
command = matrix(c('rds','r',1,"character",
                   'type','t',2,"character",
                   'downcluster','D',0,"logical",
                   'maxcells','m',2,"numeric",
                   'cluster','c',2,"character",
                   'ptsize','p',2,"numeric",
                   'BEAM','B',0,"logical",
                   'output','o',1,"character")
                 ,byrow=TRUE, ncol=4)

args=getopt(command)

pbmc <- readRDS('./epi_trans.rds')
Idents(pbmc)
table(pbmc$cell_type)
#pbmc <- subset(pbmc,ident=c('Epithelial_0', 'Epithelial_1', 'Epithelial_2', 'Epithelial_3', 'Epithelial_5'))
clustercolname <- 'cell_type'

sample_cells_per_cluster <- function(seurat_object,colnames="leiden", num_cells = 1000) {
  # 获取cluster信息
  Idents(seurat_object)=colnames
  clusters <- Idents(seurat_object)
  cluster_ids <- levels(clusters)  # 获取所有cluster的名称
  
  sampled_cell_ids <- c() # 用于存储所有抽样的细胞ID
  
  for (cluster_id in cluster_ids) {
    # 获取当前cluster的所有细胞ID
    cell_ids_in_cluster <- names(clusters[clusters == cluster_id])
    
    # 确定抽样数量
    num_cells_to_sample <- min(num_cells, length(cell_ids_in_cluster))  # 如果cluster细胞数小于1000，则抽取所有
    
    # 随机抽样
    sampled_ids <- sample(cell_ids_in_cluster, size = num_cells_to_sample, replace = FALSE)
    sampled_cell_ids <- c(sampled_cell_ids, sampled_ids) # 将抽样的细胞ID添加到列表中
  }
  subseurat = subset(seurat_object,cells=sampled_cell_ids)
  return(subseurat)
}

downcluster = TRUE
if(downcluster){
  if(!is.null(args$maxcells)){maxcell = args$maxcells}else{maxcell = 1000}
  pbmc = sample_cells_per_cluster(pbmc,num_cells = maxcell,colnames = clustercolname)
  print("Cell number of each cluster:")
  print(table(pbmc@meta.data[,clustercolname]))
}

p_data <- pbmc@meta.data
print("Meta data:")
head(p_data)

f_data <- data.frame(gene_short_name=row.names(pbmc),row.names=row.names(pbmc))
pd <- new('AnnotatedDataFrame',data = p_data)
fd <- new('AnnotatedDataFrame',data = f_data)
cds <- newCellDataSet(pbmc@assays$RNA@counts,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())#count适用于negbinomial.size,data适用于tobit且后续无需运行estimateDispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr=0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed>10))


diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr ="~cell_type" )
deg <- subset(diff,qval<0.01)
deg <- deg %>% arrange(qval,desc(num_cells_expressed))
deg <- deg[order(deg$qval,decreasing = F),]
head(deg)
deg_output <- data.frame(Genes = rownames(deg),deg[,c("pval","qval","num_cells_expressed")])
write.table(deg_output,file="./monocle/monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
length(ordergene)
cds <- setOrderingFilter(cds,ordergene[1:min(length(ordergene),2000)])
#cds <- setOrderingFilter(cds,ordergene)


pdf("./monocle/monocle.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
png("./monocle/monocle.ordergenes.png")
plot_ordering_genes(cds)
dev.off()

cds <- reduceDimension(cds,max_components = 2,reduction_method='DDRTree')#,residualModelFormulaStr = "~sample"
cds <- orderCells(cds)#重设起始点使用参数root_state=1
#cds <- orderCells(cds,root_state = 2)#重设起始点使用参数root_state=1


monocleplot <- function(cds,colnames,ptsize=1){
  
  p = plot_cell_trajectory(cds,color_by=colnames,size=ptsize,show_backbone=TRUE)
  ggsave(paste0("monocle.trajectory.",colnames,".png"),p,width = 7,height = 7.5)
  ggsave(paste0("monocle.trajectory.",colnames,".pdf"),p,width = 7,height = 7.5)
  
  ncols = ceiling(sqrt(length(unique(cds@phenoData@data[,colnames]))))
  nrows = ceiling(length(unique(cds@phenoData@data[,colnames]))/ncols)
  p = plot_cell_trajectory(cds,color_by=colnames,size=ptsize,show_backbone=TRUE)+facet_wrap(vars(!!sym(colnames)),nrow=nrows)
  ggsave(paste0("monocle.trajectory.",colnames,".faceted.png"),p,width = 4*ncols,height = 4.2*nrows)
  ggsave(paste0("monocle.trajectory.",colnames,".faceted.pdf"),p,width = 4*ncols,height = 4.2*nrows)
  
}

plotcol = c("group","State","cell_type")
#if(!is.null(args$type)){plotcol = c(plotcol,typeuse)}
#if("cell_type" %in% colnames(pbmc@meta.data)){plotcol = c(plotcol,"cell_type")}
for(i in plotcol){
  if(i %in% colnames(cds@phenoData@data)){
    monocleplot(cds,colnames = i,ptsize=1)
  }else{
    print(paste0(i," is not exist in colnames"))
  }
  
}

clustercolname <- 'group'
p = plot_cell_trajectory(cds,color_by="Pseudotime",size=0.01,show_backbone=TRUE)
p
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)

p = plot_complex_cell_trajectory(cds,x=1,y=2,color_by = clustercolname)
p
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)

ncols = ceiling(sqrt(length(unique(cds@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(cds@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(cds,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
p
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)


pseudotime_df = cds@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State",'group')]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"monocle.metadata.xls",sep = "\t",quote = F,row.names = F)


outpath = '/mnt/NexData_file/home/xzx/project/SC2507294017/shouhou4/Epithelial_all/monocle/'
for(branch_point in 1:((length(levels(cds@phenoData@data[,"State"]))-1)/2)){
  outbeam = paste0(outpath,"/BEAM_analysis.Node",branch_point)
  print(paste0("Analysis Beam Node ",branch_point))
  if(!dir.exists(outbeam)){
    dir.create(outbeam)
  }
  
  BEAM_res = BEAM(cds,branch_point = branch_point,cores = 4, progenitor_method = "duplicate")
  BEAM_res = BEAM_res[,c("gene_short_name","pval","qval")]
  BEAM_res = BEAM_res[order(BEAM_res$qval),]
  write.csv(BEAM_res,paste0(outbeam,"/BEAM_res.csv"),row.names = F)
  selectgenes = row.names(subset(BEAM_res,qval<0.05))
  p = plot_genes_branched_heatmap(cds[selectgenes[1:min(100,length(selectgenes))],],
                                  branch_point = branch_point,
                                  num_clusters = 4,
                                  cores=1,
                                  use_gene_short_name = T,
                                  show_rownames = T,return_heatmap = T)
  ggsave(paste0(outbeam,"/BEAM_heatmap.png"),p$ph_res,width=7,height = 8)
  ggsave(paste0(outbeam,"/BEAM_heatmap.pdf"),p$ph_res,width=7,height = 8)
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "State",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.pdf"),p,width=15,height = 12)
  
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "group",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_group_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_group_top20.pdf"),p,width=15,height = 12)
  
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "cell_type",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_celltype_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_celltype_top20.pdf"),p,width=15,height = 12)
}


outpath = '/mnt/NexData_file/home/xzx/project/SC2508064183/shouhou2/Epithelium/monocle/'
#结果报告
file.copy("/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/pseudotime/report.rmd",paste0(outpath,"/report.rmd"))
rmarkdown::render(paste0(outpath,"/report.rmd"),output_file=paste0(outpath,"/report.html"))
file.remove(paste0(outpath,"/report.rmd"))


saveRDS(cds,"./monocle_1/monocle.rds")
save.image("./monocle_1/mono.Rdata")


for(branch_point in 1:((length(levels(cds@phenoData@data[,"State"]))-1)/2)){
  outbeam = paste0(outpath,"/BEAM_analysis.Node",branch_point)
  print(paste0("Analysis Beam Node ",branch_point))
  if(!dir.exists(outbeam)){
    dir.create(outbeam)
  }
  
  BEAM_res = BEAM(cds,branch_point = branch_point,cores = 4, progenitor_method = "duplicate")
  BEAM_res = BEAM_res[,c("gene_short_name","pval","qval")]
  BEAM_res = BEAM_res[order(BEAM_res$qval),]
  write.csv(BEAM_res,paste0(outbeam,"/BEAM_res.csv"),row.names = F)
  selectgenes = row.names(subset(BEAM_res,qval<0.05))
  p = plot_genes_branched_heatmap(cds[selectgenes[1:min(100,length(selectgenes))],],
                                  branch_point = branch_point,
                                  num_clusters = 4,
                                  cores=1,
                                  use_gene_short_name = T,
                                  show_rownames = T,return_heatmap = T)
  ggsave(paste0(outbeam,"/BEAM_heatmap.png"),p$ph_res,width=7,height = 8)
  ggsave(paste0(outbeam,"/BEAM_heatmap.pdf"),p$ph_res,width=7,height = 8)
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "group",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.pdf"),p,width=15,height = 12)
}

p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "group",ncol = 4)
p
ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.png"),p,width=15,height = 12)
ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.pdf"),p,width=15,height = 12)

s_gene <- as.data.frame(selectgenes)



# 1. 定义你要绘制的指定基因（替换为你的目标基因名，可任意数量）
specified_genes <- c("SHMT1","SHMT2","MTHFD2","PHGDH","PSAT1","PSPH","SFXN3","SLC1A4","SLC1A5")  # 示例基因，按需修改
valid_genes <- intersect(specified_genes, row.names(cds))
cat("有效基因数量：", length(valid_genes), "\n")
cat("有效基因列表：", valid_genes, "\n")
cat("当前使用的分支点：", branch_point, "\n")
branch_point <- 3
# 2. 绘制分支伪时间图（仅展示指定基因）
p <- plot_genes_branched_pseudotime(
  cds = cds[specified_genes, ],  # 核心修改：用指定基因筛选cds对象
  branch_point = branch_point,   # 保持原分支点参数
  color_by = "group",            # 保持原分组着色参数
  ncol = 4                       # 保持原子图列数（可按需调整，如3、5等）
)

p
ggsave(paste0(outbeam,"/BEAM_pseudotime_selectgene.png"),p,width=15,height = 12)
ggsave(paste0(outbeam,"/BEAM_pseudotime_selectgene.pdf"),p,width=15,height = 12)

# 3. 显示图片（可选，若需直接输出）
print(p)