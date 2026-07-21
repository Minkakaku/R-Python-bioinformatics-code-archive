setwd("/mnt/NexData_file/home/xzx/project/SC2508064183/shouhou2/")

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

pbmc <- readRDS('./Epithelium/h5_trans.rds')
Idents(pbmc)
pbmc <- subset(pbmc,ident=c('Epithelial_0', 'Epithelial_1', 'Epithelial_2', 'Epithelial_3', 'Epithelial_5'))
clustercolname <- 'cell_type'
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
write.table(deg_output,file="monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
length(ordergene)
#cds <- setOrderingFilter(cds,ordergene[1:min(length(ordergene),2000)])
cds <- setOrderingFilter(cds,ordergene)


pdf("monocle.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
png("monocle.ordergenes.png")
plot_ordering_genes(cds)
dev.off()

cds <- reduceDimension(cds,max_components = 2,reduction_method='DDRTree')#,residualModelFormulaStr = "~sample"
cds <- orderCells(cds)#重设起始点使用参数root_state=1


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

plotcol = c("sample","group","State","cell_type")
#if(!is.null(args$type)){plotcol = c(plotcol,typeuse)}
#if("cell_type" %in% colnames(pbmc@meta.data)){plotcol = c(plotcol,"cell_type")}
for(i in plotcol){
  if(i %in% colnames(cds@phenoData@data)){
    monocleplot(cds,colnames = i,ptsize=1)
  }else{
    print(paste0(i," is not exist in colnames"))
  }
  
}

clustercolname <- 'cell_type'
p = plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE)
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


pseudotime_df = cds@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"monocle.metadata.xls",sep = "\t",quote = F,row.names = F)


outpath = '/mnt/NexData_file/home/xzx/project/SC2508064183/shouhou2/'
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
}


outpath = '/mnt/NexData_file/home/xzx/project/SC2508064183/shouhou2/Epithelium/monocle/'
#结果报告
file.copy("/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/pseudotime/report.rmd",paste0(outpath,"/report.rmd"))
rmarkdown::render(paste0(outpath,"/report.rmd"),output_file=paste0(outpath,"/report.html"))
file.remove(paste0(outpath,"/report.rmd"))


saveRDS(cds,"monocle.rds")

