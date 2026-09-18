# rm(list = ls())
.libPaths("/home/tt/R/x86_64-conda-linux-gnu-library/4.3")
library(monocle)
library(Seurat)
library(ggplot2)
library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)
library(getopt)

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

#trace(monocle:::project2MST,edit=T)
#trace(monocle::buildBranchCellDataSet,edit = T)

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
monocleplot <- function(cds,colnames,ptsize=1,colour){
	p = plot_cell_trajectory(cds,color_by=colnames,size=ptsize,show_backbone=TRUE)
	if(length(colour)!=1){
		p = p + scale_color_manual(values = colour)
	}
  ggsave(paste0("monocle.trajectory.",colnames,".png"),p,width = 7,height = 7.5)
  ggsave(paste0("monocle.trajectory.",colnames,".pdf"),p,width = 7,height = 7.5)

  ncols = ceiling(sqrt(length(unique(cds@phenoData@data[,colnames]))))
  nrows = ceiling(length(unique(cds@phenoData@data[,colnames]))/ncols)
	p = plot_cell_trajectory(cds,color_by=colnames,size=ptsize,show_backbone=TRUE)+ facet_wrap(vars(!!sym(colnames)),nrow=nrows)
	if(length(colour)!=1){
		p = p + scale_color_manual(values = colour)
	}
  ggsave(paste0("monocle.trajectory.",colnames,".faceted.png"),p,width = 4*ncols,height = 4.2*nrows)
  ggsave(paste0("monocle.trajectory.",colnames,".faceted.pdf"),p,width = 4*ncols,height = 4.2*nrows)
  
}

#读取参数和文件
rdsfile = args$rds
if(!is.null(args$cluster)){clustercolname = args$cluster}else{clustercolname = "leiden"}
if(is.null(args$downcluster)){args$downcluster = FALSE}
if(!is.null(args$downcluster)){downcluster = args$downcluster}else{downcluster=T}
if(is.null(args$BEAM)){args$BEAM = FALSE}
if(!is.null(args$BEAM)){BEAM = args$BEAM}else{BEAM=T}

pbmc <- readRDS(rdsfile)
pbmc<-subset(pbmc,cell_type %in% c("Classical Monocytes CD302","Classical Monocytes IFI44","Classical Monocytes NFKBIA","Classical Monocytes RETN ","Macrophage","Non-classical Monocytes"))

if(!is.null(args$type)){
  typeuse = str_split(args$type,pattern = ",")[[1]]
}

outpath = args$output
if(!dir.exists(outpath)){dir.create(outpath)}
setwd(outpath)

#随机抽取细胞
if(downcluster){
  if(!is.null(args$maxcells)){maxcell = args$maxcells}else{maxcell = 1000}
  pbmc = sample_cells_per_cluster(pbmc,num_cells = maxcell,colnames = clustercolname)
  print("Cell number of each cluster:")
  print(table(pbmc@meta.data[,clustercolname]))
}

#colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))]="Type"
#colnames(pbmc@meta.data)[length(colnames(pbmc@meta.data))-1]="leiden"
#head(pbmc@meta.data)
# expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts),'sparseMatrix')
p_data <- pbmc@meta.data
print("Meta data:")
head(p_data)

f_data <- data.frame(gene_short_name=row.names(pbmc),row.names=row.names(pbmc))
pd <- new('AnnotatedDataFrame',data = p_data)
fd <- new('AnnotatedDataFrame',data = f_data)
cds <- newCellDataSet(pbmc@assays$RNA$counts,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())#count适用于negbinomial.size,data适用于tobit且后续无需运行estimateDispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr=0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed>10))

#monocle选择的高变基因
# disp_table <- dispersionTable(cds)
# disp.genes <- subset(disp_table,mean_expression >=0.1 &dispersion_empirical >=1*dispersion_fit)$gene_id
# cds <- setOrderingFilter(cds,disp.genes)
# plot_ordering_genes(cds)

#dpFeature
print("Analysis dpFeature:")
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

print("ReducDimension and order:")
#ordergene <- row.names(deg)[order(deg$qval)][1:400]
cds <- reduceDimension(cds,max_components = 2,reduction_method='DDRTree', num_dim = 6)#,residualModelFormulaStr = "~sample"
cds <- orderCells(cds)#重设起始点使用参数root_state=1
cds <- orderCells(cds,root_state=5)
ggdim = plot_pc_variance_explained(cds,
    return_all = F,
    max_components = 30)
ggsave("Pseudotime_dim.png", ggdim, width = 8, height = 5, limitsize = FALSE, bg = "white")
ggsave("Pseudotime_dim.pdf", ggdim, width = 8, height = 5, limitsize = FALSE, bg = "white")
	
saveRDS(cds,"monocle.rds")
cds<-readRDS("monocle.rds")
plotcol = c("group",clustercolname)
if(!is.null(args$type)){plotcol = c(plotcol,typeuse)}
if("cell_type" %in% colnames(pbmc@meta.data)){plotcol = c(plotcol,"cell_type")}
for(i in plotcol){
  if(i %in% colnames(cds@phenoData@data)){
		colorlist<-c('Classical Monocytes CD302'="#1F77B4",'Classical Monocytes IFI44'="#4E9BD4",'Classical Monocytes NFKBIA'="#7FB6E1",'Classical Monocytes RETN '="#A5C8ED",'Macrophage'="#FF7F0E",'Non-classical Monocytes'="#9467BD")
		if(i=='group'){colorlist<-c('control'='#2F6BFF','case'='#E74C3C')}
    monocleplot(cds,colnames = i,ptsize=1,colour=colorlist)
  }else{
    print(paste0(i," is not exist in colnames"))
  }
  
}
plotcol = c("sample","State")
for(i in plotcol){
  if(i %in% colnames(cds@phenoData@data)){
    monocleplot(cds,colnames = i,ptsize=1,colour='F')
  }else{
    print(paste0(i," is not exist in colnames"))
  }
  
}

p = plot_cell_trajectory(cds,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)
p = plot_complex_cell_trajectory(cds,x=1,y=2,color_by = clustercolname)
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)
ncols = ceiling(sqrt(length(unique(cds@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(cds@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(cds,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)


pseudotime_df = cds@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"monocle.metadata.xls",sep = "\t",quote = F,row.names = F)

# outpath = "D:\\单细胞项目\\14.FWSC20241735\\3.重新细分\\Pericyte_Harmony\\monocle"
# cds = readRDS(paste0(outpath,"/pseudotime.rds"))
if(BEAM){
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
}

#结果报告
file.copy("/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/pseudotime/report.rmd",paste0(outpath,"/report.rmd"))
rmarkdown::render(paste0(outpath,"/report.rmd"),output_file=paste0(outpath,"/report.html"))
file.remove(paste0(outpath,"/report.rmd"))
