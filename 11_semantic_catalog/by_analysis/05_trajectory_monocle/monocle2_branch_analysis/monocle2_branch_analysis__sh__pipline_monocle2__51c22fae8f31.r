####内皮拟时序####
library(Seurat)
library(monocle)
library(dplyr)
library(reshape2)
library(cowplot)

S <- readRDS("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/Neurons/subcluster.rds")
S<-subset(S,idents=c("Neuroblast"))

# 通过seurat结果导入数据
data <- as(as.matrix(S@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = S@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#构建S4对象，CellDataSet
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size() # 单细胞稀疏矩阵的话用negbinomial.size()，如果是UMI的话不要标准化；FPKM值用tobit()；logFPKM值用gaussianff()
) 

# 估计size factors 和 dispersions
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# 过滤低质量细胞
HSMM <- detectGenes(HSMM, min_expr = 3)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
head(pData(HSMM))


##选择发育特征基因
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~group") #~后是对谁做差异分析的变量，可以为pd的任意列名

# 找出那些在不同细胞类型之间有显著差异表达的基因，以用于细胞排序
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# ordering_genes_top <- row.names(diff_test_res)[order(diff_test_res$qval)][1:400]  #  选择top2000的基因 没用上

HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM_top <- setOrderingFilter(HSMM, ordering_genes)

plot_ordering_genes(HSMM_top)


# 降维 & 排序
HSMM_top <- reduceDimension(HSMM_top, max_components = 2,
#                            residualModelFormulaStr = as.formula("~sample"),
                            reduction_method = 'DDRTree') # DDRTree方式

HSMM_top <- orderCells(HSMM_top)

saveRDS(HSMM_top,"/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_monocle2/Neurons_monocle2.RDS")


# 函数名: PsePlot
# 参数: CDS - monocle_cds对象
# 返回值: 画图
# 功能: 画各种图
PsePlot = function(HSMM_top){
  gg =  plot_cell_trajectory(HSMM_top, color_by = "State")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.State.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.State.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(HSMM_top, color_by = "cell_type")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.cell_type.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.cell_type.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(HSMM_top, color_by = "Pseudotime")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.Pseudotime.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.Pseudotime.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(HSMM_top, color_by = "cell_type") + 
    facet_wrap(~cell_type, nrow = floor(sqrt(length(unique(HSMM_top$cell_type)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.cell_typeEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.cell_typeEach.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(HSMM_top, color_by = "State") +
    facet_wrap(~State, nrow = floor(sqrt(length(unique(HSMM_top$State)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.StateEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.StateEach.pdf"),limitsize = FALSE)
  gg = plot_complex_cell_trajectory(HSMM_top, color_by = 'cell_type', 
                                    show_branch_points = T, cell_size = 0.5, 
                                    cell_link_size = 0.3)+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.tree.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.tree.pdf"),limitsize = FALSE)
}
# 函数名: PsePlotSample
# 参数: CDS - monocle_cds对象
# 返回值: 画样本图
# 功能: 画样本图
PsePlotSample = function(HSMM_top){
  if(("sample"%in% colnames(HSMM_top@phenoData@data))|("orig.ident"%in% colnames(HSMM_top@phenoData@data))){
    if("orig.ident"%in% colnames(HSMM_top@phenoData@data)){
      HSMM_top$Sample = as.factor(HSMM_top$orig.ident)
    }else {
      HSMM_top$Sample = as.factor(HSMM_top$sample)     
    }
  gg =  plot_cell_trajectory(HSMM_top, color_by = "Sample")+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.Sample.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.Sample.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(HSMM_top@phenoData@data)){
  gg =  plot_cell_trajectory(HSMM_top, color_by = "group")+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.group.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.group.pdf"),limitsize = FALSE)
  }
  if("Sample"%in% colnames(HSMM_top@phenoData@data)){
  gg =  plot_cell_trajectory(HSMM_top, color_by = "cell_type") + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(HSMM_top$Sample)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.SampleEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.SampleEach.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(HSMM_top@phenoData@data)){
  gg =  plot_cell_trajectory(HSMM_top, color_by = "cell_type") + facet_wrap(~group, nrow = floor(sqrt(length(unique(HSMM_top$group)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.groupEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.groupEach.pdf"),limitsize = FALSE)
  }
  if("Sample"%in% colnames(HSMM_top@phenoData@data)){
  gg = plot_complex_cell_trajectory(HSMM_top, color_by = 'cell_type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(HSMM_top$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(legend.position="right",axis.text.x = element_text(angle = 30, hjust = 1))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("SampleEach.tree.png"),limitsize = FALSE)
  ggsave(paste0("SampleEach.tree.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(HSMM_top@phenoData@data)){
  gg = plot_complex_cell_trajectory(HSMM_top, color_by = 'cell_type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~group, nrow = floor(sqrt(length(unique(HSMM_top$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(legend.position="right",axis.text.x = element_text(angle = 30, hjust = 1))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("groupEach.tree.png"),limitsize = FALSE)
  ggsave(paste0("groupEach.tree.pdf"),limitsize = FALSE)
  }
}
PseDenPlot = function(forplot){
  if("Sample" %in%colnames(forplot)){
    if(length(unique(forplot$Sample))>1){
    nrow = floor(sqrt(length(unique(forplot$Sample))))
    ncol = ceiling(length(unique(forplot$Sample))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=Sample))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~Sample,nrow = nrow)
    ggsave("PseudotimeDensity.SampleEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    ggsave("PseudotimeDensity.SampleEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=Sample))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.Sample.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.Sample.pdf",gg,width= 8,height=5,limitsize = FALSE)
    }else{
    # plotSampledata=ggplot_build(gg)$data[[1]]
    # write.table(plotSampledata,file="../DCPlotForVis/orig.ident.DCPlot.txt",sep="\t",row.names=F,quot=F)
    nrow = floor(sqrt(length(unique(forplot$cell_type))))
    ncol = ceiling(length(unique(forplot$cell_type))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=cell_type))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~cell_type,nrow = nrow)
    ggsave("PseudotimeDensity.cell_typeEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    ggsave("PseudotimeDensity.cell_typeEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=cell_type))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.cell_type.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.cell_type.pdf",gg,width= 8,height=5,limitsize = FALSE)
    # plotdata=ggplot_build(gg)$data[[1]]
    # write.table(plotdata,file="../DCPlotForVis/Cluster.DCPlot.txt",sep="\t",row.names=F,quot=F)   
    }
  }
  if("group" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$group))))
    ncol = ceiling(length(unique(forplot$group))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=group))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~group,nrow = nrow)
    ggsave("PseudotimeDensity.groupEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    ggsave("PseudotimeDensity.groupEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=group))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.group.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.group.pdf",gg,width= 8,height=5,limitsize = FALSE)
    
    gg = ggplot(forplot,aes(x=Pseudotime))+geom_density(alpha=0.5,fill="blue")+theme_cowplot()
    ggsave("PseudotimeDensity.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  
}


table = table(HSMM_top$cell_type,HSMM_top$State) %>% as.data.frame() %>% 
    dcast(Var1~Var2)
  colnames(table)[1] = "cell_type"
  table %>% write.table("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_monocle2/Pseudotime.Statistics.txt",sep = "\t",quote = F)
  HSMM_top@phenoData@data %>% data.frame(CellName=row.names(HSMM_top@phenoData@data)) %>%
    write.table("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_monocle2/Pseudotime.Summary.txt",sep = "\t",quote = F,row.names = F)
  PsePlot(HSMM_top)
  if (length(unique(HSMM_top$sample))>1|length(unique(HSMM_top$orig.ident))>1) {
    PsePlotSample(HSMM_top)
  }
  if("sample"%in% colnames(HSMM_top@phenoData@data)){
  forplot = data.frame(Sample=HSMM_top$sample,Cluster=HSMM_top$cell_type,Pseudotime=HSMM_top$Pseudotime)
  PseDenPlot(forplot)
  }
  if("orig.ident"%in% colnames(HSMM_top@phenoData@data)){
  forplot = data.frame(Sample=HSMM_top$orig.ident,Cluster=HSMM_top$cell_type,Pseudotime=HSMM_top$Pseudotime)
  PseDenPlot(forplot)
  }
  if("group"%in% colnames(HSMM_top@phenoData@data)){
  forplot = data.frame(group=HSMM_top$group,Cluster=HSMM_top$cell_type,Pseudotime=HSMM_top$Pseudotime)
  PseDenPlot(forplot)

HSMM_top$cell_type<-factor(HSMM_top$cell_type,levels=c("Interneuron","Bergmann glial cell", "Excitatory neuron", "Inhibitory neuron" , "Neuroblast"))
gg =  plot_cell_trajectory(HSMM_top, color_by = "cell_type") +scale_color_manual(values=c("Bergmann glial cell" = "grey", "Excitatory neuron" = "grey", "Inhibitory neuron" = "grey" , "Interneuron" = "grey", "Neuroblast" = "red"))+ 
facet_wrap(~group, nrow = floor(sqrt(length(unique(HSMM_top$group)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.groupEach_Neuroblast.png"),limitsize = FALSE)
ggsave(paste0("Pseudotime.groupEach_Neuroblast.pdf"),limitsize = FALSE)



#可视化
plot_cell_trajectory(HSMM_top, color_by = "Pseudotime")

plot_cell_trajectory(HSMM_top, color_by = "group") + scale_color_manual(values = c("#729ECE","#66CC66","#ED665D"))

plot_cell_trajectory(HSMM_top, color_by = "sample")+ scale_color_manual(values = c("lightgrey","lightgrey","lightgrey",
                                                                                   "#CDCC5D","#FF9E4A","#ED665D",
                                                                                   "lightgrey","lightgrey","lightgrey"))


# 沿时间轴的细胞密度图
# colour = c("#729ECE","#66CC66","#ED665D")
library(ggpubr)
df <- pData(HSMM_top)
ggplot(df, aes(Pseudotime, colour = group, fill = group)) +
  geom_density(bw = 0.5, linewidth = 1,alpha = 0.5) + 
  theme_classic2() +
  scale_fill_manual(name = "", values = colour) +
  scale_color_manual(name= "", values = colour )


# 寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图。
Time_diff <- differentialGeneTest(HSMM_top[ordering_genes,],
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <-Time_diff[, c(5, 2, 3, 4, 1, 6, 7)] # 把gene放前面，也可以不改
#write.csv(Time_diff, "Time_T_cell_all.csv", row.names = F)

Time_genes <- Time_diff %>% arrange(qval) %>% head(50) %>% pull(gene_short_name) %>% as.character() # 提取Top50gene
# gg=plot_pseudotime_heatmap(HSMM_top[Time_genes,], 
#                         cores = 6, 
#                         num_clusters = 3, 
#                         show_rownames = T, return_heatmap = T,
#                         hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
# ggsave(paste0("pseudotime_heatmap.png"),limitsize = FALSE)
# ggsave(paste0("pseudotime_heatmap.pdf"),limitsize = FALSE)


# 直接查看一些基因随细胞状态等的表达变化,选择前4个top基因并将其对象取出
# keygenes <- head(Time_genes, 20)
HSMM_subset <- HSMM_top[Time_genes,]

# 可视化：以state/cluster/pseudotime进行
# gg=plot_genes_in_pseudotime(HSMM_subset, color_by = "sample")
gg=plot_genes_in_pseudotime(HSMM_subset, color_by = "group", ncol = 6)
ggsave(paste0("pseudotime_genes_group.png"),limitsize = FALSE,width= 24,height= 20)
ggsave(paste0("pseudotime_genes_group.pdf"),limitsize = FALSE,width= 24,height= 20)
gg=plot_genes_in_pseudotime(HSMM_subset, color_by = "Pseudotime")
ggsave(paste0("pseudotime_genes.png"),limitsize = FALSE)
ggsave(paste0("pseudotime_genes.pdf"),limitsize = FALSE)


# 也可以自己指定基因
s.genes <-c("Adrb2","Ecm1","Ccn2","Acadl","Ptgs1","Cd24a")

plot_genes_jitter(HSMM_top[s.genes,],grouping = "group", color_by = "group")
plot_genes_violin(HSMM_top[s.genes,],grouping = "State", color_by = "State")
plot_genes_in_pseudotime(HSMM_top[s.genes,], color_by = "group", ncol = 3) + scale_color_manual(values = c("#729ECE","#66CC66","#ED665D"))