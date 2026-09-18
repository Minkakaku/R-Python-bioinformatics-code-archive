library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(clusterProfiler)
library(purrr)
library(grid)
# library(ggplotify)
library(patchwork)

source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20250070/02.Mission/mapping_recluster_Score/"
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20250070/01.Mission/mapping_recluster/mapping_recluster.rds")
# outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20250070/02.Mission/PTC cells_Score/"
# adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20250070/01.Mission/PTC cells/PTC cells.rds")
if(!dir.exists(outdir)) {dir.create(outdir, recursive = T)}
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
Idents(adata)<-factor(Idents(adata),levels=c("T cells","NK cells","Neutrophils","Macrophages","B cells","Proximal tubule cells","Podocyte","Myofibroblast","Loop of henle","Distal tuble cells","Collecting duct_intercalated cells","Collecting duct_principal cells","EC_1","EC_2"))

kegg_geneset <- clusterProfiler::read.gmt("/home/sh/sc_pipline/geneset_score/cell death.gmt")
kegg_geneset_list <- split(kegg_geneset$gene, kegg_geneset$term)

adata<-AddModuleScore(
    object=adata,
    features=kegg_geneset_list,
    ctr1=100, #默认值是100
    name='Score'
)
colnames(adata@meta.data)[(ncol(adata@meta.data) - length(kegg_geneset_list) + 1):ncol(adata@meta.data)] <- names(kegg_geneset_list)
mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2",names(kegg_geneset_list),"group"))
groups <- levels(mydata$group)
titles <- setNames(groups, groups)
write.table(adata@meta.data,paste0(outdir,"meatdata_AddModuleScore.csv"),sep=",",quote=F,row.names=F,col.names=T)

create_group_plot <- function(data, group_name, color_var) {
  ggplot(data[data$group == group_name, ], aes(x = UMAP_1, y = UMAP_2, colour = !!sym(color_var))) +
    geom_point(size = 1) +
    scale_color_gradientn(colours = c('lightgrey', '#8cd0e1', '#f1681e', '#eb2d4a', '#fab6f3')) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
		legend.position = 'none') +
    labs(title = group_name)
}

for(m in 1:length(kegg_geneset_list)){
prefix<-names(kegg_geneset_list[m])
p <- VlnPlot(adata,features=prefix,pt.size = 0) + geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
ggsave(paste0(outdir,prefix,".VlnPlot.pdf"),width = 12,height = 10)
ggsave(paste0(outdir,prefix,".VlnPlot.png"),width = 12,height = 10)

p <- VlnPlot(adata,features=prefix,split.by = 'group',pt.size = 0) + geom_boxplot(width=0.2,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,prefix,"_group.VlnPlot.pdf"),width = 18,height = 10)
ggsave(paste0(outdir,prefix,"_group.VlnPlot.png"),width = 18,height = 10)
# head(ggplot_build(p)$data[[1]])
p <- get_sig_ana(p,as.formula(paste0(prefix,'~split')),group.by='ident',prefix=paste0(prefix,"_group"),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(aki='#F8766D',cdn='#078992',control='#E38900'))
get_sig_ana(p$data,as.formula(paste0(prefix,'~split')),group.by='ident',prefix=paste0(prefix,"_group"),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,prefix,"_group_sig_ana.VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,prefix,"_group_sig_ana.VlnPlot.png"),width = 15,height = 10)


p <- ggplot(mydata, aes(x = UMAP_1, y = UMAP_2, colour = !!sym(prefix))) +
    geom_point(size = 1) +
    scale_color_gradientn(colours = c('lightgrey', '#8cd0e1', '#f1681e', '#eb2d4a', '#fab6f3')) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
		legend.position = 'right',plot.margin=margin(1,2,1,1,unit='cm')) +
    ggtitle(prefix)
ggsave(paste0(outdir,prefix,".UMAP.png"), plot = p, width = 10, height = 6)
ggsave(paste0(outdir,prefix,".UMAP.pdf"), plot = p, width = 10, height = 6)

unique_groups <- unique(mydata$group)
plots <- lapply(unique_groups, function(group_name) create_group_plot(mydata, group_name, prefix))
combined_plot <- do.call(plot_grid, c(plots, nrow = 1))
ggsave(paste0(outdir,prefix,"_group.UMAP.png"), plot = combined_plot, width = 18, height = 6)
ggsave(paste0(outdir,prefix,"_group.UMAP.pdf"), plot = combined_plot, width = 18, height = 6)
}

