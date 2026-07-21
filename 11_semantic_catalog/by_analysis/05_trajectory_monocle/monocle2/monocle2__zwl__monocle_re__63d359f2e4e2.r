################monocle###########################
library(monocle)
rm(list = ls())
gc()
HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2604092832-01/require2/monocle/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2604092832-01/require2/monocle/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 2)

# 重新画图
setwd('/mnt/NexData_file/home/zwl/Project/SC2604092832-01/require2/monocle')
getwd()
clustercolname = 'cell_type'

p = plot_cell_trajectory(HSMM_myo,color_by="Pseudotime",size=1,show_backbone=TRUE)
ggsave("monocle.trajectory.pseudotime.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.pseudotime.png",p,width = 7,height = 7.5)

p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)
ggsave("monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("monocle.trajectory.tree.png",p,width = 7,height = 7.5)
ncols = ceiling(sqrt(length(unique(HSMM_myo@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(HSMM_myo@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(HSMM_myo,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
ggsave("monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)

pseudotime_df = HSMM_myo@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"3.monocle.metadata_reorder.xls",sep = "\t",quote = F,row.names = F)

saveRDS(file = '3.obj_reorder.rds',HSMM_myo,compress = F)


# 寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图。
diff_test_res <- differentialGeneTest(HSMM_myo[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = '3.obj_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "3.monocle.top50gene_cortime.png",p,h=8,w=6,dpi = 150)
ggsave(filename = "3.monocle.top50gene_cortime.pdf",p,h=8,w=6,dpi = 150)

Time_genes <- diff_test_res %>% head(10) %>% pull(gene_short_name) %>% as.character() # 提取Top10gene
# 直接查看一些基因随细胞状态等的表达变化,选择前10个top基因并将其对象取出
HSMM_myo_subset <- HSMM_myo[Time_genes,]
# 可视化：以state/cluster/pseudotime进行
gg=plot_genes_in_pseudotime(HSMM_myo_subset, color_by = "group", ncol = 5)
ggsave(filename = "3.pseudotime_genes_group.png",gg,width= 20,height= 10)
ggsave(filename = "3.pseudotime_genes_group.pdf",gg,width= 20,height= 10)
