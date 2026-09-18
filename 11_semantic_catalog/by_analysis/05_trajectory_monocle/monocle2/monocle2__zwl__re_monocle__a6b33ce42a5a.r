setwd('/mnt/NexData_file/home/zwl/Project/SC2511146601/require5/1.monocle')
library(monocle)
HSMM = readRDS(paste0('/mnt/NexData_file/home/zwl/Project/SC2511146601/require5/1.monocle/monocle.rds'))
deg <- read.csv(paste0('/mnt/NexData_file/home/zwl/Project/SC2511146601/require5/1.monocle/monocle.DEG.xls'),sep="\t")
rownames(deg) = deg$Genes
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
# 重新排序
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 4)

# 重新画图

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
Time_diff <- differentialGeneTest(HSMM_myo[ordergene,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <-Time_diff[, c(5, 3, 4)] # 把gene放前面，也可以不改
Time_genes <- Time_diff %>% filter(pval < 0.05) %>% arrange(qval) %>% distinct(gene_short_name, .keep_all = TRUE)
write.csv(Time_genes, paste0("./3.Pseudotime_diff_all.csv"), row.names = F)
Time_genes <- Time_genes %>% head(10) %>% pull(gene_short_name) %>% as.character() # 提取Top10gene
# 直接查看一些基因随细胞状态等的表达变化,选择前10个top基因并将其对象取出
seuCDS_subset <- HSMM_myo[Time_genes,]
# 可视化：以state/cluster/pseudotime进行
gg=plot_genes_in_pseudotime(seuCDS_subset, color_by = "group", ncol = 5)
ggsave(paste0("./3.pseudotime_genes_group.png"),limitsize = FALSE,width= 20,height= 10)
ggsave(paste0("./3.pseudotime_genes_group.pdf"),limitsize = FALSE,width= 20,height= 10)
