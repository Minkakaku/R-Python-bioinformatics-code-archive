setwd('/mnt/NexData_file/home/zwl/Project/FW2023-179/require4')
rm(list = ls())
gc()

library(sceasy)
library(Seurat)
library(reticulate)

library(scales)


#############h5ad转变成seurat对象####################
sceasy::convertFormat('../B_celltype/Bcell_raw.h5ad', from="anndata", to="seurat",
                      outFile='Bcell_raw.rds')

################导入注释信息####################

obj <- readRDS(file = './adata_BP.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('group','cell_type'))
p = DotPlot(obj,features = c('CCR10','CCR9','CCR7','CCR6','CXCR5','CXCR4','CXCR3','CX3CR1'),
            group.by = 'group',cols = 'RdBu')+
  ggplot2::coord_flip()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        plot.background = element_rect(fill = 'white'))
ggsave(filename = "../B_CXC_CCR_dotplot.png",p,w=5,h=6,dpi = 150)
ggsave(filename = "../B_CXC_CCR_dotplot.pdf",p,w=5,h=6,dpi = 150)
###################split###########################
# cellchat_sp = SplitObject(obj,split.by = "group")

################1 构建对象#################
library(dplyr)
library(monocle)
f_data = data.frame(gene_short_name = rownames(obj),
                    row.names = rownames(obj))
count = as.matrix(LayerData(obj, assay="RNA", layer='counts'))

pd <- new("AnnotatedDataFrame", data = obj@meta.data)
fd <- new("AnnotatedDataFrame", data = f_data)
HSMM <- newCellDataSet(count, phenoData = pd, featureData = fd,
                       lowerDetectionLimit = 0.5,expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(Biobase::fData(HSMM),num_cells_expressed >= 5))

################2挑选celltype差异gene作为特征#####################
##2000左右##
diff = differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~cell_type")
ordering_genes <- row.names(subset(diff, qval < 0.05))

HSMM <- setOrderingFilter(HSMM, ordering_genes[1:2000])
plot_ordering_genes(HSMM)


##降维##
HSMM <- reduceDimension(HSMM, max_components = 2,
                        reduction_method = 'DDRTree')
#residualModelFormulaStr = "~sample"
save(file = 'obj_reduce2.RData',HSMM,ordering_genes,diff)
##排序##
HSMM <- orderCells(HSMM)

##重新排序##
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo = HSMM
#HSMM_myo <- orderCells(HSMM, root_state = 3)
saveRDS(file = 'obj_order.rds',HSMM_myo,compress = F)

##颜色##
HB_col = c('#ea7070', '#fdc4b6', '#e59572', '#2694ab', '#96ceb4', '#ffeead')
names(HB_col) = unique(obj$cell_type)

#################################可视化###################################
p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
  scale_color_manual(values = HB_col)
ggsave(filename = "2.obj_trajectory_cell_type.png",p,w=6,h=6,dpi = 100)
ggsave(filename = "2.obj_trajectory_cell_type.pdf",p,w=6,h=6,dpi = 100)

p1 = plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
ggsave(filename = "2.obj_trajectory_Pseudotime.png",p1,w=4,h=4,dpi = 100)
ggsave(filename = "2.obj_trajectory_Pseudotime.pdf",p1,w=4,h=4,dpi = 100)

p0 = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
  scale_color_manual(values = HB_col)+facet_wrap(~group)
ggsave(filename = "2.obj_trajectory_cell_type_group.png",p0,w=6,h=5,dpi = 100)
ggsave(filename = "2.obj_trajectory_cell_type_group.pdf",p0,w=6,h=5,dpi = 100)

p = plot_cell_trajectory(HSMM_myo, color_by = "State")
ggsave(filename = "2.obj_trajectory_state.png",p,w=6,h=5,dpi = 100)
ggsave(filename = "2.obj_trajectory_state.pdf",p,w=6,h=5,dpi = 100)
##cell_type split##
#HC_bar = colnames(HSMM_myo)[HSMM_myo$group=='NC']
p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
  scale_color_manual(values = HB_col)+facet_grid(group~cell_type)
ggsave(filename = "2.obj_trajectory_cell_type_group2.png",p,w=9,h=9,dpi = 100)
ggsave(filename = "2.obj_trajectory_cell_type_group2.pdf",p,w=9,h=9,dpi = 100)

p = plot_cell_trajectory(HSMM_myo, color_by = "group")
ggsave(filename = "2.obj_trajectory_group.png",p,w=6,h=5,dpi = 100)
ggsave(filename = "2.obj_trajectory_group.pdf",p,w=6,h=5,dpi = 100)
##########################4 拟时相关的差异gene############################
##将ordergene提出来做回归，找与时间是否显著相关##
##否则会计算全部gene##
diff_test_res <- differentialGeneTest(HSMM_myo[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = './2.obj_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 3,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "2.obj_top50gene_cortime.png",p,h=8,w=6,dpi = 150)
ggsave(filename = "2.obj_top50gene_cortime.pdf",p,h=8,w=6,dpi = 150)

pie_p = data.frame(table(HSMM_myo$State,HSMM_myo$cell_type,HSMM_myo$group))
pie_p1 = pie_p[apply(pie_p!= 0 , 1 , all),]
pie_p1$Var1 = factor(as.character(pie_p1$Var1),levels = c('5','4','3','6','2','7','1'))
p = ggplot(pie_p1, aes(x=Var3,y=Freq,fill=Var2)) + geom_col(position = 'fill') +
  facet_wrap(~ Var1,nrow = 2)+
  scale_fill_manual(values = HB_col)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  labs(x=NULL,y=NULL,title = 'Percentage of celltype in different State')
ggsave(filename = './monocle_obj_Mega_State_pro.png',w=10,h=6,p)
ggsave(filename = './monocle_obj_Mega_State_pro.pdf',w=10,h=6,p)



############################BEAM分支变化#############################
HSMM_myo <- readRDS('./obj_order.rds')

BEAM_res1 <- BEAM(HSMM_myo, branch_point = 3, cores = 1)
BEAM_res1 <- BEAM_res1[order(BEAM_res1$qval),]
write.csv(file = 'BEAM_Megakaryocyte_Erythrocyte.csv',BEAM_res1)
BEAM_res1 = read.csv(file = './BEAM_Megakaryocyte_Erythrocyte.csv')
BEAM_res <- BEAM_res1 %>% dplyr::filter(num_cells_expressed >10,qval <1e-4) %>%
  dplyr::select(gene_short_name,pval,qval)

#state3:巨核细胞
#state4：红细胞
A = plot_genes_branched_heatmap(HSMM_myo[BEAM_res$gene_short_name,],
                                branch_point = 3,
                                num_clusters = 4,
                                cores = 1,
                                return_heatmap=T,
                                use_gene_short_name = T,
                                show_rownames = F)

png(filename = 'branch_analysis.png',w=1000,h=1800,res = 200)
plot_genes_branched_heatmap(HSMM_myo[BEAM_res$gene_short_name[1:100],],
                            branch_point = 3,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()
lung_genes <- row.names(subset(fData(HSMM_myo),
                               gene_short_name %in% c("IFNGR1", "IFNAR1")))
plot_genes_branched_pseudotime(HSMM_myo[s.genes,],
                               branch_point = 1,
                               color_by = "cell_type",
                               ncol = 1)
