setwd('/mnt/NexData_file/home/zwl/Project/FWSC20250843/require5/T_cell/')
rm(list = ls())
gc()

library(sceasy)
library(Seurat)
library(reticulate)
#############颜色#################
library(scales)
pal= ggsci::pal_d3('category20')(16)
#names(pal) <- 
show_col(pal)

#############h5ad转变成seurat对象####################

sceasy::convertFormat('./CD8T.h5ad', from="anndata", to="seurat",
                      outFile='CD8T_raw.rds')

################导入注释信息####################

obj <- readRDS(file = './CD8T_raw.rds')
obj <- FindVariableFeatures(obj,nfeatures = 2000)
obj$cell_type = gsub('8','2',gsub('1','2',obj$leiden_res_0.8))
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('leiden_res_0.8','cell_type'))

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
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 5))

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
save(file = 'obj_reduce2.RData',HSMM,ordering_genes)
##排序##
HSMM <- orderCells(HSMM)

##重新排序##
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)
saveRDS(file = 'obj_order.rds',HSMM_myo,compress = F)

##颜色##
HB_col = c('#FF7F0EFF','#2CA02CFF','#9467BDFF','#E377C2FF')
names(HB_col) = unique(HSMM_myo$cell_type)[c(1,2,4,3)]

#################################可视化###################################
p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
  scale_color_manual(values = HB_col)
ggsave(filename = "2.obj_trajectory_cell_type.png",p,w=6,h=6,dpi = 100)

p1 = plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
ggsave(filename = "2.obj_trajectory_Pseudotime.png",p1,w=4,h=4,dpi = 100)

p0 = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
  scale_color_manual(values = HB_col)+facet_wrap(~group)
ggsave(filename = "2.obj_trajectory_cell_type_group.png",p0,w=6,h=5,dpi = 100)

p = plot_cell_trajectory(HSMM_myo, color_by = "State")
ggsave(filename = "2.obj_trajectory_state.png",p,w=6,h=5,dpi = 100)
##cell_type split##
#HC_bar = colnames(HSMM_myo)[HSMM_myo$group=='NC']
p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
  scale_color_manual(values = HB_col)+facet_grid(group~cell_type)
ggsave(filename = "2.obj_trajectory_cell_type_group2.png",p,w=15,h=5,dpi = 100)


##########################4 拟时相关的差异gene############################
##将ordergene提出来做回归，找与时间是否显著相关##
##否则会计算全部gene##
diff_test_res <- differentialGeneTest(HSMM_myo[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = './obj_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "2.obj_top50gene_cortime.png",p,h=8,w=6,dpi = 150)

pie_p = data.frame(table(HSMM_myo$State,HSMM_myo$cell_type,HSMM_myo$group))
pie_p1 = pie_p[apply(pie_p!= 0 , 1 , all),]
pie_p1$Var1 = factor(as.character(pie_p1$Var1),levels = c('1','2','4','5'))
p = ggplot(pie_p1, aes(x=Var3,y=Freq,fill=Var2)) + geom_col(position = 'fill') +
  facet_wrap(~ Var1,nrow = 2)+
  scale_fill_manual(values = HB_col)
ggsave(filename = './monocle_obj_Mega_State_pro.png',w=8,h=6,p)
# ggplot(data=pie_p1, aes(x=" ", y=Freq, group=Var2, colour=Var2, fill=Var2)) +
#   geom_bar(width = 1, stat = "identity") +
#   coord_polar("y", start=0) + 
#   facet_grid(.~ Var1) +theme_void()
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
##########################与时间相关的基因曲线图
s.genes <- c('Pdcd1','Entpd1','Tox','Havcr2','Lag3','Tigit',
             'Tcf7','Slamf6','Cd28','Icos','Mki67')
p1 <- plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "State")
p2 = plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "cell_type")
p1|p2
ggsave(filename = paste0('./gene_in_pseudotime.png'),p1|p2,w=8,h=8)

###########################seurat###############
obj@meta.data = HSMM_myo@phenoData@data[rownames(obj@meta.data),]
Seurat::DimPlot(obj,reduction = 'umap',group.by = c('State','cell_type'))
Idents(obj) = paste0(obj$State,'-',obj$sample)
VlnPlot(object = obj,features = 'Pdcd1',group.by = 'sample',split.by = 'cell_type')
obj$PDCD1_expression <- as.numeric(GetAssayData(obj, layer = "data")["Pdcd1", ])

forplot = obj@meta.data[,c('PDCD1_expression','group','State')] %>% dplyr::filter(State=='4') 
forplot$score = 'Pdcd1'

p_value = forplot%>%
  group_by(score) %>% 
  rstatix::t_test(formula = PDCD1_expression~group)%>%
  add_significance() %>% add_xy_position(x = 'score')%>%
  dplyr::mutate(p.signif = ifelse(p.signif == '','ns',p.signif))

gg = ggplot(forplot, aes(score, PDCD1_expression, colour = group)) +
  geom_boxplot(aes(fill = group),position = position_dodge(width = 0.8), outlier.size = 0, width = 0.2, show.legend = FALSE) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", alpha = 0.5, width = 0.6, size = 0.8) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size=0.5,stroke=0.5, shape=16) +
  scale_color_manual(values = c('CKO'="#ff7f0e", 'WT'="#1f77b4"))+
  scale_fill_manual(values = c('CKO'="#ff7f0e", 'WT'="#1f77b4"))+
  stat_pvalue_manual(p_value, label = "p.signif",hide.ns = T)+
  theme_classic()+labs(y=i,x=NULL)
gg
