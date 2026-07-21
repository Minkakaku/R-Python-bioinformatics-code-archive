setwd("D:/0_project/FW2024-225-G/require5/")
rm(list = ls())
gc()

library(sceasy)
library(Seurat)
library(reticulate)
use_condaenv('myname')
#############颜色#################
library(scales)
pal= ggsci::pal_d3('category20')(16)
#names(pal) <- 
show_col(pal)

#############h5ad转变成seurat对象####################

sceasy::convertFormat('D:/0_project/FW2024-225-G/require5/HSPC_redsc.h5ad', from="anndata", to="seurat",
                      outFile='HSPC_redsc.rds')

################导入注释信息####################

HSPC_red <- readRDS("D:/0_project/FW2024-225-G/require5/HSPC_redsc.rds")
HSPC_red = Seurat::FindVariableFeatures(HSPC_red,nfeatures = 2000)
metadata <-read.csv(file = 'D:/0_project/FW2024-225-G/subtype/3.ALL_subtype_metadata.csv',row.names = 1)

HSPC_red$celltype = metadata[colnames(HSPC_red),'sub_celltype2']
HSPC_red$group = HSPC_red$disease
HSPC_red = subset(HSPC_red,celltype %in% c("Erythroid progenitor","HSC/MPP","GMP/LMPP",#"Megakaryocyte","Erythrocyte",
                                           "MkP","MEP"))
table(HSPC_red$celltype)

###################split###########################
cellchat_sp = SplitObject(HSPC_red,split.by = "group")

################1 构建对象#################
library(dplyr)
library(monocle)
f_data = data.frame(gene_short_name = rownames(HSPC_red),
                    row.names = rownames(HSPC_red))
count = as.matrix(LayerData(HSPC_red, assay="RNA", layer='counts'))

pd <- new("AnnotatedDataFrame", data = HSPC_red@meta.data)
fd <- new("AnnotatedDataFrame", data = f_data)
HSMM <- newCellDataSet(count, phenoData = pd, featureData = fd,
                       lowerDetectionLimit = 0.5,expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

################2挑选celltype差异gene作为特征#####################
##2000左右##
diff = differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~celltype")
ordering_genes <- row.names(subset(diff, qval < 0.01))

HSMM <- setOrderingFilter(HSMM, ordering_genes[1:2000])
plot_ordering_genes(HSMM)


##降维##
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree',
                        residualModelFormulaStr = "~sample")
save(file = 'HSPC_red_reduce2.RData',HSMM,ordering_genes)
##排序##
HSMM <- orderCells(HSMM)

##重新排序##
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~celltype,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 1)
saveRDS(file = 'HSPC_red_order.rds',HSMM_myo,compress = F)

##颜色##
HB_col = pal[1:length(unique(HSMM_myo$celltype))]
names(HB_col) = unique(HSMM_myo$celltype)

#################################可视化###################################
p = plot_cell_trajectory(HSMM_myo, color_by = "celltype")+
  scale_color_manual(values = HB_col)
ggsave(filename = "2.HSPC_red_trajectory_celltype.png",p,w=6,h=6,dpi = 100)

p1 = plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
ggsave(filename = "2.HSPC_red_trajectory_Pseudotime.png",p1,w=4,h=4,dpi = 100)

p0 = plot_cell_trajectory(HSMM_myo, color_by = "celltype")+
  scale_color_manual(values = HB_col)+facet_wrap(~disease)
ggsave(filename = "2.HSPC_red_trajectory_celltype_group.png",p0,w=6,h=5,dpi = 100)

p = plot_cell_trajectory(HSMM_myo, color_by = "State")
ggsave(filename = "2.HSPC_red_trajectory_state.png",p,w=6,h=5,dpi = 100)
##celltype split##
#HC_bar = colnames(HSMM_myo)[HSMM_myo$disease=='NC']
p = plot_cell_trajectory(HSMM_myo, color_by = "celltype")+
  scale_color_manual(values = HB_col)+facet_grid(disease~celltype)
ggsave(filename = "2.HSPC_red_trajectory_celltype_group2.png",p,w=15,h=5,dpi = 100)


pie <- ggplot(df,aes(x = '', fill = factor(celltype))) +
  facet_wrap(~State) + 
  geom_bar(width = 1) + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


##########################4 拟时相关的差异gene############################
##将ordergene提出来做回归，找与时间是否显著相关##
##否则会计算全部gene##
setwd('D:/0_project/FW2024-225-G/require5/ALL/')
diff_test_res <- differentialGeneTest(HSMM_myo[ordering_genes[1:2000],],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = './HSPC_red_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "2.HSPC_red_top50gene_cortime.png",p,h=8,w=6,dpi = 150)

pie_p = data.frame(table(HSMM_myo$State,HSMM_myo$celltype,HSMM_myo$disease))
pie_p1 = pie_p[apply(pie_p!= 0 , 1 , all),]
pie_p1$Var1 = factor(as.character(pie_p1$Var1),levels = c('7','1','2','3','4','6','5'))
p = ggplot(pie_p1, aes(x=Var3,y=Freq,fill=Var2)) + geom_col(position = 'fill') +
  facet_wrap(~ Var1,nrow = 2)+
  scale_fill_manual(values = HB_col)
ggsave(filename = './monocle_HSPC_red_Mega_State_pro.png',w=8,h=6,p)
# ggplot(data=pie_p1, aes(x=" ", y=Freq, group=Var2, colour=Var2, fill=Var2)) +
#   geom_bar(width = 1, stat = "identity") +
#   coord_polar("y", start=0) + 
#   facet_grid(.~ Var1) +theme_void()
############################BEAM分支变化#############################
HSMM_myo <- readRDS('./HSPC_red_order.rds')

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
                               color_by = "celltype",
                               ncol = 1)
##########################与时间相关的基因曲线图
s.genes <- c("SELL","CCR7","IL7R", "CD84","CCL5","CD9","S100A4",
             'KLF2','TNFAIP3','EREG','CTSW','MIF','HBG2')
p1 <- plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "State")
p2 = plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "celltype")
p1|p2
ggsave(filename = paste0('./gene_in_pseudotime.png'),p1|p2,w=8,h=8)

s.genes <- c("CD9","CD34")
p1 <- plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "State")
p2 = plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "celltype")
p1|p2

s.genes <- c("KLF2","TNFAIP3",'EREG','CTSW','MIF','HBG2')
p1 <- plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "State")
p2 = plot_genes_in_pseudotime(HSMM_myo[s.genes,],  color_by = "celltype")
p1|p2
ggsave(filename = paste0('./gene_in_pseudotime1.png'),p1|p2,w=8,h=8)
##########################filter后的gene 富集分析##################################
library(clusterProfiler)
library(org.Hs.eg.db)

df = data.frame(A[["annotation_row"]])
res = split(rownames(df), df$Cluster)
for (i in 1:length(res)) {
  ego_CC <- enrichGO(gene          = res[[i]],
                     OrgDb         = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  p = barplot(ego_CC) + theme(panel.grid=element_blank()) #富集柱形图
  ggsave(filename = paste0('./Cluster',i,'GO_CC_enrichment.png'),p,w=8,h=6)
  write.csv(file = paste0('./Cluster',i,'GO_CC_enrichment.csv'),ego_CC@result)
  
  ego_MF <- enrichGO(gene          = res[[i]],
                     OrgDb         = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  p = barplot(ego_MF) + theme(panel.grid=element_blank())  #富集柱形图
  ggsave(filename = paste0('./Cluster',i,'GO_MF_enrichment.png'),p,w=8,h=6)
  write.csv(file = paste0('./Cluster',i,'GO_MF_enrichment.csv'),ego_MF@result)
  
  ego_BP <- enrichGO(gene          = res[[i]],
                     OrgDb         = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  p = barplot(ego_BP) + theme(panel.grid=element_blank())  #富集柱形图
  ggsave(filename = paste0('./Cluster',i,'GO_BP_enrichment.png'),p,w=8,h=6)
  write.csv(file = paste0('./Cluster',i,'GO_BP_enrichment.csv'),ego_BP@result)
}



