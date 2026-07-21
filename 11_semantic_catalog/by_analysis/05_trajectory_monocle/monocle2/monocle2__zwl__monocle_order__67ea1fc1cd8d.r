setwd('/home/zwl/Project/FW2024-428-02/7.cytotrace/')
library(monocle)
load('./Epithelial_reduce.RData')

##order monocle
HSMM <- orderCells(HSMM)

##plot trace and reassign the root##
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~celltype2,ncol = 4)
#HSMM_myo <- orderCells(HSMM, root_state = 3)
saveRDS(file = 'HSPC_B_unre_order.rds',HSMM_myo,compress = F)
HSMM_myo <- readRDS('./HSPC_B_unre_order.rds')

##use ggsci to chose color,because the library have many categorical color
library(ggsci)
pal= ggsci::pal_d3('category20')(16)
show_col(pal)
col = pal[1:length(unique(HSMM_myo$celltype2))]
names(col) = unique(HSMM_myo$celltype2)

#################################visualization###################################
p = plot_cell_trajectory(HSMM_myo, color_by = "celltype2")
  #scale_color_manual(values = col)
ggsave(filename = "2.Epithelium_trajectory_celltype.png",p,w=6,h=6,dpi = 100)

p1 = plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
ggsave(filename = "2.Epithelium_trajectory_Pseudotime.png",p1,w=4,h=4,dpi = 100)

p0 = plot_cell_trajectory(HSMM_myo, color_by = "celltype2")+
  facet_wrap(~group)
ggsave(filename = "2.Epithelium_trajectory_celltype_group.png",p0,w=6,h=5,dpi = 100)

##celltype split##
p = plot_cell_trajectory(HSMM_myo, color_by = "celltype2")+
  facet_grid(group~celltype2)
ggsave(filename = "2.Epithelium_trajectory_celltype_group2.png",p,w=15,h=5,dpi = 100)



##########################4 Differential genes associated with Pseudotime############################
##select order gene to regression，find them wither associate with time##
##if not select ,it will calculate all gene##
library(dplyr)
diff_test_res <- differentialGeneTest(HSMM_myo[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::arrange(qval) %>% head(50) %>% row.names() 
write.csv(file = './Epithelium_Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                            num_clusters = 6,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = "2.Epithelium_top50gene_cortime.png",p,h=8,w=6,dpi = 150)


