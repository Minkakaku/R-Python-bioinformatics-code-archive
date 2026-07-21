setwd('/mnt/NexData_file/home/zwl/Project/FW2024-428-02/require6')
library(ggsci)
# 根据cytotrace寻找起始点
cyto = read.csv('/mnt/NexData_file/home/zwl/Project/FW2024-428-02/require1_3/7.cytotrace/cytotrace_metadata.csv',row.names = 1)
pal= ggsci::pal_d3('category20')(16);       show_col(pal)
col = pal[1:length(unique(cyto$cell_type))]
group = 'cell_type'
potencyBoxplot_byPheno <- ggplot(cyto, aes(x = reorder(!!sym(group), -CytoTRACE2_Score), y = CytoTRACE2_Score)) +
  geom_boxplot(aes(fill = !!sym(group)), width = 0.8, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = col) +
  geom_jitter(width = 0.1, size = 0.05) +  # 添加散点
  theme_classic() +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    sec.axis = sec_axis(
      transform = ~.,
      breaks = seq(0, 1, by = 1/12),
      labels = c("", "Differentiated", "", "Unipotent", "", "Oligopotent", 
                 "", "Multipotent", "", "Pluripotent", "", "Totipotent", "")
    )
  ) +
  labs(x = group, y = "Potency score") +
  ggtitle(sprintf("Developmental potential by %s", group)) +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),  # 添加灰色网格
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # x轴标签旋转90度
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.ticks.y.right = element_line(color = c("black", NA, "black", NA, "black", NA, 
                                                "black", NA, "black", NA, "black", NA, "black")),
    aspect.ratio = 0.8,
    #axis.ticks.length.y.right = unit(0.3, "cm")
  )
ggsave(filename = './1.cytotrace2_potency_score_by_phenotype2.png',potencyBoxplot_byPheno,w = 8,h=6,dpi = 200)
ggsave(filename = './1.cytotrace2_potency_score_by_phenotype2.pdf',potencyBoxplot_byPheno,w = 8,h=6,dpi = 200)


library(monocle)
rm(list = ls())
gc()
load('./Epithelial_reduce.RData')
##order monocle
HSMM <- orderCells(HSMM)

##plot trace and reassign the root##
plot_cell_trajectory(HSMM, color_by = "cell_type")
plot_cell_trajectory(HSMM, color_by = "State")+
  facet_wrap(~cell_type,ncol = 4)
HSMM_myo <- orderCells(HSMM, root_state = 4)#SOX2，LYZ root4；HSPA5 root7

saveRDS(file = 'Epi_order4.rds',HSMM_myo,compress = F)
HSMM_myo <- readRDS('./Epi_order7.rds')

HSMM_myo@phenoData@data[["celltype"]] = gsub('LYZ_Epi','Stem Epi',
                                             gsub('SOX2_Epi','Stem Epi',
                                              gsub('Epithelial cells','Epi',
                                                   HSMM_myo@phenoData@data[["cell_type"]])))
# plot_cell_trajectory(HSMM_myo, color_by = "State")+
#   facet_wrap(~celltype,ncol = 4)
# HSMM_myo <- orderCells(HSMM_myo,root=3)
# saveRDS(file = 'Epi_re_order_LYZ.rds',HSMM_myo,compress = F)
##use ggsci to chose color,because the library have many categorical color
for (ci in c('HSPA5_Epithelial','SOX2_Epithelial','LYZ_Epithelial')) {
  ci = 'SOX2-LYZ_Epithelial'#'HSPA5_Epithelial'
  library(ggsci)
  library(scales)
  pal= ggsci::pal_d3('category20')(16)
  show_col(pal)
  col = pal[1:length(unique(HSMM_myo$cell_type))]
  names(col) = unique(HSMM_myo$cell_type)
  
  #################################visualization###################################
  p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")
  #scale_color_manual(values = col)
  ggsave(filename = paste0("2.Epithelium_trajectory_celltype_",ci,".png"),p,w=6,h=6,dpi = 100)
  
  p1 = plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
  ggsave(filename = paste0("2.Epithelium_trajectory_Pseudotime_",ci,".png"),p1,w=4,h=4,dpi = 100)
  
  p0 = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
    facet_wrap(~group)
  ggsave(filename = paste0("2.Epithelium_trajectory_celltype_group_",ci,".png"),p0,w=6,h=5,dpi = 100)
  
  ##celltype split##
  p = plot_cell_trajectory(HSMM_myo, color_by = "cell_type")+
    facet_grid(group~cell_type)
  ggsave(filename = paste0("2.Epithelium_trajectory_celltype_group2_",ci,".png"),p,w=15,h=5,dpi = 100)
  
  
  
  ##########################4 Differential genes associated with Pseudotime############################
  ##select order gene to regression，find them wither associate with time##
  ##if not select ,it will calculate all gene##
  library(dplyr)
  diff_test_res <- differentialGeneTest(HSMM_myo[ordering_genes,],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_gene_names <- diff_test_res %>% dplyr::arrange(qval) %>% head(50) %>% row.names() 
  write.csv(file = paste0('./Epithelium_Pseudotime_ordergene_',ci,'.csv'),diff_test_res)
  p = plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                              num_clusters = 6,
                              cores = 1,
                              show_rownames = T,return_heatmap = TRUE)
  ggsave(filename = paste0("2.Epithelium_top50gene_cortime_",ci,".png"),p,h=8,w=6,dpi = 150)
  
}


