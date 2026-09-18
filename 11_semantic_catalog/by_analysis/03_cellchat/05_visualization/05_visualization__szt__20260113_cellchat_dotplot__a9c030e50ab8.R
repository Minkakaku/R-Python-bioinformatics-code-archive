library(tidyverse)
library(Seurat)
library(CellChat)


out_dir_dotplot <- "20260113 dotplot of selected ligand and receptor pairs"
dir.create(out_dir_dotplot)


all_lr_dat <-
  list.files('cellchat_dir_annov5_v2 3grps - remove pku/', pattern = "net_lr", full.names = TRUE) %>%
  setNames(., str_remove(str_remove(basename(.), "^.*\\s"), "\\..*$") ) %>%
  map(read.table, sep = "\t", header = TRUE) %>%
  enframe(name = "cluster") %>%
  unnest("value")




## 1. epi -> macro ------------------
lr_epi_mac_1 <- c(
  'APP_CD74',
  
  'MIF_CD74_CXCR4',
  "MIF_CD74_CD44",
  "MIF_ACKR3",
  
  "CD99_PILRA",
  "CD99_CD99"
)
  
# 
# "caf"
# "EC" "Pericyte"
# "gial_cell"
# "Macrophage" "Monocyte" "DC"    
# "Mast Cells"
# "T" "B_plasma"
# "normal_epi"
# "pmc_infla"
# "OLFM4+tumor_cells"
# "pre_tumor"
# "tumor_cells" 


p_ct_cellchat_1 <-
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  dplyr::filter(interaction_name   %in% lr_epi_mac_1) %>%
  dplyr::filter(source %in% c("normal_epi", "pmc_infla", "pre_tumor", "OLFM4+tumor_cells", "tumor_cells")) %>%
  dplyr::filter(target %in% c("Macrophage")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "1. bubble plot - epi - macropahge.pdf"),
  height = 7,
  width  = 6
)
print(p_ct_cellchat_1)
dev.off()


##2. macro -> epi ------------------
lr_mac_epi_2 <- c(
  "LGALS9_CD45", 
  "LGALS9_HAVCR2",
  "LGALS9_CD44",
  "LGALS9_P4HB",
  
  "TNF_TNFRSF1A",
  "TNF_TNFRSF1B",
  
  "IL1B_IL1R1_IL1RAP",
  "IL1B_IL1R2",
  
  "HBEGF_EGFR", 
  "HBEGF_EGFR_ERBB2",
  "AREG_EGFR",
  "AREG_EGFR_ERBB2"
)


p_ct_cellchat_2 <-
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  dplyr::filter(interaction_name   %in% lr_mac_epi_2) %>%
  dplyr::filter(source %in% c("Macrophage")) %>%
  dplyr::filter(target %in% c("normal_epi", "pmc_infla", "pre_tumor", "OLFM4+tumor_cells", "tumor_cells")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "2. bubble plot - macrophage - epi.pdf"),
  height = 10,
  width  = 6
)
print(p_ct_cellchat_2)
dev.off()



##3. epi -> caf ------------------
lr_epi_caf_3 <- c(
  "MDK_SDC1",
  "MDK_SDC2",
  "MDK_SDC4",
  "MDK_PTPRZ1",
  "MDK_ITGA4_ITGB1",
  "MDK_ITGA6_ITGB1",
  "MDK_LRP1",
  "MDK_NCL",
  "MDK_ALK",
  
  "AREG_EGFR",
  "AREG_EGFR_ERBB2",
  "BTC_EGFR",
  "BTC_EGFR_ERBB2",
  "HBEGF_EGFR",
  "HBEGF_EGFR_ERBB2",
  "EREG_EGFR",
  "EREG_EGFR_ERBB2",
  "VEGFA_VEGFR1",
  "VEGFA_VEGFR2",
  "VEGFC_VEGFR3",
  "VEGFC_VEGFR2",
  "PGF_VEGFR1",
  "VEGFA_VEGFR1R2",
  "VEGFC_VEGFR2R3",
  "VEGFB_VEGFR1",
  
  'GDF15_TGFBR2',
  
  "PRSS3_F2R",
  "PRSS3_F2RL1",
  "PRSS3_F2RL3",
  "PRSS3_PARD3"
)


p_ct_cellchat_3 <-
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  dplyr::filter(interaction_name   %in% lr_epi_caf_3) %>%
  dplyr::filter(source %in% c("normal_epi", "pmc_infla", "pre_tumor", "OLFM4+tumor_cells", "tumor_cells")) %>%
  dplyr::filter(target %in% c("caf")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "3. bubble plot - epi - caf.pdf"),
  height = 10,
  width  = 6
)
print(p_ct_cellchat_3)
dev.off()




##4. caf -> epi ------------------
lr_caf_epi_4 <- c(
  "MDK_SDC1",
  "MDK_SDC2",
  "MDK_SDC4",
  "MDK_PTPRZ1",
  "MDK_ITGA4_ITGB1",
  "MDK_ITGA6_ITGB1",
  "MDK_LRP1",
  "MDK_NCL",
  "MDK_ALK",
  
  "NRG1_ERBB3",
  "NRG1_ERBB2_ERBB3",
  "NRG1_INTEGRIN_A6B4_COMPLEX",
  
  "PLAU_PLAUR",
  
  "COL6A2_ITGA1_ITGB1",
  "COL6A2_ITGA2_ITGB1",
  "COL6A2_ITGA3_ITGB1",
  "COL6A2_ITGA9_ITGB1",
  "COL6A2_ITGA11_ITGB1",
  "COL6A2_ITGAV_ITGB8",
  "COL6A2_CD44",
  "COL6A2_SDC1",
  "COL6A2_SDC4",
  
  "THBS2_ITGA3_ITGB1",
  "THBS2_SDC1",
  "THBS2_SDC4",
  "THBS2_CD36",
  "THBS2_CD47"
)


p_ct_cellchat_4 <-
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  dplyr::filter(interaction_name   %in% lr_caf_epi_4) %>%
  dplyr::filter(source %in% c("caf")) %>%
  dplyr::filter(target %in% c("normal_epi", "pmc_infla", "pre_tumor", "OLFM4+tumor_cells", "tumor_cells")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "4. bubble plot - caf - epi.pdf"),
  height = 15,
  width  = 6
)
print(p_ct_cellchat_4)
dev.off()


##5. epi -> immune ------------------
lr_epi_immune_5 <- c(
  "NECTIN3_NECTIN2",
  "NECTIN3_PVR",
  "NECTIN3_TIGIT",
  "NECTIN2_TIGIT",
  "CADM1_NECTIN3",
  "CADM3_NECTIN3",
  "NECTIN2_CD226"
)


p_ct_cellchat_5 <-
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  dplyr::filter(interaction_name   %in% lr_epi_immune_5) %>%
  dplyr::filter(source %in% c("caf")) %>%
  dplyr::filter(target %in% c("Macrophage", "Monocyte", "DC", "Mast Cells", "T", "B_plasma")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "5. bubble plot - caf - immune.pdf"),
  height = 6,
  width  = 6
)
print(p_ct_cellchat_5)
dev.off()



## T -> macro ------------
p_ct_cellchat_6 <-
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  # dplyr::filter(interaction_name   %in% lr_epi_immune_5) %>%
  dplyr::filter(source %in% c("T")) %>%
  dplyr::filter(target %in% c("Macrophage", "Monocyte", "DC", "Mast Cells", "T")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "6. bubble plot - t - macro.pdf"),
  height = 15,
  width  = 6
)
print(p_ct_cellchat_6)
dev.off()

## macro -> T ------------
p_ct_cellchat_7 <- 
  all_lr_dat %>%
  mutate(cluster = factor(cluster, levels = c("NAG", "pre_cancer", "cancer"))) %>%
  # dplyr::filter(interaction_name   %in% lr_epi_immune_5) %>%
  dplyr::filter(source %in% c("Macrophage", "Monocyte", "DC", "Mast Cells", "T")) %>%
  dplyr::filter(target %in% c("T")) %>%
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, ncol = 1) +
  scale_size(transform = "reverse") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

pdf(
  file.path(out_dir_dotplot, "7. bubble plot - macro - t.pdf"),
  height = 18,
  width  = 6
)
print(p_ct_cellchat_7)
dev.off()



