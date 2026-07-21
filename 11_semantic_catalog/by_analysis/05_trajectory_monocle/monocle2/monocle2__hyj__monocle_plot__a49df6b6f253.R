
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)
library(tibble)



monocle_outdir <- "/mnt/NexData_file/home/liuran/projects/3_Xenium/20260113_FWSP20252695_Thyroid_Xenium/monocle_results_260407"


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


cds <- readRDS(file.path(monocle_outdir, "cds_deg-3.rds"))
pData(cds)$RCTD_subtype2 <- factor(pData(cds)$RCTD_subtype2, levels = c(paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)))

pdf(file.path(monocle_outdir, "cell trajectory.pdf"), width = 8, height = 6)
plot_cell_trajectory(cds,color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds,color_by = "RCTD_subtype2", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - - split by RCTD_subtype2.pdf"), width = 25, height = 25)
plot_cell_trajectory(cds, color_by = "RCTD_subtype2") + facet_wrap("~RCTD_subtype2", nrow = 5) + theme(legend.position = "right") + scale_color_manual(values = plan) 
dev.off()

pdf(file.path(monocle_outdir, "cell trajectory - - split by state.pdf"), width = 15, height = 10)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 2) + scale_color_manual(values = plan)
dev.off()

pdf(file.path(monocle_outdir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 5)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("PA", "PC"))
  plot_cell_trajectory(cds, color_by = "RCTD_subtype2", size = 1, show_backbone = TRUE) + 
    scale_color_manual(values = plan) + 
    facet_wrap(~ group) + 
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(monocle_outdir, "cell trajectory - tree.pdf"), width = 30, height = 30)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "RCTD_subtype2") + theme(legend.position = "right") + 
  scale_color_manual(values = plan) +
  facet_wrap(~ RCTD_subtype2)
dev.off()

# cell barplot 
p_col_cell_state <- 
  cds %>% 
  pData() %>% 
  dplyr::select(RCTD_subtype2, State) %>% 
  group_by(RCTD_subtype2, State) %>% 
  mutate(RCTD_subtype2 = as.character(RCTD_subtype2)) %>% 
  mutate(RCTD_subtype2 = factor(RCTD_subtype2, levels = c(
    paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
  ))) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = RCTD_subtype2, y = freq)) + 
  geom_col(aes(fill = State)) + 
  scale_fill_manual(values = plan) + 
  scale_y_continuous(labels = scales::percent) + 
  coord_flip() +
  labs(x = NULL, y = NULL) + 
  ggpubr::theme_pubr() + 
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_outdir, "cell stat barplot.pdf"), height = 12, width = 8)






# 沿时间轴的细胞密度图 -----------
colour = c('#C6DADF', '#E26EB8', '#918FBE', '#A2C5A6', '#EED784', '#D4B796', 
           '#F5BDC8', '#F0BDD8', '#77A1EC', '#F8D39C', '#F17591', '#CCC9E6')

df <- pData(cds)[, c("Pseudotime", "RCTD_subtype2")]
p_density <- ggplot(df, aes(Pseudotime, colour = RCTD_subtype2, fill = RCTD_subtype2)) + 
  geom_density(bw = 0.5, size = 0.5, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = plan) + 
  scale_color_manual(name = "", values = plan)
ggsave(p_density, filename = file.path(monocle_outdir, "cell density plot.pdf"), height = 6, width = 12)



