
# spatial plot ---------------------
sp_output_filt <- "ROI_analysis/spatial plot"
dir.create(sp_output_filt, recursive = T)

all_data <- readRDS(file.path("data_for_analysis", "all_data_ROI_annotation.rds"))

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

cols_plan <- c(
  "Cortex"      = "#ea7070", 
  "Hippocampus" = "#600473", 
  "Other"       = "gray80"
)


sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = cols_plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - ROI.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

cols_plan <- c(
  "Glutamatergic neurons" = "#ea7070",           # 红色
  "GABAergic neurons" = "#fdc4b6",          # 浅橙色
  "Cholinergic neurons" = "#e59572",    # 橙色
  "Astrocytes" = "#2694ab", # 青色
  "Oligo" = "#96ceb4",           # 浅青色
  "Endothelial" = "#600473",                # 紫色
  "Other" = "gray90",                  # 蓝色
  "Microglial" = "#7F95D1",                 # 灰色
  "OPC" = "#ef5285",                # 浅蓝色
  "Epen" = "#a5dff9"
)


sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_Cortex, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = cols_plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })



# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_Cortex.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })




sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()

    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_Hippocampus, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = cols_plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })



# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - annotation_Hippocampus.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })




# annotation barplot ---------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

barplot_outdir1 <- file.path("ROI_analysis", "annotation ratio")
dir.create(barplot_outdir1)

data <- all_data[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 7, 6, 42)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

# data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(barplot_outdir1, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(barplot_outdir1, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(barplot_outdir1, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(barplot_outdir1, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

# orders <- read.table("order.txt", header = F, sep = "\t", quote = "", check.names = F)
# data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
# data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])



plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Sample", y = "Percent(%)", title = "") + 
  scale_fill_manual(values = cols_plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2)) + NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Group", y = "", title = "") + 
  scale_fill_manual(values = cols_plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 15, colour = "black"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2) + 
  plot_layout(widths = c(4, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(barplot_outdir1, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(barplot_outdir1, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)








# annotation barplot ---------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

barplot_outdir2 <- file.path("ROI_analysis", "annotation Cortex ratio")
dir.create(barplot_outdir2)

data <- all_data[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 7, 6, 45)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(barplot_outdir2, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(barplot_outdir2, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(barplot_outdir2, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(barplot_outdir2, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

# orders <- read.table("order.txt", header = F, sep = "\t", quote = "", check.names = F)
# data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
# data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])



plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Sample", y = "Percent(%)", title = "") + 
  scale_fill_manual(values = cols_plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2)) + NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Group", y = "", title = "") + 
  scale_fill_manual(values = cols_plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 15, colour = "black"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2) + 
  plot_layout(widths = c(4, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(barplot_outdir2, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(barplot_outdir2, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)






# annotation barplot ---------
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

barplot_outdir3 <- file.path("ROI_analysis", "annotation Hippocampus ratio")
dir.create(barplot_outdir3)

data <- all_data[[]] %>% rownames_to_column("Cell")

colnames(data)
data <- data[, c(1, 7, 6, 46)]
colnames(data) <- c("Cells", "Sample", "Group", "Annotation")

data <- data %>% filter(Annotation != "Other")

mydata <- data %>% count(Sample, Annotation)
mymtrix <- dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(barplot_outdir3, "Cell_type_stat.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Sample"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind,lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot1 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = file.path(barplot_outdir3, "Cell_type_stat_percent.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1) <- c("Sample", "Celltype", "Abund")

# samplemap <- read.table("sample_map_group.txt", header=F, sep="\t", check.names=F)
# colnames(samplemap)<-c("Sample","Group")
# final_data <- data %>% left_join(samplemap, by = "Sample")

final_data <- data

# write.table(final_data,file="2. dimension reduction and cluster/res.0.5/Ratio/test.txt",quote=F,row.names=F,col.names=T,sep="\t")
mydata <- final_data %>% count(Group, Annotation)
mymtrix <- dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)] <- 0
mymtrix_out <- t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = file.path(barplot_outdir3, "Cell_type_stat_group.txt"), quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix) <- mymtrix[, "Group"]
mymtrix <- mymtrix[, -1]
mymtrix <- as.matrix(mymtrix)
percent <- do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i,] <- round(mymtrix[i,]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent) <- rownames(mymtrix)
data4plot2 <- melt(percent)
percent <- t(percent)
DATA <- cbind(rownames(percent), percent)
DATA <- rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = file.path(barplot_outdir3, "Cell_type_stat_percent_group.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2) <- c("Group", "Celltype", "Abund")

# 限制样本顺序和分组顺序
sample_map <- read.table("sample_map_group.txt", header = F, sep = "\t", quote = "", check.names = F)
data4plot1$Sample <- factor(data4plot1$Sample, levels = sample_map[, 1])
data4plot2$Group <- factor(data4plot2$Group, levels = unique(sample_map[, 2]))

# orders <- read.table("order.txt", header = F, sep = "\t", quote = "", check.names = F)
# data4plot1$Celltype <- factor(data4plot1$Celltype, levels = orders[, 1])
# data4plot2$Celltype <- factor(data4plot2$Celltype, levels = orders[, 1])



plot_list <- list()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Sample", y = "Percent(%)", title = "") + 
  scale_fill_manual(values = cols_plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2)) + NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.7) + 
  labs(x = "Group", y = "", title = "") + 
  scale_fill_manual(values = cols_plan) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1), 
    axis.text.y = element_text(size = 15, colour = "black"), 
    legend.text = element_text(size = 15, colour = "black"), 
    legend.title = element_text(size = 15, colour = "black", face = "bold"), 
    axis.title = element_text(size = 15, colour = "black", face = "bold"), 
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5), 
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2) + 
  plot_layout(widths = c(4, 1)) # 两张图各自的宽度
ggsave(p, file = file.path(barplot_outdir3, "ratio_barplot.pdf"), width = 15, height = 6)
ggsave(p, file = file.path(barplot_outdir3, "ratio_barplot.png"), width = 15, height = 6, dpi = 300)













# 读取数据
all_data <- readRDS(file.path("data_for_analysis", "all_data_ROI_annotation.rds"))


# FindMarkers 海马差异分析 --------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 

all_data_Hippocampus <- all_data[, all_data@meta.data$ROI == "Hippocampus"]

unique(all_data_Hippocampus$annotation)
unique(all_data_Hippocampus$group)

table(all_data_Hippocampus$annotation, all_data_Hippocampus$group)

# Date <- Sys.Date()
# deg_outdir <- str_glue("FindMarkers-{Date}")
# dir.create(deg_outdir)


All_type <- as.character(unique(all_data_Hippocampus$annotation))
print(All_type)
length(All_type)


All_type %>% 
  setNames(.,rep('all_data_Hippocampus', 9)) %>% 
  imap(function(celltype,name){ # browser()
    tryCatch({
      object <- get(name)
      Idents(object)<-'annotation'
      # print(celltype)
      # print(name)
      # return()
      # 验证亚型隔离
      cat("Processing:", celltype, "| Cells:", sum(Idents(object) == celltype), "\n")
      
      # 使用子集对象提高效率
      subset_obj <- subset(JoinLayers(object), idents = celltype)
      
      diffmarkers <-
        FindMarkers(
          subset_obj,
          ident.1 = 'CTE',
          ident.2 = 'NC',
          test.use = "wilcox",
          logfc.threshold = 0.25,
          min.pct = 0.1,
          group.by = 'group'
        )
      # ==== 新增：检查空结果 ====
      if (nrow(diffmarkers) == 0) {
        cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
        write.csv(data.frame(), str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/NO_DEGs.csv'))
        return(NULL)
      }
      
      diffmarkers_up <- diffmarkers %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
      diffmarkers_down <- diffmarkers %>% filter(avg_log2FC < 0) %>% arrange(desc(avg_log2FC))
      
      dir.create(file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/')), recursive = TRUE)
      write.csv(diffmarkers,file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} group diffmarker.csv')))
      write.csv(diffmarkers_up,file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} group diffmarker up.csv')))
      write.csv(diffmarkers_down,file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} group diffmarker down.csv')))
      
      # GO up
      enrichGO_up_ls <- enrichGO(
        rownames(diffmarkers_up),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_up_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} Go up enriment table.csv')),row.names = FALSE)
      
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_up_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} Go up dotplot.pdf')),
        height = 6,
        width = 8)
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_up_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} Go up barplot.pdf')),
        height = 6,
        width = 8)
      
      
      # GO down
      enrichGO_down_ls <- enrichGO(
        rownames(diffmarkers_down),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_down_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} Go down enriment table.csv')),row.names = FALSE)
      
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_down_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} Go down dotplot.pdf')),
        height = 6,
        width = 8)
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_down_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} Go down barplot.pdf')),
        height = 6,
        width = 8)
      
      # KEGG up
      de_up_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers_up),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_up_ls <- enrichKEGG(
        de_up_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        # use_internal_data = FALSE
        use_internal_data = TRUE
      )
      enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_up_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} KEGG up enriment table.csv')),row.names = FALSE)
      
      enrichKEGG_up_dot <- enrichplot::dotplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50)
      enrichKEGG_up_bar <- barplot(
        enrichKEGG_up_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichKEGG_up_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} KEGG up dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichKEGG_up_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} KEGG up barplot.pdf')),
        height = 6,
        width = 8)
      
      
      # KEGG down
      de_down_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers_down),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_down_ls <- enrichKEGG(
        de_down_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        # use_internal_data = FALSE
        use_internal_data = TRUE
      )
      enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_down_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} KEGG down enriment table.csv')),row.names = FALSE)
      
      enrichKEGG_down_dot <- enrichplot::dotplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50)
      enrichKEGG_down_bar <- barplot(
        enrichKEGG_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichKEGG_down_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} KEGG down dotplot.pdf')),
        height = 6,
        width = 8)
      ggsave(
        enrichKEGG_down_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/{celltype} KEGG down barplot.pdf')),
        height = 6,
        width = 8)
      
      
      top <- rbind(diffmarkers %>% top_n(n = 25, wt = avg_log2FC),diffmarkers %>% top_n(n = -25, wt = avg_log2FC))
      p<-DoHeatmap(object, label = F,
                   features = rownames(top),
                   group.by = "group",
                   disp.max= 1,
                   disp.min= -1,
                   assay = 'RNA')+ #设置组别颜色
        scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
      pdf(file.path(str_glue("ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/diffgene_heatmap.pdf")), height = 8, width =10)
      print(p)
      dev.off()
      
      p_volcano <- EnhancedVolcano(diffmarkers,
                                   lab = rownames(diffmarkers),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   title = str_glue('{celltype} volcano plot'),
                                   pCutoff = 10e-2,
                                   FCcutoff = 0,
                                   pointSize = 3.0,
                                   labSize = 6.0)
      ggsave(p_volcano,
             filename = file.path(str_glue("ROI_analysis/FindMarkers_Hippocampus/{name}/{celltype}/diffgene valcano.pdf")),
             height = 6,
             width = 8)
    }, error = function(e) {
      # ===== 错误处理代码 ===== #
      message("\nERROR in ", celltype, ": ", e$message)
      # 创建专门存放错误对象的目录
      error_dir <- "ROI_analysis/FindMarkers_Hippocampus/all_data_error_objects/"
      dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
      # 保存引发错误的数据对象
      saveRDS(
        list(celltype = celltype, 
             error_msg = e$message, 
             subset_obj = if(exists("subset_obj")) subset_obj else NULL),
        file.path(error_dir, paste0(celltype, "_error.rds"))
      )
      return(NULL)
    })
  })












# CTE vs NC ---------
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr) 



all_data_filt <- all_data[, all_data@meta.data$ROI != "Other"]


All_type <- as.character(unique(all_data_filt$ROI))
print(All_type)
length(All_type)


All_type %>% 
  setNames(.,rep('all_data_filt', 2)) %>% 
  imap(function(celltype, name){ # browser()
    tryCatch({
      object <- get(name)
      Idents(object)<-'ROI'
      # print(celltype)
      # print(name)
      # return()
      # 验证亚型隔离
      cat("Processing:", celltype, "| Cells:", sum(Idents(object) == celltype), "\n")
      
      # 使用子集对象提高效率
      subset_obj <- subset(JoinLayers(object), idents = celltype)
      
      diffmarkers <- 
        FindMarkers(
          subset_obj, 
          ident.1 = 'CTE', 
          ident.2 = 'NC', 
          test.use = "wilcox", 
          logfc.threshold = 0.25, 
          min.pct = 0.1, 
          group.by = 'group'
        )
      # ==== 新增：检查空结果 ====
      if (nrow(diffmarkers) == 0) {
        cat("No DEGs found for", celltype, ". Skipping enrichment analysis.\n")
        write.csv(data.frame(), str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/NO_DEGs.csv'))
        return(NULL)
      }
      
      diffmarkers_up <- diffmarkers %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC))
      diffmarkers_down <- diffmarkers %>% filter(avg_log2FC < 0) %>% arrange(desc(avg_log2FC))
      
      dir.create(file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/')), recursive = TRUE)
      write.csv(diffmarkers, file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} group diffmarker.csv')))
      write.csv(diffmarkers_up, file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} group diffmarker up.csv')))
      write.csv(diffmarkers_down, file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} group diffmarker down.csv')))
      
      # GO up
      enrichGO_up_ls <- enrichGO(
        rownames(diffmarkers_up), 
        OrgDb = org.Mm.eg.db, 
        ont='ALL', 
        pAdjustMethod = 'BH', 
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05, 
        keyType = 'SYMBOL')
      write.csv(enrichGO_up_ls@result, file = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} Go up enriment table.csv')), row.names = FALSE)
      
      enrichGO_up_dot <- enrichplot::dotplot(
        enrichGO_up_ls, 
        showCategory = 10, 
        label_format = 50)
      ggsave(
        enrichGO_up_dot, 
        filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} Go up dotplot.pdf')),
        height = 6, 
        width = 8)
      enrichGO_up_bar <- barplot(
        enrichGO_up_ls, 
        showCategory = 10, 
        label_format = 50)
      ggsave(
        enrichGO_up_bar, 
        filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} Go up barplot.pdf')),
        height = 6,
        width = 8)
      
      
      # GO down
      enrichGO_down_ls <- enrichGO(
        rownames(diffmarkers_down),
        OrgDb = org.Mm.eg.db, 
        ont='ALL',
        pAdjustMethod = 'BH',
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        keyType = 'SYMBOL')
      write.csv(enrichGO_down_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} Go down enriment table.csv')),row.names = FALSE)
      
      enrichGO_down_dot <- enrichplot::dotplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_down_dot,
        filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} Go down dotplot.pdf')),
        height = 6,
        width = 8)
      enrichGO_down_bar <- barplot(
        enrichGO_down_ls,
        showCategory = 10,
        label_format = 50)
      ggsave(
        enrichGO_down_bar,
        filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} Go down barplot.pdf')),
        height = 6,
        width = 8)
      
      # KEGG up
      de_up_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers_up),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_up_ls <- enrichKEGG(
        de_up_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        # use_internal_data = FALSE
        use_internal_data = TRUE
      )
      enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_up_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} KEGG up enriment table.csv')),row.names = FALSE)
      
      # enrichKEGG_up_dot <- enrichplot::dotplot(
      #   enrichKEGG_up_ls,
      #   showCategory = 10,
      #   label_format = 50)
      # enrichKEGG_up_bar <- barplot(
      #   enrichKEGG_up_ls,
      #   showCategory = 10,
      #   label_format = 50)
      # ggsave(
      #   enrichKEGG_up_dot,
      #   filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} KEGG up dotplot.pdf')),
      #   height = 6,
      #   width = 8)
      # ggsave(
      #   enrichKEGG_up_bar,
      #   filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} KEGG up barplot.pdf')),
      #   height = 6,
      #   width = 8)
      
      
      # KEGG down
      de_down_genes_id_for_enrichmenet <- 
        clusterProfiler::bitr(
          rownames(diffmarkers_down),
          OrgDb = "org.Mm.eg.db",
          fromType = "SYMBOL",
          toType = "ENTREZID"
        ) %>%
        .$ENTREZID
      enrichKEGG_down_ls <- enrichKEGG(
        de_down_genes_id_for_enrichmenet,
        organism = "mmu", 
        keyType = "ncbi-geneid",
        # keyType = "kegg", 
        pvalueCutoff = 0.05, 
        pAdjustMethod = "BH", 
        minGSSize = 10,
        maxGSSize = 500, 
        qvalueCutoff = 0.2,
        # use_internal_data = FALSE
        use_internal_data = TRUE
      )
      enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      write.csv(enrichKEGG_down_ls@result,file = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} KEGG down enriment table.csv')),row.names = FALSE)
      
      # enrichKEGG_down_dot <- enrichplot::dotplot(
      #   enrichKEGG_down_ls,
      #   showCategory = 10,
      #   label_format = 50)
      # enrichKEGG_down_bar <- barplot(
      #   enrichKEGG_down_ls,
      #   showCategory = 10,
      #   label_format = 50)
      # ggsave(
      #   enrichKEGG_down_dot,
      #   filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} KEGG down dotplot.pdf')),
      #   height = 6,
      #   width = 8)
      # ggsave(
      #   enrichKEGG_down_bar,
      #   filename = file.path(str_glue('ROI_analysis/FindMarkers/{name}/{celltype}/{celltype} KEGG down barplot.pdf')),
      #   height = 6,
      #   width = 8)
      
      
      top <- rbind(diffmarkers %>% top_n(n = 25, wt = avg_log2FC),diffmarkers %>% top_n(n = -25, wt = avg_log2FC))
      p<-DoHeatmap(object, label = F,
                   features = rownames(top),
                   group.by = "group",
                   disp.max= 1,
                   disp.min= -1,
                   assay = 'RNA')+ #设置组别颜色
        scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
      pdf(file.path(str_glue("ROI_analysis/FindMarkers/{name}/{celltype}/diffgene_heatmap.pdf")), height = 8, width =10)
      print(p)
      dev.off()
      
      p_volcano <- EnhancedVolcano(diffmarkers,
                                   lab = rownames(diffmarkers),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj',
                                   title = str_glue('{celltype} volcano plot'),
                                   pCutoff = 10e-2,
                                   FCcutoff = 0,
                                   pointSize = 3.0,
                                   labSize = 6.0)
      ggsave(p_volcano,
             filename = file.path(str_glue("ROI_analysis/FindMarkers/{name}/{celltype}/diffgene valcano.pdf")),
             height = 6,
             width = 8)
    }, error = function(e) {
      # ===== 错误处理代码 ===== #
      message("\nERROR in ", celltype, ": ", e$message)
      # 创建专门存放错误对象的目录
      error_dir <- "ROI_analysis/FindMarkers/all_data_error_objects/"
      dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
      # 保存引发错误的数据对象
      saveRDS(
        list(celltype = celltype, 
             error_msg = e$message, 
             subset_obj = if(exists("subset_obj")) subset_obj else NULL),
        file.path(error_dir, paste0(celltype, "_error.rds"))
      )
      return(NULL)
    })
  })













