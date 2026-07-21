library(clusterProfiler)
library(RColorBrewer)
library(tidyverse)

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20253251/08.Mission/Enrich/"

setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/Bcells/")
groups <- c("Activate B","Memory B","Naive B","pDC")
all_top5 <- list()
for(group in groups) {
  file_path <- paste0("./cluster", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  deg_top5 <- deg[1:5, ]
  deg_top5$group <- group
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
# 绘图
p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = brewer.pal(7, "RdYlBu")) +
  scale_size_continuous(
    range = c(0, 25),
    breaks = c(0, 10, 20, 30, 40, 50),
    # labels = c("50", "100", "150", "200", "250", "300"),
    name = "Count"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),      # X和Y轴标题字体大小
    axis.text = element_text(size = 15)       # X和Y轴刻度标签字体大小
  ) +
  labs(x = "group", fill = "Adjusted.P.value")

ggsave(paste0(outdir, "Bcells_GOBP_dotplot.pdf"),width = 15,height = 15,dpi=300,p)
ggsave(paste0(outdir, "Bcells_GOBP_dotplot.png"),width = 15,height = 15,dpi=300,p)


setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/DC/")
groups <- c("cDC1-Ccr7","cDC1-Clec9a","cDC2-Cd209a")
all_top5 <- list()
for(group in groups) {
  file_path <- paste0("./cluster", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  deg_top5 <- deg[1:5, ]
  deg_top5$group <- group
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
# 绘图
p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = brewer.pal(7, "RdYlBu")) +
  scale_size_continuous(
    range = c(0, 25),
    breaks = c(0, 10, 20, 30, 40, 50),
    # labels = c("50", "100", "150", "200", "250", "300"),
    name = "Count"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),      # X和Y轴标题字体大小
    axis.text = element_text(size = 15)       # X和Y轴刻度标签字体大小
  ) +
  labs(x = "group", fill = "Adjusted.P.value")

ggsave(paste0(outdir, "DC_GOBP_dotplot.pdf"),width = 15,height = 15,dpi=300,p)
ggsave(paste0(outdir, "DC_GOBP_dotplot.png"),width = 15,height = 15,dpi=300,p)


setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Tcells/")
groups <- c("Activated T", "Cytotoxic T", "NK", "Naive T", "Proliferative T", "Tfh", "Th1", "Th2", "Treg")
all_top5 <- list()
for(group in groups) {
  file_path <- paste0("./cluster", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  deg_top5 <- deg[1:5, ]
  deg_top5$group <- group
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
# 绘图
p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = brewer.pal(7, "RdYlBu")) +
  scale_size_continuous(
    range = c(0, 25),
    breaks = c(0, 10, 20, 30, 40, 50),
    # labels = c("50", "100", "150", "200", "250", "300"),
    name = "Count"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 25),      # X和Y轴标题字体大小
    axis.text = element_text(size = 25)       # X和Y轴刻度标签字体大小
  ) +
  labs(x = "group", fill = "Adjusted.P.value")

ggsave(paste0(outdir, "Tcells_GOBP_dotplot.pdf"),width = 30,height = 32,dpi=300,p)
ggsave(paste0(outdir, "Tcells_GOBP_dotplot.png"),width = 30,height = 32,dpi=300,p)


setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Neutrophils/")
groups <- c("Cmss1 Neutrophils","Il1rn Neutrophils","Pglyrp1 Neutrophils")
all_top5 <- list()
for(group in groups) {
  file_path <- paste0("./cluster", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  deg_top5 <- deg[1:5, ]
  deg_top5$group <- group
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
# 绘图
p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = brewer.pal(7, "RdYlBu")) +
  scale_size_continuous(
    range = c(0, 25),
    breaks = c(0, 10, 20, 30, 40, 50),
    # labels = c("50", "100", "150", "200", "250", "300"),
    name = "Count"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),      # X和Y轴标题字体大小
    axis.text = element_text(size = 15)       # X和Y轴刻度标签字体大小
  ) +
  labs(x = "group", fill = "Adjusted.P.value")

ggsave(paste0(outdir, "Neutrophils_GOBP_dotplot.pdf"),width = 15,height = 15,dpi=300,p)
ggsave(paste0(outdir, "Neutrophils_GOBP_dotplot.png"),width = 15,height = 15,dpi=300,p)


setwd("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Macrophages/")
groups <- c("Aoah Macrophage","Gda Macrophage","Gpx3 Macrophage","Inpp5d Macrophage","Jun Macrophage")
all_top5 <- list()
for(group in groups) {
  file_path <- paste0("./cluster", group, "/GO_Biological_Process_2023.Mouse.enrichr.reports.txt")
  deg <- read.csv(file = file_path, header = TRUE, sep = '\t')
  deg_top5 <- deg[1:5, ]
  deg_top5$group <- group
  all_top5[[group]] <- deg_top5
}
# 合并所有数据
combined_top5 <- do.call(rbind, all_top5)
# 查看结果
print(combined_top5)
combined_top5$Description <- gsub("\\s*\\([^)]*\\)", "", combined_top5$Term)
combined_top5$Count <- sub("/.*", "", combined_top5$Overlap)
combined_top5$Count <- as.numeric(combined_top5$Count)
combined_top5$Description <- factor(combined_top5$Description, levels = unique(combined_top5$Description))
combined_top5$group <- factor(combined_top5$group, levels = groups)
# 绘图
p = ggplot(combined_top5, aes(x = group, y = Description, fill = Adjusted.P.value, size = Count)) +
  geom_point(shape = 21, colour = "black") +  # 黑圈
  scale_fill_gradientn(colors = brewer.pal(7, "RdYlBu")) +
  scale_size_continuous(
    range = c(0, 25),
    breaks = c(0, 10, 20, 30, 40, 50),
    # labels = c("50", "100", "150", "200", "250", "300"),
    name = "Count"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 15),      # X和Y轴标题字体大小
    axis.text = element_text(size = 15)       # X和Y轴刻度标签字体大小
  ) +
  labs(x = "group", fill = "Adjusted.P.value")

ggsave(paste0(outdir, "Macrophages_GOBP_dotplot.pdf"),width = 20,height = 20,dpi=300,p)
ggsave(paste0(outdir, "Macrophages_GOBP_dotplot.png"),width = 20,height = 20,dpi=300,p)