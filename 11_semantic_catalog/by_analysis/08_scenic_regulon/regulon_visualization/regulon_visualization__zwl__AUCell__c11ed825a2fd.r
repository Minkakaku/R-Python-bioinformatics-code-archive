setwd('/mnt/NexData_file/home/zwl/Project/SC2602040915/require1')
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/SC2602040915/require1/Epithelium/res_raw.h5ad', from="anndata", to="seurat",
                      outFile='/mnt/NexData_file/home/zwl/Project/SC2602040915/require1/Epithelium/res_raw.rds')

library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(GSEABase)
sce<-readRDS('/mnt/NexData_file/home/zwl/Project/SC2602040915/require1/Epithelium/res_raw.rds')
cells_rankings <- AUCell_buildRankings(sce@assays$RNA@counts,nCores=1, plotStats=TRUE)
geneSets <- c('SFTPC', 'SFTPB', 'SFTPA1', 'SFTPA2', 'SFTPD', 'ABCA3', 'LAMP3', 'LPCAT1', 'SLC34A2',
              'ETV5', 'MUC1', 'MUC5', 'NAPSA', 'NKX2-1', 'HT2-280', 'GPRC5A', 'GRAMD2', 'RTKN2')
extraGeneSets <- c(GeneSet(sample(geneSets),setName="AT2 score"))
geneSets1 <- GeneSetCollection(extraGeneSets)
## 计算AUC
cells_AUC <- AUCell_calcAUC(geneSets1, cells_rankings)
cells_AUC

pdf('./Epithelium/celltype/cells_assignment.pdf',width = 6,height = 5)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
dev.off()

##set gene set of interest here for plotting
geneSet <- "AT2 score"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sce$AUC <- aucs
df<- data.frame(sce@meta.data, sce@reductions[["umap"]]@cell.embeddings)
colnames(df)

class_avg <- df %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

p1 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = AUC)) +
  viridis::scale_color_viridis(option="H") +
  ggrepel::geom_text_repel(aes(label = cell_type), data = class_avg, size = 6, segment.color = NA) +
  theme_classic() +
  theme(legend.position = "none")
library(ggpubr)
p.auc <- ggviolin(df, x = "cell_type", y = "AUC", color = "cell_type", add = 'mean_sd', fill = 'cell_type',
                  add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  scale_color_manual(values = my9color) +
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  NoLegend() +
  labs(x = '')

# 上下拼图
p_final <- p1 / p.auc

# 保存拼图
ggsave('./Epithelium/celltype/umap_AT2_score.png', p1, width = 8, height = 10, dpi = 300)
