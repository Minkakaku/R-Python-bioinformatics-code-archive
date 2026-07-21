# 拟时序 ------------
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


all_data <- readRDS("/mnt/NexData_file/home/hyj/Rstudio/XE2511036121-01_mengliqiang/result/Fibro_monocle_result/Fibro_mono.rds")
table(all_data@meta.data$subtype)

allCells = names(Idents(all_data))
allType = levels(Idents(all_data))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(all_data)== x ]
  num = round(length(cgCells)*0.1)
  cg=sample(cgCells,num)
  cg
}))

all_data_sub <- subset(all_data, cells = choose_Cells)
table(all_data_sub@meta.data$subtype)



monocle_outdir <- "/mnt/NexData_file/home/hyj/shell/XE2511036121-01_mengliqiang/Fibro_monocle_result"
dir.create(monocle_outdir)
# expr_dat_monocle <- as(LayerData(JoinLayers(sc_fib), "counts"), "sparseMatrix")
# all_data_subll<-subset(all_data_sub,subtype == 'M2'|subtype == 'M1'|subtype == 'Cycling M2')

count <- GetAssayData(object = JoinLayers(all_data_sub)[["RNA"]], layer = "counts")
fdata <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))
pdata <- all_data_sub@meta.data
#构建monocle对象
cds <- newCellDataSet(count,
                      featureData = new("AnnotatedDataFrame", data =fdata),
                      phenoData = new("AnnotatedDataFrame", data = pdata),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#标准化sacle及筛选，即使原先的Seurat对象已经归一化过，官方推荐在monocle中重新进行归一化处理
#标准化细胞之间的mRNA的差异
cds <- estimateSizeFactors(cds)
#离散度值可以帮助我们进行后续的差异分析
cds <- estimateDispersions(cds)

#在fData(cds)中添加“num_cell_expressed”列，可用于过滤细胞
cds <- detectGenes(cds, min_expr = 0.1)

head(pData(cds))
head(fData(cds))

expressed_genes<-row.names(subset(fData(cds),num_cells_expressed >= 10))
length(expressed_genes)

#向Epidata(cds)中添加“UMI”列
pdata$UMI <- Matrix::colSums(exprs(cds))

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

### 筛选高变基因(离散度高的gene) ----------
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit) 
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-1.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

# cds=orderCells(cds,root_state =3) 
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-1.pdf"))
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-1.pdf"), height = 10, width = 25)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

saveRDS(cds, file.path(monocle_outdir, "cds_disp-1.rds"))

gc()

### 高变基因 -----------
# group1_data<- group1_data[['RNA']] %>% as('Assay5')
# group1_data <- group1_data%>%NormalizeData(verbose=FALSE)%>%FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# group1_data <- FindVariableFeatures(group1_data, selection.method = "vst", nfeatures = 2000)
ordergene <- VariableFeatures(all_data_sub)
cds <- setOrderingFilter(cds, ordergene)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-2.pdf"), height = 10, width = 25)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-2.pdf"), height = 10, width = 10)
print(p)
dev.off()

saveRDS(cds, file.path(monocle_outdir, "cds_VariableFeatures-2.rds"))
gc()

### 差异基因 ----------
deg.cluster <- FindAllMarkers(JoinLayers(all_data_sub))
express_genes <- subset(deg.cluster,p_val_adj <0.05)$gene
express_genes <- unique(express_genes)
cds <- setOrderingFilter(cds, express_genes)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds<- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "subtype")
pdf(file.path(monocle_outdir, "monocle_cluster-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)
pdf(file.path(monocle_outdir, "monocle_facet-subtype-3.pdf"), height = 10, width = 25)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

saveRDS(cds, file.path(monocle_outdir, "cds_deg-3.rds"))
gc()
