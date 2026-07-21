# install.packages("BiocManager")
# BiocManager::install(c("copykat"))
# install.packages(c("Seurat", "dplyr"))

library(Seurat)
library(copykat)

obj  = readRDS("copykat.RDS")

# downsample
 downsample_seurat <- function(seu, n.cells = 1000) {
   if (ncol(seu) <= n.cells) {
     message("Cell number <= target, skip downsampling")
     return(seu)
   }
   cells.use <- sample(colnames(seu), n.cells)
   subset(seu, cells = cells.use)
 }
# obj <- downsample_seurat(obj, n.cells = 15000)

DefaultAssay(obj) <- "RNA"

# Seurat v5: 优先取 counts layer
obj = JoinLayers(obj)
mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
mat <- as.matrix(mat)

# 建议过滤：去掉全0基因
mat <- mat[rowSums(mat) > 0, ]
table(obj@meta.data$cell_type)
# 例：用免疫细胞当 reference（按你自己的列名改）
ref_cells <- colnames(obj)[
  obj$cell_type %in% c("DC")
]

set.seed(123)

ck <- copykat(
  rawmat = mat,
  id.type = "S", # "S"=SYMBOL; 如果是 Ensembl 用 "E"
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1, # 常用 0.1~0.15，越大越“保守”
  sam.name = "sample1",
  distance = "euclidean",
  norm.cell.names = ref_cells,
  n.cores = 20,
  genome = "mm10"
)

# ck$prediction 里是每个细胞的类别：aneuploid / diploid 等
head(ck$prediction)
