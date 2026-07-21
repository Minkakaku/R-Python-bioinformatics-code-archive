library(argparser)
library(Seurat)
library(tidyverse)

# usage
# Rscript sc_analysis.R --dims 15 --rds all_data_dim8.rds --save --k.param 20


argv <- arg_parser('')
# argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--dims", help="dimenstion")
argv <- add_argument(argv,"--rds", help="Seurat rds")
argv <- add_argument(argv,"--save", help="save Seurat rds", flag = TRUE)
argv <- add_argument(argv,"--k.param", help="k.param to find neighbors", default = 20)

argv <- parse_args(argv)

# print all arguments
print(argv)

dim  <- argv$dims
rds  <- argv$rds
save <- argv$save
k    <- argv$k.param


# analysis -----------------------
all_data <- readRDS(rds)

dir.out <- str_glue("dim{dim}")
dir.create(dir.out)

pdf(str_glue("{dir.out}/Dimplot dim{dim} project.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% RunUMAP(dims = 1:dim)

all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dim,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:dim, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_"
)


pdf(str_glue("{dir.out}/Dimplot dim{dim} project - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### cluster
all_data <-
  FindNeighbors(
    all_data,
    reduction = "harmony",
    dims = 1:dim,
    k.param = k
  )
all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(str_glue("{dir.out}/Dimplot dim{dim} cluster - harmony.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# dotplot
marker_ls <- c(
  "Cardiomyocytes" = "Myh6",
  "Cardiomyocytes" = "Ttn",
  "Cardiomyocytes" = "Tnnt2",
  "Cardiomyocytes" = "Actn2",
  "Cardiomyocytes" = "Gja1",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Mylk",
  "Smooth muscle cells" = "Cald1",
  "Smooth muscle cells" = "Pdgfrb",
  "Smooth muscle cells" = "Rgs5",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  "Fibroblasts" = "Col13a1",
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Dcn",
  "Immune cells" = "Ctsc",
  "Immune cells" = "Mrc1",
  "Immune cells" = "Cd74",
  "Immune cells" = "F13a1",
  "Immune cells" = "Csf1r",
  "Endothelial cells" = "Fabp4",
  "Endothelial cells" = "Egfl7",
  "Endothelial cells" = "Flt1",
  "Endothelial cells" = "Cdh5",
  "Endothelial cells" = "Emcn"
)


pdf(str_glue("{dir.out}/dotplot dims{dim}.pdf"), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls)
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls)
dev.off()


if(save){
  saveRDS(all_data, str_glue("all_data_dim{dim}.rds"))
}


