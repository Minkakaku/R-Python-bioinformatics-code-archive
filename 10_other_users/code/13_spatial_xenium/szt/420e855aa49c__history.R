scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 5,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
)
deconProp_for_plot %>%
dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 20,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
)
deconProp_for_plot %>%
dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 20,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
deconProp_for_plot %>%
dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 10,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
deconProp_for_plot %>%
# dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 10,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot, # %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
deconProp_for_plot %>%
dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 10,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
deconProp_for_plot %>%
dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 15,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
deconProp_for_plot %>%
# dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 15,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot， # %>% dplyr::filter(x < 15000, y < 15000),
deconProp_for_plot %>%
# dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 15,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot, # %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
View(deconProp_for_plot)
p_deconProp_for_plot <-
deconProp_for_plot %>%
# dplyr::filter(x < 15000, y < 15000) %>%
ggplot(aes(x = x, y = y)) +
# coord_equal(
#   xlim = c(0, deconProp_for_plot$x %>% range() %>% diff),
#   ylim = c(0, deconProp_for_plot$y %>% range() %>% diff)
# ) +
coord_fixed(ratio = 1/1) +
scatterpie::geom_scatterpie(
aes(x = x, y = y, group = Row.names, r = 15,
# color = Pixel.Groups
),
pie_scale = 1,
# lwd = lwd,
data = deconProp_for_plot, # %>% dplyr::filter(x < 15000, y < 15000),
cols = as.character(1:4), # topicColumns,
legend_name = "CellTypes"
) +
theme_void()
ggsave(
p_deconProp_for_plot,
file = file.path(dr_outdir, "STdeconvolve spatial plot.pdf"),
height = 32,
width  = 32
)
## cluster ------
sc_n <-
FindNeighbors(
sc_n,
reduction = "pca",
dims = 1:dims,
k.param = 20
)
sc_n <-
FindClusters(
sc_n,
resolution = c(0.2, 0.4, 0.6, 0.8)
)
dims
# dimplot
pdf(file.path(dr_outdir, str_glue("sc_n Dimplot dim{dims} cluster - pca.pdf")), height = 6, width = 8)
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_n, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()
# spatial-dimplot
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} cluster - pca.pdf")), height = 12, width = 12)
str_c('Spatial.008um_snn_res.', seq(0.2, 0.8, by = 0.2) ) %>%
map(function(res){
SpatialDimPlot(
sc_n,
group.by = res,
label = T,
repel = T,
alpha = 0.5,
shape = 22,
label.size = 4,
pt.size.factor = 6,
image.scale = "hires"
# images = name
) +
ggsci::scale_fill_igv()
})
dev.off()
## dotplot -------
major_marker_ls <- c(
'Myoid cell' = 'ACTA2',
'Myoid cell' = 'MYH11',
'Myoid cell' = 'PTCH1',
'Myoid cell' = 'PTCH2',
'Sertoli cell' = 'WT1',
'Sertoli cell' = 'SOX9',
'Sertoli cell' = 'BEX2',
'Sertoli cell' = 'GATA4',
'Macrophage cell' = 'CD68',
'Macrophage cell' = 'MS4A6A',
'Macrophage cell' = 'LYZ',
'Macrophage cell' = 'IL3RA',
'Leydig cell' = 'CYP11A1',
'Leydig cell' = 'STAR',
'Leydig cell' = 'IGF1',
'Leydig cell' = 'IGFBP5',
'Endothelial cell' = 'VWF',
'Endothelial cell' = 'PECAM1',
'Endothelial cell' = 'CD34',
'Endothelial cell' = 'TGFBR2',
'T cell' = 'CD3D',
'T cell' = 'CD3E',
'T cell' = 'CD3G',
'T cell' = 'CD2',
Spermatogonia = 'FGFR3',
Spermatogonia = 'STRA8',
Spermatogonia = 'DMRT1',
Spermatogonia = 'UCHL1',
Spermatogonia = 'TCF3',
Spermatocytes = 'SYCP1',
Spermatocytes = 'SMC1B',
Spermatocytes = 'SYCP3',
Spermatocytes = 'OVOL2',
Spermatids = 'TNP1',
Spermatids = 'PRM1',
Spermatids = 'TXNDC2',
Spermatids = 'PRM3'
)
theme_90angle <-
theme(
axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)
pdf(file.path(dr_outdir, str_glue("sc_n dotplot dims{dims} pca main_marker.pdf")), width = 25, height = 6)
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sc_n, group.by = "Spatial.008um_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()
major_marker_ls2 <-
major_marker_ls %>%
enframe("pathway", "gene") %>%
group_by(pathway) %>%
nest() %>%
deframe() %>%
map(unlist, use.names = FALSE)
pathway_score_ls <-
major_marker_ls2 %>%
imap(function(x, name){# browser()
sc_n %>%
AddModuleScore(features = list(x) ) %>%
FeaturePlot(
features = "Cluster1"
) +
ggtitle(name)
})
pdf(
file.path(dr_outdir, str_glue("cell annotation pathway scores dims{dims}.pdf")),
height = 18,
width  = 24
)
pathway_score_ls %>% patchwork::wrap_plots()
dev.off()
## find markers -------
Idents(sc_n) <- "Spatial.008um_snn_res.0.6"
DefaultAssay(sc_n) <- "Spatial.008um"
all_markers <-
FindAllMarkers(
JoinLayers(sc_n),
test.use = "wilcox",
logfc.threshold = 0.25,
min.pct = 0.1,
only.pos = TRUE
)
all_markers %>%
group_by(cluster) %>%
dplyr::filter(p_val_adj < 0.05) %>%
dplyr::filter(pct.1 > 0.25) %>%
arrange(desc(avg_log2FC)) %>%
slice_head(n = 10) %>%
nest() %>%
deframe() %>%
map(pull, 'gene') %>%
imap(~{
xx <- str_c(.x, collapse = ",")
str_glue("cluster_{.y}:{xx}")
}) %>%
walk(cat, "\n")
major_marker_ls %>% names %>% unique
## annotation ------
cell_cluster <- c(
'0'   = "Spermatids",
'1'   = "Myoid cells", ##
'2'   = "Spermatocytes",
'3'   = "Sertoli cell", #
'4'   = "Leydig cell", #
'5'   = "Spermatids",
'6'   = "Myoid cell",
'7'   = "Leydig cell",
'8'   = "Spermatids",
'9'   = "Leydig cell", ##
'10'  = "Myoid cell",
'11'  = "Myoid cells", ##
'12'  = "Myoid cell",
'13'  = "Myoid cell",
'14'  = "Myoid cell"
)
sc_n <- RenameIdents(sc_n, cell_cluster)
sc_n$annotation <- Idents(sc_n)
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} annotation - pca.pdf")), height = 12, width = 12)
str_c('Spatial.008um_snn_res.', seq(0.2, 0.8, by = 0.2) ) %>%
map(function(res){
SpatialDimPlot(
sc_n,
group.by = res,
label = T,
repel = T,
alpha = 0.5,
shape = 22,
label.size = 4,
pt.size.factor = 6,
image.scale = "hires"
# images = name
) +
ggsci::scale_fill_igv()
})
# spatial-dimplot
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} annotation - pca.pdf")), height = 12, width = 12)
SpatialDimPlot(
sc_n,
group.by = "annotation",
label = T,
repel = T,
alpha = 0.5,
shape = 22,
label.size = 4,
pt.size.factor = 6,
image.scale = "hires"
# images = name
) +
ggsci::scale_fill_igv()
dev.off()
## annotation ------
cell_cluster <- c(
'0'   = "Spermatids",
'1'   = "Myoid cells", ##
'2'   = "Spermatocytes",
'3'   = "Sertoli cell", #
'4'   = "Leydig cell", #
'5'   = "Spermatids",
'6'   = "Myoid cell",
'7'   = "Leydig cell",
'8'   = "Spermatids",
'9'   = "Leydig cell", ##
'10'  = "Myoid cell",
'11'  = "Myoid cell", ##
'12'  = "Myoid cell",
'13'  = "Myoid cell",
'14'  = "Myoid cell"
)
sc_n <- RenameIdents(sc_n, cell_cluster)
## annotation ------
Idents(sc_n) <- "Spatial.008um_snn_res.0.6"
cell_cluster <- c(
'0'   = "Spermatids",
'1'   = "Myoid cells", ##
'2'   = "Spermatocytes",
'3'   = "Sertoli cell", #
'4'   = "Leydig cell", #
'5'   = "Spermatids",
'6'   = "Myoid cell",
'7'   = "Leydig cell",
'8'   = "Spermatids",
'9'   = "Leydig cell", ##
'10'  = "Myoid cell",
'11'  = "Myoid cell", ##
'12'  = "Myoid cell",
'13'  = "Myoid cell",
'14'  = "Myoid cell"
)
sc_n <- RenameIdents(sc_n, cell_cluster)
sc_n$annotation <- Idents(sc_n)
# spatial-dimplot
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} annotation - pca.pdf")), height = 12, width = 12)
SpatialDimPlot(
sc_n,
group.by = "annotation",
label = T,
repel = T,
alpha = 0.5,
shape = 22,
label.size = 4,
pt.size.factor = 6,
image.scale = "hires"
# images = name
) +
ggsci::scale_fill_igv()
dev.off()
## annotation ------
Idents(sc_n) <- "Spatial.008um_snn_res.0.6"
cell_cluster <- c(
'0'   = "Spermatids",
'1'   = "Myoid cell", ##
'2'   = "Spermatocytes",
'3'   = "Sertoli cell", #
'4'   = "Leydig cell", #
'5'   = "Spermatids",
'6'   = "Myoid cell",
'7'   = "Leydig cell",
'8'   = "Spermatids",
'9'   = "Leydig cell", ##
'10'  = "Myoid cell",
'11'  = "Myoid cell", ##
'12'  = "Myoid cell",
'13'  = "Myoid cell",
'14'  = "Myoid cell"
)
sc_n <- RenameIdents(sc_n, cell_cluster)
sc_n$annotation <- Idents(sc_n)
# spatial-dimplot
pdf(file.path(dr_outdir, str_glue("sc_n SpatialDimplot dim{dims} annotation - pca.pdf")), height = 12, width = 12)
SpatialDimPlot(
sc_n,
group.by = "annotation",
label = T,
repel = T,
alpha = 0.5,
shape = 22,
label.size = 4,
pt.size.factor = 6,
image.scale = "hires"
# images = name
) +
ggsci::scale_fill_igv()
dev.off()
sc_n %>% ElbowPlot()
sc_n$nCount_Spatial.008um
sc_n$nCount_Spatial.008um %>% density()
sc_n$nCount_Spatial.008um %>% density() %>% as.array.default()
sc_n$nCount_Spatial.008um %>% density() %>% as.data.frame()
sc_n$nCount_Spatial.008um %>% density() %>% .$x
sc_n$nCount_Spatial.008um %>% density() %>% {data.frame(x = .$x, y = .%y)}
sc_n$nCount_Spatial.008um %>% density() %>% {data.frame(x = .$x, y = .$y)}
sc_n$nCount_Spatial.008um %>% density() %>% {data.frame(x = .$x, y = .$y)} %>% head
sc_n$nCount_Spatial.008um %>% density() %>% {data.frame(x = .$x, y = .$y)} -> ll
ll
ll %>% head
den_df_delta <-
data.frame(
delda_x = den_df$x[-1] - den_df$x[-nrow(den_df)],
delda_y = den_df$y[-1] - den_df$y[-nrow(den_df)]
)
den_df <- sc_n$nCount_Spatial.008um %>% density() %>% {data.frame(x = .$x, y = .$y)}
den_df_delta
den_df$x[-1] - den_df$x[-nrow(den_df)]
den_df$y[-1] - den_df$y[-nrow(den_df)]
data.frame(
delda_x = den_df$x[-1] - den_df$x[-nrow(den_df)],
delda_y = den_df$y[-1] - den_df$y[-nrow(den_df)]
)
den_df_delta <-
data.frame(
delda_x = den_df$x[-1] - den_df$x[-nrow(den_df)],
delda_y = den_df$y[-1] - den_df$y[-nrow(den_df)]
)
den_df_delta
den_df_delta %>% head
den_df_delta %>% plot
den_df_delta %>% min
d_df <- den_df_delta$delda_y/den_df_delta$delda_x
d_df %>% plot
sc_n$nCount_Spatial.008um %>% density()  %>% plot
plot(den_df$x[-1], d_df)
library(argparser)
library(Seurat)
library(tidyverse)
history()
