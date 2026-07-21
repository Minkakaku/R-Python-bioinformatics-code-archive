
dir.create("Gene Expression - 20251201")

pdf("Gene Expression - 20251201/Zbp1 vlnplot.pdf", width = 12, height = 6)
VlnPlot(all_data, features = "Zbp1", group.by = "cell_type_res05", pt.size = 0, raster = FALSE)
VlnPlot(all_data, features = "Zbp1", group.by = "cell_type_res05", pt.size = 0.1, raster = FALSE)
VlnPlot(all_data, features = "Zbp1", group.by = "group", pt.size = 0, raster = FALSE)
VlnPlot(all_data, features = "Zbp1", group.by = "group", pt.size = 0.1, raster = FALSE)
dev.off()


pdf("Gene Expression - 20251201/Zbp1 featureplot.pdf", width = 14, height = 10)
p_list <- FeaturePlot(all_data, features = "Zbp1", reduction = "full.harmonyumap", raster = FALSE, split.by = "group", combine = FALSE)
combined_plot <- wrap_plots(p_list, ncol = 2)
print(combined_plot)
dev.off()


pdf("Gene Expression - 20251201/Zbp1 dotplot.pdf", width = 8, height = 6)
DotPlot(all_data, features = "Zbp1", group.by = "cell_type_res05", dot.scale = 10)
dev.off()




pdf("Gene Expression - 20251201/Zbp1 Mac vlnplot.pdf", width = 12, height = 6)
VlnPlot(sc_mac_1, features = "Zbp1", group.by = "RNA_snn_res.0.4", pt.size = 0, raster = FALSE)
VlnPlot(sc_mac_1, features = "Zbp1", group.by = "RNA_snn_res.0.4", pt.size = 0.1, raster = FALSE)
VlnPlot(sc_mac_1, features = "Zbp1", group.by = "group", pt.size = 0, raster = FALSE)
VlnPlot(sc_mac_1, features = "Zbp1", group.by = "group", pt.size = 0.1, raster = FALSE)
dev.off()


pdf("Gene Expression - 20251201/Zbp1 Mac featureplot.pdf", width = 14, height = 10)
p_list <- FeaturePlot(sc_mac_1, features = "Zbp1", reduction = "harmonyumap", raster = FALSE, split.by = "group", combine = FALSE)
combined_plot <- wrap_plots(p_list, ncol = 2)
print(combined_plot)
dev.off()


pdf("Gene Expression - 20251201/Zbp1 Mac dotplot.pdf", width = 8, height = 6)
DotPlot(sc_mac_1, features = "Zbp1", group.by = "RNA_snn_res.0.4", dot.scale = 10)
dev.off()
