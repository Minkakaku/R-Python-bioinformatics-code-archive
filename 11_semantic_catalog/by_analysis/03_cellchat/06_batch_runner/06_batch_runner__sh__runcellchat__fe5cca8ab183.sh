Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20250829/04.Mission/" \
--species  human \
--comparison "case:control"

# bdata<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/all_cells_features.xls",sep="\t",header=TRUE,row.names=1)
# adata_cell_names <- rownames(bdata)
# bdata_cell_names <- rownames(adata@meta.data)
# adata_cell_types <- as.character(bdata$cell_type)
# bdata_cell_types <- as.character(adata@meta.data$cell_type)
# mapped_indices <- match(bdata_cell_names, adata_cell_names)
# matched_cell_types <- adata_cell_types[mapped_indices]
# bdata_cell_types[!is.na(mapped_indices)] <- matched_cell_types[!is.na(mapped_indices)]
# adata@meta.data$cell_type <- as.factor(bdata_cell_types)
# Idents(adata) <- adata@meta.data$sub_cell_type

# p <- DimPlot(adata, reduction = 'umap',group.by='cell_type', label = FALSE, pt.size = 0.5)
# pdf('/mnt/NexData_file/home/sh/Projects/FWSC20250829/04.Mission/cell_type.pdf',width=11,height=7)
# print(p)
# dev.off()
# png('/mnt/NexData_file/home/sh/Projects/FWSC20250829/04.Mission/cell_type.png',width=11*80,height=7*80)
# print(p)
# dev.off()
# adata@meta.data$raw_cell_type <- adata@meta.data$cell_type
# adata@meta.data$cell_type <- adata@meta.data$sub_cell_type
# adata@meta.data <- adata@meta.data[, !names(adata@meta.data) %in% c("sub_cell_type","nFeature_RNA","n_genes")]