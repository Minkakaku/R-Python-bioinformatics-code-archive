suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
})

pkgs <- c(
  "Seurat", "qs", "GSVA", "GSEABase", "pheatmap", "ComplexHeatmap",
  "msigdbr", "limma", "openxlsx", "ggplot2", "dplyr", "tidyr", "tibble",
  "Matrix", "clusterProfiler", "enrichplot"
)

cat("R:", R.version.string, "\n")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(p)) else "NA"
  cat(p, ok, ver, "\n")
}

obj_path <- "/mnt/NexData_file/home/lyh/sc_RNAseq_aftersale/Quanming/20260615_Epi/Epi_recluster_4subgroups/Epi_recluster.qs"
obj <- qread(obj_path)
cat("object dim:", paste(dim(obj), collapse = " x "), "\n")
cat("assays:", paste(Assays(obj), collapse = ","), "\n")
cat("default_assay:", DefaultAssay(obj), "\n")
cat("reductions:", paste(Reductions(obj), collapse = ","), "\n")
cat("metadata columns:\n")
print(colnames(obj@meta.data))
if ("Epi_subcluster" %in% colnames(obj@meta.data)) {
  cat("Epi_subcluster counts:\n")
  print(table(obj$Epi_subcluster, useNA = "ifany"))
}
if ("Group" %in% colnames(obj@meta.data)) {
  cat("Group counts:\n")
  print(table(obj$Group, useNA = "ifany"))
}
