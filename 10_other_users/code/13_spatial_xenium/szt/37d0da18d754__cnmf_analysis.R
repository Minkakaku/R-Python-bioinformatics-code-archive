# cNMF ----------------
library(tidyverse)
library(Seurat)

cnmf_out_dir <- "cnmf"
dir.create(cnmf_out_dir)

# all epi cells arge tumor
# sc_epi 
sc_epi$RNA_snn_res.0.05 %>% unique

sc_epi_subset_meta <-
  sc_epi[['RNA_snn_res.0.05']] %>%
  rownames_to_column("cell_id") %>%
  group_by(RNA_snn_res.0.05) %>%
  slice_sample(n = 50000)


# sc_epi_expr_counts <- sc_epi %>% FetchData(layer = "counts", vars = rownames(sc_epi))
sc_epi_expr_counts <- 
  sc_epi[, sc_epi_subset_meta$cell_id] %>%
  JoinLayers() %>%
  LayerData(layer = "counts")
sc_epi_expr_counts <- t(sc_epi_expr_counts)

# colSums
data.frame(count = rowSums(sc_epi_expr_counts)) %>% 
  ggplot(aes(x = 1, y = count)) +
  geom_violin() +
  scale_y_log10()

# hvg_gene
epi_hvg_gene <- sc_epi %>% FindVariableFeatures(assay = "RNA", nfeatures = 4000) %>% VariableFeatures()

# if error occurs: x cells have zero counts of overdispersed genes.
# remove them, rerun cnmf.
sc_epi_expr_counts_subset <- 
  sc_epi_expr_counts[
    rowSums(sc_epi_expr_counts) > 100,
    epi_hvg_gene
  ]

write.table(
  sc_epi_expr_counts_subset,
  file = file.path(cnmf_out_dir, "counts.txt"),
  quote = F,
  sep = "\t",
  row.names = T,
  col.names = T
)

##########run cnmf in linux command window ############

mkdir cnmf && cd cnmf
cnmf prepare --output-dir . --name cNMF -c counts.txt -k 5 6 7 8 9 10 11 12 13 --n-iter 100 --seed 14

cnmf factorize --output-dir . --name cNMF --worker-index 0 --total-workers 1 # 2.5h

cnmf combine --output-dir . --name cNMF

cnmf k_selection_plot --output-dir . --name cNMF

cnmf consensus --output-dir . --name cNMF --components 10 --local-density-threshold 0.01 --show-clustering
cnmf consensus --output-dir . --name cNMF --components 10 --local-density-threshold 0.1 --show-clustering
cnmf consensus --output-dir . --name cNMF --components 8 --local-density-threshold 0.1 --show-clustering
cnmf consensus --output-dir . --name cNMF --components 13 --local-density-threshold 0.1 --show-clustering














