#adata_raw = sc.read_h5ad("./data/adata_filtered.h5ad")
#adata = sc.read_h5ad("./data/adata_celltype.h5ad")
#adata_raw.obs["cell_id_FOV"] = adata_raw.obs["fov"].astype(str) + '_'+adata_raw.obs["cell_id"].astype(str)
#adata_raw.obs = adata_raw.obs.set_index('cell_id_FOV')
#adata_raw.X = adata_raw.X.astype(int).copy()
#adata_raw.X = adata_raw.X.astype(int).copy()
#adata.obs[["umap1",'umap2']] = adata.obsm['X_umap'].copy()
#adata.obs["cell_id_FOV"] = adata.obs["fov"].astype(str) + '_'+adata.obs["cell_id"].astype(str)
#adata.obs = adata.obs.set_index('cell_id_FOV')
#adata.obs.to_parquet('meta.parquet')
pacman::p_load(anndataR,tidyverse,arrow,qs,Seurat,data.table)
## ----------------------------
## 0. User settings
## ----------------------------
work_dir <- "/mnt/NexData_file/home/lyh/Xenium_aftersale/Sundantong/20260320_report/R_OUT"
Xenium_orig_path <- "/data/xenium/XE2601190441-01"
anndata_raw_path <- "/mnt/NexData_file/home/lyh/Xenium_aftersale/Sundantong/20260320_report/output/data/adata_forR.h5ad"
metadata_path <- "/mnt/NexData_file/home/lyh/Xenium_aftersale/Sundantong/20260320_report/meta.parquet"
sample_info_path <- "/mnt/NexData_file/home/lyh/Xenium_aftersale/Sundantong/20260320_report/sample_info.txt"

## ----------------------------
## setwd
## ----------------------------
if (!dir.exists(work_dir)) {
  dir.create(work_dir)
}
setwd(work_dir)
getwd()

obj <- read_h5ad(anndata_raw_path,as = "Seurat")

LayerData(obj,layer = "counts") <- LayerData(obj,layer = "X")
LayerData(obj,layer = "X") <- NULL

meta <- read_parquet(metadata_path)  %>% as.data.frame()  %>% column_to_rownames("cell_id_FOV")
tmp_col <- setdiff(colnames(meta),colnames(obj@meta.data))
obj <- AddMetaData(obj,meta[,tmp_col])

#set umap reduciton
umap_mat <- as.matrix(meta[, c("umap1", "umap2")])
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
umap_mat <- umap_mat[colnames(obj), ]
obj[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_",assay = DefaultAssay(obj))

obj@meta.data[,c('X','Y')]  <- obj@reductions$spatial@cell.embeddings

obj <- NormalizeData(obj)

qsave(obj,"reprot_cell_anno.qs")

#for boundreys
proj_name <- read.table(sample_info_path, header = TRUE, sep = "\t", check.names = FALSE)
proj_map <- unique(proj_name[, c("fov","raw")])
boundaries_df <- rbindlist(
    lapply(seq_len(nrow(proj_map)), function(k) {
    fov_i <- proj_map$fov[k]
    raw_i <- proj_map$raw[k]
    tmp <- file.path(Xenium_orig_path, raw_i, "cell_boundaries.parquet")
    dt <- as.data.table(read_parquet(tmp))
    dt[, fov := paste0(fov_i, "_", cell_id)]
  }),
  use.names = TRUE,
  fill = TRUE
)
setDF(boundaries_df)
boundaries_df <- boundaries_df  %>%  left_join(meta %>% select(sample) %>% rownames_to_column('fov'))  %>% filter(!is.na(sample))  %>% select(-fov)
boundary_list <- split(boundaries_df, boundaries_df$sample)

qsave(boundary_list,'boundary_list.qs')
