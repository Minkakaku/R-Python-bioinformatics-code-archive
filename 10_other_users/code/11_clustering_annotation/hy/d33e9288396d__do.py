"""
Virtual Knockout Workflow for Specified Cell Type
Compatible with scTenifoldpy >=0.2.0 + scanpy
"""
import scanpy as sc
import anndata as ad
import pandas as pd
from scTenifold import scTenifoldKnk
from scTenifold.plotting import plot_network_heatmap, plot_qqplot
import matplotlib.pyplot as plt

# ==========================
# 1. Data Loading and Preprocessing
# ==========================
# Load h5ad data
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15.h5ad")

# Filter for target cell type
target_celltype = "WNT+ ISC"
sub_adata = adata[adata.obs["celltype"] == target_celltype].copy()
sub_adata.obs_names_make_unique()  # Resolve duplicate observation name warnings

# Limit maximum number of cells
if sub_adata.n_obs > 5000:
    print(f"Limiting cells from {sub_adata.n_obs} to 5000")
    sub_adata = sub_adata[:5000].copy()

# ---------- Highly Variable Gene Selection ----------
sc.pp.filter_genes(sub_adata, min_cells=5)  # Minimum expressing cells threshold
sc.pp.highly_variable_genes(
    sub_adata,
    n_top_genes=3000,
    flavor="seurat_v3"
)
sub_adata = sub_adata[:, sub_adata.var["highly_variable"]]

print(f"Data preparation complete | Cells: {sub_adata.n_obs} | Genes: {sub_adata.n_vars}")

# ==========================
# 2. Virtual Knockout Configuration
# ==========================
ko_genes = ["CD44"]

# Check if KO genes exist in the filtered dataset
available_ko_genes = []
for gene in ko_genes:
    if gene in sub_adata.var_names:
        available_ko_genes.append(gene)
    else:
        print(f"Warning: Gene '{gene}' not found in highly variable genes. Skipping.")

if not available_ko_genes:
    raise ValueError("No specified KO genes found in the dataset. Please check gene names.")

knk = scTenifoldKnk(
    data=sub_adata,
    ko_genes=available_ko_genes,  # Use only available genes
    # Most permissive QC parameters (no cell filtering)
    qc_kws={
        "min_lib_size": 0,
        "max_mito_ratio": 1.0,
        "min_percent": 0,
        "min_exp_avg": 0,
        "min_exp_sum": 0,
        "plot": False
    },
    nc_kws={
        "n_nets": 100,
        "n_samp_cells": 500,
        "q": 0,
        "backend": "serial"
    },
    td_kws={"K": 2, "max_iter": 200, "init": "random"},
    ma_kws={"d": 2},
    ko_method="default"
)

# Execute workflow steps
for step in ["qc", "nc", "td", "ko", "ma", "dr"]:
    knk.run_step(step)

# ==========================
# 3. Results Extraction and Filtering (FDR < 0.1)
# ==========================
results = knk.d_regulation.copy()
results["significant"] = results["adjusted p-value"] < 0.1

print("\nSignificantly differentially regulated genes (FDR < 0.1):")
print(results[results["significant"]].sort_values("p-value").head(10))

# ==========================
# 4. Save Results
# ==========================
save_prefix = "ko_test"

# Save affected genes table
results.to_csv(f"{save_prefix}_differential_genes.csv", index=False)
print(f"Affected genes table saved as: {save_prefix}_differential_genes.csv")

# Save visualizations
plt.figure(figsize=(8, 6))
plot_network_heatmap(knk.tensor_dict["WT"].to_numpy())
plt.title(f"Wild-Type Network ({target_celltype})")
plt.savefig(f"{save_prefix}_wt_network.png", dpi=300, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8, 6))
plot_network_heatmap(knk.tensor_dict["KO"].to_numpy())
plt.title(f"Knockout Network ({target_celltype}, KO: {available_ko_genes})")
plt.savefig(f"{save_prefix}_ko_network.png", dpi=300, bbox_inches="tight")
plt.close()

plt.figure(figsize=(6, 4))
plot_qqplot(results)
plt.title("Differential Regulation QQ Plot (FDR < 0.1)")
plt.savefig(f"{save_prefix}_qqplot.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"\nAll outputs saved with prefix: {save_prefix}_*")
