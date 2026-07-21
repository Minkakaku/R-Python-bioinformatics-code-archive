import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

adata # merged xenium object

# Dimension reduce
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
adata.raw = adata.copy()

# HVG
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    batch_key="sample"
)
adata = adata[:, adata.var['highly_variable']]

sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)

sc.pp.pca(adata)
with rc_context({'figure.figsize': (10, 6)}):
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)

# harmony
sc.external.pp.harmony_integrate(adata, key="sample")

# umap
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=20, metric='cosine', use_rep='X_pca_harmony', random_state=10)
sc.tl.umap(adata, min_dist=0.3, spread=1, random_state=10)

# cluster
res = [0.2, 0.4, 0.6, 0.8]
for r in res:
    sc.tl.leiden(
        adata, 
        key_added=f"res_{r}", 
        resolution=res, 
        flavor="igraph"
    )

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(
        adata, 
        color=['sample'] + [ f'res_{r}' for r in res ],
        legend_loc='on data', 
        legend_fontsize=9, 
        legend_fontoutline=3,
        frameon=True,
        # size=size,
        show=False
    )