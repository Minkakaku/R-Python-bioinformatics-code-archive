from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import os, sys
import loompy as lp
import numpy as np
import anndata as ad
output = "/home/css/scRNA"
adata = ad.read_h5ad('/data/NAS01/css/scRNAseq/FW2024-428_organoids/FW2024-428-01/result.h5ad')
row_attrs = {"Gene": np.array(adata.raw.var_names),}
col_attrs = {"CellID": np.array(adata.raw.obs_names)}
lp.create(os.path.join(output,"sce.loom"),adata.raw.X.transpose(),row_attrs,col_attrs)

os.system(f'sh /data/NAS01/css/scRNAseq/FW2024-428_organoids/FW2024-428-01/SCENIC.sh {os.path.join(output,"sce.loom")} {os.path.join(output,"grn_result.csv")} {os.path.join(output,"ctx_result.csv")} {os.path.join(output,"aue_result.loom")}')