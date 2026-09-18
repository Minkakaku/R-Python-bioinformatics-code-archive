import os, sys
import loompy as lp
import numpy as np
import anndata as ad
output = "/mnt/NexData_file/home/zwl/Project/FW2024-428-02/require4/scenicpy"
adata = ad.read_h5ad('/home/zwl/Project/FW2024-428-02/result.h5ad')
row_attrs = {"Gene": np.array(adata.raw.var_names),}
col_attrs = {"CellID": np.array(adata.raw.obs_names)}
lp.create(os.path.join(output,"sce.loom"),adata.raw.X.transpose(),row_attrs,col_attrs)

os.system(f'sh /mnt/NexData_file/home/zwl/Project/FW2024-428-02/require4/scenicpy/SCENIC/pyscenic.sh {os.path.join(output,"sce.loom")} {os.path.join(output,"grn_result.csv")} {os.path.join(output,"ctx_result.csv")} {os.path.join(output,"aue_result.loom")}')
