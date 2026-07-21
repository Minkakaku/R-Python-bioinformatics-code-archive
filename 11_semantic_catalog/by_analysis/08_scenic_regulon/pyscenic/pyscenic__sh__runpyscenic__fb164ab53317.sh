python "/mnt/NexData_file/share/scRNA-tools/SCENIC/new_code/scenic.py" \
--adata_path "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/Myeloid/M1_2_Monocyte.h5ad" \
--output_dir  "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/scenic/" \
--species human

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/script/scenic.r" \
-o "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/scenic/" \
-u "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/scenic/umap_coordinates.csv" \
-m "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/scenic/metadata.csv" \
-c 'cell_type' \
-l "/mnt/NexData_file/home/sh/Projects/FWSC20252616/06.Mission/scenic/sce.loom"