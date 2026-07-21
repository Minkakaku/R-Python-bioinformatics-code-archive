python "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/scenic.py" \
--adata_path "/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/Myeloid.h5ad" \
--output_dir "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/SENIC/" \
--species 'human' \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/scenic.r" \
-o "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/SENIC/" \
-u "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/SENIC/umap_coordinates.csv" \
-m "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/SENIC/metadata.csv" \
-c 'group' \
-l "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/SENIC/aue_result.loom"

python "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/scenic.py" \
--adata_path "/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/Myeloid.h5ad" \
--output_dir "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/SENIC_cell_type/" \
--species 'human' \