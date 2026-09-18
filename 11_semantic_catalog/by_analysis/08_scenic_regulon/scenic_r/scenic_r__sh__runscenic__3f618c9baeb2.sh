python "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/script/scenic.py" \
--adata_path "/mnt/NexData_file/scRNA/scVision_output/SC2507163761_fin/result.h5ad" \
--output_dir "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/SENIC/" \
--species 'mouse' \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/scenic.r" \
-o "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/SENIC/" \
-u "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/SENIC/umap_coordinates.csv" \
-m "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/SENIC/metadata.csv" \
-c 'group' \
-l "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/SENIC/sce.loom"

# python "/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/scenic.py" \
# --adata_path "/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/Myeloid.h5ad" \
# --output_dir "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/SENIC_cell_type/" \
# --species 'human' \