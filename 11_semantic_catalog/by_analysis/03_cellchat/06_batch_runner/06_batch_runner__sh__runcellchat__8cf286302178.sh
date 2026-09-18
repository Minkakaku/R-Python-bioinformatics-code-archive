# conda activate "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/"
Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/all_data_sc_doublet_clean.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/" \
--splitby group \
--species human \
--comparison "LAHS:NLAHS,LAHS:NC,NLAHS:NC"