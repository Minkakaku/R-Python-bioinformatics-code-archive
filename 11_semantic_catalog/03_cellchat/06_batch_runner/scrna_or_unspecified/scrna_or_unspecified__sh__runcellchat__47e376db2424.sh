Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/data/NAS02/sh/Projects/FWSC2024-192/08.Mission/annotation2/all_data_sc_doublet.rds" \
--outpath /data/NAS02/sh/Projects/FWSC2024-192/08.Mission/ \
--species  human \
--comparison "LAHS:NLAHS" \
--subcelltype "CD14 Monocyte,FCGR3A Monocyte,Dendritic cell,TPX2 Monocyte,CD8_Trm,CD8_Tex,CD8_MAIT,CD8_TNaive,CD8_Tem"