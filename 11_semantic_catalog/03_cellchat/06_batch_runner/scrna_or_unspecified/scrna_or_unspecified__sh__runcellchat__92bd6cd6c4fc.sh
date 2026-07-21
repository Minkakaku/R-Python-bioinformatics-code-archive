Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
--rds "/mnt/NexData_file/home/sh/Projects/AE2511106372-01/05.Mission/Main_mapping.rds" \
--outpath /mnt/NexData_file/home/sh/Projects/AE2511106372-01/05.Mission/ \
--splitby group \
--species mouse \
--comparison "MI_DZ:MI"

Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
--rds "/mnt/NexData_file/home/sh/Projects/AE2511106372-01/05.Mission/Main_mapping_raw.rds" \
--outpath /mnt/NexData_file/home/sh/Projects/AE2511106372-01/05.Mission/ \
--splitby group \
--species mouse \
--comparison "MI_DZ:Sham,MI:Sham"