Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
--rds "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/Main_mapping/Main_mapping.rds" \
--outpath /mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/ \
--splitby group1 \
--species mouse \
--comparison "OVX_bone-CNS:sham_bone-CNS"

Rscript "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/Main_mapping/Main_mapping.rds" \
--outpath /mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/ \
--splitby sample1 \
--species mouse \
--comparison "OVX_bone:sham_bone,OVX_CNS:sham_CNS"