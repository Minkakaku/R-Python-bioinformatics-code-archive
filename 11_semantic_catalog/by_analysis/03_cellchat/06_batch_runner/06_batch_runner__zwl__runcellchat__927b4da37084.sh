Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
--rds "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/Main_mapping.rds" \
--outpath /mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/ \
--splitby group2 \
--species mouse \
--comparison "OVX_CNS-defect:sham_CNS-defect"

Rscript "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/Main_mapping.rds" \
--outpath /mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/ \
--splitby sample1 \
--species mouse \
--comparison "OVX_defect:sham_defect"