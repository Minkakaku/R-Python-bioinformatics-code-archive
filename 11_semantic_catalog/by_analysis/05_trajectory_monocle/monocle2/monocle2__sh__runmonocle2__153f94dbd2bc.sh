Rscript "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/Monocle2_pseudotime.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/monocle/"

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/Monocle2_pseudotime.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/monocle_sub/" \
--RootState 4 \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/CytoTRACE_cytotrace.R" \
-r "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds" \
-o "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/CytoTRACE_sub/"

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/Monocle2_pseudotime.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/sc_mono.rds" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/monocle_bp/" \