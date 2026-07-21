Rscript "/home/sh/sc_pipline/monocle2/Monocle2_pseudotime.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/monocle2/" \
--clusters "Non-classical monocyte,PTTG1 classical monocyte,ZFP36L1 classical monocyte,S100A8 classical monocyte,Classical monocyte"

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/script/Monocle2_pseudotime.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/monocle2/" \
--clusters "Non-classical monocyte,PTTG1 classical monocyte,ZFP36L1 classical monocyte,S100A8 classical monocyte,Classical monocyte" \
--RootState 4