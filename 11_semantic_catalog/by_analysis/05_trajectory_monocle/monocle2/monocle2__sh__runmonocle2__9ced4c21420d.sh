Rscript "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/script/monocle2.analysis.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/Myeloid cells/Myeloid.rds" \
--output "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/monocle2_group/" \
--cluster cell_type \

Rscript "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/script/Monocle2_pseudotime.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/Myeloid cells/Myeloid.rds" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/monocle2/" \
--clusters "VCAN Classical monocyte","ISG15 Classical monocyte","DUSP1 Classical monocyte","Non-classical monocyte"