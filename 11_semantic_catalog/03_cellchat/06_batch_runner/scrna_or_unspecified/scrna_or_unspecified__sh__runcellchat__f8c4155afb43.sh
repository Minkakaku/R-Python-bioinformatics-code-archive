# conda activate "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/"
Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20251028/06.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20251028/06.Mission/cellchat_all.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20251028/06.Mission/" \
--splitby group \
--species mouse \
--comparison "HGF:ELA,HGF:Control,ELA:Control"