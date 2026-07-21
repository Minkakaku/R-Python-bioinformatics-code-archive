Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_cell_type.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_all/" \
--species  mouse \
--comparison "F" \

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_sub_cell_type.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_sub/" \
--species  mouse \
--comparison "F" \

Rscript "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_cell_type.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/" \
--species  mouse \
--comparison "aday1:sham,aday3:sham,aday7:sham,aday14:sham" \