Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
 --rds "/mnt/NexData_file/home/wmx/Project/SC2510175738-05/data/All_fixed.rds" \
 --outpath /mnt/NexData_file/home/wmx/Project/SC2510175738-05/cellchat/Original/ \
 --splitby group \
 --species mouse \
 --comparison "poor:control,stress:control"

Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R \
 --rds "/mnt/NexData_file/home/wmx/Project/SC2510175738-05/data/All_mapping_fixed.rds" \
 --outpath /mnt/NexData_file/home/wmx/Project/SC2510175738-05/cellchat/sub_mapping/ \
 --splitby group \
 --species mouse \
 --comparison "poor:control,stress:control"
