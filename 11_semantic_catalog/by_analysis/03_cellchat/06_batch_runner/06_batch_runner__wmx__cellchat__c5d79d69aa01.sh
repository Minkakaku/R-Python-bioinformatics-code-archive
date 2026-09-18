Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R  \
--rds "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/data/mapping_fixed.rds" \
--outpath /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/cellchat/yaqun/ \
--splitby group \
--species human \
--comparison "SLE:HC" 

Rscript /mnt/NexData_file/share/scRNA-tools/pip_cellchat.R  \
--rds "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/data/mapping_fixed_daqun.rds" \
--outpath /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/cellchat/daqun/ \
--splitby group \
--species human \
--comparison "SLE:HC" 
