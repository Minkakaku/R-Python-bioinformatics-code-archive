# # 1.CD4T 
# Rscript /mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/Cd4T_fixed.rds" \
#   -l 8 \
#   -o //mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle_result/CD4T/ \
#   -P cell_type,state,pseudotime,group

# # 2.CD8T 
# Rscript /mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/Cd8T_fixed.rds" \
#   -l 8 \
#   -o //mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle_result/CD8T/ \
#   -P cell_type,state,pseudotime,group

# # 3.Mac ,这个轨迹太丑了，有一些离群点，用mac_clean脚本去除了
# Rscript /mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/mac_fixed.rds" \
#   -l 8 \
#   -o //mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle_result/mac/ \
#   -P cell_type,state,pseudotime,group

# # 4.DC
# Rscript /mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/DC_fixed.rds" \
#   -l 8 \
#   -o //mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle_result/DC/ \
#   -P cell_type,state,pseudotime,group

Rscript /mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/mac_fixed_clean.rds" \
  -l 8 \
  -o //mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle_result/mac_clean/ \
  -P cell_type,state,pseudotime,group
