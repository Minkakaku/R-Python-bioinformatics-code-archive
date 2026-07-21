# Rscript /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/data/CD8T.rds" \
#   -o /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/CD8T/ \
#   -R TRUE \
#   -p 0.8 \
#   -M 1000 \
#   -P cell_type,state,pseudotime,sample,group 

# Rscript /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/data/CD4T.rds" \
#   -R TRUE \
#   -p 0.7 \
#   -M 1000 \
#   -o /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/CD4T/ \
#   -P cell_type,state,pseudotime,sample,group 

# Rscript /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/data/B_fixed.rds" \
#   -o /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/B/ \
#   -P cell_type,state,pseudotime,sample,group \

# 查看了 结果之后确定的root state
Rscript /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle_root.R \
  -r "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/CD8T/Pseudotime/Pseudotime.monocle.rds" \
  -o /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/CD8T/ \
  -P cell_type,state,pseudotime,sample,group \
  -S 5 

# 查看了 结果之后确定的root state
Rscript /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle_root.R \
  -r "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/CD4T/Pseudotime/Pseudotime.monocle.rds"  \
  -o /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/CD4T/ \
  -P cell_type,state,pseudotime,sample,group \
  -S 5

# 查看了 结果之后确定的root state
Rscript /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle_root.R \
  -r "/mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/B/Pseudotime/Pseudotime.monocle.rds" \
  -o /mnt/NexData_file/home/wmx/Project/SC2603312640-01/shouhou11/Monocle/B/ \
  -P cell_type,state,pseudotime,sample,group \
  -S 5