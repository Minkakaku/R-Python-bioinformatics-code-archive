# 根据常规知识修改root节点

# 1.CD4T , 2个group 分别跑
# 以state 2 为起点，CD4 Naive T
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_root.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD4T/Pseudotime/Pseudotime.monocle.rds" \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD4T/ \
#   -S 2 \
#   -P sub_cell_type,state,pseudotime,sample

# # 以state 6 为起点，CD4 Naive T
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_root.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD4T/Pseudotime/Pseudotime.monocle.rds" \
#   -S 6 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD4T/ \
#   -P sub_cell_type,state,pseudotime,sample

# # 2.CD8T , 2个group 分别跑
# # state 5 为起点，CD8 Naive T
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_root.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD8T/Pseudotime/Pseudotime.monocle.rds" \
#   -S 5 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD8T/ \
#   -P sub_cell_type,state,pseudotime,sample

# state 5 为起点，CD8 Naive T
Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_root.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD8T/Pseudotime/Pseudotime.monocle.rds" \
  -S 5 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD8T/ \
  -P sub_cell_type,state,pseudotime,sample

# 3.BandPlasma , 2个group 分别跑
# 以state 5 为起点，Naive B
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_root.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_BandPlasma/Pseudotime/Pseudotime.monocle.rds" \
#   -S 5 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_BandPlasma/ \
#   -P sub_cell_type,state,pseudotime,sample

# # 以state 4 为起点，Naive B
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_root.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_BandPlasma/Pseudotime/Pseudotime.monocle.rds" \
#   -S 4 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_BandPlasma/ \
#   -P sub_cell_type,state,pseudotime,sample


