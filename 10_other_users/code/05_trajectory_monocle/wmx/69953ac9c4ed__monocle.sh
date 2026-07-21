# 1.CD4T , 2个group 分别跑
Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/LARC_CD4T.rds" \
  -l 8 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD4T/ \
  -R TRUE \
  -p 0.4 \
  -M 1000 \
  -P sub_cell_type,state,pseudotime,sample

Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/N_CD4T.rds" \
  -l 8 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD4T/ \
  -P sub_cell_type,state,pseudotime,sample

# 2.CD8T , 2个group 分别跑
Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/LARC_CD8T.rds" \
  -l 8 \
  -R TRUE \
  -p 0.6 \
  -M 1000 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_CD8T/ \
  -P sub_cell_type,state,pseudotime,sample

Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/N_CD8T.rds" \
  -l 8 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_CD8T/ \
  -P sub_cell_type,state,pseudotime,sample

# # 3.BandPlasma , 2个group 分别跑
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/LARC_BandPlasma.rds" \
#   -l 8 \
#   -R TRUE \
#   -p 0.6 \
#   -M 1000 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_BandPlasma/ \
#   -P sub_cell_type,state,pseudotime,sample

# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/N_BandPlasma.rds" \
#   -l 8 \
#   -R TRUE \
#   -p 0.4 \
#   -M 1000 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_BandPlasma/ \
#   -P sub_cell_type,state,pseudotime,sample

# # 4.Endo , 2个group 分别跑
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/LARC_Endothelial_cells.rds" \
#   -l 8 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_Endothelial_cells/ \
#   -P sub_cell_type,state,pseudotime,sample

# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/N_Endothelial_cells.rds" \
#   -l 8 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_Endothelial_cells/ \
#   -P sub_cell_type,state,pseudotime,sample

# # 5.Fibro , 2个group 分别跑
# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/LARC_Fibroblasts.rds" \
#   -l 8 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/LARC_Fibroblasts/ \
#   -P sub_cell_type,state,pseudotime,sample

# Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
#   -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_rds/N_Fibroblasts.rds" \
#   -l 8 \
#   -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle_result/N_Fibroblasts/ \
#   -P sub_cell_type,state,pseudotime,sample

