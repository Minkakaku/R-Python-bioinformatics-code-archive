Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/MP.rds" \
  -l 8 \
  -R TRUE \
  -p 0.6 \
  -M 1000 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/Monocle_result/MP/ \
  -P MP_group,state,pseudotime,sample

  Rscript /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/MP_remove_unsolved.rds" \
  -l 8 \
  -R TRUE \
  -p 0.6 \
  -M 1000 \
  -o /mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/Monocle_result/MP_remove_unsolved/ \
  -P MP_group,state,pseudotime,sample