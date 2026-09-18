Rscript /mnt/NexData_file/home/wmx/Project/SC2510175738-05/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/SC2510175738-05/data/Spermatids_fixed.rds" \
  -l 2 \
  -o "/mnt/NexData_file/home/wmx/Project/SC2510175738-05/Monocle/" \
  -R TRUE \
  -p 0.6 \
  -M 1000 \
  -P cell_type,state,pseudotime,group
