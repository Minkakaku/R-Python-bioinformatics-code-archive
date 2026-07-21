Rscript /mnt/NexData_file/home/wmx/Project/SC2601290784/shouhou8/Monocle.R \
  -r "/mnt/NexData_file/home/wmx/Project/SC2601290784/shouhou8/data/Macro_mono_fixed.rds" \
  -l 5 \
  -o //mnt/NexData_file/home/wmx/Project/SC2601290784/shouhou8/Monocle_result/Macro_mono/ \
  -P cell_type,state,pseudotime,sample,group \
  -R TRUE \
  -p 0.5 \
  -M 1000 \
  -S 2 \
  -G "Ly6c2,Ccr2,Pparg,Marco,Mertk,Axl,Anxa1,Fpr2,Il1b,Ccl2,Spp1,Areg,Igf1"
