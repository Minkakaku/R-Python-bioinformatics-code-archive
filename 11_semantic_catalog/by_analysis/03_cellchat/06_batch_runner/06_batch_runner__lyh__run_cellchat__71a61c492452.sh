nohup singularity exec --bind /home/lyh/ /mnt/sif/my_sandbox/  Rscript /home/lyh/Script/Cellchat_single.R  -i /home/lyh/sc_RNA_seq/Xiao/T_cellsanno/anno_1.qs -o /home/lyh/sc_RNA_seq/Xiao/Cellchat/MI_PD1  -s mouse -c sub --subset_col Group  --subset_group  MI_PD1   --python_path /usr/local/sc_project/.venv/bin/python --cores 1  > work.log &

nohup singularity exec --bind /home/lyh/ /mnt/sif/my_sandbox/  Rscript /home/lyh/Script/Cellchat_single.R  -i /home/lyh/sc_RNA_seq/Xiao/T_cellsanno/anno_1.qs -o /home/lyh/sc_RNA_seq/Xiao/Cellchat/WT_PD1  -s mouse -c sub --subset_col Group  --subset_group  WT_PD1   --python_path /usr/local/sc_project/.venv/bin/python --cores 1  > work1.log &


# compare
nohup singularity exec --bind /home/lyh/ /mnt/sif/my_sandbox/  Rscript /home/lyh/Script/Cellchat_compare.R --input_ctrl  /home/lyh/sc_RNA_seq/Xiao/Cellchat/WT_PD1/cellchat_object.rds   --input_treat  /home/lyh/sc_RNA_seq/Xiao/Cellchat/MI_PD1/cellchat_object.rds  --name_ctrl  WT_PD1 --name_treat  MI_PD1 -o /home/lyh/sc_RNA_seq/Xiao/Cellchat/Compare --python_path /usr/local/sc_project/.venv/bin/python  > work_compare.log  &