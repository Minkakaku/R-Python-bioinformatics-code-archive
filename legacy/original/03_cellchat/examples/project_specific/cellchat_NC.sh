source activate /mnt/NexData_file/home/wyt/software/minicoda/envs/R_envs
Rscript /home/wyt/test/protocol/cellchat_Spatial/muti_group_cellchat_Spatial.R \
--gro NC \
--infle /mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/total_anno/xenium_1_human_skin_anno_drop_low.rds \
--path /mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/cellchat/ \
--assay Xenium \
--seg_type xenium \
--species human \
--cell_type celltype_fac \
--sample_label P1-NC P2-NC P3-NC \
--fov_label  P1-NC P2-NC P3-NC