Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/cellchat_H0/" \
--species  mouse \
--comparison "F" \
--subgroup "H0" \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/cellchat_H2/" \
--species  mouse \
--comparison "F" \
--subgroup "H2" \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/cellchat_H8/" \
--species  mouse \
--comparison "F" \
--subgroup "H8" \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/cellchat_H24/" \
--species  mouse \
--comparison "F" \
--subgroup "H24" \

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/cellchat_H72/" \
--species  mouse \
--comparison "F" \
--subgroup "H72" \
# levels(adata$cell_type)<-c("B cells","Endothelial cells","FABP5 Epi","Fibroblasts","KRT5 Epi","MUC16 Epi","Macrophages","Mast cells","Fibroblasts","Fibroblasts","Neutrophil","Plasma cells","SPINK1 Epi","SPRR3 Epi","Smooth muscle cells","T cells","Fibroblasts","WFDC2 Epi")