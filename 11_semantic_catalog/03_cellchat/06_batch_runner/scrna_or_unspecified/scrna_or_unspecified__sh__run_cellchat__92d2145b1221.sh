Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/Main_raw.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/cellchat_sub/" \
--species  human \
--comparison "ADC:SCC"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/Main.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/cellchat/" \
--species  human \
--comparison "ADC:SCC"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/Main_raw.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/cellchat_sub_F/" \
--species  human \
--comparison "F"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/Main.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/cellchat_F/" \
--species  human \
--comparison "F"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/Main_sub.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2507033447-07/05.Mission/cellchat_sub/" \
--species  human \
--comparison "F"
# levels(adata$cell_type)<-c("B cells","Endothelial cells","FABP5 Epi","Fibroblasts","KRT5 Epi","MUC16 Epi","Macrophages","Mast cells","Fibroblasts","Fibroblasts","Neutrophil","Plasma cells","SPINK1 Epi","SPRR3 Epi","Smooth muscle cells","T cells","Fibroblasts","WFDC2 Epi")