# Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20241287/06.Mission/script/pip_cellchat.R" \
# --rds "/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/recluster.rds" \
# --outpath "/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/" \
# --species  mouse \
# --comparison "MODEL:CONTROL,PF11-H:MODEL,PF11-H:CONTROL" \
# --subcelltype  "B cells,DC,Distal tubular cells,Endothelial cells,Fibroblasts,Immature Mac C1,Immature Mac C2,Intercalated cells,Mesangial cells,Podocytes,Principal cells,Pro-inflammatory Mac,Pro-repair Mac,Proximal tubule cells,Remodeling Mac-Svep1,T cells"

Rscript "/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/script/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/recluster.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/newcellchat/" \
--species  mouse \
--comparison "MODEL:CONTROL,PF11-H:MODEL,PF11-H:CONTROL" \
--subcelltype  "B cells,DC,Distal tubular cells,Endothelial cells,Fibroblasts,Immature Mac C1,Immature Mac C2,Intercalated cells,Mesangial cells,Podocytes,Principal cells,Pro-inflammatory Mac,Pro-repair Mac,Proximal tubule cells,Remodeling Mac-Svep1,T cells"