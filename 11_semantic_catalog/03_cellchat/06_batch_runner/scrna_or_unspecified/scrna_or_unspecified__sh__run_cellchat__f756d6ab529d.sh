Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/adata_mapped_celltype.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2512227671-01/02.Mission/" \
--subcelltype "PT_Camk1d,PT_Cmss1,PT_Gramd1b,PT_Slc5a12,PT_Slc13a3,DC,Mac_Cmss1,Mac_Glis1,Mac_Mrc1,Mac_Tgfbr1,Mac_lsaamp,Mono,pDC" \
--species  mouse \
--comparison "JS:K"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/adata_mapped_celltype.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2512227671-01/02.Mission/" \
--subcelltype "PT_Camk1d,PT_Cmss1,PT_Gramd1b,PT_Slc5a12,PT_Slc13a3,DC,Mac_Cmss1,Mac_Glis1,Mac_Mrc1,Mac_Tgfbr1,Mac_lsaamp,Mono,pDC" \
--species  mouse \
--comparison "YS:K"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/adata_mapped_celltype.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2512227671-01/02.Mission/" \
--subcelltype "PT_Camk1d,PT_Cmss1,PT_Gramd1b,PT_Slc5a12,PT_Slc13a3,Fib_Adgrd1,Fib_Cmksr2,Fib_Cmss1,Fib_Csmd1,Fib_Mbnl3,Fib_Mug1,Fib_Nav3,Fib_Ndst3,Fib_Plin1" \
--species  mouse \
--comparison "SQM:K"

Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/adata_mapped_celltype.rds" \
--outpath "/mnt/NexData_file/home/sh/Projects/AE2512227671-01/02.Mission/" \
--subcelltype "PT_Camk1d,PT_Cmss1,PT_Gramd1b,PT_Slc5a12,PT_Slc13a3,Fib_Adgrd1,Fib_Cmksr2,Fib_Cmss1,Fib_Csmd1,Fib_Mbnl3,Fib_Mug1,Fib_Nav3,Fib_Ndst3,Fib_Plin1" \
--species  mouse \
--comparison "ZC:K"

# levels(adata$cell_type)<-c("B cells","Endothelial cells","FABP5 Epi","Fibroblasts","KRT5 Epi","MUC16 Epi","Macrophages","Mast cells","Fibroblasts","Fibroblasts","Neutrophil","Plasma cells","SPINK1 Epi","SPRR3 Epi","Smooth muscle cells","T cells","Fibroblasts","WFDC2 Epi")