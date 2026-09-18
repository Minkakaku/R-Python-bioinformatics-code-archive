Rscript "/home/sh/sc_pipline/cellchat/cellchat.R" \
--rds "/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.rds" \
--prefix Main \
--outpath "/data/NAS02/sh/Projects/FWSC20241287/03.Mission/cellchat/" \
--species  mouse \
--comparison "CONTROL:MODEL,CONTROL:PF11-H,MODEL:PF11-H"

# argv$rds="/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.rds"
# argv$prefix="Main"
# argv$outpath="/data/NAS02/sh/Projects/FWSC20241287/03.Mission/cellchat/"
# argv$species="mouse"
# argv$comparison="MODEL:PF11-H"