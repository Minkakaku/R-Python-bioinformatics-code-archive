Rscript "/home/sh/sc_pipline/cellchat/cellchat.R" \
--rds "/data/NAS02/sh/Projects/FWSC20250070/01.Mission/mapping_recluster/mapping_recluster.rds" \
--prefix Main \
--outpath "/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Main_cellchat/" \
--species  mouse \
--comparison "control:aki,control:cdn"

# argv$rds="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/mapping_recluster/mapping_recluster.rds"
# argv$prefix="Main"
# argv$outpath="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Main_cellchat/"
# argv$species="mouse"
# argv$comparison="control:cdn"