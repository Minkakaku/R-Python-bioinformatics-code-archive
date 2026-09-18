Rscript "/mnt/NexData_file/share/scRNA-tools/pip_cellchat.R" \
--rds "/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/mapping_recluster.rds" \
--outpath "/data/NAS02/sh/Projects/FWSC20250579/02.Mission/" \
--species  mouse \
--comparison "KO:CON,CON:KO"

# argv$rds="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/mapping_recluster.rds"
# argv$outpath="/data/NAS02/sh/Projects/FWSC20250579/02.Mission/"
# argv$species="mouse"
# argv$comparison="KO:CON,CON:KO"