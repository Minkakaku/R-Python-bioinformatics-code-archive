Rscript "/home/sh/sc_pipline/monocle2/Monocle2_pseudotime.R" \
--rds "/data/NAS02/sh/Projects/FWSC20241287/02.Mission/Main_Fibroblasts_mapping/recluster.rds" \
--outdir "/data/NAS02/sh/Projects/FWSC20241287/02.Mission/monocle2/" \
--clusters "Fibroblasts_Cfh,Fibroblasts_Slc6a6,Mesangial cell" \
--CellsPercent 1

# args$rds="/data/NAS02/sh/Projects/FWSC20241287/02.Mission/Main_Fibroblasts_mapping/recluster.rds"
# args$outdir="/data/NAS02/sh/Projects/FWSC20241287/02.Mission/monocle2/"
# args$clusters="Fibroblasts_Cfh,Fibroblasts_Slc6a6,Mesangial cell"
# args$CellsPercent=1

Rscript "/home/sh/sc_pipline/monocle2/Monocle2_Pseudotime_BEAM.R" \
--rds "/data/NAS02/sh/Projects/FWSC20241287/02.Mission/monocle2/Pseudotime.monocle.rds" \
--outdir "/data/NAS02/sh/Projects/FWSC20241287/02.Mission/beam/" \
--branchPoint 1 \
--genelist "/data/NAS02/sh/Projects/FWSC20241287/02.Mission/monocle2/Pseudotime_diff_all.csv" \
--numClusters 3
