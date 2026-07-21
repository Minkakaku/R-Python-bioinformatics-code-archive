Rscript "/data/NAS02/sh/Projects/FWAE20242338/01.Mission/script/Monocle2_pseudotime.R" \
--rds "/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/Neurons/subcluster.rds" \
--lambda  '8'  \
--outpath "/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_monocle2" \
--CellsPercent  '1' \
--clusters "Bergmann glial cell,Excitatory neuron,Inhibitory neuron,Interneuron,Neuroblast"
