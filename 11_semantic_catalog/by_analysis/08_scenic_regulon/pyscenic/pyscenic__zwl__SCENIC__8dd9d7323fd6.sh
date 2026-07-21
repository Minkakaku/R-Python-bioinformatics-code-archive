#!/bin/bash
docker run --name SCENIC_wzb_SMC --rm \
-v /home/css:/home/css \
-v /data/NAS01/css:/data/NAS01/css \
aertslab/pyscenic:0.12.1 pyscenic grn --num_workers 64 \
--sparse \
--method grnboost2 \
--output $2 \
$1 \
/data/NAS01/css/scRNAseq/FW2024-428_organoids/allTFs_hg38.txt

docker run --name SCENIC_wzb_SMC --rm \
-v /home/css:/home/css \
-v /data/NAS01/css:/data/NAS01/css \
aertslab/pyscenic:0.12.1 pyscenic ctx --num_workers 64 \
--output $3 \
--expression_mtx_fname $1 \
--mode "custom_multiprocessing" \
--annotations_fname /data/NAS01/css/scRNAseq/FW2024-428_organoids/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
$2 \
/data/NAS01/css/scRNAseq/FW2024-428_organoids/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

docker run --name SCENIC_wzb_SMC --rm \
-v /home/css:/home/css \
-v /data/NAS01/css:/data/NAS01/css \
aertslab/pyscenic:0.12.1 pyscenic aucell --num_workers 64 \
--output $4	\
$1	\
$3