#!/bin/bash
docker run --name SCENIC_wzb_SMC --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-v /mnt/NexData_file/home/zwl:/mnt/NexData_file/home/zwl \
-v /home/zwl:/home/zwl \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic grn --num_workers 64 \
--sparse \
--method grnboost2 \
--output $2 \
$1 \
/mnt/NexData_file/home/zwl/Project/FW2024-428-02/require4/scenicpy/SCENIC/allTFs_hg38.txt

docker run --name SCENIC_wzb_SMC --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-v /mnt/NexData_file/home/zwl:/mnt/NexData_file/home/zwl \
-v /home/zwl:/home/zwl \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic ctx --num_workers 64 \
--output $3 \
--expression_mtx_fname $1 \
--mode "custom_multiprocessing" \
--annotations_fname /mnt/NexData_file/home/zwl/Project/FW2024-428-02/require4/scenicpy/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
$2 \
/mnt/NexData_file/home/zwl/Project/FW2024-428-02/require4/scenicpy/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

docker run --name SCENIC_wzb_SMC --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-v /mnt/NexData_file/home/zwl:/mnt/NexData_file/home/zwl \
-v /home/zwl:/home/zwl \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic aucell --num_workers 64 \
--output $4	\
$1	\
$3
