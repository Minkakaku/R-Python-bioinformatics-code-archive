#!/bin/bash
docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic grn --num_workers 32 \
--sparse \
--method grnboost2 \
--output $2 \
$1 \
/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/

docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic ctx --num_workers 32 \
--output $3 \
--expression_mtx_fname $1 \
--mode "custom_multiprocessing" \
--annotations_fname /mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
$2 \
/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic aucell --num_workers 32 \
--output $4	\
$1	\
$3

#可视化
/mnt/NexData_file/share/scRNA-tools/envs/pip_SCENICR/bin/Rscript \
/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/scenic.r \
-o $5 \
-u './umap_coordinates.csv' \
-m './metadata.csv' \
-c 'cell_type' \
-l $4

#报告
cp /mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools/SENIC/new_code/report_scenic.rmd $5
cd $5
/mnt/NexData_file/share/scRNA-tools/envs/pip_SCENICR/bin/Rscript  -e "rmarkdown::render('./report_scenic.rmd',output_file = 'report.html')"