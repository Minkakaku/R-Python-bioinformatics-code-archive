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
/lvdata/wzb/pipline/function/SCENIC/mm_mgi_tfs.txt

docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic ctx --num_workers 32 \
--output $3 \
--expression_mtx_fname $1 \
--mode "custom_multiprocessing" \
--annotations_fname /lvdata/wzb/pipline/function/SCENIC/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
$2 \
/lvdata/wzb/pipline/function/SCENIC/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

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
docker run --name scenic_zwl --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-v /mnt/NexData_file/home/wzb/backup:/mnt/NexData_file/home/wzb/backup \
-v /mnt/NexData_file/:/mnt/NexData_file/ \
zwl_rstudio:latest Rscript /mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC/new_code/scenic.r \
-o $5 \
-u './umap_coordinates.csv' \
-m './metadata.csv' \
-c 'cell_type' \
-l $4

#报告
docker run --name scenic_zwl --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-v /mnt/NexData_file/home/wzb/backup:/mnt/NexData_file/home/wzb/backup \
-v /mnt/NexData_file/:/mnt/NexData_file/ \
zwl_rstudio:latest Rscript -e "rmarkdown::render(
  '/mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC/new_code/report_scenic.rmd',
  output_file = 'report.html',
  output_dir = '${5}')
  )" $5  