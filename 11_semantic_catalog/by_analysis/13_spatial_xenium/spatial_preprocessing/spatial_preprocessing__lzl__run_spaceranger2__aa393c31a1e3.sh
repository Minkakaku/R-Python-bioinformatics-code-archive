#!/bin/bash

/mnt/NexData_file/home/ztt/software/spaceranger-3.1.2/bin/spaceranger count \
	--id=WT2 \
	--transcriptome=/home/ztt/database/refdata-gex-mm10-2020-A \
	--fastqs=/mnt/NexData_file/nas/rawdata/rawdata-sp/VH2604022693-01/WT2 \
	--sample=WT2 \
	--image=/mnt/NexData_file/home/lzl/VisiumHD_data/VH2604022693-01/images/WT2/WT-2-2.tif \
	--slide=H1-26CX8Y3 \
	--area=D1 \
	--localcores=64 \
	--localmem=256 \
	--create-bam false \
	--probe-set /mnt/NexData_file/home/ztt/software/spaceranger-3.1.2/external/tenx_feature_references/targeted_panels/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv  \
	--cytaimage=/mnt/NexData_file/home/lzl/VisiumHD_data/VH2604022693-01/images/WT2/CAVG10719_2026-04-22_14-06-01_2026-04-22_13-49-24_H1-26CX8Y3_D1_A253066WT2.tif \
	--loupe-alignment=/mnt/NexData_file/home/lzl/VisiumHD_data/VH2604022693-01/images/WT2/H1-26CX8Y3-D1-fiducials-image-registration.json
