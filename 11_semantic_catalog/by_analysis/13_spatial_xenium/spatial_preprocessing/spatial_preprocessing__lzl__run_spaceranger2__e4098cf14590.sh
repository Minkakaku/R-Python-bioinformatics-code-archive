#!/bin/bash

/mnt/NexData_file/home/ztt/software/spaceranger-3.1.2/bin/spaceranger count \
	--id=KB4_1-KB4_2-KB5_1 \
	--transcriptome=/home/ztt/database/refdata-gex-GRCh38-2020-A \
	--fastqs=/mnt/NexData_file/nas/rawdata/rawdata-zc-2025/VH2508184396-04/KB4_1-KB4_2-KB5_1 \
	--sample=VH2508184396-KB4_1-KB4_2-KB5_1 \
	--image=/mnt/NexData_file/home/lzl/VisiumHD_data/VH2508184396-04/images/KB4_1-KB4_2-KB5_1/HD2-2.tif \
	--slide=H1-7DMZFXH \
	--area=D1 \
	--localcores=64 \
	--localmem=256 \
	--create-bam false \
	--probe-set /mnt/NexData_file/home/ztt/software/spaceranger-3.1.2/external/tenx_feature_references/targeted_panels/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv  \
	--cytaimage=/mnt/NexData_file/home/lzl/VisiumHD_data/VH2508184396-04/images/KB4_1-KB4_2-KB5_1/CAVG10719_2026-05-22_11-21-02_2026-05-22_10-28-14_H1-7DMZFXH_D1_hd-2.tif \
	--loupe-alignment=/mnt/NexData_file/home/lzl/VisiumHD_data/VH2508184396-04/images/KB4_1-KB4_2-KB5_1/H1-7DMZFXH-D1-fiducials-image-registration.json
