#!/bin/bash

/mnt/NexData_file/home/ztt/software/spaceranger-3.1.2/bin/spaceranger count \
	--id=24M-16 \
	--transcriptome=/home/ztt/database/refdata-gex-mm10-2020-A \
	--fastqs=/mnt/NexData_file/nas/rawdata/rawdata-sp/FWSP20250661/FWSP20250661-24M-16 \
	--sample=FWSP20250661-24M-16 \
	--image=/mnt/NexData_file/home/ztt/spatial_projects/visium_hd/FWSP20250661_mouse_knee/images/24M_16/24M-16.tif \
	--slide=H1-K22YHYH \
	--area=D1 \
	--localcores=64 \
	--localmem=256 \
	--create-bam false \
	--probe-set /mnt/NexData_file/home/ztt/software/spaceranger-3.1.2/external/tenx_feature_references/targeted_panels/Visium_Mouse_Transcriptome_Probe_Set_v2.0_mm10-2020-A.csv  \
	--cytaimage=/mnt/NexData_file/home/ztt/spatial_projects/visium_hd/FWSP20250661_mouse_knee/images/24M_16/CAVG10719_2025-02-25_11-15-22_2025-02-25_10-38-48_H1-K22YHYH_D1_24m-16.tif \
	--loupe-alignment=/mnt/NexData_file/home/lzl/VisiumHD_data/FWSP20250661_mouse_knee_spaceranger/images/24M_16/H1-K22YHYH-D1-fiducials-image-registration.json
