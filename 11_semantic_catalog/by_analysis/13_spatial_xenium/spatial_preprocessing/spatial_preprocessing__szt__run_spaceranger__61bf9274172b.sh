#!/bin/bash

spaceranger count \
	--id=control2-disease1 \
	--transcriptome=/home/ztt/database/refdata-gex-GRCh38-2020-A/ \
	--fastqs=/data/rawdata/rawdata-sp/FWSP20240200/rawdata \
	--sample=control2-disease1 \
	--image=/home/ztt/spatial/visium_HD/FWSP20240200/image_realigment/control2-disease1_lpf-sj.tif \
	--slide=H1-9K98NK7 \
	--area=D1 \
	--localcores=64 \
	--localmem=256 \
	--create-bam false \
	--probe-set /home/ztt/software/spaceranger-3.0.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
	--cytaimage=/home/ztt/spatial/visium_HD/FWSP20240200/image_realigment/CAVG10719_2024-07-04_14-15-50_2024-07-04_14-00-35_H1-9K98NK7_D1_lpf-sj.tif


spaceranger count \
	--id=disease2-disease3 \
	--transcriptome=/home/ztt/database/refdata-gex-GRCh38-2020-A/ \
	--fastqs=/data/rawdata/rawdata-sp/FWSP20240200/rawdata \
	--sample=disease2-disease3 \
	--image=/home/ztt/spatial/visium_HD/FWSP20240200/image_realigment/disease2-disease3_gh-zlb.tif \
	--slide=H1-9K98NK7 \
	--area=A1 \
	--localcores=64 \
	--localmem=256 \
	--create-bam false \
	--probe-set /home/ztt/software/spaceranger-3.0.0/probe_sets/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
	--cytaimage=/home/ztt/spatial/visium_HD/FWSP20240200/image_realigment/CAVG10719_2024-07-04_14-15-50_2024-07-04_14-00-35_H1-9K98NK7_A1_gh-zlb.tif \
	--loupe-alignment=/home/ztt/spatial/visium_HD/FWSP20240200/image_realigment/H1-9K98NK7-A1-fiducials-image-registration.json