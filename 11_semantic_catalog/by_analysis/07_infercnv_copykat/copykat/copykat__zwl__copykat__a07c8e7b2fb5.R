setwd("/home/zwl/Project/FWSC20241362/09.require3/")
rm(list = ls())
gc()

library(sceasy)
library(Seurat)
library(reticulate)
use_condaenv('scirpy')
#############颜色#################
library(scales)
pal= ggsci::pal_d3('category20')(16)
#names(pal) <- 
show_col(pal)

#############h5ad转变成seurat对象####################
#测试sceasy和diopy哪个好用：都一样
sceasy::convertFormat('/home/zwl/Project/FWSC20241362/09.require3/allgene_anndata.h5ad',
                      from="anndata", to="seurat",outFile='sceasy_test.rds')
obj = readRDS('./sceasy_test.rds')
obj = FindVariableFeatures(obj,assay = 'RNA')
#diopy
obj1 = readRDS('../test.rds')
rm(obj1)

#############