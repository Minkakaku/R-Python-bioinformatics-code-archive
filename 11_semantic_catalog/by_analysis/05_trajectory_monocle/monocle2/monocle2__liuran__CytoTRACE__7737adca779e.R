# singularity exec --bind /data,/mnt:/mnt/NexData_file/home/liuran/Container/Single_cell_mRNA/ /mnt/NexData_file/home/liuran/Container/Personal_Sif/ACGT101_10X_mRNA_3.4_yzh/ACGT101_10X_mRNA_3.4_yzh.sif R

library(reticulate)
use_condaenv(condaenv = "/home/liuran/software/miniconda3/envs/rpy/bin/python")

reticulate::py_config()
py_available()
reticulate::py_module_available("scanoramaCT")
reticulate::py_module_available("numpy")

library(CytoTRACE)
library(Seurat)
library(SCOPfunctions)
library(dplyr)

