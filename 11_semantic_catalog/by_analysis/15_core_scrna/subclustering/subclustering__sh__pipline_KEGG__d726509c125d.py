import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_find_marker import find_marker_gene

# kegg富集分析
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC cells/PTC cells.h5ad")
species="mouse"
output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC_Enrich/"

find_marker_gene(adata=adata_raw,group='cell_type',species=species,out_dir=output)