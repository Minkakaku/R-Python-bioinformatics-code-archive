#!/usr/bin/env scRNA_p(python 3.8.11)
# -*- coding: utf-8 -*-


# export OPENBLAS_NUM_THREADS=1
# export MKL_NUM_THREADS=1
# export OMP_NUM_THREADS=1
# export NUMEXPR_NUM_THREADS=1

import os
import argparse
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use("Agg")


# force single-thread
import os
# 强制所有计算库使用单线程
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"


# load functions
import sys
sys.path.append('/mnt/NexData_file/share/scRNA-tools/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion




def parse_args():
    parser = argparse.ArgumentParser(description="Subcluster each cell type.")
    parser.add_argument('--adata', required=True, help='Input h5ad file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--labels', nargs='+', required=True, help='Cell types to iterate')
    parser.add_argument('--resolution', type=float, required=True, help='Resolution for recluster')
    parser.add_argument('--n_pcs', type=int, required=True, help='Number of PCs')
    parser.add_argument('--batch_effect', type=int, required=True, help='Batch effect parameter')
    parser.add_argument('--species', required=True, help='human or mouse')
    return parser.parse_args()




def main():
    args = parse_args()


    os.makedirs(args.outdir, exist_ok=True)


    adata = sc.read_h5ad(args.adata)


    for label in args.labels:
        outdir1 = os.path.join(args.outdir, label)
        os.makedirs(outdir1, exist_ok=True)


        adata_tmp = adata[adata.obs['cell_type'] == label].copy()


        adata_tmp = recluster(adata_tmp, species=args.species,
                              neighbors=10,
                              resolution=args.resolution,
                              n_pcs=args.n_pcs,
                              batch_effect=args.batch_effect,output=outdir1)


        adata_tmp = adata_tmp.sub_cluster()




if __name__ == '__main__':
    main()

## use
## conda deactivate ;conda activate /mnt/NexData_file/share/scRNA-rawdata/envs/scRNA_p
#   python /mnt/NexData_file/home/xzx/scripts/recluster.py \
#   --adata /mnt/NexData_file/home/xzx/project/FWSC20252744/shouhou1/adata_new.h5ad \
#   --outdir /mnt/NexData_file/home/xzx/project/FWSC20252744/shouhou1/ \
#   --labels ["Epithelial cells","stromal"] \
#   --resolution 0.5 \
#   --n_pcs 10 \
#   --batch_effect 2 \
#   --species "human"