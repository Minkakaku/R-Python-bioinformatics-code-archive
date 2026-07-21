import os
import shutil
import random
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import streamlit as st
from PIL import Image
import celltypist
from celltypist import models
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

# import pandas as pd
import scipy.io
import anndata
from scipy.sparse import csr_matrix


import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA_platform/scVision5/tmp/.scanpy_cached"
sc.settings.set_figure_params(figsize=(5,5), facecolor='white')


@st.fragment
def make_report_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for i in ['01-Dnbc4tools', '02-Data_QC', '03-Cell_filter', '04-Dimension_reduction_cluster', '05-Marker_gene_analysis', '06-Cell_type_annotation']:
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)

def read_anndata_C4(path):
    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat ,dtype="float32")
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[0] .values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    adata.obs_names = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.var_names_make_unique()
    return adata



@st.fragment
def get_C4_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    for k,v in sample_dict.items():
        for i in v:
            name =i
            if os.path.exists(os.path.join(data_path, name)):
                path = os.path.join(data_path, name, 'output', 'filter_matrix')
                data = read_anndata_C4(path)
                # data = sc.read_10x_mtx(os.path.join(data_path, name, 'outs', 'filtered_feature_bc_matrix'), var_names='gene_symbols', cache=True)
                data.var_names_make_unique()
                data.obs['sample'] = name
                data.obs['group'] = k
                adatas_list.append(data)
                os.makedirs(os.path.join(out_dir, "01-Dnbc4tools", name), exist_ok=True)
                if not os.path.exists(os.path.join(out_dir,"01-Dnbc4tools", name, "filter_matrix")):
                    os.system(f'cp {os.path.join(data_path, name, "output", "metrics_summary.xls")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                    # os.system(f'cp {os.path.join(data_path, name, "output", "cloupe.cloupe")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                    os.system(f'cp {os.path.join(data_path, name, "output", name+"_scRNA_report.html")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                    os.system(f'cp -r {os.path.join(data_path, name, "output", "filter_matrix")} {os.path.join(out_dir, "01-Dnbc4tools", name)}')
                else:
                    print("File exists")
    adata = adatas_list[0].concatenate(adatas_list[1:], index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    del adata.obs['batch']

    return adata

