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
    for i in ['01-seeksoultools', '02-Data_QC', '03-Cell_filter', '04-Dimension_reduction_cluster', '05-Marker_gene_analysis', '06-Cell_type_annotation']:
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)
    return out_dir


@st.fragment
def get_xunyin_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    for i in os.listdir(data_path):
        name = i.split('.')[0]
        for k,v in sample_dict.items():
            print(k,v)
            if name in v:
                os.makedirs(os.path.join(out_dir, "01-seeksoultools", name), exist_ok=True)
                if not os.path.exists(os.path.join(out_dir, "01-seeksoultools", name, 'filtered_feature_bc_matrix')):
                    os.system(f'cp {os.path.join(data_path, name, name, name+"_summary.csv")} {os.path.join(out_dir, "01-seeksoultools", name)}')
                # os.system(f'cp {os.path.join(data_path, name, "output", "cloupe.cloupe")} {os.path.join(out_dir, "01-seeksoultools", name)}')
                    os.system(f'cp {os.path.join(data_path, name, name, name+"_report.html")} {os.path.join(out_dir, "01-seeksoultools", name)}')
                    os.system(f'cp -r {os.path.join(data_path,name,name, "step3", "filtered_feature_bc_matrix")} {os.path.join(out_dir, "01-seeksoultools", name)}')
                else:
                    print('File exist')
                path = os.path.join(out_dir, "01-seeksoultools", name, 'filtered_feature_bc_matrix')
                #data = read_anndata_xunyin(path)
                data = sc.read_10x_mtx(path, var_names='gene_symbols', cache=False)
                data.var_names_make_unique()
                data.obs['sample'] = name
                data.obs['group'] = k
                adatas_list.append(data)
                
    adata = adatas_list[0].concatenate(adatas_list[1:], index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    del adata.obs['batch']
    #adata.write_h5ad(os.path.join(out_dir,'adata.h5ad'))

    return adata

