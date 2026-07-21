import os
import re
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
from matplotlib import font_manager as fm
import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA_platform/scVision5/tmp/.scanpy_cached"
sc.settings.set_figure_params(figsize=(5,5), facecolor='white')

plt.tight_layout()
@st.fragment
def make_report_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for i in ['01-CellRanger', '02-Data_QC', '03-Cell_filter', '04-Dimension_reduction_cluster', '05-Marker_gene_analysis', '06-Cell_type_annotation']:
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)
    return out_dir


@st.fragment
def get_10x_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    for k,v in sample_dict.items():
        for i in v:
            name =i
            if os.path.exists(os.path.join(data_path, name)):
                os.makedirs(os.path.join(out_dir, "01-CellRanger", name), exist_ok=True)
                if not os.path.exists(os.path.join(out_dir, "01-CellRanger", name, 'filtered_feature_bc_matrix')):
                    os.system(f'cp {os.path.join(data_path, name, "outs", "metrics_summary.csv")} {os.path.join(out_dir, "01-CellRanger", name)}')
                    os.system(f'cp {os.path.join(data_path, name, "outs", "cloupe.cloupe")} {os.path.join(out_dir, "01-CellRanger", name)}')
                    os.system(f'cp {os.path.join(data_path, name, "outs", "web_summary.html")} {os.path.join(out_dir, "01-CellRanger", name)}')
                    os.system(f'cp -r {os.path.join(data_path, name, "outs", "filtered_feature_bc_matrix")} {os.path.join(out_dir, "01-CellRanger", name)}')
                else:
                    print("File exists")
                data = sc.read_10x_mtx(os.path.join(out_dir, "01-CellRanger",name,'filtered_feature_bc_matrix'), var_names='gene_symbols', cache=False)

                data.var_names_make_unique()
                data.obs['sample'] = name
                data.obs['group'] = k
                adatas_list.append(data)
                
    adata = adatas_list[0].concatenate(adatas_list[1:], index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    del adata.obs['batch']

    return adata

