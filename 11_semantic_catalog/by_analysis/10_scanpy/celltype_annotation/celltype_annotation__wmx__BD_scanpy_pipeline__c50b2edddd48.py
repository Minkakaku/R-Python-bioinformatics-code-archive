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
import csv
import zipfile


import warnings
warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA_platform/scVision5/tmp/.scanpy_cached"
sc.settings.set_figure_params(figsize=(5,5), facecolor='white')


@st.fragment
def make_report_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for i in ['01-Rhapsody', '02-Data_QC', '03-Cell_filter', '04-Dimension_reduction_cluster', '05-Marker_gene_analysis', '06-Cell_type_annotation']:
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)

def extract_metrics(input_file, output_file):
    """
    从BD Rhapsody分析报告中提取关键指标并保存到新文件
    
    参数:
    input_file -- 输入CSV文件路径
    output_file -- 输出CSV文件路径
    """
    # 定义需要提取的关键词列表
    keywords = [
        'Library', 'Putative_Cell_Count', 'Pct_Reads_from_Putative_Cells',
        'Mean_Reads_per_Cell', 'Median_Reads_per_Cell', 'Mean_Molecules_per_Cell',
        'Median_Molecules_per_Cell', 'Mean_Bioproducts_per_Cell',
        'Median_Bioproducts_per_Cell', 'Total_Bioproducts_Detected',
        'Total_CellLabel_UMI_Reads', 'Annotated_Transcriptome_Pct', 'Introns_Pct',
        'Intergenic_Regions_Pct', 'Antisense_Pct', 'Not_Unique_Pct', 'Ambiguous_Pct',
        'Unaligned_Pct', 'Total_Filtered_Reads', 'Pct_Q30_Bases_in_Filtered_R2',
        'Pct_CellLabel_UMI', 'Pct_CellLabel_UMI_Aligned_Uniquely', 'Pct_Reads_Useful',
        'Sequencing_Saturation'
    ]
    
    # 初始化结果字典
    results = {key: "NA" for key in keywords}  # 用"NA"表示未找到的值
    
    # 定义需要处理的节及其对应的数据行特征
    sections = {
        "#Library Quality#": "Total_Filtered_Reads",
        "#Alignment Categories#": "Total_CellLabel_UMI_Reads",
        "#Reads and Molecules#": "Aligned_Reads_By_Type",
        "#Cells#": "Putative_Cell_Count"
    }
    
    try:
        # 读取文件内容
        with open(input_file, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
        
        current_section = None
        headers = []
        
        for line in lines:
            # 跳过注释行
            if line.startswith("##"):
                continue
                
            # 检测节标题
            if line.startswith('#') and line.endswith('#'):
                current_section = line
                continue
            
            # 处理目标节的数据
            if current_section in sections:
                # 识别列标题行
                if line.startswith(sections[current_section]):
                    headers = line.split(',')
                    continue
                
                # 识别数据行
                if headers and line[0].isdigit():
                    data = line.split(',')
                    
                    # 将列名与数据配对
                    for header, value in zip(headers, data):
                        if header in keywords:
                            results[header] = value
                    
                    # 重置当前节处理完成
                    current_section = None
                    headers = []
        
        # 写入结果到新CSV文件
        with open(output_file, 'w', newline='') as f_out:
            writer = csv.writer(f_out)
            # 写入标题行
            writer.writerow(keywords)
            # 写入数据行（按keywords顺序）
            writer.writerow([results[key] for key in keywords])
            
        print(f"成功提取指标到 {output_file}")
        return True
    
    except Exception as e:
        print(f"处理文件时出错: {str(e)}")
        return False

def read_anndata_BD(path):
    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = csr_matrix(mat)
    adata = anndata.AnnData(mat ,dtype="float32")
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[1] .values
    var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[1].values

    barcode = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.obs_names = np.char.add('CellLabel_', barcode.astype(str))
    
    adata.var_names_make_unique()
    return adata

def unzip_mex(path,sample):
    # zip 文件路径
    zip_file = os.path.join(path,sample,f"{sample}_RSEC_MolsPerCell_MEX.zip")
    
    # 目标解压目录
    out_dir = os.path.join(path,sample,"filtered_feature_bc_matrix")
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(zip_file):
        print(f"未找到文件: {zip_file}")
        return
    
    # 解压
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        zip_ref.extractall(out_dir)
        print(f"已解压到: {out_dir}")

@st.fragment
def get_BD_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    for k,v in sample_dict.items():
        for i in v:
            name =i
            if os.path.exists(os.path.join(data_path, name)):
                path = os.path.join(data_path, name, 'filtered_feature_bc_matrix')
                if not os.path.exists(os.path.join(path, "matrix.mtx.gz")):
                    print(f"⚠️ 未找到 {path}/matrix.mtx.gz，尝试解压 zip 文件...")
                    unzip_mex(data_path, name)
                
                # 再次检查是否解压成功
                if not os.path.exists(os.path.join(path, "matrix.mtx.gz")):
                    raise FileNotFoundError(f"解压后仍然找不到 {path}/matrix.mtx.gz")
                data = read_anndata_BD(path)
                data.var_names_make_unique()
                data.obs['sample'] = name
                data.obs['group'] = k
                adatas_list.append(data)
                if not os.path.exists(os.path.join(out_dir, "01-Rhapsody", name,"filtered_feature_bc_matrix")):
                    os.makedirs(os.path.join(out_dir, "01-Rhapsody", name), exist_ok=True)
                    os.system(f'cp {os.path.join(data_path, name,  name+"_Metrics_Summary.csv")} {os.path.join(out_dir, "01-Rhapsody", name)}')
                    os.system(f'cp {os.path.join(data_path, name, name+"_Pipeline_Report.html")} {os.path.join(out_dir, "01-Rhapsody", name)}')
                    os.system(f'cp -r {os.path.join(data_path, name, "filtered_feature_bc_matrix")} {os.path.join(out_dir, "01-Rhapsody", name)}')
                else:
                    print("File exists")
    adata = adatas_list[0].concatenate(adatas_list[1:], index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    del adata.obs['batch']
    return adata


