import os
import shutil
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from collections import OrderedDict
import argparse
import scipy

parser = argparse.ArgumentParser()
parser.add_argument("--data_path", required=True)
parser.add_argument("--sample_file", required=True)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--platform", required=True)
# parser.add_argument("--output_file", required=True)
args = parser.parse_args()

data_path = args.data_path
sample_file = args.sample_file
out_dir = args.out_dir
platform = args.platform
# output_file = args.output_file




# data_path = snakemake.input["data_path"]
# sample_file = snakemake.input["sample_file"]
# out_dir = snakemake.params["out_dir"]
# platform = snakemake.params["platform"]
# output_file = snakemake.output[0]

out_cellranger=None
adata=None

def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()
def read_sample_dict(sample_file):
    df = pd.read_csv(sample_file, sep="\t")
    sample_dict = OrderedDict()
    for _, row in df.iterrows():
        sample_dict.setdefault(row["group"], []).append(row["sample"])
    return sample_dict

def get_10x_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    out_cellranger = os.path.join(out_dir, "01-CellRanger")
    os.makedirs(out_cellranger, exist_ok=True)
    for group, samples in sample_dict.items():
        for sample in samples:
            sample_in = os.path.join(data_path, sample)
            sample_out = os.path.join(out_cellranger, sample)
            os.makedirs(sample_out, exist_ok=True)
            target_matrix = os.path.join(sample_out, "filtered_feature_bc_matrix")
            src_matrix = os.path.join(sample_in, "outs", "filtered_feature_bc_matrix")
            if os.path.exists(target_matrix):
                raise FileExistsError(f"❌ 目录已存在，终止运行：{target_matrix}")
            elif os.path.exists(src_matrix):
                shutil.copytree(src_matrix, target_matrix)
                for f in ["metrics_summary.csv", "web_summary.html"]:
                    src = os.path.join(sample_in, "outs", f)
                    if os.path.exists(src):
                        shutil.copy(src, sample_out)
            data = sc.read_10x_mtx(target_matrix, var_names="gene_symbols", cache=False)
            data.var_names_make_unique()
            data.obs["sample"] = sample
            data.obs["group"] = group
            adatas_list.append(data)

    adata = ad.concat(adatas_list, join="outer", label=None, index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    if "batch" in adata.obs:
        del adata.obs["batch"]

    # 按 sample.txt 顺序排列
    sample_order = [s for samples in sample_dict.values() for s in samples]
    group_order = list(sample_dict.keys())
    adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=sample_order, ordered=True)
    adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=group_order, ordered=True)

    return adata, out_cellranger

# 创建函数读取C4 filter matrix
def read_anndata_C4(path):
    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = scipy.sparse.csr_matrix(mat)
    adata = ad.AnnData(mat ,dtype="float32")
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[0] .values
    var_names = ad.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    adata.obs_names = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.var_names_make_unique()
    return adata

def get_C4_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    out_cellranger = os.path.join(out_dir, "01-Dnbc4tools")
    os.makedirs(out_cellranger, exist_ok=True)
    for group, samples in sample_dict.items():
        for sample in samples:
            sample_in = os.path.join(data_path, sample)
            sample_out = os.path.join(out_cellranger, sample)
            os.makedirs(sample_out, exist_ok=True)
            target_matrix = os.path.join(sample_out, "filter_matrix")
            src_matrix = os.path.join(sample_in, "output", "filter_matrix")
            if os.path.exists(target_matrix):
                raise FileExistsError(f"❌ 目录已存在，终止运行：{target_matrix}")
            elif os.path.exists(src_matrix):
                shutil.copytree(src_matrix, target_matrix)
                for f in ["metrics_summary.xls", f"{sample}_scRNA_report.html"]:
                    src = os.path.join(sample_in, "output", f)
                    if os.path.exists(src):
                        shutil.copy(src, sample_out)
            data = read_anndata_C4(target_matrix)
            data.var_names_make_unique()
            data.obs["sample"] = sample
            data.obs["group"] = group
            adatas_list.append(data)

    adata = ad.concat(adatas_list, join="outer", label=None, index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    if "batch" in adata.obs:
        del adata.obs["batch"]

    # 按 sample.txt 顺序排列
    sample_order = [s for samples in sample_dict.values() for s in samples]
    group_order = list(sample_dict.keys())
    adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=sample_order, ordered=True)
    adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=group_order, ordered=True)

    return adata, out_cellranger


# scanpy直接读取xunyin filter matrix
def get_xunyin_mtx(data_path, sample_dict, out_dir):
    adatas_list = []
    out_cellranger = os.path.join(out_dir, "01-seeksoultools")
    os.makedirs(out_cellranger, exist_ok=True)
    for group, samples in sample_dict.items():
        for sample in samples:
            sample_in = os.path.join(data_path, sample)
            sample_out = os.path.join(out_cellranger, sample)
            os.makedirs(sample_out, exist_ok=True)
            target_matrix = os.path.join(sample_out, "filtered_feature_bc_matrix")
            src_matrix = ''
            src_matrix1 = os.path.join(sample_in, "Outs", "filtered_feature_bc_matrix")
            src_matrix2 = os.path.join(sample_in, sample, "Outs", "filtered_feature_bc_matrix")
            src_matrix = src_matrix1 if os.path.exists(src_matrix1) else src_matrix2
            
            if os.path.exists(target_matrix):
                raise FileExistsError(f"❌ 目录已存在，终止运行：{target_matrix}")
            elif os.path.exists(src_matrix):
                shutil.copytree(src_matrix, target_matrix)
                for f in [f"{sample}_summary.csv", f"{sample}_report.html"]:
                    src = os.path.join(os.path.dirname(src_matrix), f)
                    if os.path.exists(src):
                        shutil.copy(src, sample_out)
            data = sc.read_10x_mtx(target_matrix, var_names="gene_symbols", cache=False)
            data.var_names_make_unique()
            data.obs["sample"] = sample
            data.obs["group"] = group
            adatas_list.append(data)

    adata = ad.concat(adatas_list, join="outer", label=None, index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    if "batch" in adata.obs:
        del adata.obs["batch"]

    # 按 sample.txt 顺序排列
    sample_order = [s for samples in sample_dict.values() for s in samples]
    group_order = list(sample_dict.keys())
    adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=sample_order, ordered=True)
    adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=group_order, ordered=True)

    return adata, out_cellranger


# 创建函数读取BD数据：数字index添加CellLabel_ 标签
def read_anndata_BD(path):

    mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz")).astype("float32")
    mat = mat.transpose()
    mat = scipy.sparse.csr_matrix(mat)
    adata = ad.AnnData(mat, dtype="float32")
    
    genes = pd.read_csv(os.path.join(path, "features.tsv.gz"), header=None, sep='\t')
    var_names = genes[0].values
    var_names = ad.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_symbols'] = genes[0].values
    
    # 读取barcode并添加CellLabel_前缀
    barcode = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"), header=None)[0].values
    adata.obs_names = np.char.add('CellLabel_', barcode.astype(str))
    adata.var_names_make_unique()
    return adata

def get_BD_mtx(data_path, sample_dict, out_dir):

    adatas_list = []
    out_cellranger = os.path.join(out_dir, "01-Rhapsody")
    os.makedirs(out_cellranger, exist_ok=True)
    
    for group, samples in sample_dict.items():
        for sample in samples:
            sample_in = os.path.join(data_path, sample)
            sample_out = os.path.join(out_cellranger, sample)
            os.makedirs(sample_out, exist_ok=True)
            target_matrix = os.path.join(sample_out, "filtered_feature_bc_matrix")
            src_matrix = os.path.join(sample_in,"filtered_feature_bc_matrix")
            if os.path.exists(target_matrix):
                raise FileExistsError(f"❌ 目录已存在，终止运行：{target_matrix}")
            elif os.path.exists(src_matrix):
                shutil.copytree(src_matrix, target_matrix)
                for f in [f"{sample}_Metrics_Summary.csv", f"{sample}_scRNA_report.html"]:
                    src = os.path.join(sample_in, f)
                    if os.path.exists(src):
                        shutil.copy(src, sample_out)
            data = read_anndata_BD(target_matrix)
            data.var_names_make_unique()
            data.obs["sample"] = sample
            data.obs["group"] = group
            adatas_list.append(data)

    adata = ad.concat(adatas_list, join="outer", label=None, index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    if "batch" in adata.obs:
        del adata.obs["batch"]
    
    sample_order = [s for samples in sample_dict.values() for s in samples]
    group_order = list(sample_dict.keys())
    adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=sample_order, ordered=True)
    adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=group_order, ordered=True)
    
    return adata, out_cellranger

def merge_cellranger_results(out_cellranger):
    """
    将 out_cellranger 下各子文件夹中的 .xls 或 .csv 文件合并为一个表，
    并添加 Sample Name 列（即子文件夹名称）。
    输出到 out_cellranger/merged_results.csv。
    """
    all_data = []

    # 遍历所有子文件夹
    for sample_name in os.listdir(out_cellranger):
        sample_path = os.path.join(out_cellranger, sample_name)
        if not os.path.isdir(sample_path):
            continue  # 跳过非文件夹

        # 遍历每个子文件夹中的文件
        for file in os.listdir(sample_path):
            file_path = os.path.join(sample_path, file)
            if file.endswith((".csv", ".xls")):
                try:
                    if file.endswith(".csv"):
                        df = pd.read_csv(file_path)
                    else:
                        df = pd.read_csv(file_path, sep='\t')
                    df["Sample Name"] = sample_name
                    cols = list(df)
                    cols = [cols[-1]] + cols[:-1]
                    df = df[cols]
                    all_data.append(df)
                except Exception as e:
                    print(f"⚠️ 读取失败: {file_path} ({e})")

    if all_data==[]:
        print(f"❌ 没有在 {out_cellranger} 下找到任何 .csv 或 .xls 文件！")
    else:
    # 合并所有表
        merged_df = pd.concat(all_data, ignore_index=True)

        # 输出文件路径
        output_path = os.path.join(out_cellranger, "metrics_summary.csv")
        merged_df.to_csv(output_path, index=False)

def data_qc(adata, out_dir):   
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-') 
    print(adata.var_names)
    adata.var['ribo'] = adata.var_names.str.upper().str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=90, multi_panel=True, show=False)
    axes = plt.gcf().axes
    for ax in axes:
        ax.set_xticklabels(ax.get_xticklabels(), horizontalalignment='right')
    save_fig(plt, out_dir, file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
    save_fig(fig, out_dir, file_name='orig_nCount_nFeature_scatter')
    
    pd.DataFrame(adata.obs).to_csv(os.path.join(out_dir, 'all_cell_features_raw.xls'), sep='\t')

    for s in adata.obs['sample'].drop_duplicates().tolist():
        os.makedirs(os.path.join(out_dir, s), exist_ok=True)
        sc.pl.violin(adata[adata.obs['sample']==s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby='sample', rotation=90, multi_panel=True, show=False)
        save_fig(plt, os.path.join(out_dir, s), file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(adata[adata.obs['sample']==s], x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(adata[adata.obs['sample']==s], x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, os.path.join(out_dir, s), file_name='orig_nCount_nFeature_scatter')
        pd.DataFrame(adata[adata.obs['sample']==s].obs).to_csv(os.path.join(out_dir, s, 'cell_features_raw.xls'), sep='\t')
    
    return adata

os.makedirs(out_dir, exist_ok=True)
sample_dict = read_sample_dict(sample_file)

if platform == "10X":
    adata, out_cellranger = get_10x_mtx(data_path, sample_dict, out_dir)
    print(f"✅ 合并完成，结果保存至: {os.path.join(out_dir,'01-CellRanger','merged.h5ad')}")
elif platform == "华大":
    adata, out_cellranger = get_C4_mtx(data_path, sample_dict, out_dir)
    print(f"✅ 合并完成，结果保存至: {os.path.join(out_dir,'01-CellRanger','merged.h5ad')}")
elif platform == "寻因":
    adata, out_cellranger = get_xunyin_mtx(data_path, sample_dict, out_dir)
    print(f"✅ 合并完成，结果保存至: {os.path.join(out_dir,'01-CellRanger','merged.h5ad')}")
elif platform == "BD":
    adata, out_cellranger = get_BD_mtx(data_path, sample_dict, out_dir)
    print(f"✅ 合并完成，结果保存至: {os.path.join(out_dir,'01-CellRanger','merged.h5ad')}")
else:
    pass

merge_cellranger_results(out_cellranger)
out_qc=os.path.join(out_dir, "02-Data_QC")
os.makedirs(out_qc, exist_ok=True)
adata_q= data_qc(adata, out_qc)
sample_num = adata_q.obs['sample'].drop_duplicates()
cell_num_raw = adata_q.X.shape[0]
Parameter_dict={}
Parameter_dict['sample_num']=sample_num
Parameter_dict['cell_num_raw']=cell_num_raw
csv_path = os.path.join(out_dir, 'Parameter_dict.csv')
pd.DataFrame(Parameter_dict.items(), columns=['Key', 'Value']).to_csv(csv_path, index=False)
adata_q.write_h5ad(os.path.join(out_qc, "merged.h5ad"))
adata_q.write_h5ad(os.path.join(out_qc, "merged.h5ad"))