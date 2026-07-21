import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import random
import shutil
import seaborn as sns
from scipy.stats import ttest_ind
from scipy.stats import ranksums
import argparse  # 导入argparse库

# 创建临时目录
def create_output_dir(output):
    out_dir = os.path.join(output)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

# 保存图像的函数
def save_fig(plt, out_dir=os.getcwd(), file_name=None, fig_size={'w': 5, 'h': 5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches='tight', pad_inches=0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches='tight', pad_inches=0.1)
    plt.clf()

def sc_score(h5ad_file_path, gmt_file_path, out_dir='out_dir'):
    """
    Perform gene set scoring and UMAP plotting for the given H5AD and GMT files.
    
    Args:
    - h5ad_file_path (str): Path to the input H5AD file (scRNA-seq data).
    - gmt_file_path (str): Path to the input GMT file (gene sets).
    - out_dir (str): Directory to save results and plots. Defaults to 'out_dir'.
    """
    # Load H5AD data
    adata = sc.read_h5ad(h5ad_file_path)
    adata.raw = adata
    adata = adata.raw.to_adata()
    adata.raw = adata
    adata.obs_names_make_unique()
    
    print(f"Successfully loaded {h5ad_file_path}")
    print(f"Data shape: {adata.shape}")
    group = adata.obs['group'].cat.categories.to_list()
    
    # Load GMT file for gene sets
    gene_sets = pd.read_csv(gmt_file_path, sep="\t", header=None, index_col=0)
    print(f"Successfully loaded {gmt_file_path}")
    print(f"Gene sets:\n{gene_sets.head()}")
    
    # Ensure the output directory exists
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'umap'), exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'boxplot'), exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'violin'), exist_ok=True)

    output_data = pd.DataFrame()
    # For each gene set term, perform scoring and generate UMAP plots
    num_rows = len(gene_sets)
    for i in range(num_rows):
        term = gene_sets.index[i]
        genes = gene_sets.iloc[i, 2:].dropna().to_list()  # Get gene list for the gene set
        
        # Score genes using Scanpy
        sc.tl.score_genes(adata, gene_list=genes, score_name=term)
        # Plot UMAP based on the current gene set score
        fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group) * 6, 6)) if len(group) > 1 else plt.subplots(figsize=(6, 6))
        
        # if len(group) > 1:
        #     # Generate and save UMAP for each group
        #     for i in range(len(group)):
        #         tmp = adata[adata.obs["group"] == group[i]]
        #         sc.pl.umap(tmp, size=20, color=term, title=group[i], ax=axs[i], show=False)
        #     save_fig(plt, out_dir=os.path.join(out_dir, 'umap'), file_name=f'umap_{term}_split_groups', fig_size={'w': len(group) * 6, 'h': 6})

        # else:
            # Generate UMAP for the entire dataset
        sc.pl.umap(adata, color=term, size=25, show=False)
        save_fig(plt, out_dir=os.path.join(out_dir, 'umap'), file_name=f'umap_{term}', fig_size={'w': 6, 'h': 6})
    output_data = adata.obs[['cell_type', 'group']]  # 获取cell_type和group列
    for i in range(len(gene_sets)):
        term = gene_sets.index[i]
        output_data[term] = adata.obs[term]     # 将每个基因集得分列加入到结果中        
    output_data.to_csv(os.path.join(out_dir, 'sc_score.xls'), sep='\t', index=True)
    print("Gene set scoring and UMAP plotting complete.")

    os.system(f"Rscript /mnt/NexData_file/home/xzx/scripts/scscore/script/sc_score_boxplot.R --dat {os.path.join(out_dir, 'sc_score.xls')} --output {os.path.join(out_dir, 'boxplot')}")
    os.system(f"Rscript /mnt/NexData_file/home/xzx/scripts/scscore/script/sc_score_violin.R --dat {os.path.join(out_dir, 'sc_score.xls')} --output {os.path.join(out_dir, 'violin')}")

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Perform gene set scoring and UMAP plotting for single-cell RNA-seq data.')
    
    # Add arguments for input files and output directory
    parser.add_argument('--h5ad', type=str, required=True, help='Path to the input H5AD file (scRNA-seq data).')
    parser.add_argument('--gmt', type=str, required=True, help='Path to the input GMT file (gene sets).')
    parser.add_argument('--output', type=str, required=True, help='Path to the output directory.')
    
    # Parse the arguments
    args = parser.parse_args()
    out_dir=create_output_dir(args.output)
    # Perform gene set scoring and UMAP plotting
    sc_score(args.h5ad, args.gmt, out_dir)

if __name__ == "__main__":
    main()
