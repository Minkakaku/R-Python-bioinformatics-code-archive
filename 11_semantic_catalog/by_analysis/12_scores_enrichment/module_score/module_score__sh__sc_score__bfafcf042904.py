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
import csv

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
    # adata.raw = adata
    adata = adata.raw.to_adata() #需要全基因，不是top
    # adata.raw = adata
    adata.obs_names_make_unique()
    
    print(f"Successfully loaded {h5ad_file_path}")
    print(f"Data shape: {adata.shape}")
    group = adata.obs['group'].cat.categories.to_list()
    
    # Load GMT file for gene sets
    # gene_sets = pd.read_csv(gmt_file_path, sep="\t", header=None, index_col=0)
    rows = []
    with open(gmt_file_path, newline='', encoding='utf-8') as f:
        rdr = csv.reader(f, delimiter='\t')   # 或 '\t'，按实际分隔符改
        for row in rdr:
            rows.append(row)
    # 找到最长行的长度
    max_cols = max(len(r) for r in rows)
    # 统一填充成矩形表（NaN 填充缺失）
    gene_sets = pd.DataFrame(
        [r + [None]*(max_cols - len(r)) for r in rows],
                columns=[f'col{i}' for i in range(max_cols)])
    gene_sets.index = gene_sets['col0']
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
        # term = gene_sets.index[i]
        term = gene_sets.index[i]
        genes = gene_sets.iloc[i, 1:].dropna().to_list()  # Get gene list for the gene set !!!这里应该是从1开始
        
        # Score genes using Scanpy
        sc.tl.score_genes(adata, gene_list=genes, score_name=term)
        # Plot UMAP based on the current gene set score
        fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group) * 6, 6)) if len(group) > 1 else plt.subplots(figsize=(6, 6))
        
        if len(group) > 1:
            axs = tuple(map(lambda x: f'ax{x+1}', range(len(group)))) #加上这两句没问题了
            fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
            # Generate and save UMAP for each group
            for j in range(len(group)):
                tmp = adata[adata.obs["group"] == group[j]]
                sc.pl.umap(tmp, size=20, color=term, title=group[j], ax=axs[j], show=False) #
            save_fig(plt, out_dir=os.path.join(out_dir, 'umap'), file_name=f'umap_{term}_split_groups', fig_size={'w': len(group) * 6, 'h': 6})

        else:
            # Generate UMAP for the entire dataset
            sc.pl.umap(adata, color=term, size=25, show=False)
            save_fig(plt, out_dir=os.path.join(out_dir, 'umap'), file_name=f'umap_{term}', fig_size={'w': 6, 'h': 6})
    # output_data = adata.obs[['cell_type', 'group']]  # 获取cell_type和group列
    output_data = adata.obs[['cell_type', 'group']]  # 获取cell_type和group列
    for i in range(len(gene_sets)):
        term = gene_sets.index[i]
        output_data[term] = adata.obs[term]     # 将每个基因集得分列加入到结果中        
    output_data.to_csv(os.path.join(out_dir, 'sc_score.xls'), sep='\t', index=True)
    print("Gene set scoring and UMAP plotting complete.")
    if len(group) > 1:
        os.system(f"Rscript /mnt/NexData_file/home/sh/Projects/FWSC20252616/01.Mission/script/sc_score_boxplot.R --dat {os.path.join(out_dir, 'sc_score.xls')} --output {os.path.join(out_dir, 'boxplot')}") #函数识别不到term
        os.system(f"Rscript /mnt/NexData_file/home/sh/Projects/FWSC20252616/01.Mission/script/sc_score_violin.R --dat {os.path.join(out_dir, 'sc_score.xls')} --output {os.path.join(out_dir, 'violin')}")
    generate_rmd_report(gene_sets, out_dir)
    os.system(f"Rscript -e \"rmarkdown::render('{os.path.join(out_dir, 'report.rmd')}', output_file='{os.path.join(out_dir, 'report.html')}')\"")
    os.remove(os.path.join(out_dir, 'report.rmd'))  # 可选：渲染后删除中间文件


def generate_rmd_report(gene_sets, out_dir):
    report_path = os.path.join(out_dir, "report.rmd")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("""---
            title: "<center> 基因集打分结果报告 </center>"
            author:
              - name: "<center> 杭州广科安德生物科技有限公司 </center>"
            runtime: static
            output:
              html_document:
                keep_md: false
                toc: yes
                toc_float: yes
                theme: paper
                number_sections: true
            ---\n
            ```{r global_options, include=FALSE}
            knitr::opts_chunk$set(echo = FALSE, out.width = "90%", message = FALSE, warning = FALSE)
            ```\n
            """)
        f.write('<style>body .main-container {max-width: 1500px;margin-left: 100px;margin-right: 100px;}</style>\n\n')
        f.write("## Gene Set Scoring Overview\n\n")

        for i in range(len(gene_sets)):
            term = gene_sets.index[i]

            f.write(f"### {term} 基因集打分结果 {{.tabset .tabset-fade .tabset-pills}}\n\n")
            
            f.write(f"#### UMAP图\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('umap/umap_{term}_split_groups.png')\n")
            f.write("```\n\n")

            f.write(f"#### 组间比较箱型图\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('boxplot/{term}_Group_Comparison_Boxplot.png')\n")
            f.write("```\n\n")

            f.write(f"#### 组间比较小提琴图\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('violin/{term}_Group_Comparison_ViolinPlot.png')\n")
            f.write("```\n\n")


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Perform gene set scoring and UMAP plotting for single-cell RNA-seq data.')
    
    # Add arguments for input files and output directory
    parser.add_argument('--h5ad', type=str, required=True, help='Path to the input H5AD file (scRNA-seq data).')
    parser.add_argument('--gmt', type=str, required=True, help='Path to the input GMT file (gene sets).')
    parser.add_argument('--output', type=str, required=True, help='Path to the output directory.')
    
    # Parse the arguments
    args = parser.parse_args()
    # Perform gene set scoring and UMAP plotting
    sc_score(args.h5ad, args.gmt, args.output)

if __name__ == "__main__":
    main()
