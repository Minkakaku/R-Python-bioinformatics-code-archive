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

def read_gmt(gmt_file_path):
    data = []
    with open(gmt_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # 跳过格式不正确的行
            pathway_name = parts[0]
            pathway_id = parts[1]
            genes = parts[2:]
            data.append([pathway_name, pathway_id] + genes)
    return pd.DataFrame(data)

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
    # adata = adata[adata.obs.cell_type.isin(['Proliferating Chondrocyte','Hypertrophic Chondrocyte']), :]
    # adata.raw = adata
    adata = adata.raw.to_adata()
    adata.raw = adata
    adata.obs_names_make_unique()
    
    print(f"Successfully loaded {h5ad_file_path}")
    print(f"Data shape: {adata.shape}")
    group = adata.obs['group'].cat.categories.to_list()
    
    # Load GMT file for gene sets
    # gene_sets = pd.read_csv(gmt_file_path, sep="\t", header=None, index_col=0)
    gene_sets = read_gmt(gmt_file_path)
    gene_sets.columns = ['Pathway Name', 'Pathway ID'] + [f'Gene {i+1}' for i in range(len(gene_sets.columns) - 2)]
    gene_sets.set_index('Pathway Name', inplace=True)
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
        
        if len(group) > 1:
            # Generate and save UMAP for each group
            for i in range(len(group)):
                tmp = adata[adata.obs["group"] == group[i]]
                sc.pl.umap(tmp, size=20, color=term, title=group[i], ax=axs[i], show=False)
            save_fig(plt, out_dir=os.path.join(out_dir, 'umap'), file_name=f'umap_{term}_split_groups', fig_size={'w': len(group) * 6, 'h': 6})

        else:
            # Generate UMAP for the entire dataset
            sc.pl.umap(adata, color=term, size=25, show=False)
            save_fig(plt, out_dir=os.path.join(out_dir, 'umap'), file_name=f'umap_{term}', fig_size={'w': 6, 'h': 6})
    output_data = adata.obs[['cell_type', 'group']]  # 获取cell_type和group列
    for i in range(len(gene_sets)):
        term = gene_sets.index[i]
        output_data[term] = adata.obs[term]     # 将每个基因集得分列加入到结果中        
    output_data.to_csv(os.path.join(out_dir, 'sc_score.xls'), sep='\t', index=True)
    print("Gene set scoring and UMAP plotting complete.")
    if len(group) > 1:
        os.system(f"Rscript /mnt/NexData_file/home/sh/Projects/CX2510296062-01/01.Mission/script/sc_score_boxplot.R --dat {os.path.join(out_dir, 'sc_score.xls')} --output {os.path.join(out_dir, 'boxplot')}")
        os.system(f"Rscript /mnt/NexData_file/home/sh/Projects/CX2510296062-01/01.Mission/script/sc_score_violin.R --dat {os.path.join(out_dir, 'sc_score.xls')} --output {os.path.join(out_dir, 'violin')}")
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
                editor_options: 
                  markdown: 
                    wrap: 72
                ---\n
                ```{r global_options, include=FALSE}
                knitr::opts_chunk$set(echo = FALSE, out.width = "90%", message = FALSE, warning = FALSE)
                ```\n
                ```{=html}
                <style type ="text/css">body .main-container {    max-width: 1500px;    margin-left: 100px;    margin-right: 100px;} td { font-size: 15px; text-align: center; font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";} th { font-size: 15px; text-align: center; font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif"; font-weight: 600;} h1 { font-size: 25px; font-weight: 600; color: #377EB8; font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";} h2 { font-size: 22px; font-weight: 600; color: #377EB8; font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";} h3 { font-size: 20px; font-weight: 600; color: #377EB8; font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";} h4 { font-size: 18px; font-weight: 600; color: #377EB8; font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";} h5 {  line-height: 200%;  text-align: center;  font-weight: 600;  font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";} img {    margin-left: auto;    margin-right:auto;    display:block;}p {    line-height: 180%;    font-family: "Times","Microsoft YaHei","黑体","SimHei","宋体","SimSun","sans-serif";    text-align: justify;    font-size: 15px;}
                </style>
                ```\n
                # 基因集打分\n
                ## 目的\n
                  在单细胞RNA测序（scRNA-seq）分析中，基因集打分（Gene Set Scoring）是一种常用的技术，旨在通过评估特定基因集在每个细胞中的活跃程度，进而分析细胞的生物学功能、状态或通路的表达情况。\n
                ## 不同分组间比较" + "{.tabset .tabset-fade .tabset-pills}\n
                  首先，从不同分组基因集打分的UMAP图中可以看到不同分组中该基因集表达的大体差异。
            """)
        for i in range(len(gene_sets)):
            term = gene_sets.index[i]

            f.write(f"### {term} 基因集UMAP结果 " + "" + "{.tabset .tabset-fade .tabset-pills}}\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('umap/umap_{term}_split_groups.png')\n")
            f.write("```\n\n")

        f.write("## 全部细胞中基因集打分组间比较统计检验 " + "{.tabset .tabset-fade .tabset-pills}\n\n")
        for i in range(len(gene_sets)):
            term = gene_sets.index[i]
            f.write(f"### {term}基因集全部细胞中组间比较箱型图 " + "{.tabset .tabset-fade .tabset-pills}\n")
            f.write(f"  {term}基因集全部细胞中不同分组（group）箱型图，使用wilcox检验组间差异。\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('boxplot/{term}_Group_Comparison_Boxplot.png')\n")
            f.write("```\n\n")
            f.write("~注 **横轴:** 不同分组；**纵轴:** 基因集分值；**填充颜色:** 组别；~\n\n")

            f.write(f"### {term}基因集全部细胞中组间比较小提琴图 " + "{.tabset .tabset-fade .tabset-pills}\n")
            f.write(f"  {term}基因集全部细胞中不同分组（group）小提琴图，使用wilcox检验组间差异。\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('violin/{term}_Group_Comparison_ViolinPlot.png')\n")
            f.write("```\n\n")
            f.write("~注 **横轴:** 不同分组；**纵轴:** 基因集分值；**填充颜色:** 组别；~\n\n")
        f.write("## 不同细胞类型中基因集打分组间比较统计检验 " + "{.tabset .tabset-fade .tabset-pills}\n\n")
        for i in range(len(gene_sets)):
            term = gene_sets.index[i]
            f.write(f"### {term}基因集不同细胞类型中组间比较箱型图 " + "{.tabset .tabset-fade .tabset-pills}\n")
            f.write(f"  {term}基因集不同细胞类型中不同分组（group）箱型图，使用wilcox检验组间差异。\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('boxplot/{term}_Scores_Boxplot.png')\n")
            f.write("```\n\n")
            f.write("~注 **横轴:** 不同分组；**纵轴:** 基因集分值；**填充颜色:** 组别；**分面:** 不同细胞类型；~\n\n")

            f.write(f"### {term}基因集不同细胞类型中组间比较小提琴图 " + "{.tabset .tabset-fade .tabset-pills}\n")
            f.write(f"  {term}基因集不同细胞类型中不同分组（group）小提琴图，使用wilcox检验组间差异。\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('violin/{term}_Scores_ViolinPlot.png')\n")
            f.write("```\n\n")
            f.write("~注 **横轴:** 不同分组；**纵轴:** 基因集分值；**填充颜色:** 组别；**分面:** 不同细胞类型；~\n\n")
        f.write("## 全部细胞中基因集打分细胞类型间比较统计检验 " + "{.tabset .tabset-fade .tabset-pills}\n\n")
        for i in range(len(gene_sets)):
            term = gene_sets.index[i]
            f.write(f"### {term}基因集细胞类型间比较箱型图 " + "{.tabset .tabset-fade .tabset-pills}\n")
            f.write(f"  {term}基因集全部细胞中不同细胞类型（cell_type）箱型图，使用wilcox检验细胞类型间差异。\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('boxplot/{term}_CellType_Comparison_Boxplot.png')\n")
            f.write("```\n\n")
            f.write("~注 **横轴:** 不同分组；**纵轴:** 基因集分值；**填充颜色:** 组别；~\n\n")

            f.write(f"### {term}基因集全部细胞中细胞类型间比较小提琴图 " + "{.tabset .tabset-fade .tabset-pills}\n")
            f.write(f"  {term}基因集全部细胞中不同细胞类型（cell_type）小提琴图，使用wilcox检验细胞类型间差异。\n")
            f.write(f"```{{r out.width='70%'}}\n")
            f.write(f"knitr::include_graphics('violin/{term}_CellType_Comparison_ViolinPlot.png')\n")
            f.write("```\n\n")
            f.write("~注 **横轴:** 不同细胞类型；**纵轴:** 基因集分值；**填充颜色:** 细胞类型；~\n\n")
        f.write("## 各组内基因集打分细胞类型间比较统计检验 " + "{.tabset .tabset-fade .tabset-pills}\n\n")
        for i in range(len(gene_sets)):
            term = gene_sets.index[i]
            for j in group:
                f.write(f"### {term}基因集{j}组内细胞类型间比较箱型图 " + "{.tabset .tabset-fade .tabset-pills}\n")
                f.write(f"  {term}基因集{j}组内不同细胞类型（cell_type）箱型图，使用wilcox检验细胞类型间差异。\n")
                f.write(f"```{{r out.width='70%'}}\n")
                f.write(f"knitr::include_graphics('boxplot/{term}_{j}_CellType_Comparison_Boxplot.png')\n")
                f.write("```\n\n")
                f.write("~注 **横轴:** 不同细胞类型；**纵轴:** 基因集分值；**填充颜色:** 细胞类型；~\n\n")
                f.write(f"### {term}基因集{j}组内细胞类型间比较小提琴图 " + "{.tabset .tabset-fade .tabset-pills}\n")
                f.write(f"  {term}基因集{j}组内不同细胞类型（cell_type）小提琴图，使用wilcox检验细胞类型间差异。\n")
                f.write(f"```{{r out.width='70%'}}\n")
                f.write(f"knitr::include_graphics('violin/{term}_{j}_CellType_Comparison_ViolinPlot.png')\n")
                f.write("```\n\n")
                f.write("~注 **横轴:** 不同细胞类型；**纵轴:** 基因集分值；**填充颜色:** 细胞类型；~\n\n")
    


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
