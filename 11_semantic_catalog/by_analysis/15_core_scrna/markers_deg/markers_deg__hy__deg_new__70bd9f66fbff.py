import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# 创建输出目录
os.makedirs("DEG/sample_violins", exist_ok=True)

def deg_analysis(adata, 
                target_genes,
                group_key='group',
                sample_key='sample',  # 新增样本列名参数
                control_group='normal',
                case_group='disease',
                output_dir='DEG'):
    """
    差异表达分析主函数（支持按样本绘制）
    """
    # 修复观测名重复问题
    adata.obs_names_make_unique()
    
    # 验证分组信息
    print("各组样本数量：")
    print(adata.obs[group_key].value_counts())
    
    # 初始化结果存储
    results = []
    
    # 差异分析流程保持不变...
    # [此处保持之前的差异分析代码不变]

    # 绘制按样本的小提琴图
    plot_sample_violins(adata, 
                       target_genes=target_genes,
                       group_key=group_key,
                       sample_key=sample_key,
                       control_group=control_group,
                       case_group=case_group,
                       output_dir=output_dir)

def plot_sample_violins(adata, target_genes, group_key, sample_key, control_group, case_group, output_dir):
    """
    按样本绘制小提琴图的新函数
    """
    # 准备颜色映射
    sample_groups = adata.obs[[sample_key, group_key]].drop_duplicates()
    group_palette = {control_group: '#1f77b4', case_group: '#ff7f0e'}
    
    for gene in target_genes:
        if gene not in adata.var_names:
            continue
            
        plt.figure(figsize=(12, 6))
        
        # 准备绘图数据
        plot_data = adata.obs[[sample_key, group_key]].join(adata[:, gene].to_df())
        
        # 绘制小提琴图
        ax = sns.violinplot(
            x=sample_key,
            y=gene,
            hue=group_key,
            data=plot_data,
            palette=group_palette,
            inner=None,
            cut=0,
            dodge=False  # 关闭自动分组
        )
        
        # 优化图形显示
        plt.title(f"{gene} Expression by Sample")
        plt.ylabel("Normalized expression")
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Sample ID")
        
        # 添加分组指示线
        ymin, ymax = plot_data[gene].agg(['min', 'max'])
        group_boundaries = _find_group_boundaries(plot_data, sample_key, group_key)
        for boundary in group_boundaries:
            ax.axvline(boundary-0.5, color='gray', linestyle='--', alpha=0.6)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/sample_violins/{gene}_by_sample.pdf", dpi=300)
        plt.close()

def _find_group_boundaries(plot_data, sample_key, group_key):
    """
    辅助函数：找到组间分界位置
    """
    sorted_samples = plot_data.sort_values(group_key)[sample_key].unique()
    boundaries = []
    prev_group = None
    
    for i, sample in enumerate(sorted_samples):
        current_group = plot_data[plot_data[sample_key]==sample][group_key].iloc[0]
        if current_group != prev_group and prev_group is not None:
            boundaries.append(i)
        prev_group = current_group
    
    return boundaries

# 其余函数保持不变...
# [保持之前的差异分析和其他绘图函数]

if __name__ == '__main__':
    # 加载数据
    adata = sc.read("adata_celltype.h5ad")
    
    # 目标基因列表
    target_genes = ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS6","IL2","IL4","IL7","IL9","IL15","IL21","TSLP","IL2RA","IL2RB","IL2RG","IL4R","IL7R","IL9R","IL15RA","IL21R","TSLPR","IL3","IL5","GMCSF","IL27A","IL31","OSM","LIF","CNTF","CT1","CT2","CLCF1","IL6","IL11","IL13","IL10","IL19","IL22","IL20","IL24","IL26","IL12A","IL12B","IL23A","EPO","GH","PRL","THPO","GCSF","LEPTIN","Upd","IFNA","IFNB","IFNE","IFNK","IFNL1","IFNL2","IFNL3","IFNW","IFNG","EGF","PDGFA","PDGFB","IL6R","IL6ST","IL11RA","IL13RA1","IL13RA2","IL27RA","IL31RA","OSMR","LIFR","CNTFR","IL10RA","IL10RB","IL20RA","IL20RB","IL22RA1","IL22RA2","IL12RB1","IL12RB2","IL23R","EPOR","GHR","PRLR","TPOR","CSF3R","LEPR","IFNAR1","IFNAR2","IFNLR1","IFNGR1","IFNGR2","EGFR","PDGFRA","PDGFRB","IRF9","CISH","BCL2","MCL1","BCL2L1","PIM1","MYC","CCND1","CCND2","CCND3","CDKN1A","AOX","GFAP","STAM","PTPN2","PTPN6","EP300","PIAS1","PIAS2","PIAS3","PIAS4","FHL1","PTPN11","GRB2","SOS","HRAS","RAF1","PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","PIK3R3","AKT","MTO"]  # 你的基因列表
    
    # 运行分析
    deg_analysis(adata=adata,
                target_genes=target_genes,
                group_key='group',
                sample_key='sample',  # 指定样本列名
                control_group='normal',
                case_group='disease')
