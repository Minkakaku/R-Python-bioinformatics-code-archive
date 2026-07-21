import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from matplotlib.colors import Normalize
import warnings
import os
from matplotlib.colors import LinearSegmentedColormap

# 1. 读取h5ad数据并加载聚类信息
adata = sc.read_h5ad("adata_celltype.h5ad")  # 替换为你的h5ad文件路径
#print(adata)

# 配置输出目录
output_dir = "optimized_dot_plots"
os.makedirs(output_dir, exist_ok=True)

# 自定义颜色映射
colors = [(0.94, 0.94, 0.94), (1, 0, 0)]
cmap = LinearSegmentedColormap.from_list('custom_red', colors, N=256)

# 需要可视化的基因列表
genes = ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS6","IL2","IL4","IL7","IL9","IL15","IL21","TSLP","IL2RA","IL2RB","IL2RG","IL4R","IL7R","IL9R","IL15RA","IL21R","TSLPR","IL3","IL5","GMCSF","IL27A","IL31","OSM","LIF","CNTF","CT1","CT2","CLCF1","IL6","IL11","IL13","IL10","IL19","IL22","IL20","IL24","IL26","IL12A","IL12B","IL23A","EPO","GH","PRL","THPO","GCSF","LEPTIN","Upd","IFNA","IFNB","IFNE","IFNK","IFNL1","IFNL2","IFNL3","IFNW","IFNG","EGF","PDGFA","PDGFB","IL6R","IL6ST","IL11RA","IL13RA1","IL13RA2","IL27RA","IL31RA","OSMR","LIFR","CNTFR","IL10RA","IL10RB","IL20RA","IL20RB","IL22RA1","IL22RA2","IL12RB1","IL12RB2","IL23R","EPOR","GHR","PRLR","TPOR","CSF3R","LEPR","IFNAR1","IFNAR2","IFNLR1","IFNGR1","IFNGR2","EGFR","PDGFRA","PDGFRB","IRF9","CISH","BCL2","MCL1","BCL2L1","PIM1","MYC","CCND1","CCND2","CCND3","CDKN1A","AOX","GFAP","STAM","PTPN2","PTPN6","EP300","PIAS1","PIAS2","PIAS3","PIAS4","FHL1","PTPN11","GRB2","SOS","HRAS","RAF1","PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","PIK3R3","AKT","MTOR"]

# 过滤有效基因
valid_genes = [g for g in genes if g in adata.var_names]
print(f"可绘制基因数量：{len(valid_genes)}")

# 全局绘图参数
plt.rcParams['figure.dpi'] = 600
sc.set_figure_params(dpi=600, frameon=False)

# 获取分组样本信息
disease_samples = sorted(adata.obs[adata.obs['group'] == 'disease']['sample'].unique())
normal_samples = sorted(adata.obs[adata.obs['group'] == 'normal']['sample'].unique())
max_cols = max(len(disease_samples), len(normal_samples)) or 1

for gene in valid_genes:
    # 动态计算画布尺寸
    fig_width = 5 * max_cols+2  # 增加列宽
    fig_height = 14  # 增加总高度
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # 优化布局参数
    gs = fig.add_gridspec(2, max_cols, 
                         hspace=0.4, wspace=0.25,
                         left=0.07, right=0.85,
                         top=0.93, bottom=0.1)
    
    # 计算颜色标准化
    vmax = np.percentile(adata[:, gene].X.toarray(), 98)
    
    # 动态点大小计算函数
    def calculate_dot_size(adata_sub):
        n_cells = adata_sub.n_obs
        spatial_range = np.ptp(adata_sub.obsm['spatial'], axis=0)
        area = spatial_range[0] * spatial_range[1]
        
        # 自适应公式：基础大小与细胞密度成反比
        base_size = 1000 / (np.sqrt(n_cells) + 10)
        # 限制点大小范围
        return np.clip(base_size, 1, 5)  # 最小1px，最大5px

    # 绘制疾病组样本
    for i, sample in enumerate(disease_samples):
        ax = fig.add_subplot(gs[0, i])
        adata_sample = adata[adata.obs['sample'] == sample]
        spatial = adata_sample.obsm['spatial']
        
        # 动态计算点大小
        dot_size = calculate_dot_size(adata_sample)
        
        scat = ax.scatter(
            spatial[:,0], spatial[:,1],
            c=adata_sample[:, gene].X.toarray().flatten(),
            cmap=cmap, 
            s=1,  # 应用动态大小
            vmin=0, 
            vmax=vmax,
            alpha=1,   # 增加透明度
            edgecolors='none',
            rasterized=False  # 启用光栅化加速渲染
        )
        ax.set_title(f"Disease: {sample}\n(n={adata_sample.n_obs} cells)", fontsize=10, pad=8)
        ax.axis('off')

    # 绘制正常组样本
    for i, sample in enumerate(normal_samples):
        ax = fig.add_subplot(gs[1, i])
        adata_sample = adata[adata.obs['sample'] == sample]
        spatial = adata_sample.obsm['spatial']
        
        dot_size = calculate_dot_size(adata_sample)
        
        scat = ax.scatter(
            spatial[:,0], spatial[:,1],
            c=adata_sample[:, gene].X.toarray().flatten(),
            cmap=cmap, 
            s=1,
            vmin=0, 
            vmax=vmax,
            alpha=1,
            edgecolors='none',
            rasterized=False
        )
        ax.set_title(f"Normal: {sample}\n(n={adata_sample.n_obs} cells)", fontsize=10, pad=8)
        ax.axis('off')

    # 添加公共colorbar
    cax = fig.add_axes([0.88, 0.15, 0.015, 0.7])
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(0, vmax), cmap=cmap)
    fig.colorbar(sm, cax=cax, label='Expression Level')

    # 保存图片
    plt.suptitle(f"Gene: {gene} (Dot size adapts to cell density)", y=0.97, fontsize=14)
    plt.savefig(
        os.path.join(output_dir, f"{gene}.png"),
        bbox_inches='tight',
        dpi=300
    )
    plt.close()

print(f"优化后的点图已保存至：{os.path.abspath(output_dir)}")
