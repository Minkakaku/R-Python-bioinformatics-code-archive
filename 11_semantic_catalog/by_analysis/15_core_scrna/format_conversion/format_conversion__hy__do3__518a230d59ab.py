import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import seaborn as sns

# 1. 加载你的h5ad数据
adata = sc.read_h5ad('result.h5ad')  # 请替换为你的文件路径

# 2. 定义你感兴趣的基因集
# 例如，定义一个T细胞毒性相关基因集
gene_set = {"Cytotoxicity":[
    "ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1","BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","H1-0","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6","IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1","TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP"]
}

# 3. 确保使用标准化后的数据（例如log1p后的数据）
# 如果尚未标准化，请取消以下注释
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# 4. 计算AUCell评分
# 注意：AUCell在Python中没有直接的Scanpy集成，但可以手动计算或使用其他包如scampy
# 这里提供一个基于Scanpy内置函数或自定义计算的替代方案

# 由于AUCell在Python中安装可能复杂，我们可以使用Scanpy的score_genes方法
# 使用Scanpy的score_genes方法计算基因集评分[1](@ref)
sc.tl.score_genes(adata, gene_list=gene_set["Cytotoxicity"], score_name='cytotoxicity_score')

# 5. 可视化：按样本分组展示评分
plt.figure(figsize=(10, 6))
# 小提琴图
sc.pl.violin(adata, keys='cytotoxicity_score', groupby='sample', rotation=45)
plt.title('Cytotoxicity Score by Sample')
plt.tight_layout()
plt.savefig('cytotoxicity_score_violin.png', dpi=300)
plt.show()

# 可选：UMAP图上可视化评分
sc.pl.umap(adata, color=['cytotoxicity_score', 'sample'], frameon=False)
plt.savefig('cytotoxicity_score_umap.png', dpi=300)
plt.show()

# 可选：热图可视化（展示基因集中基因的表达）
# 首先提取基因集中存在于数据集中的基因
available_genes = [gene for gene in gene_set["Cytotoxicity"] if gene in adata.var_names]
# 然后可以绘制热图
sc.pl.heatmap(adata, var_names=available_genes, groupby='sample', figsize=(12, 8))
plt.savefig('gene_expression_heatmap.png', dpi=300)
plt.show()
