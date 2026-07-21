import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. 加载您的h5ad数据
adata = sc.read_h5ad('result.h5ad')  # 请替换为您的文件路径

# 2. 定义您的基因集（请替换为您图片中实际的基因列表）
# 示例基因集 - 请用您图片中的实际基因替换这些示例基因
gene_set = ["ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1","BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","H1-0","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6","IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1","TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP"
    # 请继续添加图片中的所有基因名称...
    # 确保基因名称与您的单细胞数据中的命名一致
]

# 3. 过滤掉数据中不存在的基因
available_genes = [gene for gene in gene_set if gene in adata.var_names]
missing_genes = [gene for gene in gene_set if gene not in adata.var_names]

print(f"找到 {len(available_genes)} 个基因在数据集中")
print(f"缺失 {len(missing_genes)} 个基因: {missing_genes}")

# 4. 使用scanpy的score_genes方法计算基因集评分
sc.tl.score_genes(adata, gene_list=available_genes, score_name='geneset_score')

# 5. 可视化 - 创建子图展示不同样本的评分分布
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# UMAP图（按样本分组）
sc.pl.umap(adata, color='sample', ax=axes[0], show=False, frameon=False)
axes[0].set_title('UMAP by Sample')

# UMAP图（按基因集评分）
sc.pl.umap(adata, color='geneset_score', ax=axes[1], show=False, frameon=False)
axes[1].set_title('Gene Set Score on UMAP')

plt.tight_layout()
plt.savefig('geneset_score_umap.png', dpi=300, bbox_inches='tight')
plt.show()

# 6. 可选：按样本分组展示评分
# 创建分组 violin plot
sc.pl.violin(adata, keys='geneset_score', groupby='sample', rotation=45)
plt.tight_layout()
plt.savefig('geneset_score_by_sample.png', dpi=300, bbox_inches='tight')
plt.show()

# 7. 保存结果
# 保存添加了评分的数据
adata.write('adata_with_geneset_score.h5ad')

# 保存评分结果
score_df = pd.DataFrame({
    'cell_barcode': adata.obs_names,
    'sample': adata.obs['sample'],
    'geneset_score': adata.obs['geneset_score']
})
score_df.to_csv('geneset_scores_per_cell.csv', index=False)

print("分析完成！结果已保存")
