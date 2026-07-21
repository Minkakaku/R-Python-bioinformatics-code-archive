import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. 加载您的h5ad数据
adata = sc.read_h5ad('result.h5ad')  # 请替换为您的文件路径

# 2. 从图片中提取的基因集
gene_set = ["ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1","BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","H1-0","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6","IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9","NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA","RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1","TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP"]

# 3. 确保使用raw矩阵
if adata.raw is not None:
    print("使用raw矩阵进行基因评分计算")
    # 创建临时对象使用raw矩阵
    raw_adata = adata.raw.to_adata()
else:
    print("没有raw矩阵，使用当前矩阵")
    raw_adata = adata.copy()

# 4. 过滤掉数据中不存在的基因
available_genes = [gene for gene in gene_set if gene in raw_adata.var_names]
missing_genes = [gene for gene in gene_set if gene not in raw_adata.var_names]

print(f"找到 {len(available_genes)} 个基因在数据集中")
print(f"缺失 {len(missing_genes)} 个基因: {missing_genes}")

# 5. 使用scanpy的score_genes方法计算基因集评分
sc.tl.score_genes(raw_adata, gene_list=available_genes, score_name='geneset_score')

# 6. 将评分结果添加回原始adata对象
adata.obs['geneset_score'] = raw_adata.obs['geneset_score']

# 7. 获取样本列表（假设样本信息存储在adata.obs['sample']中）
samples = adata.obs['sample'].unique()
print(f"找到 {len(samples)} 个样本: {samples}")

# 8. 创建4个子图（按样本拆分UMAP）
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
axes = axes.ravel()  # 将2x2的轴数组展平

# 为每个样本创建一个UMAP图
for i, sample in enumerate(samples[:4]):  # 只取前4个样本
    # 提取当前样本的数据
    sample_mask = adata.obs['sample'] == sample
    sample_adata = adata[sample_mask].copy()
    
    # 绘制当前样本的UMAP图
    sc.pl.umap(sample_adata, color='geneset_score', 
               ax=axes[i], show=False, frameon=False,
               title=f'Sample: {sample}', 
               color_map='viridis')
    
    # 设置坐标轴标签
    axes[i].set_xlabel('UMAP1')
    axes[i].set_ylabel('UMAP2')

# 调整布局并保存
plt.tight_layout()
plt.savefig('geneset_score_umap_by_sample.png', dpi=300, bbox_inches='tight')
plt.show()

# 9. 可选：绘制所有样本在一起的UMAP图
plt.figure(figsize=(10, 8))
sc.pl.umap(adata, color='geneset_score', show=False, frameon=False,
           title='Gene Set Score on UMAP (All Samples)', 
           color_map='viridis')
plt.savefig('geneset_score_umap_all_samples.png', dpi=300, bbox_inches='tight')
plt.show()

# 10. 可选：按样本分组展示评分分布
plt.figure(figsize=(12, 6))
sc.pl.violin(adata, keys='geneset_score', groupby='sample', rotation=45)
plt.tight_layout()
plt.savefig('geneset_score_by_sample.png', dpi=300, bbox_inches='tight')
plt.show()

# 11. 保存结果
# 保存添加了评分的数据
adata.write('adata_with_geneset_score.h5ad')

# 保存评分结果
score_df = pd.DataFrame({
    'cell_barcode': adata.obs_names,
    'sample': adata.obs['sample'],
    'geneset_score': adata.obs['geneset_score']
})
score_df.to_csv('geneset_scores_per_cell.csv', index=False)

print("分析完成！UMAP图已按样本拆分")
