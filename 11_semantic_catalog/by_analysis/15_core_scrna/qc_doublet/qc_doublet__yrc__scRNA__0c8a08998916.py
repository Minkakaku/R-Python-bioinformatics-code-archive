# 载入包
import scanpy as sc
import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt
import scipy

# 载入文件
adata = sc.read_10x_mtx(
    './filtered_gene_bc_matrices/hg19/',  # mtx 文件目录
    var_names='gene_symbols',             # 使用 gene_symbols 作为变量名
    cache=True)                           # 写入缓存，可以更快的读取文件

# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# mitochondrial genes, "MT-" for human, "Mt-" for mouse 线粒体
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes 核糖体
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes 血红蛋白
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo','pct_counts_hb'],
             jitter=0.4, multi_panel=True)
plt.savefig("./filtered_gene_bc_matrices/QC.png")


# 获取线粒体基因占比 5% 以下；核糖体基因占比 5% 以上的细胞样本
adata = adata[adata.obs.pct_counts_mt < 5, :]
adata = adata[adata.obs.pct_counts_ribo > 5, :]

#删除数据集中的线粒体基因，血红蛋白基因（慎重）
#mt_genes = adata.var_names.str.startswith('MT-')
#hb_genes = adata.var_names.str.startswith('^HB[^(P)]')
#remove = np.add(mt_genes, hb_genes)
#keep = np.invert(remove)
#adata = adata[:,keep]

# 获取表达基因数在 200-2500 的细胞样本
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.n_genes_by_counts > 200, :]

#原始矩阵
mat = pd.DataFrame(data=adata.X.todense(), index=adata.obs_names, columns=adata.var_names)
mat 

# 标准化，使得不同细胞样本间可比
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 计算
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
 
#保留原始数据
adata.raw = adata

# 获取只有特异性基因的数据集
adata = adata[:, adata.var.highly_variable]

# 回归每个细胞的总计数和表达的线粒体基因的百分比的影响。
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# 将每个基因缩放到单位方差。阈值超过标准偏差 10。
sc.pp.scale(adata, max_value=10)

# PCA 图
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca(adata, color='CST3')
#plt.savefig("./filtered_gene_bc_matrices/CST3.png")
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig("./filtered_gene_bc_matrices/PCA.png")

#降维
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# 聚类
sc.tl.leiden(adata, resolution=1)
# 绘制
sc.pl.umap(adata, color=['leiden'], add_outline=True)
plt.savefig("./filtered_gene_bc_matrices/leiden_umap.png")
adata.write("./filtered_gene_bc_matrices/umap.h5ad")

sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("./filtered_gene_bc_matrices/deg_ttest.png")

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("./filtered_gene_bc_matrices/deg_wilcoxon.png")


marker_genes_dict = {
    'T cells': ['CD3D', 'CD8A', 'CD4'], #T细胞
    'NK cells': ['NKG7', 'KLRB1'], #NK细胞
    'B cells': ['MS4A1', 'CD79A', 'CD79B'], #B细胞
    'Monocytes': ['CD68', 'LYZ'], #单核细胞
    'Dendritic cells': ['FCER1A'], #树突状细胞
}

#stacked-violin plot
sc.pl.stacked_violin(adata, marker_genes_dict, groupby="leiden", swap_axes=False, dendrogram=False)
plt.savefig("./filtered_gene_bc_matrices/marker_stackplot.png",bbox_inches = 'tight')

#matrixplot
sc.pl.matrixplot(
    adata,
    marker_genes_dict,
    "leiden",
    dendrogram=False,
    cmap="Blues",
    standard_scale="var",
    colorbar_title="column scaled\nexpression",
)
plt.savefig("./filtered_gene_bc_matrices/marker_matrixplot.png",bbox_inches = 'tight')

#dotplot
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden', dendrogram=False)
plt.savefig("./filtered_gene_bc_matrices/marker_dotplot.png",bbox_inches = 'tight')

#heatmap
sc.pl.heatmap(adata, marker_genes_dict, groupby="leiden", cmap="viridis", dendrogram=False)
plt.savefig("./filtered_gene_bc_matrices/marker_heatmap1.png",bbox_inches = 'tight')

#violinplot
sc.pl.violin(adata, ['CD3D', 'CD8A', 'CD4', 'NKG7', 'CD79A', 'MS4A1', 'LYZ', 'FCER1A' ], groupby='leiden')
plt.savefig("./filtered_gene_bc_matrices/marker_vlnplot.png",bbox_inches = 'tight')

#featureplot
sc.pl.umap(adata, color=['CD3D', 'CD8A', 'CD4', 'NKG7', 'CD79A', 'MS4A1', 'LYZ', 'FCER1A'], frameon=False, ncols=4)
plt.savefig("./filtered_gene_bc_matrices/marker_featureplot.png",bbox_inches = 'tight')

annotation = {
    '0': 'CD4 T',
    '1': 'CD14+ Mono',
    '2': 'B',
    '3': 'CD8 T',
    '4': 'FCGR3A+ Mono',
    '5': 'NK',
    '6': 'DC'
}
adata.obs['celltype'] = adata.obs['leiden'].map(annotation).astype('category')


sc.pl.umap(adata, color='celltype', legend_loc='on data',title='celltype', add_outline=True,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
plt.savefig("./filtered_gene_bc_matrices/anno_umap.png")

sc.pl.embedding(
    adata,
    basis="X_umap",
    color='celltype',
    title='RNA-seq',
    frameon=False,
    ncols=3,
    return_fig=True,
    show=False,
)
plt.savefig("./filtered_gene_bc_matrices/celltype.png")

# marker气泡图
sc.tl.dendrogram(adata,groupby='celltype')
sc.pl.dotplot(
    adata,
    marker_genes_dict,
    groupby='celltype',
    dendrogram=False,
    color_map="Blues",
    swap_axes=False,
    use_raw=True,
    standard_scale="var",
)
plt.savefig("./filtered_gene_bc_matrices/marker_dotplot.png")

sc.pl.heatmap(adata, marker_genes_dict, groupby='celltype', cmap='viridis', dendrogram=True)
plt.savefig("./filtered_gene_bc_matrices/marker_heatmap.png")

sc.pl.tracksplot(adata, marker_genes_dict, groupby='celltype', dendrogram=True)
plt.savefig("./filtered_gene_bc_matrices/marker_trackplot.png")

# 每个簇前3个高变基因表达情况
#Visualize marker genes using matrixplot
sc.pl.rank_genes_groups_matrixplot(adata, n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled')
#Visualize marker genes using stacked violin plots
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=3, cmap='viridis_r')
#Visualize marker genes using heatmap
sc.pl.rank_genes_groups_heatmap(adata, n_genes=3, use_raw=False, swap_axes=True, vmin=-3, vmax=3, cmap='bwr', layer='scaled', figsize=(10,7), show=False);
#Visualize marker genes using tracksplot
sc.pl.rank_genes_groups_tracksplot(adata, n_genes=3)

adata.raw.to_adata().write("./write/pbmc3k_withoutX.h5ad")
adata.write("./write/pbmc3k.h5ad")

#归一化后data.X数据格式变化，需要转换后才能看矩阵情况
adata = sc.read_h5ad("./write/pbmc3k.h5ad")
adata.X = scipy.sparse.csr_matrix(adata.X)
mat = pd.DataFrame(adata.X.todense(),index = adata.obs_names,columns = adata.var_names)
mat