import scanpy as sc

# 1. 读取 h5ad 文件
adata = sc.read_h5ad('combined_data.h5ad')

# 2. 确认数据中包含 UMAP 坐标（通常存储在 .obsm['X_umap']）
if 'X_umap' not in adata.obsm:
    print("数据中未找到 UMAP 坐标，正在计算 UMAP...")
    sc.pp.neighbors(adata)          # 先计算邻居图
    sc.tl.umap(adata)               # 再计算 UMAP
else:
    print("UMAP 坐标已存在，直接使用。")

# 3. 检查 'celltype' 列是否存在
if 'celltype' not in adata.obs:
    raise KeyError("adata.obs 中未找到 'celltype' 列，请确认列名是否正确。")

# 4. 绘制 UMAP，按 celltype 着色
sc.pl.umap(adata, color='celltype', title='UMAP colored by celltype',save='umap_celltype.pdf')


# 可选：如果你想同时查看 res_0.8 的聚类结果，可以将 color 设为 ['celltype', 'res_0.8']
# sc.pl.umap(adata, color=['celltype', 'res_0.8'], ncols=2)
