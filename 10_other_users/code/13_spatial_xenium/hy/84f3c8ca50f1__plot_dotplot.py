import scanpy as sc
import matplotlib.pyplot as plt
import os

# ==============================
# 1. 读取数据
# ==============================
adata = sc.read_h5ad("data_with_group_final.h5ad")  # 请根据实际路径调整

# ==============================
# 2. 定义标记基因字典（已根据数据过滤）
# ==============================
marker_genes_dict = {
    "B cell": ["CD19", "CD79A", "MS4A1"],
    "Plasma": ["TNFRSF17", "MZB1", "IGJ", "JCHAIN", "XBP1"],
    "T cell": ['CD3G', "CD3E", 'CD3D', 'TRAC', "CD4","FOXP3", "CD8A", "GZMB"],
    "Proliferating": ["MKI67", "TOP2A"],
    "Macrophage": ["CD68", "GNLY", "NKG7", "CD163", 'CSF1R'],
    'Neutrophi': ["CSF3R", "FCGR3B", "S100A8", "S100A9", "G0S2", "RETNLG"],
    "Endothelial": ["PECAM1","AQP1","VWF"],
    "Epithelial": ["EPCAM", "KRT5", "KRT20", 'KRT14'],
    "Fibroblasts": ["PDGFRA", "CCDC80", "SFRP4" ,"ACTA2","COL1A1","DCN", "COL1A2"],
}

# 仅保留数据中实际存在的基因
valid_genes = set(adata.var_names)
marker_genes_dict_filtered = {
    ct: [gene for gene in genes if gene in valid_genes]
    for ct, genes in marker_genes_dict.items()
    if any(gene in valid_genes for gene in genes)
}

# ==============================
# 3. 绘制 DotPlot（按 celltype 分组）
# ==============================
groupby_col = "celltype"  # 确保 adata.obs 中存在该列，否则请修改为实际列名

dp = sc.pl.dotplot(
    adata,
    var_names=marker_genes_dict_filtered,
    groupby=groupby_col,
    use_raw=True,
    dendrogram=False,
    return_fig=True
)

# 美化样式
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)

# ==============================
# 4. 保存图片为 PNG 和 PDF
# ==============================
output_dir = "dotplot_results"
os.makedirs(output_dir, exist_ok=True)

png_path = os.path.join(output_dir, "dotplot_celltype.png")
pdf_path = os.path.join(output_dir, "dotplot_celltype.pdf")

dp.savefig(png_path, dpi=300, bbox_inches='tight')
dp.savefig(pdf_path, bbox_inches='tight')

plt.close('all')
print(f"点图已保存至 {output_dir} 文件夹。")
