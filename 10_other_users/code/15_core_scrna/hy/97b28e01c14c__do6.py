import scanpy as sc
import matplotlib.pyplot as plt
import os
from pathlib import Path

# 1. 加载你的h5ad数据
adata = sc.read_h5ad('result.h5ad')  # 请替换为你的文件路径

# 2. 定义你感兴趣的基因集
gene_list = ['EPCAM','CDH1','CK20','CDX2']  # 请替换为你的目标基因集

# 3. 检查基因是否存在
available_genes = [gene for gene in gene_list if gene in adata.var_names]
missing_genes = [gene for gene in gene_list if gene not in adata.var_names]
print(f"以下基因在数据集中未找到: {missing_genes}")
print(f"将绘制以下基因: {available_genes}")

# 4. 获取样本列表
sample_list = adata.obs['sample'].unique()
print(f"样本列表: {sample_list}")

# 5. 创建输出目录
output_dir = Path("./feature_plots_by_sample")
output_dir.mkdir(parents=True, exist_ok=True)

# 6. 设置美观的绘图参数
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, frameon=False)
# 建议使用viridis等感知均匀的色阶
plt.rcParams['image.cmap'] = 'viridis'

# 7. 循环为每个基因绘图
for gene in available_genes:
    # 创建一个带有两个子图的大画布
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1行，2列
    fig.suptitle(f'Gene Expression: {gene}', fontsize=16, fontweight='bold')  # 总标题

    # 遍历每个样本
    for i, sample in enumerate(sample_list):
        # 提取当前样本的数据
        adata_sample = adata[adata.obs['sample'] == sample].copy()
        
        # 在当前样本的子图上绘制特征图
        sc.pl.umap(adata_sample, color=gene, show=False, ax=axes[i], title=f'{sample}', frameon=False, color_map='viridis')
        # 可以移除坐标轴以更专注于颜色表达
        axes[i].set_xlabel('UMAP1')
        axes[i].set_ylabel('UMAP2')
    
    # 调整布局
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # 为总标题留出空间
    
    # 保存当前基因的图片
    plt.savefig(output_dir / f"{gene}_expression_by_sample.png", bbox_inches='tight')
    plt.close()  # 关闭图形，避免内存泄漏

print(f"所有图像已保存至 {output_dir.absolute()}")
