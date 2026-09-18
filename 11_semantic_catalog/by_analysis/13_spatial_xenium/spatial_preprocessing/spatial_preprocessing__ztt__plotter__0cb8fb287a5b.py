import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np

class Plotter:
    def __init__(self):
        plt.style.use('default')
        plt.rcParams['figure.facecolor'] = 'white'
        
    def plot_qc_violin(self, adata):
        fig, axes = plt.subplots(1, 4, figsize=(15, 4))      
        keys = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'log10GenesPerUMI']
        for ax, key in zip(axes, keys):
            sc.pl.violin(adata, 
                         key,
                         groupby='sample',
                         jitter=False,
                         rotation=90,
                         ax=ax,
                         show=False)
        plt.tight_layout()
        return fig


    # def plot_expression_distribution(self, adata: sc.AnnData) -> plt.Figure:
    #     """绘制基因表达分布图"""
    #     fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
    #     # 每个细胞的基因数分布
    #     axes[0].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
    #     axes[0].set_xlabel('基因数')
    #     axes[0].set_ylabel('细胞数')
    #     axes[0].set_title('每个细胞的基因数分布')
    #     axes[0].axvline(adata.obs['n_genes_by_counts'].median(), color='red', linestyle='--', label='中位数')
    #     axes[0].legend()
        
    #     # 每个细胞的UMI数分布
    #     axes[1].hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7)
    #     axes[1].set_xlabel('UMI计数')
    #     axes[1].set_ylabel('细胞数')
    #     axes[1].set_title('每个细胞的UMI计数分布')
    #     axes[1].axvline(adata.obs['total_counts'].median(), color='red', linestyle='--', label='中位数')
    #     axes[1].legend()
        
    #     plt.tight_layout()
    #     return fig
    
    # def plot_cell_statistics(self, adata: sc.AnnData) -> plt.Figure:
    #     """绘制细胞统计图"""
    #     fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
    #     # 基因数 vs UMI数
    #     axes[0, 0].scatter(adata.obs['n_genes_by_counts'], adata.obs['total_counts'], alpha=0.5)
    #     axes[0, 0].set_xlabel('基因数')
    #     axes[0, 0].set_ylabel('UMI计数')
    #     axes[0, 0].set_title('基因数 vs UMI计数')
        
    #     # 线粒体基因比例 vs UMI数
    #     axes[0, 1].scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'], alpha=0.5)
    #     axes[0, 1].set_xlabel('UMI计数')
    #     axes[0, 1].set_ylabel('线粒体基因比例 (%)')
    #     axes[0, 1].set_title('UMI计数 vs 线粒体基因比例')
        
    #     # 基因数 vs 线粒体基因比例
    #     axes[1, 0].scatter(adata.obs['n_genes_by_counts'], adata.obs['pct_counts_mt'], alpha=0.5)
    #     axes[1, 0].set_xlabel('基因数')
    #     axes[1, 0].set_ylabel('线粒体基因比例 (%)')
    #     axes[1, 0].set_title('基因数 vs 线粒体基因比例')
        
    #     # 综合散点图
    #     scatter = axes[1, 1].scatter(adata.obs['n_genes_by_counts'], adata.obs['total_counts'], 
    #                                  c=adata.obs['pct_counts_mt'], cmap='viridis', alpha=0.7)
    #     axes[1, 1].set_xlabel('基因数')
    #     axes[1, 1].set_ylabel('UMI计数')
    #     axes[1, 1].set_title('基因数 vs UMI计数 (颜色=线粒体比例)')
    #     plt.colorbar(scatter, ax=axes[1, 1])
        
    #     plt.tight_layout()
    #     return fig
    
    # def plot_qc_violin(self, adata_raw: sc.AnnData, adata_filtered: sc.AnnData, 
    #                   show_before_after: bool = True) -> plt.Figure:
    #     """绘制质控小提琴图"""
    #     if show_before_after:
    #         fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            
    #         # 过滤前
    #         sc.pl.violin(adata_raw, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
    #                     show=False, ax=axes[0, :])
    #         axes[0, 0].set_title('过滤前 - 基因数')
    #         axes[0, 1].set_title('过滤前 - UMI计数')
    #         axes[0, 2].set_title('过滤前 - 线粒体比例')
            
    #         # 过滤后
    #         sc.pl.violin(adata_filtered, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
    #                     show=False, ax=axes[1, :])
    #         axes[1, 0].set_title('过滤后 - 基因数')
    #         axes[1, 1].set_title('过滤后 - UMI计数')
    #         axes[1, 2].set_title('过滤后 - 线粒体比例')
    #     else:
    #         fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    #         sc.pl.violin(adata_filtered, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
    #                     show=False, ax=axes)
    #         axes[0].set_title('基因数分布')
    #         axes[1].set_title('UMI计数分布')
    #         axes[2].set_title('线粒体比例分布')
        
    #     plt.tight_layout()
    #     return fig
    
    # def plot_pca_variance(self, adata: sc.AnnData) -> plt.Figure:
    #     """绘制PCA方差解释率"""
    #     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
    #     # 方差解释率
    #     ax1.bar(range(1, len(adata.uns['pca']['variance_ratio']) + 1), 
    #             adata.uns['pca']['variance_ratio'])
    #     ax1.set_xlabel('主成分')
    #     ax1.set_ylabel('方差解释率')
    #     ax1.set_title('PCA方差解释率')
        
    #     # 累积方差解释率
    #     cumvar = np.cumsum(adata.uns['pca']['variance_ratio'])
    #     ax2.plot(range(1, len(cumvar) + 1), cumvar, 'bo-')
    #     ax2.axhline(0.9, color='red', linestyle='--', label='90%方差')
    #     ax2.set_xlabel('主成分')
    #     ax2.set_ylabel('累积方差解释率')
    #     ax2.set_title('累积方差解释率')
    #     ax2.legend()
        
    #     plt.tight_layout()
    #     return fig
    
    # def plot_umap(self, adata: sc.AnnData, color_by: str = 'leiden') -> plt.Figure:
    #     """绘制UMAP图"""
    #     fig, ax = plt.subplots(figsize=(10, 8))
    #     sc.pl.umap(adata, color=color_by, show=False, ax=ax)
    #     ax.set_title(f'UMAP - 按{color_by}着色')
    #     return fig
    
    # def plot_clusters(self, adata: sc.AnnData) -> plt.Figure:
    #     """绘制聚类结果"""
    #     fig, ax = plt.subplots(figsize=(10, 8))
    #     sc.pl.umap(adata, color='leiden', show=False, ax=ax, legend_loc='on data')
    #     ax.set_title('聚类结果')
    #     return fig
    
    # def plot_annotated_umap(self, adata: sc.AnnData, show_confidence: bool = True) -> plt.Figure:
    #     """绘制带注释的UMAP"""
    #     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        
    #     # 细胞类型注释
    #     sc.pl.umap(adata, color='cell_type', show=False, ax=ax1, legend_loc='right margin')
       