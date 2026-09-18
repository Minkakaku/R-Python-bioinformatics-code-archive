import scanpy as sc
import pandas as pd
import numpy as np
from typing import Optional, Dict, Any

class DataHandler:
    def __init__(self):
        self.adata = None
        self.qc_metrics = {}
        
    def load_data(self, file_path: str, platform: str = None) -> sc.AnnData:
        """加载单细胞数据"""
        if file_path.endswith('.h5ad'):
            adata = sc.read_h5ad(file_path)
        elif file_path.endswith('.csv'):
            adata = sc.read_csv(file_path)
        elif file_path.endswith('.txt') or file_path.endswith('.tsv'):
            adata = sc.read_text(file_path)
        else:
            raise ValueError(f"不支持的文件格式: {file_path}")
        
        self.adata = adata
        return adata
    
    def quality_control(self, adata: sc.AnnData, min_genes: int = 200, min_cells: int = 3,
                       max_mt_percent: float = 0.05) -> sc.AnnData:
        """质量控制"""
        # 计算质控指标
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        
        # 过滤细胞和基因
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        
        # 过滤线粒体基因比例过高的细胞
        adata = adata[adata.obs.pct_counts_mt < max_mt_percent * 100, :]
        
        return adata
    
    def filter_cells_and_genes(self, adata: sc.AnnData, min_genes: int = 200, max_genes: int = 2500,
                              min_counts: int = 500, max_counts: int = 35000,
                              max_mt_percent: float = 0.05) -> sc.AnnData:
        """过滤细胞和基因"""
        adata = adata.copy()
        
        # 过滤细胞
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_cells(adata, max_genes=max_genes)
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, max_counts=max_counts)
        
        # 过滤线粒体基因比例
        adata = adata[adata.obs.pct_counts_mt < max_mt_percent * 100, :]
        
        # 过滤基因
        sc.pp.filter_genes(adata, min_cells=3)
        
        return adata
    
    def normalize_data(self, adata: sc.AnnData, target_sum: float = 1e4) -> sc.AnnData:
        """数据标准化"""
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=target_sum)
        sc.pp.log1p(adata_norm)
        return adata_norm
    
    def get_qc_metrics(self, adata: sc.AnnData) -> Dict[str, Any]:
        """获取质控指标"""
        metrics = {
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'total_counts': adata.X.sum(),
            'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean(),
            'mean_counts_per_cell': adata.obs['total_counts'].mean(),
            'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
            'median_counts_per_cell': adata.obs['total_counts'].median(),
            'pct_mt': adata.obs['pct_counts_mt'].mean()
        }
        return metrics