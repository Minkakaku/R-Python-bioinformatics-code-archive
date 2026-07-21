import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


class Plotter:
    def __init__(self):
        plt.style.use('default')
        plt.rcParams['figure.facecolor'] = 'white'
        
    def plot_qc_violin(self, adata):
        keys = ['n_genes_by_counts', 'total_counts', 'log10GenesPerUMI']
        if 'pct_counts_mt' in adata.obs.columns:
            keys = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'log10GenesPerUMI']
        fig, axes = plt.subplots(1, len(keys), figsize=(15, 4))
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