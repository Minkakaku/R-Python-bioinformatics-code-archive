import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

def plot_stacked_bar_chart(adata, sample_col='sample', celltype_col='celltype', 
                          figsize=(18, 12), color_palette='tab20'):
    """
    Plot stacked bar chart of cell types across samples
    
    Parameters:
    ----------
    adata : AnnData
        Single-cell data object
    sample_col : str
        Column name for sample information, default 'sample'
    celltype_col : str
        Column name for cell type information, default 'celltype'
    figsize : tuple
        Figure size, default (18, 12)
    color_palette : str
        Color palette name, default 'tab20'
    """
    
    # 1. Check if required columns exist
    if sample_col not in adata.obs.columns:
        print(f"Error: Sample column '{sample_col}' not found")
        print(f"Available columns: {list(adata.obs.columns)}")
        return None
    
    if celltype_col not in adata.obs.columns:
        print(f"Error: Cell type column '{celltype_col}' not found")
        print(f"Available columns: {list(adata.obs.columns)}")
        return None
    
    # 2. Get sample and cell type information
    samples = adata.obs[sample_col].astype(str)
    celltypes = adata.obs[celltype_col].astype(str)
    
    # 3. Create cross-tabulation table
    cross_tab = pd.crosstab(samples, celltypes)
    
    # Sort samples by total cell count
    cross_tab = cross_tab.loc[cross_tab.sum(axis=1).sort_values(ascending=False).index]
    
    print(f"Number of samples: {cross_tab.shape[0]}")
    print(f"Number of cell types: {cross_tab.shape[1]}")
    print(f"Total cells: {cross_tab.sum().sum():,}")
    
    # 4. Calculate percentages
    cross_tab_percent = cross_tab.div(cross_tab.sum(axis=1), axis=0) * 100
    
    # 5. Create figure
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # 5.1 Stacked bar chart (absolute counts)
    ax1 = axes[0, 0]
    
    # Get colors
    n_celltypes = cross_tab.shape[1]
    if n_celltypes <= 20:
        colors = plt.cm.get_cmap(color_palette)(np.linspace(0, 1, n_celltypes))
    else:
        colors = plt.cm.viridis(np.linspace(0, 1, n_celltypes))
    
    # Plot stacked bar chart
    bottom = np.zeros(len(cross_tab))
    for i, (celltype, counts) in enumerate(cross_tab.items()):
        ax1.bar(cross_tab.index, counts, bottom=bottom, 
                label=celltype, color=colors[i], edgecolor='white', linewidth=0.5)
        bottom += counts
    
    ax1.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Cell Count', fontsize=12, fontweight='bold')
    ax1.set_title('Cell Type Composition by Sample (Absolute Counts)', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # 5.2 Stacked bar chart (percentages)
    ax2 = axes[0, 1]
    
    bottom = np.zeros(len(cross_tab_percent))
    for i, (celltype, percents) in enumerate(cross_tab_percent.items()):
        ax2.bar(cross_tab_percent.index, percents, bottom=bottom, 
                label=celltype, color=colors[i], edgecolor='white', linewidth=0.5)
        bottom += percents
    
    ax2.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cell Percentage (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Cell Type Composition by Sample (Percentages)', fontsize=14, fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 5.3 Bar chart of total cells per sample
    ax3 = axes[1, 0]
    sample_totals = cross_tab.sum(axis=1)
    bars = ax3.bar(range(len(sample_totals)), sample_totals.values, color='steelblue', alpha=0.7)
    ax3.set_xlabel('Sample', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Total Cells', fontsize=12, fontweight='bold')
    ax3.set_title('Total Cell Count per Sample', fontsize=14, fontweight='bold')
    ax3.set_xticks(range(len(sample_totals)))
    ax3.set_xticklabels(sample_totals.index, rotation=45)
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for i, v in enumerate(sample_totals.values):
        ax3.text(i, v + max(sample_totals)*0.01, f'{v:,}', 
                ha='center', va='bottom', fontsize=8, rotation=0)
    
    # 5.4 Bar chart of cell type frequencies
    ax4 = axes[1, 1]
    celltype_totals = cross_tab.sum(axis=0).sort_values(ascending=False)
    celltype_colors = colors[:len(celltype_totals)] if n_celltypes <= 20 else colors
    
    bars = ax4.bar(range(len(celltype_totals)), celltype_totals.values, 
                  color=celltype_colors, alpha=0.7)
    ax4.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Total Cells', fontsize=12, fontweight='bold')
    ax4.set_title('Total Cells per Cell Type (All Samples)', fontsize=14, fontweight='bold')
    ax4.set_xticks(range(len(celltype_totals)))
    ax4.set_xticklabels(celltype_totals.index, rotation=45, ha='right')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Add percentage labels on bars
    total_cells = celltype_totals.sum()
    for i, v in enumerate(celltype_totals.values):
        percentage = v / total_cells * 100
        ax4.text(i, v + max(celltype_totals)*0.01, f'{percentage:.1f}%', 
                ha='center', va='bottom', fontsize=8, rotation=0)
    
    # 6. Add legend
    if n_celltypes <= 15:
        ax1.legend(title='Cell Types', bbox_to_anchor=(1.05, 1), loc='upper left', 
                  fontsize=8, title_fontsize=10)
    else:
        # Create separate legend figure
        fig_legend, ax_legend = plt.subplots(figsize=(8, min(20, n_celltypes*0.5)))
        ax_legend.axis('off')
        
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=colors[i], edgecolor='black', 
                               label=f'{celltype} ({celltype_totals[celltype]:,})')
                         for i, celltype in enumerate(celltype_totals.index)]
        
        ax_legend.legend(handles=legend_elements, loc='center', 
                        fontsize=9, title='Cell Types (Counts)', 
                        title_fontsize=12, ncol=2 if n_celltypes > 30 else 1)
        
        plt.tight_layout()
        fig_legend.savefig('celltype_legend.pdf', dpi=300, bbox_inches='tight')
        fig_legend.savefig('celltype_legend.png', dpi=300, bbox_inches='tight')
        plt.close(fig_legend)
        print(f"Cell type legend saved as: celltype_legend.pdf and celltype_legend.png")
    
    plt.suptitle(f'Cell Type Distribution Across Samples (Total Cells: {total_cells:,})', 
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    return fig, cross_tab, cross_tab_percent

def main():
    # 1. Read h5ad file
    print("Reading combined_data.h5ad file...")
    try:
        adata = sc.read_h5ad("combined_data.h5ad")
        print(f"Successfully loaded data, shape: {adata.shape} (cells: {adata.shape[0]}, genes: {adata.shape[1]})")
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # 2. Check data
    print(f"\nObservation variables (obs) columns:")
    print(list(adata.obs.columns))
    
    # 3. Confirm column names
    sample_col = 'sample'
    celltype_col = 'celltype'
    
    # If column names are not default, let user specify
    if sample_col not in adata.obs.columns:
        print(f"\nWarning: Column '{sample_col}' not found")
        print("Please select sample column from the following:")
        for i, col in enumerate(adata.obs.columns):
            unique_vals = adata.obs[col].nunique()
            print(f"  {i}: {col} ({unique_vals} unique values)")
        
        try:
            sample_idx = int(input("Enter sample column number: "))
            sample_col = adata.obs.columns[sample_idx]
            print(f"Using '{sample_col}' as sample column")
        except:
            print("Invalid input, trying to find possible sample columns...")
            possible_cols = [col for col in adata.obs.columns 
                           if any(keyword in col.lower() 
                                  for keyword in ['sample', 'patient', 'donor', 'batch', 'library'])]
            if possible_cols:
                sample_col = possible_cols[0]
                print(f"Auto-selected '{sample_col}' as sample column")
            else:
                print("Error: Cannot determine sample column")
                return
    
    if celltype_col not in adata.obs.columns:
        print(f"\nWarning: Column '{celltype_col}' not found")
        print("Please select cell type column from the following:")
        for i, col in enumerate(adata.obs.columns):
            unique_vals = adata.obs[col].nunique()
            print(f"  {i}: {col} ({unique_vals} unique values)")
        
        try:
            celltype_idx = int(input("Enter cell type column number: "))
            celltype_col = adata.obs.columns[celltype_idx]
            print(f"Using '{celltype_col}' as cell type column")
        except:
            print("Invalid input, trying to find possible cell type columns...")
            possible_cols = [col for col in adata.obs.columns 
                           if any(keyword in col.lower() 
                                  for keyword in ['celltype', 'cell_type', 'cluster', 'annotation', 'label'])]
            if possible_cols:
                celltype_col = possible_cols[0]
                print(f"Auto-selected '{celltype_col}' as cell type column")
            else:
                print("Error: Cannot determine cell type column")
                return
    
    # 4. Handle missing values
    if adata.obs[sample_col].isnull().any():
        print(f"Warning: '{sample_col}' column has {adata.obs[sample_col].isnull().sum()} missing values")
        adata = adata[~adata.obs[sample_col].isnull()].copy()
        print(f"Removed cells with missing sample information, remaining cells: {adata.shape[0]}")
    
    if adata.obs[celltype_col].isnull().any():
        print(f"Warning: '{celltype_col}' column has {adata.obs[celltype_col].isnull().sum()} missing values")
        adata = adata[~adata.obs[celltype_col].isnull()].copy()
        print(f"Removed cells with missing cell type information, remaining cells: {adata.shape[0]}")
    
    # 5. Display basic information
    print(f"\nSample column: '{sample_col}'")
    print(f"Number of samples: {adata.obs[sample_col].nunique()}")
    
    print(f"\nCell type column: '{celltype_col}'")
    print(f"Number of cell types: {adata.obs[celltype_col].nunique()}")
    
    # 6. Generate stacked bar chart
    print("\nGenerating stacked bar chart...")
    result = plot_stacked_bar_chart(adata, sample_col=sample_col, celltype_col=celltype_col)
    
    if result is not None:
        fig, cross_tab, cross_tab_percent = result
        
        # 7. Save figures in both PNG and PDF formats
        output_png = "sample_celltype_stacked_bar.png"
        output_pdf = "sample_celltype_stacked_bar.pdf"
        
        fig.savefig(output_png, dpi=300, bbox_inches='tight')
        fig.savefig(output_pdf, bbox_inches='tight')
        print(f"\nFigures saved as:")
        print(f"  - {output_png}")
        print(f"  - {output_pdf}")
        
        # 8. Show figure
        plt.show()
        
    else:
        print("Failed to generate figure, please check your data.")

if __name__ == "__main__":
    main()
    print("\nScript execution completed!")
