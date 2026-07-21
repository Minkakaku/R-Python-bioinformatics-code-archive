import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import numpy as np
import scanpy as sc
from pathlib import Path
import matplotlib.cm as cm
from matplotlib.colors import to_hex

def load_and_process_data(h5ad_path, boundaries_path):
    """
    Load and process h5ad file and boundary data, safely handle duplicate cell_id
    """
    print("Loading h5ad file...")
    adata = sc.read_h5ad(h5ad_path)
    
    # Check necessary columns
    required_cols = ['annotation', 'CNs_10', 'sample']
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"Column {col} not found, available columns: " + ", ".join(adata.obs.columns))
    
    # Handle duplicate cell_id
    if 'cell_id' in adata.obs.columns:
        duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep=False)
        duplicate_count = duplicate_mask.sum()
        
        if duplicate_count > 0:
            print(f"Found {duplicate_count} duplicate cell_id, processing...")
            keep_mask = ~adata.obs.duplicated(subset='cell_id', keep='first')
            adata = adata[keep_mask, :].copy()
            print(f"De-duplication completed: {adata.n_obs} cells")
        
        adata.obs.set_index('cell_id', inplace=True)
        print("Using cell_id column as index")
    else:
        print("Warning: cell_id column not found, using default index")
    
    print(f"Cell count: {adata.n_obs}, Gene count: {adata.n_vars}")
    print("CNs_10 distribution:")
    print(adata.obs['CNs_10'].value_counts(dropna=False))
    print("Sample count:", adata.obs['sample'].nunique())
    
    # Read boundary file
    print("Loading cell boundary file...")
    boundaries = pd.read_csv(boundaries_path)
    
    if 'cell_id' not in boundaries.columns:
        raise ValueError("Boundary file must contain cell_id column")
    
    # Determine coordinate columns
    x_col, y_col = None, None
    coordinate_columns = boundaries.columns.tolist()
    coordinate_columns.remove('cell_id')
    
    numeric_columns = boundaries[coordinate_columns].select_dtypes(include=[np.number]).columns
    if len(numeric_columns) >= 2:
        x_col, y_col = numeric_columns[0], numeric_columns[1]
        print(f"Auto-detected coordinate columns: {x_col} and {y_col}")
    else:
        possible_x = ['x', 'X', 'vertex_x', 'coord_x', 'Coordinate_X']
        possible_y = ['y', 'Y', 'vertex_y', 'coord_y', 'Coordinate_Y']
        
        for col in possible_x:
            if col in boundaries.columns:
                x_col = col
                break
        for col in possible_y:
            if col in boundaries.columns:
                y_col = col
                break
    
    if x_col is None or y_col is None:
        raise ValueError("Coordinate columns not found, ensure file contains x, y coordinates")
    
    print(f"Using coordinate columns: {x_col} and {y_col}")
    print(f"Boundary file contains {boundaries['cell_id'].nunique()} unique cell IDs")
    
    return adata, boundaries, x_col, y_col

def get_target_annotations():
    """
    Return target cell type list
    """
    target_annotations = ["Proliferative_Epithelium","Luminal_Epithelium","EMT_Epithelium","ECM_CAF","FGFR2+_Basal_Epithelium","Stress-responsive_glandular_Epithelium","Epi-CAF","Basal_Epithelium","myCAF","EMT_M2"]
    return target_annotations

def create_colors_for_annotations(target_annotations, exclude_colors=None):
    """
    Create color mapping for target cell types, excluding specified colors
    """
    if exclude_colors is None:
        exclude_colors = ["#ADD8E6", "#D3D3D3"]  # Default exclude light blue and light gray
    
    # Use specified color scheme
    color_scheme = [
        "#D32F2F",  # Dark red
        "#FF8F00",  # Dark orange
        "#FFD600",  # Bright yellow
        "#388E3C",  # Dark green
        "#795548",  # Brown
        "#E91E63",  # Pink
        "#1976D2",  # Standard blue
        "#9E9D24",  # Yellow-green
        "#C2185B",  # Rose red
        "#00E5FF",  # Bright cyan
    ]
    
    # Create color dictionary
    color_dict = {}
    for i, annotation in enumerate(target_annotations):
        if i < len(color_scheme):
            color_dict[annotation] = color_scheme[i]
        else:
            # If more annotations than colors, use tab20c
            tab20c = cm.get_cmap('tab20c')
            color_idx = (i - len(color_scheme)) % 20
            color_dict[annotation] = to_hex(tab20c(color_idx))
    
    # Add colors for special categories
    color_dict['others-CN8'] = "#D3D3D3"  # Light gray
    color_dict['others'] = "#ADD8E6"  # Light blue
    
    return color_dict

def process_cell_annotations(adata, cell_ids, sample_name, target_annotations):
    """
    Process cell annotation labels
    """
    processed_annotations = []
    
    for cell_id in cell_ids:
        if cell_id in adata.obs.index:
            # Get original annotation
            orig_annotation = adata.obs.loc[cell_id, 'annotation']
            if isinstance(orig_annotation, pd.Series):
                orig_annotation = orig_annotation.iloc[0]
            
            # Get CNs_10 value
            cns_value = adata.obs.loc[cell_id, 'CNs_10']
            if isinstance(cns_value, pd.Series):
                cns_value = cns_value.iloc[0]
            
            # Process annotation
            if str(cns_value) == 'CN8':
                if orig_annotation in target_annotations:
                    processed_annotations.append(orig_annotation)
                else:
                    processed_annotations.append('others-CN8')
            else:
                processed_annotations.append('others')
        else:
            processed_annotations.append('others')
    
    return processed_annotations

def create_cell_polygons(boundaries, adata, x_col, y_col, sample_name, target_annotations):
    """
    Create polygon objects for cells in specified sample
    """
    print(f"Creating cell polygons for sample {sample_name}...")
    
    # Get cell IDs for current sample
    sample_cell_ids = set(adata.obs[adata.obs['sample'] == sample_name].index)
    
    # Filter boundary data to only include cells from current sample
    sample_boundaries = boundaries[boundaries['cell_id'].isin(sample_cell_ids)]
    
    if len(sample_boundaries) == 0:
        print(f"Warning: No matching boundary data found for sample {sample_name}")
        return [], [], []
    
    # Group boundary points by cell ID
    cell_groups = sample_boundaries.groupby('cell_id')
    
    polygons = []
    poly_cellids = []
    
    matched_count = 0
    for cell_id, group in cell_groups:
        # Check if cell has information
        if cell_id in adata.obs.index:
            # Get coordinate points
            points = group[[x_col, y_col]].values
            
            if len(points) >= 3:  # At least 3 points needed to create polygon
                # Create polygon object
                polygon = patches.Polygon(points, closed=True)
                polygons.append(polygon)
                poly_cellids.append(cell_id)
                matched_count += 1
    
    print(f"Sample {sample_name}: Successfully matched {matched_count} cells")
    
    # Process cell annotation labels
    processed_annotations = process_cell_annotations(adata, poly_cellids, sample_name, target_annotations)
    
    return polygons, processed_annotations, poly_cellids

def plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir):
    """
    Draw spatial in-situ maps for each sample
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Get all samples
    if 'sample' not in adata.obs.columns:
        print("Warning: sample column not found, treating all cells as one sample")
        samples = ['all_samples']
        adata.obs['sample'] = 'all_samples'
    else:
        samples = sorted(adata.obs['sample'].unique())
    
    print(f"Starting to process {len(samples)} samples")
    
    # Get target cell type list
    target_annotations = get_target_annotations()
    print(f"Target cell type list ({len(target_annotations)} types):")
    for i, ann in enumerate(target_annotations, 1):
        print(f"  {i}. {ann}")
    
    # Create color mapping
    color_dict = create_colors_for_annotations(target_annotations)
    
    # Sample-level statistics
    sample_stats = {}
    
    for sample in samples:
        print(f"\n=== Processing sample: {sample} ===")
        
        # Create cell polygons
        polygons, processed_annotations, poly_cellids = create_cell_polygons(
            boundaries, adata, x_col, y_col, sample, target_annotations
        )
        
        if len(polygons) == 0:
            print(f"Sample {sample} has no drawable cells, skipping")
            sample_stats[sample] = {'total_cells': 0, 'cn2_cells': 0}
            continue
        
        # Count each category
        annotation_counts_in_plot = pd.Series(processed_annotations).value_counts()
        
        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(15, 12))
        
        # Create color list
        colors = [color_dict.get(ann, "#808080") for ann in processed_annotations]
        
        # Use PatchCollection to add polygons in batch
        p = PatchCollection(polygons, alpha=0.8, edgecolor='white', linewidth=0.3)
        p.set_facecolor(colors)
        ax.add_collection(p)
        
        # Set figure properties
        ax.set_xlabel('X Coordinate', fontsize=12)
        ax.set_ylabel('Y Coordinate', fontsize=12)
        
        # Get statistics
        cn2_cells = sum(1 for ann in processed_annotations if ann != 'others')
        target_annotation_count = sum(annotation_counts_in_plot.get(ann, 0) for ann in target_annotations)
        others_cn2_count = annotation_counts_in_plot.get('others-CN8', 0)
        others_count = annotation_counts_in_plot.get('others', 0)
        
        title = f'CN8 Cell Type Spatial Distribution - Sample: {sample}\n'
        title += f'CN8 Cells: {cn2_cells} | Target Types: {target_annotation_count} | Others-CN8: {others_cn2_count} | Others: {others_count}'
        ax.set_title(title, fontsize=16, pad=20)
        
        # Auto-adjust axis range
        if len(polygons) > 0:
            all_vertices = np.vstack([poly.get_xy() for poly in polygons])
            x_min, x_max = all_vertices[:, 0].min(), all_vertices[:, 0].max()
            y_min, y_max = all_vertices[:, 1].min(), all_vertices[:, 1].max()
            
            # Add 5% margin
            x_margin = (x_max - x_min) * 0.05
            y_margin = (y_max - y_min) * 0.05
            
            ax.set_xlim(x_min - x_margin, x_max + x_margin)
            ax.set_ylim(y_min - y_margin, y_max + y_margin)
            
            # Critical fix: set equal aspect ratio
            ax.set_aspect('equal')
        
        # Create legend elements
        legend_elements = []
        
        # Add target annotations
        for annotation in target_annotations:
            if annotation in color_dict and annotation in annotation_counts_in_plot:
                count = annotation_counts_in_plot.get(annotation, 0)
                if count > 0:  # Only show types that actually exist
                    legend_elements.append(patches.Patch(
                        color=color_dict[annotation], 
                        label=f'{annotation} ({count})',
                        alpha=0.8
                    ))
        
        # Add others-CN8
        if 'others-CN8' in annotation_counts_in_plot:
            legend_elements.append(patches.Patch(
                color=color_dict['others-CN8'],
                label=f'others-CN8 ({annotation_counts_in_plot.get("others-CN8", 0)})',
                alpha=0.8
            ))
        
        # Add others
        if 'others' in annotation_counts_in_plot:
            legend_elements.append(patches.Patch(
                color=color_dict['others'],
                label=f'others ({annotation_counts_in_plot.get("others", 0)})',
                alpha=0.8
            ))
        
        if len(legend_elements) > 0:
            # Adjust number of columns based on number of legend items
            ncols = 1
            if len(legend_elements) > 10:
                ncols = 2
            
            ax.legend(handles=legend_elements, loc='center left', 
                     bbox_to_anchor=(1.05, 0.5), fontsize=9, framealpha=0.9,
                     ncol=ncols, title="Cell Categories")
        
        # Add statistics information
        stats_text = f'Total cells: {len(polygons)}\nCN8 cells: {cn2_cells}\nTarget types: {target_annotation_count}\nOthers-CN8: {others_cn2_count}\nOthers: {others_count}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # Save image
        sample_filename = f"cn2_selected_spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.png"
        sample_filename_pdf = f"cn2_selected_spatial_map_{str(sample).replace(' ', '_').replace('/', '_')}.pdf"
        output_path = output_dir / sample_filename
        output_path_pdf = output_dir / sample_filename_pdf
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        print(f"Image saved: {output_path}")
        plt.close(fig)
        
        # Record statistics
        sample_stats[sample] = {
            'total_cells': len(polygons),
            'cn2_cells': cn2_cells,
            'target_annotation_count': target_annotation_count,
            'others_cn2_count': others_cn2_count,
            'others_count': others_count,
            'target_annotations': target_annotations,
            'output_path': output_path
        }
    
    # Print overall statistics
    print("\n=== Processing completed ===")
    print("Sample statistics:")
    for sample, stats in sample_stats.items():
        print(f"Sample {sample}: {stats['total_cells']} cells, CN8 cells: {stats['cn2_cells']}, "
              f"Target types: {stats['target_annotation_count']}, Others-CN8: {stats['others_cn2_count']}, "
              f"Others: {stats['others_count']}")
    
    return sample_stats

def create_summary_statistics(adata, sample_stats, output_dir):
    """
    Create summary statistics and charts
    """
    output_dir = Path(output_dir)
    
    # Create summary dataframe
    summary_data = []
    for sample, stats in sample_stats.items():
        if stats['total_cells'] > 0:
            cn2_proportion = stats['cn2_cells'] / stats['total_cells'] * 100
        else:
            cn2_proportion = 0
        
        summary_data.append({
            '样本': sample,
            '总细胞数': stats['total_cells'],
            'CN8细胞数': stats['cn2_cells'],
            'CN8细胞比例': cn2_proportion,
            '目标类型CN8细胞数': stats.get('target_annotation_count', 0),
            'Others-CN8细胞数': stats.get('others_cn2_count', 0),
            'Others细胞数': stats.get('others_count', 0)
        })
    
    if not summary_data:
        print("No valid statistics data")
        return
    
    summary_df = pd.DataFrame(summary_data)
    
    # 1. CN8 cell proportion bar chart
    plt.figure(figsize=(12, 8))
    
    # Sort by CN8 proportion
    summary_df = summary_df.sort_values('CN8细胞比例', ascending=False)
    
    bars = plt.bar(range(len(summary_df)), summary_df['CN8细胞比例'], 
                   color='skyblue', edgecolor='black')
    
    plt.title('CN8细胞比例分布 (按样本, 指定细胞类型分析)', fontsize=16, pad=20)
    plt.xlabel('样本', fontsize=14)
    plt.ylabel('CN8细胞比例 (%)', fontsize=14)
    plt.xticks(range(len(summary_df)), summary_df['样本'], rotation=45, ha='right')
    
    # Add value on each bar
    for i, bar in enumerate(bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'cn2_selected_proportion_summary.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Cell category distribution stacked bar chart
    plt.figure(figsize=(12, 8))
    
    categories = ['目标类型CN8细胞数', 'Others-CN8细胞数', 'Others细胞数']
    colors = ['#2ca02c', '#D3D3D3', '#ADD8E6']  # Green, light gray, light blue
    labels = ['目标类型CN8', 'Others-CN8', 'Others']
    
    bottom = np.zeros(len(summary_df))
    for i, (category, color, label) in enumerate(zip(categories, colors, labels)):
        values = summary_df[category].values
        plt.bar(range(len(summary_df)), values, bottom=bottom, 
                color=color, edgecolor='black', label=label)
        bottom += values
    
    plt.title('细胞类别分布 (按样本, 指定细胞类型分析)', fontsize=16, pad=20)
    plt.xlabel('样本', fontsize=14)
    plt.ylabel('细胞数量', fontsize=14)
    plt.xticks(range(len(summary_df)), summary_df['样本'], rotation=45, ha='right')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'cn2_selected_category_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Target cell type distribution heatmap
    target_annotations = get_target_annotations()
    
    # Collect count of each target cell type in each sample
    heatmap_data = []
    for sample, stats in sample_stats.items():
        sample_data = {'样本': sample}
        
        # Get cell information for current sample
        sample_cell_ids = set(adata.obs[adata.obs['sample'] == sample].index)
        
        for ann in target_annotations:
            # Calculate CN8 cell count of current cell type in current sample
            count = 0
            for cell_id in sample_cell_ids:
                if cell_id in adata.obs.index:
                    cell_ann = adata.obs.loc[cell_id, 'annotation']
                    if isinstance(cell_ann, pd.Series):
                        cell_ann = cell_ann.iloc[0]
                    
                    cns_value = adata.obs.loc[cell_id, 'CNs_10']
                    if isinstance(cns_value, pd.Series):
                        cns_value = cns_value.iloc[0]
                    
                    if str(cns_value) == 'CN8' and cell_ann == ann:
                        count += 1
            
            sample_data[ann] = count
        
        heatmap_data.append(sample_data)
    
    heatmap_df = pd.DataFrame(heatmap_data)
    heatmap_df.set_index('样本', inplace=True)
    
    if not heatmap_df.empty:
        plt.figure(figsize=(12, 8))
        
        # Calculate row sums, sort by sum
        heatmap_df['总细胞数'] = heatmap_df.sum(axis=1)
        heatmap_df = heatmap_df.sort_values('总细胞数', ascending=False)
        heatmap_df = heatmap_df.drop('总细胞数', axis=1)
        
        # Draw heatmap
        im = plt.imshow(heatmap_df.values, cmap='YlOrRd', aspect='auto')
        
        plt.title('各样本中CN8目标细胞类型分布热图', fontsize=16, pad=20)
        plt.xlabel('目标细胞类型', fontsize=14)
        plt.ylabel('样本', fontsize=14)
        plt.xticks(range(len(target_annotations)), target_annotations, rotation=45, ha='right')
        plt.yticks(range(len(heatmap_df)), heatmap_df.index)
        
        # Add color bar
        plt.colorbar(im, label='细胞数量')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'cn2_selected_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save heatmap data
        heatmap_df.to_csv(output_dir / 'cn2_selected_heatmap_data.csv')
    
    # 4. Save summary data
    summary_df.to_csv(output_dir / 'cn2_selected_statistics.csv', index=False)
    
    print(f"\nSummary statistics saved to: {output_dir}")
    print(f"\nCN8 cell statistics summary (specified cell type analysis):")
    print(f"Total samples: {len(summary_df)}")
    print(f"Total CN8 cells: {summary_df['CN8细胞数'].sum()}")
    print(f"Average CN8 proportion: {summary_df['CN8细胞比例'].mean():.2f}%")
    print(f"CN8 proportion range: {summary_df['CN8细胞比例'].min():.2f}% - {summary_df['CN8细胞比例'].max():.2f}%")

def main():
    # File path configuration
    h5ad_path = "adata_with_CNs10_TMA_5.h5ad"  # Replace with your h5ad file path
    boundaries_path = "cell_boundaries.csv"  # Replace with your boundary file path
    output_dir = f"spatial_CN8_selected"
    
    print(f"=== CN8 Specified Cell Type Spatial Distribution Plot Generator ===")
    print(f"h5ad file: {h5ad_path}")
    print(f"Boundary file: {boundaries_path}")
    print(f"Output directory: {output_dir}")
    print("=" * 50)
    
    try:
        # 1. Load and process data
        adata, boundaries, x_col, y_col = load_and_process_data(h5ad_path, boundaries_path)
        
        # 2. Check cell ID matching
        h5ad_cells = set(adata.obs.index)
        boundary_cells = set(boundaries['cell_id'].unique())
        matched_cells = h5ad_cells & boundary_cells
        
        print(f"\n=== Cell ID Matching Check ===")
        print(f"Cells in h5ad file: {len(h5ad_cells)}")
        print(f"Cells in boundary file: {len(boundary_cells)}")
        print(f"Matched cells: {len(matched_cells)}")
        print(f"Matching rate: {len(matched_cells)/len(h5ad_cells)*100:.2f}%")
        
        if len(matched_cells) == 0:
            print("Error: No matching cell IDs found!")
            print("First 5 cell IDs in h5ad:", list(h5ad_cells)[:5])
            print("First 5 cell IDs in boundary file:", list(boundary_cells)[:5])
            return
        
        # 3. Check CN8 cells
        cn2_cells = adata.obs[adata.obs['CNs_10'] == 'CN8'].shape[0]
        print(f"\n=== CN8 Cell Statistics ===")
        print(f"Total CN8 cells: {cn2_cells}")
        print(f"CN8 cell proportion: {cn2_cells/adata.n_obs*100:.2f}%")
        
        if cn2_cells == 0:
            print("Warning: No CN8 cells found!")
        
        # 4. Draw CN8 cell spatial maps for each sample
        sample_stats = plot_sample_spatial_maps(adata, boundaries, x_col, y_col, output_dir)
        
        # 5. Create summary statistics
        if sample_stats:
            create_summary_statistics(adata, sample_stats, output_dir)
        
        print(f"\nAll processing completed! Images saved in: {output_dir}")
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
